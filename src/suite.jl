# suite.jl — BenchSuite, compares, run_suite.

"""Which predictions to score against which reference."""
struct CompareSpec
    pred_tool::String
    ref_tool::String
    metrics::Vector{AbstractMetric}
end

CompareSpec(pred::AbstractString, ref::AbstractString,
            metrics::AbstractVector{<:AbstractMetric} = default_metrics()) =
    CompareSpec(String(pred), String(ref), collect(AbstractMetric, metrics))

"""Runnable collection of panels, tools, and compares."""
struct BenchSuite
    name::String
    panels::Vector{AbstractPanel}
    tools::Vector{AbstractAnnotator}
    compares::Vector{CompareSpec}
    timing::TimingSpec
end

function BenchSuite(name::AbstractString;
                    panels,
                    tools,
                    compares::AbstractVector{CompareSpec},
                    timing::TimingSpec = TimingSpec())
    BenchSuite(String(name),
               AbstractPanel[p for p in panels],
               AbstractAnnotator[t for t in tools],
               collect(compares), timing)
end

"""In-memory result of [`run_suite`](@ref)."""
struct BenchResult
    name::String
    mode::String
    step::Union{Nothing,Int}
    tags::Dict{String,Any}
    metrics::Vector{Dict{String,Any}}
    timing::Vector{Dict{String,Any}}
    store_path::String
end

function resolve_reference(ref_tool::AbstractString, panel::PanelData,
                           preds::Dict{String,Vector{CallRecord}})
    if ref_tool == ":gold"
        isnothing(panel.gold) && error("panel $(panel.id) has no gold")
        return panel.gold
    else
        haskey(preds, ref_tool) || error("unknown ref tool '$ref_tool'")
        return preds[ref_tool]
    end
end

function resolve_prediction(pred_tool::AbstractString,
                            preds::Dict{String,Vector{CallRecord}})
    haskey(preds, pred_tool) || error("unknown pred tool '$pred_tool'")
    preds[pred_tool]
end

function tags_to_dict(tags)
    d = Dict{String,Any}()
    for (k, v) in pairs(tags)
        d[String(k)] = v
    end
    d
end

function calls_aligned_to_ids(ids::AbstractVector{<:AbstractString},
                              rows::AbstractVector{CallRecord})
    by = Dict{String,CallRecord}(r.sequence_id => r for r in rows)
    out = Vector{CallRecord}(undef, length(ids))
    for i in eachindex(ids)
        id = String(ids[i])
        haskey(by, id) || error("prediction missing sequence_id $id")
        out[i] = by[id]
    end
    out
end

"""
    run_suite(suite; mode, store, step, tags, cache, reuse_predictions) -> BenchResult

Library entrypoint for standalone scripts and IgFormer embedding.
Pass `store = nothing` / [`NullRunStore`](@ref) to skip disk writes.
Pass `cache = PanelCache(dir)` to freeze panel sequences across diagnostic steps.

Each tool is annotated **once** per panel; that pass is timed and scored.
Pass `reuse_predictions=true` with a [`DirectoryRunStore`](@ref) to load
existing `predictions/<panel>__<tool>.airr.tsv.gz` instead of re-annotating.
"""
function run_suite(suite::BenchSuite;
                   mode::RunMode = FullReportMode(),
                   store = nothing,
                   step = nothing,
                   tags = (;),
                   cache = nothing,
                   reuse_predictions::Bool = false)
    store_obj = isnothing(store) ? NullRunStore() : store
    tag_dict = tags_to_dict(tags)
    !isnothing(step) && (tag_dict["step"] = Int(step))

    panels_meta = Dict{String,Any}()
    metric_rows = Dict{String,Any}[]
    nested = Dict{String,Any}()
    timing_rows = Dict{String,Any}[]
    gallery_rows = Dict{String,Any}[]

    tspec = merge_timing(mode, suite.timing)
    want_gallery = include_span_gallery(mode)
    prior_timing = stored_timing_map(store_obj)

    for panel_spec in suite.panels
        panel = panel_for_mode(panel_spec, mode)
        data = load_panel_cached(panel, cache)
        panels_meta[data.id] = data.meta

        preds = Dict{String,Vector{CallRecord}}()
        for tool in suite.tools
            tname = tool_name(tool)
            pred_file = reuse_predictions ? predictions_path(store_obj, data.id, tname) : ""
            tkey = (data.id, tname)
            if !isempty(pred_file) && isfile(pred_file) && haskey(prior_timing, tkey)
                rows = calls_aligned_to_ids(data.ids, read_airr_calls(pred_file))
                preds[tname] = rows
                if tspec.repeats > 0
                    td = copy(prior_timing[tkey])
                    td["panel"] = data.id
                    td["species"] = data.species
                    merge!(td, tag_dict)
                    push!(timing_rows, td)
                end
                continue
            end
            rows, tr = annotate_timed(tool, data.sequences, data.ids, data.germline, tspec)
            preds[tname] = rows
            if store_predictions(mode)
                write_predictions!(store_obj, data.id, tname, rows)
            end
            if tspec.repeats > 0
                td = timing_to_dict(tr)
                td["panel"] = data.id
                td["species"] = data.species
                merge!(td, tag_dict)
                push!(timing_rows, td)
            end
        end

        nested[data.id] = Dict{String,Any}()
        for cmp in suite.compares
            tools_needed = String[]
            cmp.pred_tool == ":gold" || push!(tools_needed, cmp.pred_tool)
            cmp.ref_tool == ":gold" || push!(tools_needed, cmp.ref_tool)
            all(t -> haskey(preds, t), tools_needed) || continue
            if cmp.ref_tool == ":gold" && isnothing(data.gold)
                continue
            end
            if cmp.pred_tool == ":gold"
                continue
            end

            pref = resolve_prediction(cmp.pred_tool, preds)
            ref = resolve_reference(cmp.ref_tool, data, preds)
            key = "$(cmp.pred_tool)_vs_$(cmp.ref_tool)"
            nested[data.id][key] = Dict{String,Any}()
            for m in cmp.metrics
                mv = evaluate(m, pref, ref)
                row = Dict{String,Any}(
                    "panel" => data.id,
                    "pred" => cmp.pred_tool,
                    "ref" => cmp.ref_tool,
                    "metric" => metric_name(m),
                    "species" => data.species,
                )
                merge!(row, metric_dict(mv))
                merge!(row, tag_dict)
                push!(metric_rows, row)
                nested[data.id][key][metric_name(m)] = metric_dict(mv)
            end
        end
        if want_gallery && !isnothing(data.gold) && !isempty(preds)
            rng = MersenneTwister(hash(data.id))
            tool_order = [tool_name(t) for t in suite.tools]
            append!(gallery_rows, span_gallery(preds, data.gold, data.germline, data.id;
                                               rng, tool_order))
        end
    end

    payload = Dict{String,Any}(
        "schema_version" => SCHEMA_VERSION,
        "suite" => suite.name,
        "mode" => mode_name(mode),
        "step" => step,
        "tags" => tag_dict,
        "tools" => [tool_name(t) for t in suite.tools],
        "panels" => panels_meta,
        "metrics" => metric_rows,
        "timing" => timing_rows,
        "span_gallery" => gallery_rows,
    )
    write_manifest!(store_obj, Dict{String,Any}(
        "schema_version" => payload["schema_version"],
        "suite" => payload["suite"],
        "mode" => payload["mode"],
        "step" => payload["step"],
        "tags" => payload["tags"],
        "tools" => payload["tools"],
    ))
    write_panels_meta!(store_obj, panels_meta)
    write_metrics_bundle!(store_obj, metric_rows, nested)
    write_timing!(store_obj, timing_rows)
    want_gallery && write_span_gallery!(store_obj, gallery_rows)
    write_report_if_full(mode, store_obj, payload)

    BenchResult(suite.name, mode_name(mode),
                isnothing(step) ? nothing : Int(step),
                tag_dict, metric_rows, timing_rows,
                store_root(store_obj))
end
