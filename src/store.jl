# store.jl — Dashboard-ready run artifacts (JSON / JSONL / AIRR).

const SCHEMA_VERSION = 5

"""Persistence backend for a single suite run."""
abstract type AbstractRunStore end

"""Write under a filesystem directory."""
struct DirectoryRunStore <: AbstractRunStore
    root::String
end

DirectoryRunStore(root::AbstractString) = DirectoryRunStore(String(root))

store_root(s::DirectoryRunStore) = s.root

function ensure_store!(s::DirectoryRunStore)
    mkpath(s.root)
    mkpath(joinpath(s.root, "predictions"))
    s
end

"""Replace NaN/Inf with `nothing` so JSON stays portable for dashboards."""
json_sanitize(x::Float64) = isfinite(x) ? x : nothing
json_sanitize(x::AbstractDict) = Dict{String,Any}(String(k) => json_sanitize(v) for (k, v) in x)
json_sanitize(x::AbstractVector) = Any[json_sanitize(v) for v in x]
json_sanitize(x) = x

function write_json(path::AbstractString, obj)
    open(path, "w") do io
        JSON.print(io, json_sanitize(obj), 2)
    end
    path
end

function write_manifest!(s::DirectoryRunStore, obj::AbstractDict)
    ensure_store!(s)
    write_json(joinpath(s.root, "manifest.json"), obj)
end

function write_panels_meta!(s::DirectoryRunStore, obj::AbstractDict)
    ensure_store!(s)
    write_json(joinpath(s.root, "panels.json"), obj)
end

predictions_path(::AbstractRunStore, ::AbstractString, ::AbstractString) = ""

function predictions_path(s::DirectoryRunStore, panel_id::AbstractString,
                          tool::AbstractString)
    safe_panel = replace(String(panel_id), r"[^\w\.-]" => "_")
    safe_tool = replace(String(tool), r"[^\w\.-]" => "_")
    joinpath(s.root, "predictions", "$(safe_panel)__$(safe_tool).airr.tsv.gz")
end

function write_predictions!(s::DirectoryRunStore, panel_id::AbstractString,
                            tool::AbstractString, rows::AbstractVector{CallRecord})
    ensure_store!(s)
    write_airr_calls(predictions_path(s, panel_id, tool), rows)
end

function write_metrics_bundle!(s::DirectoryRunStore, rows::AbstractVector{<:AbstractDict},
                               nested::AbstractDict)
    ensure_store!(s)
    open(joinpath(s.root, "metrics.jsonl"), "w") do io
        for row in rows
            println(io, JSON.json(json_sanitize(row)))
        end
    end
    write_json(joinpath(s.root, "metrics.json"), nested)
end

function write_timing!(s::DirectoryRunStore, obj)
    ensure_store!(s)
    write_json(joinpath(s.root, "timing.json"), obj)
end

"""Prior `(panel, tool)` timing rows from `timing.json`, if present."""
function stored_timing_map(::AbstractRunStore)
    Dict{Tuple{String,String},Dict{String,Any}}()
end

function stored_timing_map(s::DirectoryRunStore)
    path = joinpath(s.root, "timing.json")
    isfile(path) || return Dict{Tuple{String,String},Dict{String,Any}}()
    rows = JSON.parsefile(path)
    out = Dict{Tuple{String,String},Dict{String,Any}}()
    for t in rows
        out[(String(t["panel"]), String(t["tool"]))] = t
    end
    out
end

function write_span_gallery!(s::DirectoryRunStore, obj)
    ensure_store!(s)
    write_json(joinpath(s.root, "span_gallery.json"), obj)
end

function write_summary!(s::DirectoryRunStore, text::AbstractString)
    ensure_store!(s)
    open(joinpath(s.root, "summary.md"), "w") do io
        print(io, text)
    end
end

"""No-op store for in-process diagnostic embedding."""
struct NullRunStore <: AbstractRunStore end
store_root(::NullRunStore) = ""
ensure_store!(s::NullRunStore) = s
write_manifest!(::NullRunStore, _) = nothing
write_panels_meta!(::NullRunStore, _) = nothing
write_predictions!(::NullRunStore, _, _, _) = nothing
write_metrics_bundle!(::NullRunStore, _, _) = nothing
write_timing!(::NullRunStore, _) = nothing
write_span_gallery!(::NullRunStore, _) = nothing
write_summary!(::NullRunStore, _) = nothing
