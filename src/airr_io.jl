# airr_io.jl — Read / write AIRR rearrangement TSV (plain or .gz).

using CodecZlib: GzipDecompressorStream, GzipCompressorStream

const AIRR_COLUMNS = (
    "sequence_id", "sequence",
    "v_call", "d_call", "j_call",
    "v_sequence_start", "v_sequence_end",
    "d_sequence_start", "d_sequence_end",
    "j_sequence_start", "j_sequence_end",
)

function open_airr_read(path::AbstractString)
    endswith(lowercase(path), ".gz") ?
        GzipDecompressorStream(open(path)) : open(path)
end

function open_airr_read(f::Function, path::AbstractString)
    io = open_airr_read(path)
    result = f(io)
    close(io)
    result
end

function open_airr_write(path::AbstractString)
    endswith(lowercase(path), ".gz") ?
        GzipCompressorStream(open(path, "w")) : open(path, "w")
end

function open_airr_write(f::Function, path::AbstractString)
    io = open_airr_write(path)
    result = f(io)
    close(io)
    result
end

function airr_span_field(start_s, end_s)
    (isempty(start_s) || isempty(end_s) || start_s == "NA" || end_s == "NA") &&
        return EMPTY_SPAN
    s = tryparse(Int, start_s)
    e = tryparse(Int, end_s)
    (isnothing(s) || isnothing(e) || e < s) && return EMPTY_SPAN
    Span(s, e)
end

span_start_str(s::Span) = isempty(s) ? "" : string(s.start)
span_stop_str(s::Span) = isempty(s) ? "" : string(s.stop)

"""
    read_airr_calls(path; max_rows, skip_nonproductive) -> Vector{CallRecord}

Stream an AIRR TSV (optionally `.gz`).
"""
function read_airr_calls(path::AbstractString;
                         max_rows::Union{Nothing,Integer} = nothing,
                         skip_nonproductive::Bool = false)
    records = CallRecord[]
    open_airr_read(path) do io
        header = split(readline(io), '\t')
        col = Dict(h => i for (i, h) in enumerate(header))
        for r in ("sequence_id", "sequence")
            haskey(col, r) || error("AIRR file missing column '$r' in $path")
        end
        has_v = haskey(col, "v_call")
        has_d = haskey(col, "d_call")
        has_j = haskey(col, "j_call")
        has_prod = haskey(col, "productive")
        vs = get(col, "v_sequence_start", 0)
        ve = get(col, "v_sequence_end", 0)
        ds = get(col, "d_sequence_start", 0)
        de = get(col, "d_sequence_end", 0)
        js = get(col, "j_sequence_start", 0)
        je = get(col, "j_sequence_end", 0)
        for line in eachline(io)
            isempty(strip(line)) && continue
            cols = split(line, '\t'; limit = length(header))
            length(cols) < length(header) && continue
            getc(name) = cols[col[name]]
            if has_prod && skip_nonproductive
                p = lowercase(getc("productive"))
                p in ("t", "true", "1", "yes") || continue
            end
            seq = getc("sequence")
            isempty(seq) && continue
            d_call = has_d ? getc("d_call") : ""
            d_call = (d_call == "NA") ? "" : d_call
            rec = CallRecord(
                getc("sequence_id"),
                seq,
                has_v ? getc("v_call") : "",
                d_call,
                has_j ? getc("j_call") : "",
                (vs == 0 || ve == 0) ? EMPTY_SPAN : airr_span_field(cols[vs], cols[ve]),
                (ds == 0 || de == 0) ? EMPTY_SPAN : airr_span_field(cols[ds], cols[de]),
                (js == 0 || je == 0) ? EMPTY_SPAN : airr_span_field(cols[js], cols[je]),
            )
            push!(records, rec)
            !isnothing(max_rows) && length(records) >= Int(max_rows) && break
        end
    end
    records
end

"""Write `CallRecord` rows as AIRR TSV (`.gz` if path ends with `.gz`)."""
function write_airr_calls(path::AbstractString, rows::AbstractVector{CallRecord})
    open_airr_write(path) do io
        println(io, join(AIRR_COLUMNS, '\t'))
        for r in rows
            println(io, join((
                r.sequence_id, r.sequence,
                r.v_call, r.d_call, r.j_call,
                span_start_str(r.v_span), span_stop_str(r.v_span),
                span_start_str(r.d_span), span_stop_str(r.d_span),
                span_start_str(r.j_span), span_stop_str(r.j_span),
            ), '\t'))
        end
    end
    path
end
