# msa.jl — Query-anchored overlap MSA for untrimmed query + germline + span clips.
#
# Star MSA: overlap-align the full query to the full germline, then project each
# tool's [start, stop] onto those columns. That is the VDJ special case of POA
# (tools are the same read cut at different V/D/J bounds), so a general graph
# aligner is unnecessary at gallery size.

const ALIGN_MATCH = 2
const ALIGN_MISMATCH = -1
const ALIGN_GAP = -2

function match_line(qaln::AbstractString, saln::AbstractString)
    n = ncodeunits(qaln)
    buf = Vector{UInt8}(undef, n)
    for i in 1:n
        qc = codeunit(qaln, i)
        sc = codeunit(saln, i)
        buf[i] = (qc == sc && qc != UInt8('-')) ? UInt8('|') :
                 (qc == UInt8('-') || sc == UInt8('-')) ? UInt8(' ') : UInt8('.')
    end
    String(buf)
end

"""
Overlap (ends-free) alignment: internal gaps are penalized, end gaps are free.

Unused prefixes/suffixes are still emitted so both sequences appear in full.
Returns `(a_aln, b_aln, match_line)`.
"""
function overlap_align(a::AbstractString, b::AbstractString)
    q = uppercase(String(a))
    s = uppercase(String(b))
    n = ncodeunits(q)
    m = ncodeunits(s)
    n == 0 && return ("-" ^ m, s, "")
    m == 0 && return (q, "-" ^ n, "")

    F = Matrix{Int}(undef, n + 1, m + 1)
    ptr = Matrix{UInt8}(undef, n + 1, m + 1)
    F[1, 1] = 0
    ptr[1, 1] = 0x03
    for j in 1:m
        F[1, j + 1] = 0
        ptr[1, j + 1] = 0x03
    end
    for i in 1:n
        F[i + 1, 1] = 0
        ptr[i + 1, 1] = 0x03
    end
    for i in 1:n
        qi = codeunit(q, i)
        for j in 1:m
            sj = codeunit(s, j)
            sub = qi == sj ? ALIGN_MATCH : ALIGN_MISMATCH
            diag = F[i, j] + sub
            up = F[i, j + 1] + ALIGN_GAP
            left = F[i + 1, j] + ALIGN_GAP
            if diag >= up && diag >= left
                F[i + 1, j + 1] = diag
                ptr[i + 1, j + 1] = 0x00
            elseif up >= left
                F[i + 1, j + 1] = up
                ptr[i + 1, j + 1] = 0x01
            else
                F[i + 1, j + 1] = left
                ptr[i + 1, j + 1] = 0x02
            end
        end
    end
    best = typemin(Int)
    bi = n
    bj = m
    for j in 0:m
        sc = F[n + 1, j + 1]
        if sc >= best
            best = sc
            bi = n
            bj = j
        end
    end
    for i in 0:(n - 1)
        sc = F[i + 1, m + 1]
        if sc >= best
            best = sc
            bi = i
            bj = m
        end
    end
    traceback_overlap(q, s, ptr, bi, bj)
end

function traceback_overlap(q::AbstractString, s::AbstractString, ptr::Matrix{UInt8},
                           i_end::Int, j_end::Int)
    n = ncodeunits(q)
    m = ncodeunits(s)
    i = i_end
    j = j_end
    qa = Char[]
    sa = Char[]
    while i > 0 && j > 0
        p = ptr[i + 1, j + 1]
        p == 0x03 && break
        if p == 0x00
            push!(qa, Char(codeunit(q, i)))
            push!(sa, Char(codeunit(s, j)))
            i -= 1
            j -= 1
        elseif p == 0x01
            push!(qa, Char(codeunit(q, i)))
            push!(sa, '-')
            i -= 1
        else
            push!(qa, '-')
            push!(sa, Char(codeunit(s, j)))
            j -= 1
        end
    end
    reverse!(qa)
    reverse!(sa)
    qleft = i == 0 ? "" : SubString(q, 1, i)
    sleft = j == 0 ? "" : SubString(s, 1, j)
    qright = i_end >= n ? "" : SubString(q, i_end + 1, n)
    sright = j_end >= m ? "" : SubString(s, j_end + 1, m)
    qstr = ("-" ^ ncodeunits(sleft)) * String(qleft) * String(qa) *
           String(qright) * ("-" ^ ncodeunits(sright))
    sstr = String(sleft) * ("-" ^ ncodeunits(qleft)) * String(sa) *
           ("-" ^ ncodeunits(qright)) * String(sright)
    qstr, sstr, match_line(qstr, sstr)
end

"""Keep query bases in `[start, stop]` (1-based, inclusive); other query bases become gaps."""
function project_query_span(query_aln::AbstractString, start::Integer, stop::Integer)
    n = ncodeunits(query_aln)
    buf = Vector{UInt8}(undef, n)
    qpos = 0
    lo = Int(start)
    hi = Int(stop)
    for i in 1:n
        c = codeunit(query_aln, i)
        if c == UInt8('-')
            buf[i] = UInt8('-')
        else
            qpos += 1
            buf[i] = (qpos >= lo && qpos <= hi) ? c : UInt8('-')
        end
    end
    String(buf)
end

"""Project a query alignment row so bases after `stop` (1-based, inclusive) become gaps."""
project_query_stop(query_aln::AbstractString, stop::Integer) =
    project_query_span(query_aln, 1, stop)

"""
MSA rows: untrimmed query, untrimmed germline, then named `[start, stop]` clips.

`trims` is `(name, start, stop, call)` in display order.
"""
function query_germline_msa(query::AbstractString, gl_name::AbstractString,
                            gl_seq::AbstractString, trims)
    qa, ga, _ = overlap_align(query, gl_seq)
    rows = Dict{String,Any}[
        Dict{String,Any}("id" => "query", "label" => "query", "kind" => "query",
                         "seq" => qa),
        Dict{String,Any}("id" => "germline", "label" => "germline", "kind" => "germline",
                         "seq" => ga, "allele" => String(gl_name)),
    ]
    for t in trims
        name = String(t[1])
        start = Int(t[2])
        stop = Int(t[3])
        call = String(t[4])
        push!(rows, Dict{String,Any}(
            "id" => name,
            "label" => name,
            "kind" => "trim",
            "seq" => project_query_span(qa, start, stop),
            "start" => start,
            "stop" => stop,
            "call" => call,
        ))
    end
    Dict{String,Any}("width" => ncodeunits(qa), "rows" => rows)
end
