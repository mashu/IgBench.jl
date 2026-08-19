# swiftig.jl — SwiftIG CLI annotator (no Julia package dependency).

"""
Annotate by shelling out to a SwiftIG binary.

Binary path: constructor argument, else `ENV["SWIFTIG_BIN"]`.
Package loads without the binary present; `annotate` errors if missing.
"""
struct SwiftIGAnnotator <: AbstractAnnotator
    name::String
    bin::String
    threads::Int
    extra_args::Vector{String}
end

function SwiftIGAnnotator(;
                          name::AbstractString = "swiftig",
                          bin::AbstractString = get(ENV, "SWIFTIG_BIN", "swiftig"),
                          threads::Integer = 1,
                          extra_args = String[])
    SwiftIGAnnotator(String(name), String(bin), Int(threads),
                     String[String(a) for a in extra_args])
end

tool_name(a::SwiftIGAnnotator) = a.name

function resolve_swiftig_bin(a::SwiftIGAnnotator)
    if isfile(a.bin) || ispath(a.bin)
        return a.bin
    end
    # allow bare name on PATH
    path_env = get(ENV, "PATH", "")
    for dir in split(path_env, ':')
        cand = joinpath(dir, a.bin)
        isfile(cand) && return cand
    end
    error("SwiftIG binary not found at '$(a.bin)' (set SWIFTIG_BIN or pass bin=)")
end

function annotate(a::SwiftIGAnnotator, sequences, ids, germline::GermlinePaths; kwargs...)
    n = length(ids)
    length(sequences) == n || error("sequences/ids length mismatch")
    isnothing(germline.d) && error("SwiftIGAnnotator requires a D germline FASTA path")
    bin = resolve_swiftig_bin(a)
    mktempdir() do dir
        query = joinpath(dir, "query.fasta")
        out_tsv = joinpath(dir, "out.tsv")
        open(query, "w") do io
            for (id, seq) in zip(ids, sequences)
                println(io, '>', id)
                println(io, seq)
            end
        end
        cmd = `$bin`
        for arg in a.extra_args
            cmd = `$cmd $arg`
        end
        cmd = `$cmd -query $query -germline_db_V $(germline.v) -germline_db_D $(germline.d) -germline_db_J $(germline.j) -out $out_tsv -outfmt 19 -num_threads $(a.threads)`
        run(cmd)
        parsed = read_airr_calls(out_tsv; max_rows = nothing, skip_nonproductive = false)
        by_id = Dict(r.sequence_id => r for r in parsed)
        out = Vector{CallRecord}(undef, n)
        for i in 1:n
            r = get(by_id, String(ids[i]), nothing)
            if isnothing(r)
                out[i] = CallRecord(ids[i], sequences[i], "", "", "")
            else
                out[i] = r
            end
        end
        out
    end
end
