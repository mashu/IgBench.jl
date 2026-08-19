# report_html.jl — Self-contained report.html for full reports.

const REPORT_TEMPLATE_PATH = joinpath(@__DIR__, "assets", "report.html")
include_dependency(REPORT_TEMPLATE_PATH)
const REPORT_DATA_TOKEN = "__IGBENCH_DATA__"

"""Serialize run payload into the HTML template (JSON inlined, `<` escaped)."""
function html_report(payload::AbstractDict)
    json = replace(JSON.json(json_sanitize(payload)), '<' => "\\u003c")
    replace(read(REPORT_TEMPLATE_PATH, String), REPORT_DATA_TOKEN => json)
end

function write_html_report!(s::DirectoryRunStore, html::AbstractString)
    ensure_store!(s)
    open(joinpath(s.root, "report.html"), "w") do io
        print(io, html)
    end
end

write_html_report!(::NullRunStore, _) = nothing

write_report_if_full(::DiagnosticMode, store, payload) = nothing

function write_report_if_full(::FullReportMode, store, payload)
    write_summary!(store, format_summary(String(payload["suite"]),
                                         payload["metrics"], payload["timing"]))
    write_html_report!(store, html_report(payload))
end
