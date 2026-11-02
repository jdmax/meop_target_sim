"""Build index.html from page_template.html and results.json.

Usage:  python3 build_page.py
The template holds the page; results.json (written by analyze.py) is embedded inline so the page is a
single static file with no fetches.
"""
import json
from pathlib import Path

here = Path(__file__).parent
tpl = (here / "page_template.html").read_text()
data = json.dumps(json.loads((here / "results.json").read_text()), separators=(",", ":"))

head, body = tpl.split("</style>", 1)
body = body.replace("__DATA__", data.replace("</", "<\\/"))
page = (
    "<!doctype html>\n<html lang=\"en\">\n<head>\n"
    "<meta charset=\"utf-8\">\n"
    "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n"
    + head + "</style>\n<style>body { margin: 0; }</style>\n</head>\n<body>\n"
    + body.lstrip() + "\n</body>\n</html>\n"
)
(here / "index.html").write_text(page)
print(f"wrote index.html ({len(page) / 1e3:.0f} kB)")
