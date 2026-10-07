#!/usr/bin/env python3
"""Typeset the Markdown review as A4 with native CJK line breaking.

uv run --python /usr/bin/python3 --with matplotlib --with markdown \
    --with weasyprint tools/render_string_document.py

Requires Times New Roman, Noto Serif CJK SC (Song), Noto Sans CJK SC (Hei).
Missing fonts are errors, not silently substituted.
"""

import base64
import html
import io
import re
import subprocess
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import mathtext
from matplotlib.font_manager import FontProperties
import markdown
from weasyprint import HTML

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "docs/string_contraction.md"
OUTPUT = ROOT / "docs/string_contraction_document.pdf"
matplotlib.rcParams["mathtext.fontset"] = "stix"

CSS = """
@page {
    size: A4; margin: 22mm 22mm 20mm;
    @bottom-center {
        content: counter(page);
        font-family: "Times New Roman"; font-size: 9pt; color: black;
    }
}
html { font-family: "Times New Roman", "Noto Serif CJK SC"; }
body { font-size: 12pt; line-height: 1.65; color: black; }
p { margin: 0 0 0.8em; orphans: 2; widows: 2; }
h1, h2, h3 {
    font-family: "Noto Sans CJK SC", sans-serif;
    font-weight: 700; color: black; line-height: 1.4; break-after: avoid;
}
h1 { font-size: 22pt; margin: 0 0 16pt; }
h2 { font-size: 16pt; margin: 19pt 0 10pt; }
h3 { font-size: 14pt; margin: 14pt 0 8pt; }
p, td, th { word-break: normal; overflow-wrap: normal; hyphens: none; }
table { width: 100%; border-collapse: collapse; margin: 10pt 0 14pt; }
thead { display: table-header-group; }
tr { break-inside: avoid; }
th, td {
    border-bottom: 0.5pt solid #ccc; padding: 6pt 7pt;
    text-align: left; vertical-align: top; font-weight: normal;
}
th { background: #f1f1f1; }
a { color: black; text-decoration: none; }
code { font-family: "Times New Roman", "Noto Serif CJK SC"; }
pre {
    font-family: "Times New Roman", "Noto Serif CJK SC";
    font-size: 10.5pt; line-height: 1.5; white-space: pre-wrap;
    padding: 10pt; background: #f5f5f5; break-inside: avoid;
}
.math-inline { vertical-align: -0.15em; max-width: 100%; }
.equation {
    text-align: center; margin: 12pt 0 15pt;
    break-inside: avoid; break-before: avoid;
}
.equation img { max-width: 100%; height: auto; }
.flow { padding: 10pt 12pt; margin: 10pt 0 14pt; background: #f5f5f5;
        font-size: 11pt; break-inside: avoid; }
.flow p { margin: 0.25em 0; }
"""


def require_font(family):
    actual = subprocess.check_output(
        ["fc-match", family, "-f", "%{family}"], text=True).strip()
    if family not in actual.split(","):
        raise RuntimeError(f"Required font {family!r} missing (matched {actual!r})")


def math_source(value):
    value = re.sub(r"\\(mathbb|overline|dot)\s+([A-Za-z])", r"\\\1{\2}", value)
    return value.replace(r"\boxed", "").strip()


def matrix_svg(prefix, body):
    """Render matrix glyphs and brackets as vectors."""
    rows = [row.strip().split("&") for row in body.strip().split(r"\\")]
    height, width = len(rows)*23+10, len(rows[0])*35+65
    fig = plt.figure(figsize=(width/72, height/72))
    ax = fig.add_axes([0, 0, 1, 1], xlim=(0, width), ylim=(0, height))
    ax.axis("off")
    ax.text(0, height/2, "$"+math_source(prefix)+"$", fontsize=13, va="center")
    left, right = 43, width-3
    for r, row in enumerate(rows):
        for c, cell in enumerate(row):
            ax.text(left+18+c*35, height-17-r*23, "$"+cell.strip()+"$",
                    fontsize=13, ha="center", va="center")
    for x, direction in ((left, 1), (right, -1)):
        ax.plot([x+direction*5, x, x, x+direction*5],
                [height-3, height-3, 3, 3], color="black", lw=.8)
    stream = io.BytesIO()
    fig.savefig(stream, format="svg", transparent=True)
    plt.close(fig)
    return stream.getvalue(), width, height


def formula(value, display):
    value = value.strip()
    matrix = re.search(r"(.*?)\\begin\{pmatrix\}(.*?)\\end\{pmatrix\}", value, re.S)
    if matrix:
        svg, width, height = matrix_svg(*matrix.groups())
    else:
        value = "$"+math_source(value.replace("\n", " "))+"$"
        prop = FontProperties(size=13 if display else 12)
        width, height, *_ = mathtext.MathTextParser("path").parse(value, prop=prop)
        stream = io.BytesIO()
        mathtext.math_to_image(value, stream, prop=prop, format="svg")
        svg = stream.getvalue()
    data = base64.b64encode(svg).decode("ascii")
    tag = (f'<img class="{"math-display" if display else "math-inline"}" '
           f'width="{width*96/72:.2f}" height="{height*96/72:.2f}" '
           f'src="data:image/svg+xml;base64,{data}"/>')
    return f'<div class="equation">{tag}</div>' if display else tag


def render_html():
    source, fragments = SOURCE.read_text(), {}

    def protect(value):
        key = f"XMVBFRAGMENT{len(fragments)}END"
        fragments[key] = value
        return key

    def code_block(match):
        language, content = match.groups()
        if language.strip() == "mermaid":
            if "旧路径" in content:
                flow = ("旧路径：共享单自旋核 → pair 组合 → 结构散射",
                        "新路径：共享单自旋核 → 系数块收缩 → 结构矩阵／向量")
            else:
                flow = ("结构向量 → string 系数 → 核作用 → image → 结构投影",)
            block = '<div class="flow">'+"".join(f"<p>{row}</p>" for row in flow)+"</div>"
        else:
            block = "<pre>"+html.escape(content.rstrip())+"</pre>"
        return "\n\n"+protect(block)+"\n\n"

    source = re.sub(r"```([^\n]*)\n(.*?)```", code_block, source, flags=re.S)
    source = re.sub(r"\$\$(.*?)\$\$", lambda m: "\n\n"+protect(
        formula(m[1], True))+"\n\n", source, flags=re.S)
    source = re.sub(r"\$([^$\n]+)\$", lambda m: protect(formula(m[1], False)), source)
    body = markdown.markdown(source, extensions=["tables"])
    for key, value in fragments.items():
        body = body.replace(f"<p>{key}</p>", value).replace(key, value)
    return ('<!doctype html><html lang="zh-CN"><head><meta charset="utf-8">'
            '<title>Unique-spin-string 直接结构收缩</title>'
            f'<style>{CSS}</style></head><body>{body}</body></html>')


def main():
    for family in ("Times New Roman", "Noto Serif CJK SC", "Noto Sans CJK SC"):
        require_font(family)
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    document = HTML(string=render_html(), base_url=str(SOURCE.parent)).render()
    # Regression for the premature Latin/CJK fragment break in the old PDF.
    if "spin-string 系数块" in SOURCE.read_text():
        lines = ["".join(getattr(child, "text", "") for child in box.descendants())
                 for page in document.pages for box in page._page_box.descendants()
                 if type(box).__name__ == "LineBox"]
        if not any("spin-string 系数块" in line for line in lines):
            raise RuntimeError("Latin/CJK line-break regression: spin-string 系数块")
    document.write_pdf(OUTPUT)
    print(f"{OUTPUT} ({len(document.pages)} pages)")


if __name__ == "__main__":
    main()
