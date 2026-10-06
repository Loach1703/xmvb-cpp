#!/usr/bin/env python3
"""Render the continuous string-contraction document as an A4 PDF.

uv run --python /usr/bin/python3 --with matplotlib --with reportlab \
    --with svglib tools/render_string_document.py
"""

import html
import io
import re
import tempfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
from matplotlib import mathtext
from matplotlib.font_manager import FontProperties
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.pagesizes import A4
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.graphics.shapes import Drawing, Line
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, KeepTogether,
)
from svglib.svglib import svg2rlg

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "docs/string_contraction.md"
OUTPUT = ROOT / "docs/string_contraction_document.pdf"
WIDTH = A4[0] - 100
INK = colors.HexColor("#243443")
BLUE = colors.HexColor("#245a80")
pdfmetrics.registerFont(TTFont("Latin", "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"))
pdfmetrics.registerFont(TTFont("CJK", "/usr/share/fonts/truetype/droid/DroidSansFallbackFull.ttf"))
STYLES = {
    "body": ParagraphStyle("body", fontName="Latin", fontSize=10.5,
                           leading=18, textColor=INK, spaceAfter=9),
    "title": ParagraphStyle("title", fontName="Latin", fontSize=23,
                            leading=33, textColor=BLUE, spaceAfter=18),
    "h2": ParagraphStyle("h2", fontName="Latin", fontSize=15,
                         leading=23, textColor=BLUE, spaceBefore=18,
                         spaceAfter=9, keepWithNext=True),
    "h3": ParagraphStyle("h3", fontName="Latin", fontSize=12,
                         leading=20, textColor=BLUE, spaceBefore=12,
                         spaceAfter=7, keepWithNext=True),
    "table": ParagraphStyle("table", fontName="Latin", fontSize=9,
                            leading=15, textColor=INK),
    "code": ParagraphStyle("code", fontName="Latin", fontSize=9,
                           leading=15, textColor=INK),
}
matplotlib.rcParams["mathtext.fontset"] = "stix"


def plain(value):
    """Use the CJK font only for CJK characters; preserve Latin coverage."""
    return re.sub(r"[\u2e80-\ua4cf\uff00-\uffef]+",
                  lambda m: '<font name="CJK">' + m[0] + '</font>',
                  html.escape(value))


def math_source(value):
    value = re.sub(r"\\(mathbb|overline|dot)\s+([A-Za-z])", r"\\\1{\2}", value)
    return value.replace(r"\boxed", "").strip()


def main():
    story = []
    with tempfile.TemporaryDirectory(prefix="string-document-") as workspace:
        inline_id = 0

        def inline(value, size=10.5):
            nonlocal inline_id
            value = re.sub(r"\[([^]]+)\]\([^)]+\)", r"\1", value)
            parts = re.split(r"(\$[^$]+\$|`[^`]+`)", value)
            output = []
            for part in parts:
                if part.startswith("$") and part.endswith("$"):
                    inline_id += 1
                    path = Path(workspace) / f"math-{inline_id}.png"
                    mathtext.math_to_image("$"+math_source(part[1:-1])+"$", path, dpi=200,
                                          prop=FontProperties(size=size), format="png")
                    from PIL import Image
                    with Image.open(path) as bitmap:
                        w, h = bitmap.size
                    output.append(f'<img src="{path}" width="{w*72/200:.2f}" '
                                  f'height="{h*72/200:.2f}" valign="middle"/>')
                elif part.startswith("`"):
                    output.append('<font color="#245a80">' + plain(part[1:-1]) + '</font>')
                else:
                    output.append(plain(part.replace("**", "")))
            return "".join(output)

        def para(value, kind="body"):
            return Paragraph(inline(value, STYLES[kind].fontSize), STYLES[kind])

        def equation(value):
            matrix = re.search(r"(.*?)\\begin\{pmatrix\}(.*?)\\end\{pmatrix\}(.*)",
                               value, re.S)
            if matrix:
                prefix, body, suffix = matrix.groups()
                rows = [[cell.strip() for cell in row.split("&")]
                        for row in body.strip().split(r"\\")]
                height, matrix_width = len(rows)*23, len(rows[0])*38
                label = svg_math(prefix.strip())
                offset = label.width + 14
                drawing = Drawing(offset+matrix_width+12, height+8)
                label.translate(0, height/2-label.height/2)
                drawing.add(label)
                for r, row in enumerate(rows):
                    for c, cell in enumerate(row):
                        glyph = svg_math(cell)
                        glyph.translate(offset+10+c*38, height-(r+1)*23+5)
                        drawing.add(glyph)
                for x, direction in ((offset, 1), (offset+matrix_width+4, -1)):
                    drawing.add(Line(x, 2, x, height+2))
                    drawing.add(Line(x, 2, x+direction*6, 2))
                    drawing.add(Line(x, height+2, x+direction*6, height+2))
                drawing.hAlign = "CENTER"
                append_equation(drawing)
                return

            drawing = svg_math(value.replace("\n", " "))
            scale = min(1, WIDTH / drawing.width)
            drawing.scale(scale, scale)
            drawing.width *= scale
            drawing.height *= scale
            drawing.hAlign = "CENTER"
            append_equation(drawing)

        def svg_math(value):
            stream = io.BytesIO()
            mathtext.math_to_image("$" + math_source(value) + "$",
                                  stream, prop=FontProperties(size=13), format="svg")
            stream.seek(0)
            return svg2rlg(stream)

        def append_equation(drawing):
            block = [Spacer(1, 6), drawing, Spacer(1, 12)]
            if story and isinstance(story[-1], Paragraph):
                block.insert(0, story.pop())
            story.append(KeepTogether(block))

        def add_table(rows, widths=None, header=True):
            table = Table(rows, colWidths=widths or [WIDTH / len(rows[0])] * len(rows[0]),
                          repeatRows=1 if header else 0, hAlign="LEFT")
            rules = [("VALIGN", (0, 0), (-1, -1), "TOP"),
                     ("LEFTPADDING", (0, 0), (-1, -1), 7),
                     ("RIGHTPADDING", (0, 0), (-1, -1), 7),
                     ("TOPPADDING", (0, 0), (-1, -1), 6),
                     ("BOTTOMPADDING", (0, 0), (-1, -1), 6),
                     ("LINEBELOW", (0, 0), (-1, -1), .3, colors.HexColor("#dce3e9"))]
            if header:
                rules.append(("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#eaf0f5")))
            table.setStyle(TableStyle(rules))
            story.extend([table, Spacer(1, 12)])

        lines = SOURCE.read_text().splitlines()
        i = 0
        while i < len(lines):
            line = lines[i].strip()
            i += 1
            if not line:
                continue
            if line == "$$":
                block = []
                while i < len(lines) and lines[i].strip() != "$$":
                    block.append(lines[i])
                    i += 1
                i += 1
                equation("\n".join(block))
            elif line.startswith("```"):
                block = []
                while i < len(lines) and not lines[i].startswith("```"):
                    block.append(lines[i])
                    i += 1
                i += 1
                if "mermaid" in line:
                    diagram = ("共享单自旋核 → pair 组合 → 结构散射\n"
                               "共享单自旋核 → 系数块收缩 → 结构矩阵／向量") if any(
                                   "旧路径" in v for v in block) else (
                               "结构向量 → string 系数 → 核作用 → image → 结构投影")
                    cells = [[para(row, "code")] for row in diagram.splitlines()]
                else:
                    cells = [[Paragraph(plain(row).replace(" ", "&#160;"), STYLES["code"])]
                             for row in block]
                panel = Table(cells, colWidths=[WIDTH])
                panel.setStyle(TableStyle([
                    ("BACKGROUND", (0, 0), (-1, -1), colors.HexColor("#f3f6f8")),
                    ("LEFTPADDING", (0, 0), (-1, -1), 12),
                    ("TOPPADDING", (0, 0), (-1, -1), 4),
                    ("BOTTOMPADDING", (0, 0), (-1, -1), 4)]))
                story.extend([KeepTogether(panel), Spacer(1, 12)])
            elif line.startswith("|"):
                rows = [line]
                while i < len(lines) and lines[i].strip().startswith("|"):
                    rows.append(lines[i].strip())
                    i += 1
                cells = [[para(v.strip(), "table") for v in row.strip("|").split("|")]
                         for row in rows if not re.match(r"^\|[\s:|\-]+\|$", row)]
                add_table(cells)
            elif line.startswith("#"):
                count = len(line) - len(line.lstrip("#"))
                story.append(para(line[count:].strip(), {1: "title", 2: "h2"}.get(count, "h3")))
            else:
                block = [line.removeprefix("> ")]
                while i < len(lines) and lines[i].strip() and not lines[i].startswith(
                        ("#", "|", "```", "$$")):
                    block.append(lines[i].strip().removeprefix("> "))
                    i += 1
                story.append(para(" ".join(block)))

        def footer(canvas, doc):
            canvas.setFont("Latin", 8)
            canvas.setFillColor(colors.HexColor("#758491"))
            canvas.drawString(50, 28, "XMVB-CPP | String-to-structure contraction | Review draft")
            canvas.drawRightString(A4[0]-50, 28, str(doc.page))

        doc = SimpleDocTemplate(str(OUTPUT), pagesize=A4, rightMargin=50,
                                leftMargin=50, topMargin=48, bottomMargin=48,
                                title="Unique-spin-string to VB structure contraction",
                                author="xmvb-cpp project")
        doc.build(story, onFirstPage=footer, onLaterPages=footer)
    print(OUTPUT)


if __name__ == "__main__":
    main()
