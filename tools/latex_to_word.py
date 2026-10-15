#!/usr/bin/env python3
"""Export the LaTeX master to an editable Word review snapshot.

Run with: uv run --with pypandoc_binary --with python-docx tools/latex_to_word.py
Choose --output for each review round; existing snapshots are never overwritten.
"""

import argparse
import re
from pathlib import Path
from zipfile import ZipFile
from xml.etree import ElementTree as ET

import pypandoc
from docx import Document
from docx.shared import Pt, RGBColor

ROOT = Path(__file__).resolve().parents[1]


def number_equations(source):
    """Preserve numbered equations and references unsupported by Pandoc's writer."""
    labels = {}
    count = 0

    def equation(match):
        nonlocal count
        count += 1
        body = match.group(1)
        for label in re.findall(r"\\label\{([^}]+)\}", body):
            labels[label] = count
        body = re.sub(r"\\label\{[^}]+\}", "", body).strip()
        return ("\\begin{equation}\n" + body +
                "\\qquad\\text{(" + str(count) + ")}\n\\end{equation}")

    source = re.sub(r"\\begin\{equation\}(.*?)\\end\{equation\}",
                    equation, source, flags=re.DOTALL)

    def reference(match):
        label = match.group(1)
        if label not in labels:
            raise ValueError(f"Unresolved equation reference: {label}")
        return f"({labels[label]})"

    return re.sub(r"\\eqref\{([^}]+)\}", reference, source)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        default=ROOT / "docs/jctc_second_order_from_latex.docx")
    args = parser.parse_args()
    if args.output.exists():
        parser.error("Choose a new output path to preserve reviewed snapshots.")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    source = (ROOT / "docs/jctc_second_order.tex").read_text()
    pypandoc.convert_text(
        number_equations(source), "docx", format="latex",
        outputfile=str(args.output),
        extra_args=["--standalone", "--number-sections",
                    "--reference-doc=" + str(ROOT / "docs/jctc_second_order_outline.docx")])
    doc = Document(args.output)
    # The reference supplies formatting, not review comments from an older draft.
    for rel_id, rel in list(doc.part.rels.items()):
        if rel.reltype.endswith("/comments"):
            del doc.part.rels[rel_id]
    # Pandoc uses these standard styles rather than the ACS custom body styles.
    for name in ("Normal", "Body Text", "First Paragraph"):
        if name not in doc.styles:
            continue
        style = doc.styles[name]
        style.font.name = "Times New Roman"
        style.font.size = Pt(12)
        style.font.color.rgb = RGBColor(0, 0, 0)
        style.paragraph_format.line_spacing = 1.3
    doc.core_properties.subject = "LaTeX-to-Word internal review snapshot"
    doc.save(args.output)
    ns = {"m": "http://schemas.openxmlformats.org/officeDocument/2006/math"}
    with ZipFile(args.output) as package:
        assert package.testzip() is None
        xml = ET.fromstring(package.read("word/document.xml"))
        math_count = len(xml.findall(".//m:oMath", ns))
        assert math_count >= 12, "Expected all twelve display equations."
    print(f"{args.output}\nNative Word math objects: {math_count}")


if __name__ == "__main__":
    main()
