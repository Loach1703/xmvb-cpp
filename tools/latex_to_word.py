#!/usr/bin/env python3
"""Export the LaTeX master to an editable Word review snapshot.

Run with: uv run --with pypandoc_binary --with python-docx tools/latex_to_word.py
Choose --output for each review round; existing snapshots are never overwritten.
"""

import argparse
from pathlib import Path
from zipfile import ZipFile
from xml.etree import ElementTree as ET

import pypandoc
from docx import Document
from docx.shared import Pt, RGBColor

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        default=ROOT / "docs/jctc_second_order_from_latex.docx")
    args = parser.parse_args()
    if args.output.exists():
        parser.error("Choose a new output path to preserve reviewed snapshots.")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    pypandoc.convert_file(
        str(ROOT / "docs/jctc_second_order.tex"), "docx", format="latex",
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
