#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Render the MinerU HTML tables inside a Markdown file as readable pipe rows.

MinerU emits tables as one long ``<table>...</table>`` blob, often on a single
very long line, which is awkward to read (and gets truncated by the editor
tool).  This helper prints every table found in the file as::

    == TABLE (line 41) ==
    | (1) | - | - | + | ... |
    | a   | + | - | - | ... |

so the content can be transcribed into a Markdown table without guesswork.

Usage::

    python html_table_txt.py FILE [--line N]
"""

import argparse
import html
import re
import sys
from pathlib import Path

TABLE_RE = re.compile(r"<table\b.*?</table>", re.S | re.I)
ROW_RE = re.compile(r"<tr\b.*?</tr>", re.S | re.I)
CELL_RE = re.compile(r"<t[dh]\b([^>]*)>(.*?)</t[dh]>", re.S | re.I)
SPAN_RE = re.compile(r"colspan\s*=\s*\"?(\d+)", re.I)
ROWSPAN_RE = re.compile(r"rowspan\s*=\s*\"?(\d+)", re.I)


def cell_text(raw: str) -> str:
    raw = re.sub(r"<img\b[^>]*>", "[image]", raw, flags=re.I)
    raw = re.sub(r"<br\s*/?>", " / ", raw, flags=re.I)
    raw = re.sub(r"<[^>]+>", "", raw)
    raw = html.unescape(raw)
    return re.sub(r"\s+", " ", raw).strip()


def render(table_html: str) -> list:
    rows = []
    for tr in ROW_RE.findall(table_html):
        cells = []
        for attrs, body in CELL_RE.findall(tr):
            text = cell_text(body)
            colspan = SPAN_RE.search(attrs)
            rowspan = ROWSPAN_RE.search(attrs)
            marks = []
            if colspan:
                marks.append(f"colspan={colspan.group(1)}")
            if rowspan:
                marks.append(f"rowspan={rowspan.group(1)}")
            if marks:
                text = f"{text} [{' '.join(marks)}]"
            cells.append(text)
        rows.append(cells)
    return rows


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("file")
    parser.add_argument("--line", type=int)
    args = parser.parse_args()

    path = Path(args.file)
    text = path.read_text(encoding="utf-8")

    for match in TABLE_RE.finditer(text):
        line = text.count("\n", 0, match.start()) + 1
        if args.line and line != args.line:
            continue
        print(f"== TABLE (line {line}) ==")
        for cells in render(match.group(0)):
            print("| " + " | ".join(cells) + " |")
        print()
    return 0


if __name__ == "__main__":
    sys.exit(main())
