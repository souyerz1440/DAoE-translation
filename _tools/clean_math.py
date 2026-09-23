#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Clean up the systematic spacing artefacts MinerU leaves inside math mode.

Rules (all applied only when they are unambiguous):

* ``9 0 + 1 0 0``  ->  ``90 + 100``      (digit, space, digit -> digits)
* ``(- 3 0)^{2}``  ->  ``(-30)^{2}``     (same rule, repeated until stable)
* ``9398. 0 0``    ->  ``9398.00``       (digit, dot, space, digit)
* ``\\frac {1}``   ->  ``\\frac{1}``     (macro name followed by space)
* ``S S _ {A}``    ->  ``SS_{A}``
* ``A B``          ->  ``AB``            (case sensitive! ``a b`` is a
                                          treatment combination)
* ``\\overline{{y}} _ {A}`` -> ``\\overline{{y}}_{A}``  (space hugging _ or ^)

Chinese prose never puts a space between two digits, and never writes
``\\frac {``, so the replacements are safe outside math as well.

Usage::

    python clean_math.py FILE [FILE ...] [--apply]
"""

import argparse
import re
import sys
from pathlib import Path

RULES = [
    (re.compile(r"(?<=\d) (?=\d)"), ""),
    (re.compile(r"(?<=\d)\. (?=\d)"), "."),
    (re.compile(r"(\\[a-zA-Z]+) \{"), r"\1{"),
    (re.compile(r"(\\[a-zA-Z]+)\["), r"\1["),
    (re.compile(r"S S _"), "SS_"),
    (re.compile(r"S S_"), "SS_"),
    (re.compile(r"\bS S\b"), "SS"),
    (re.compile(r"\bM S _"), "MS_"),
    (re.compile(r"\bM S\b"), "MS"),
    (re.compile(r"\bP S E\b"), "PSE"),
    (re.compile(r"\bA B\b"), "AB"),
    (re.compile(r"(?<=[A-Za-z0-9)\]}])\s+_\s*\{"), "_{"),
    (re.compile(r"(?<=[A-Za-z0-9)\]}])\s+\^\s*\{"), "^{"),
]


def clean(text: str) -> str:
    previous = None
    while previous != text:
        previous = text
        for pattern, replacement in RULES:
            text = pattern.sub(replacement, text)
    return text


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("files", nargs="+")
    parser.add_argument("--apply", action="store_true")
    args = parser.parse_args()

    for name in args.files:
        path = Path(name)
        before = path.read_text(encoding="utf-8")
        after = clean(before)
        if before == after:
            print(f"clean    {path}")
            continue
        if args.apply:
            path.write_text(after, encoding="utf-8")
            print(f"UPDATED  {path}")
        else:
            print(f"NEEDS CLEANUP  {path}")
            bl = before.split("\n")
            al = after.split("\n")
            for i in range(min(len(bl), len(al))):
                if bl[i] != al[i]:
                    print(f"   line {i + 1}")
                    print(f"     - {bl[i]}")
                    print(f"     + {al[i]}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
