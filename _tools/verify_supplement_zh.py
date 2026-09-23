"""校验补充材料的译文：结构完整性、公式/表格/图片是否保留、图片链接是否有效。

用法：
    python _tools/verify_supplement_zh.py            # 概要
    python _tools/verify_supplement_zh.py --detail   # 附逐文件明细
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SUPP = ROOT / "Supplemental Materials"

DISP = re.compile(r"\$\$")
INLINE = re.compile(r"(?<!\$)\$(?!\$)[^$\n]+\$(?!\$)")
MD_TABLE_ROW = re.compile(r"(?m)^\|")
HTML_TABLE = re.compile(r"<table\b", re.I)
MD_IMG = re.compile(r"!\[[^\]]*\]\(([^)]+)\)")
HTML_IMG = re.compile(r"(?i)<img\b[^>]*\bsrc\s*=\s*[\"']([^\"']+)[\"']")
CJK = re.compile(r"[\u4e00-\u9fff]")
TRUNC = re.compile(r"\[One or more long lines were truncated")


def stats(text: str) -> dict[str, int]:
    return {
        "cjk": len(CJK.findall(text)),
        "display_math": len(DISP.findall(text)) // 2,
        "inline_math": len(INLINE.findall(text)),
        "md_table_rows": len(MD_TABLE_ROW.findall(text)),
        "html_tables": len(HTML_TABLE.findall(text)),
        "images": len(MD_IMG.findall(text)) + len(HTML_IMG.findall(text)),
        "headings": len(re.findall(r"(?m)^#{1,6}\s", text)),
    }


def image_paths(text: str) -> list[str]:
    return MD_IMG.findall(text) + HTML_IMG.findall(text)


def main() -> int:
    detail = "--detail" in sys.argv
    expected = list(range(1, 9))
    problems: list[str] = []

    print(f"{'章':>4} {'源B':>7} {'译B':>7} {'汉字':>6} {'源表':>5} {'译表行':>6} "
          f"{'源公式':>6} {'译公式':>6} {'源图':>5} {'译图':>5} {'H':>3}")
    for n in expected:
        src = SUPP / f"ch{n:02d}.md"
        out = ROOT / f"Chapter {n:02d}" / "translations_zh" / f"S{n}_补充材料_zh.md"
        if not src.is_file():
            problems.append(f"ch{n:02d}: 源文件缺失 {src}")
            continue
        if not out.is_file():
            print(f"{n:>4} {src.stat().st_size:>7} {'—':>7}   ← 译文尚未生成")
            problems.append(f"ch{n:02d}: 译文缺失 {out}")
            continue

        s_text, o_text = src.read_text(encoding="utf-8"), out.read_text(encoding="utf-8")
        s, o = stats(s_text), stats(o_text)

        print(f"{n:>4} {src.stat().st_size:>7} {out.stat().st_size:>7} {o['cjk']:>6} "
              f"{s['html_tables']:>5} {o['md_table_rows']:>6} "
              f"{s['display_math']:>6} {o['display_math']:>6} "
              f"{s['images']:>5} {o['images']:>5} {o['headings']:>3}")

        if o["html_tables"]:
            problems.append(f"ch{n:02d}: 译文仍残留 {o['html_tables']} 个 <table>")
        n_trunc = len(TRUNC.findall(s_text))
        if n_trunc:
            problems.append(
                f"ch{n:02d}: 源文件有 {n_trunc} 处长行截断标记，译文需对照 PDF 核对"
            )
        if o["cjk"] < 500:
            problems.append(f"ch{n:02d}: 译文汉字仅 {o['cjk']} 个，疑似未完整翻译")
        if not o_text.lstrip().startswith("# "):
            problems.append(f"ch{n:02d}: 译文首行不是 H1")
        if s["images"] != o["images"]:
            problems.append(f"ch{n:02d}: 图片引用数不符（源 {s['images']} / 译 {o['images']}）")

        for ref in image_paths(o_text):
            if re.match(r"^(https?:|data:)", ref):
                continue
            resolved = (out.parent / ref).resolve()
            if not resolved.is_file():
                problems.append(f"ch{n:02d}: 图片链接失效 {ref}")

        if detail:
            print(f"      标题：")
            for line in o_text.splitlines():
                if line.startswith("#"):
                    print(f"        {line[:88]}")

    print()
    if problems:
        print(f"发现 {len(problems)} 个问题：")
        for p in problems:
            print(f"  - {p}")
    else:
        print("全部 8 章校验通过。")
    return 1 if problems else 0


if __name__ == "__main__":
    raise SystemExit(main())
