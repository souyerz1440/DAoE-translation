"""检查补充材料译文里可能「漏译」的地方：未译标题、连续多行纯英文段落。"""

from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
CJK = re.compile(r"[\u4e00-\u9fff]")

# 允许保持英文的行首：公式/表格/图片/标题/代码块/引用/列表符号
SKIP_PREFIX = ("$", "|", "!", "#", ">", "`", "-", "*", "  ")


def main() -> int:
    for n in range(1, 9):
        path = ROOT / f"Chapter {n:02d}" / "translations_zh" / f"S{n}_补充材料_zh.md"
        if not path.is_file():
            print(f"ch{n:02d}: 缺文件 {path.name}")
            continue
        text = path.read_text(encoding="utf-8")
        lines = text.split("\n")

        bad_headings = [
            line for line in lines
            if line.startswith("#") and not CJK.search(line)
        ]

        runs: list[tuple[int, int]] = []
        run_start, run_len = 0, 0
        in_code = False
        for i, line in enumerate(lines, 1):
            stripped = line.strip()
            if stripped.startswith("```"):
                in_code = not in_code
                continue
            english = (
                len(stripped) > 40
                and not CJK.search(stripped)
                and not stripped.startswith(SKIP_PREFIX)
            )
            if english and not in_code:
                if run_len == 0:
                    run_start = i
                run_len += 1
            else:
                if run_len >= 3:
                    runs.append((run_start, run_len))
                run_len = 0
        if run_len >= 3:
            runs.append((run_start, run_len))

        print(f"ch{n:02d}: 未译标题 {len(bad_headings)} 个；连续纯英文段 {len(runs)} 处")
        for h in bad_headings:
            print(f"    标题: {h[:88]}")
        for start, length in runs:
            print(f"    L{start} 起 {length} 行")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
