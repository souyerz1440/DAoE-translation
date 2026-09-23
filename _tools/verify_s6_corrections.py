"""核对 S6.1–S6.7 的 5 处数学校正是否已落实，并复核译文整体状态。"""

from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
CH06 = ROOT / "Chapter 06" / "translations_zh" / "S6_补充材料_zh.md"

# (说明, 应存在, 应不存在)
CHECKS: list[tuple[str, str | None, str | None]] = [
    (
        "S6.1 第 2 个正规方程首项 x_{i1}",
        r"\\hat\{\\beta\}_\{0\}\\sum_\{i=1\}\^\{4\}x_\{i1\} \+ \\hat\{\\beta\}_\{1\}",
        r"\\hat\{\\beta\}_\{0\}\\sum_\{i=1\}\^\{4\}x_\{i\} \+ \\hat\{\\beta\}_\{1\}",
    ),
    (
        "S6.2 表头 平方和 = n2^k",
        r"平方和 \$\(3\)\^\{2\}\\div n2\^\{k\}\$",
        None,
    ),
    (
        "S6.2 表头仍保留 效应估计值 = n2^{k-1}",
        r"效应估计值 \$\(3\)\\div n2\^\{k-1\}\$",
        None,
    ),
    (
        "S6.4 第 4 次重复处 1/16",
        r"\\frac\{1\}\{16\}\(1 \+ \(1\)\^\{2\} \+ \(1\)\^\{2\}\)",
        None,
    ),
    (
        "S6.4 n=3 处仍为 1/12",
        r"\\frac\{1\}\{12\}\(1 \+ \(1\)\^\{2\} \+ \(1\)\^\{2\}\)",
        None,
    ),
    (
        "S6.7 补出 beta_0",
        r"& = \\beta_\{0\} \+ \\beta_\{11\} \+ \\beta_\{22\} \+ \\dots \+ \\beta_\{kk\}",
        None,
    ),
    (
        "S6.7 假设式补出 +",
        r"H_\{0\}: \\beta_\{11\} \+ \\beta_\{22\} \+ \\dots \+ \\beta_\{kk\} = 0",
        None,
    ),
    (
        "S6.7 差值式不应含 beta_0（β0 相消）",
        r"E\(\\overline\{y\}_\{F\} - \\overline\{y\}_\{C\}\) = \\beta_\{11\}",
        None,
    ),
]


def main() -> int:
    text = CH06.read_text(encoding="utf-8")
    failures: list[str] = []

    print(f"文件：{CH06.relative_to(ROOT)}")
    for label, must_have, must_not in CHECKS:
        have = bool(re.search(must_have, text)) if must_have else True
        gone = not re.search(must_not, text) if must_not else True
        status = "OK  " if have and gone else "FAIL"
        if not (have and gone):
            failures.append(label)
        extra = ""
        if not have:
            extra = "  (未找到期望内容)"
        if not gone:
            extra = "  (仍残留旧内容)"
        print(f"  [{status}] {label}{extra}")

    marks = sorted(set(re.findall(r"\[\^(\d+)\](?!:)", text)))
    defs = sorted(set(re.findall(r"(?m)^\[\^(\d+)\]:", text)))
    print(f"  脚注标记 {marks} / 定义 {defs} -> {'OK  ' if marks == defs else 'FAIL'}")
    if marks != defs:
        failures.append("脚注标记与定义不匹配")

    print()
    if failures:
        print(f"未通过 {len(failures)} 项：")
        for f in failures:
            print(f"  - {f}")
        return 1
    print("S6 的 5 处数学校正全部就位。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
