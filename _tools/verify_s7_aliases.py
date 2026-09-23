"""用设计矩阵精确复算补充材料里「丢失/截断的运行顺序」所导致的别名结构。

背景：S7.3 说「若只实施前 8 次运行，得到的其实是下面这些量的估计」，
但原书/OCR 的表达式有残缺。本脚本按
    c_{P,Q} = (1/8) * (X_I^T X_I)_{QP}
精确算出每个被估量，从而判断原书印的是什么、正确值是什么。

用法：
    python _tools/verify_s7_aliases.py
"""

from __future__ import annotations

import itertools
import re
from pathlib import Path

# 表 2 / 表 3 的「标准顺序」（1-based），按运行顺序 1..16 列出
TABLE2_STD = [2, 12, 10, 15, 14, 4, 7, 3, 5, 8, 11, 16, 1, 9, 6, 13]
TABLE3_STD = [10, 15, 3, 6, 12, 8, 13, 1, 11, 2, 7, 14, 16, 5, 9, 4]

LETTERS = ["A", "B", "C", "D"]


def design_row(std: int) -> dict[str, int]:
    """标准顺序 -> 各因子水平（±1）。A=bit0, B=bit1, C=bit2, D=bit3 of (std-1)。"""
    n = std - 1
    return {letter: (1 if n >> i & 1 else -1) for i, letter in enumerate(LETTERS)}


def all_effects() -> list[str]:
    names: list[str] = []
    for r in range(1, 5):
        for combo in itertools.combinations(LETTERS, r):
            names.append("".join(combo))
    return names


def effect_sign(row: dict[str, int], name: str) -> int:
    value = 1
    for letter in name:
        value *= row[letter]
    return value


def alias_rows(std_order: list[int], first_n: int) -> dict[str, dict[str, float]]:
    """返回 {被估量: {真实项: 系数}}，只保留非零系数。"""
    subset = [design_row(s) for s in std_order[:first_n]]
    effects = all_effects()
    terms = ["I"] + effects
    result: dict[str, dict[str, float]] = {}
    for q in terms:
        coeffs: dict[str, float] = {}
        for p in terms:
            total = sum(
                (1 if q == "I" else effect_sign(row, q))
                * (1 if p == "I" else effect_sign(row, p))
                for row in subset
            )
            value = total / first_n
            if abs(value) > 1e-9:
                coeffs[p] = value
        result[q] = coeffs
    return result


def show(title: str, std_order: list[int], first_n: int) -> dict[str, dict[str, float]]:
    rows = alias_rows(std_order, first_n)
    print(f"===== {title}（前 {first_n} 次运行）=====")
    for q in ["I"] + all_effects():
        items = " ".join(
            f"{'+' if v > 0 else '−'}{'' if abs(v) == 1 else f'{abs(v):g}*'}{p}"
            for p, v in rows[q].items()
        )
        mnemonic = "截距" if q == "I" else q
        print(f"  [{mnemonic:>9}] = {items}")
    clean = [q for q in all_effects() if list(rows[q]) == [q] and rows[q][q] == 1]
    print(f"  -> 未被混杂（可干净估计）的效应：{clean or '无'}")
    print()
    return rows


def main() -> int:
    show("表 2：完全随机化的 2^4", TABLE2_STD, 8)
    show("表 3：四个区组、只做第 1 区组（运行 1-4）+ 第 2 区组（运行 5-8）",
         TABLE3_STD, 8)
    show("表 3：只做前三个区组（运行 1-12）", TABLE3_STD, 12)

    # 顺带核对表 3 的区组划分假设
    print("===== 表 3 各「4 次运行」分组的合法性 =====")
    for start in range(0, 16, 4):
        chunk = [design_row(s) for s in TABLE3_STD[start:start + 4]]
        signs = {name: sum(effect_sign(r, name) for r in chunk) for name in all_effects()}
        zero = [n for n, v in signs.items() if v == 0]
        four = [n for n, v in signs.items() if abs(v) == 4]
        print(f"  运行 {start + 1}-{start + 4}: 乘积为 0 的效应 {zero}；全同号效应 {four}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
