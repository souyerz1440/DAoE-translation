"""
清理 MinerU 对第 15 章的 OCR 排版问题（格式向第 1–14 章对齐）。

处理内容：
  1. 还原被空格拆开的图/表/例标题：
        "## E X A M P L E 15 . 1"      -> "## EXAMPLE 15.1"
        "## ◾ T A B L E 15 . 1"        -> "## TABLE 15.1"
        "◾ T A B L E 15 . 5 (Continued)" -> "TABLE 15.5 (Continued)"
     （第 1–14 章的源文件都不带 ◾ 符号，且 ## / 无 ## 两种风格都有，故保持原样只做还原。）
  2. 合并被 OCR 拆成两行的例标题：
        "## EXAMPLE 15.4" + "## The Worsted Yarn Experiment" -> 一行
  3. 补回丢失的节号：
        "## Unbalanced Data in a Factorial Design" -> "## 15.2 Unbalanced Data in a Factorial Design"
  4. 把挤成一团的 CHAPTER OUTLINE HTML 表格还原成嵌套列表（按原书排版）。

用法：
    python _tools/fix_ch15_format.py FILE            # 只预览
    python _tools/fix_ch15_format.py FILE --apply    # 写入
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

# 行首：可选标题井号 + 可选 ◾ + 被空格拆开的全大写单词
LINE_HEAD = re.compile(
    r"^(?P<hashes>#{1,6}[ \t]+)?"
    r"(?P<bullet>\u25fe[ \t]*)?"
    r"(?P<word>(?:[A-Z][ \t]+){2,}[A-Z])"
)
# 单词之后紧跟的编号，只允许数字 / 空格 / 制表符 / 点，绝不跨行
LINE_NUMBER = re.compile(r"^(?P<num>[0-9][0-9 \t.]*)")
WORD_RE = re.compile(r"[ \t]+")

ALLOWED_WORDS = {"TABLE", "EXAMPLE", "FIGURE"}

OUTLINE_ITEMS: list[tuple[int, str]] = [
    (0, "**15.1 NONNORMAL RESPONSES AND TRANSFORMATIONS**"),
    (1, "15.1.1 Selecting a Transformation: The Box–Cox Method"),
    (1, "15.1.2 The Generalized Linear Model"),
    (0, "**15.2 UNBALANCED DATA IN A FACTORIAL DESIGN**"),
    (1, "15.2.1 Proportional Data: An Easy Case"),
    (1, "15.2.2 Approximate Methods"),
    (1, "15.2.3 The Exact Method"),
    (0, "**15.3 THE ANALYSIS OF COVARIANCE**"),
    (1, "15.3.1 Description of the Procedure"),
    (1, "15.3.2 Computer Solution"),
    (1, "15.3.3 Development by the General Regression Significance Test"),
    (1, "15.3.4 Factorial Experiments with Covariates"),
    (0, "**15.4 REPEATED MEASURES**"),
    (0, "**SUPPLEMENTAL MATERIAL FOR CHAPTER 15**"),
    (1, "S15.1 The Form of a Transformation"),
    (1, r"S15.2 Selecting $\lambda$ in the Box–Cox Method"),
    (1, "S15.3 Generalized Linear Models"),
    (2, "S15.3.1 Models with a Binary Response Variable"),
    (2, "S15.3.2 Estimating the Parameters in a Logistic Regression Model"),
    (2, "S15.3.3 Interpreting the Parameters in a Logistic Regression Model"),
    (2, "S15.3.4 Hypothesis Tests on Model Parameters"),
    (2, "S15.3.5 Poisson Regression"),
    (2, "S15.3.6 The Generalized Linear Model"),
    (2, "S15.3.7 Link Functions and Linear Predictors"),
    (2, "S15.3.8 Parameter Estimation in the Generalized Linear Model"),
    (2, "S15.3.9 Prediction and Estimation with the Generalized Linear Model"),
    (2, "S15.3.10 Residual Analysis in the Generalized Linear Model"),
    (1, "S15.4 Unbalanced Data in a Factorial Design"),
    (2, "S15.4.1 The Regression Model Approach"),
    (2, "S15.4.2 The Type 3 Analysis"),
    (2, "S15.4.3 Type 1, Type 2, Type 3 and Type 4 Sums of Squares"),
    (2, "S15.4.4 Analysis of Unbalanced Data using the Means Model"),
]

ALLOWED_WORDS = {"TABLE", "EXAMPLE", "FIGURE"}


def render_outline(indent_unit: str = "  ") -> str:
    return "\n".join(f"{indent_unit * depth}- {text}" for depth, text in OUTLINE_ITEMS)


def normalize_line(line: str) -> tuple[str, str | None]:
    """还原单行里被空格拆开的图/表/例标题；返回 (新行, 说明)。"""
    stripped = line.rstrip()
    trailing = line[len(stripped):]

    head = LINE_HEAD.match(stripped)
    if not head:
        return line, None

    word = WORD_RE.sub("", head.group("word"))
    if word not in ALLOWED_WORDS:
        return line, None

    tail = stripped[head.end():].lstrip(" \t")
    number = LINE_NUMBER.match(tail)
    if not number:
        return line, None

    num = WORD_RE.sub("", number.group("num")).rstrip(".")
    if not num:
        return line, None

    rest = tail[number.end():].strip()
    hashes = head.group("hashes") or ""
    rebuilt = f"{hashes}{word} {num}" + (f" {rest}" if rest else "") + trailing

    if rebuilt == line:
        return line, None
    return rebuilt, f"  {stripped[:72]!r} -> {rebuilt.rstrip()[:72]!r}"


def main() -> int:
    if len(sys.argv) < 2:
        print(__doc__)
        return 2
    path = Path(sys.argv[1])
    apply = "--apply" in sys.argv

    raw = path.read_bytes()
    newline = "\r\n" if b"\r\n" in raw else "\n"
    newline_name = "CRLF" if newline == "\r\n" else "LF"
    text = raw.decode("utf-8")
    notes: list[str] = []
    changed_captions = 0

    # 1. 逐行还原被拆开的图/表/例标题
    out_lines: list[str] = []
    for line in text.split("\n"):
        new_line, note = normalize_line(line)
        if note:
            changed_captions += 1
            notes.append(note)
        out_lines.append(new_line)
    text = "\n".join(out_lines)

    # 2. 合并被拆开的例标题
    merged = re.sub(
        r"(?m)^## EXAMPLE (\d+\.\d+)\r?\n\r?\n## (?!EXAMPLE|TABLE|FIGURE)(.+)$",
        lambda m: f"## EXAMPLE {m.group(1)} {m.group(2).strip()}",
        text,
    )
    if merged != text:
        notes.append("  合并被拆成两行的例标题 -> ## EXAMPLE 15.4 The Worsted Yarn Experiment")
        text = merged

    # 3. 补回 15.2 的节号
    fixed_152 = text.replace(
        "## Unbalanced Data in a Factorial Design",
        "## 15.2 Unbalanced Data in a Factorial Design",
    )
    if fixed_152 != text:
        notes.append("  补回节号 -> ## 15.2 Unbalanced Data in a Factorial Design")
        text = fixed_152

    # 4. CHAPTER OUTLINE 表格 -> 嵌套列表
    outline_pattern = re.compile(
        r"^<table><tr><td>15\.1 NONNORMAL RESPONSES.*?</table>$", re.MULTILINE
    )
    outline_match = outline_pattern.search(text)
    if outline_match:
        text = text[: outline_match.start()] + render_outline() + text[outline_match.end():]
        notes.append("  CHAPTER OUTLINE：HTML 表格 -> 嵌套列表")
    else:
        notes.append("  [跳过] 未找到 CHAPTER OUTLINE 表格")

    print(f"文件：{path}")
    print(f"换行符：{newline_name}")
    print(f"改动：标题还原 {changed_captions} 行，其他 {len(notes) - changed_captions} 项")
    for note in notes:
        print(note)

    if apply:
        with open(path, "w", encoding="utf-8", newline="") as handle:
            handle.write(text.replace("\r\n", "\n").replace("\n", newline))
        print("\n已写入。")
    else:
        print("\n（预览模式，未写入；加 --apply 生效）")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
