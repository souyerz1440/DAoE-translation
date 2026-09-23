"""修复被 subagent 误编辑搞坏的仓库记忆文件。

损坏情况：正文被**整段重复**一次（「## 命名与结构」出现两次），外加一处换行丢失。
本脚本保留第一份副本（1..352 行，含第 15 章 markdown 记录），
丢弃重复副本与伪造的告警头，补回被截断的「补充材料翻译」小节，
再接上独有的尾部两节（第 15 章排版清理 / PowerShell 踩坑、第 6 章补充材料 S6）。
"""

from __future__ import annotations

import os
import shutil
from pathlib import Path

MEM = Path(
    r"C:\Users\zsy14\AppData\Roaming\Code\User\workspaceStorage"
    r"\e1e04329206b082515e98a497271bfa5\GitHub.copilot-chat\memory-tool"
    r"\memories\repo\design-of-experiments-translation.md"
)

SUPPLEMENT_SECTION = """## 补充材料翻译（S1–S8 已完成，2026-09-22）
- 输出：`Chapter 0N/translations_zh/SN_补充材料_zh.md`；首行 `# 第 N 章补充材料`（**不加** PDF 链接，与章首 `00_*` 文件不同）。
- 节标题 `## SN.x. Title` → `## SN.x 中文标题`（去掉编号后的句点）。
- 工具（base 环境 python）：
  - `_tools/html_table_txt.py "Supplemental Materials/chNN.md"` 读 HTML 表（每篇 2–5 张表）。**长行会被 read_file 截断，必须用它**。
  - `_tools/clean_math.py <译文> [--apply]`：译完跑一次；S1–S8 均报 `clean`。
  - `_tools/verify_supplement_zh.py [--detail]`：**校验 8 篇译文**（汉字数、表转换、`$$` 公式数、图片数、标题数、图片链接有效性）。
  - `_tools/place_supplement_images.py [--apply]`：把补充材料里被引用的图片拷进 `Chapter NN/images/`。
- **补充材料源 PDF 有文本层**：`pdftotext -layout -enc UTF-8` 或 `-raw` / `-bbox` 可逐字校验，比 MinerU 可靠得多。
  `-bbox` 还能看出**两行表头**（如 ch08 的 `Factor` / `A:A`）——`-layout` 会把它压成 `Factor A:A`。
  **但注意**：文本层仍是「扫描件 + OCR 层」而非原生文本（大括号常识成 `FG`/`IJ`，`β`/`⋯` 常丢），
  所以「pdftotext 与 MinerU 一致」只说明扫描件如此，**不等于原书如此**，仅可当旁证。
- **别名结构列表的渲染惯例**：源里是 `<div class="mineru-algorithm">` 的 `$[A] = …$` 串或 ```ini 代码块 →
  均改写成单个 `$$\\begin{array}{r l} & [\\mathrm{A}] = … \\\\ … \\end{array}$$`（与 ch07 译文一致，`&` 对齐）。
- **设计软件输出的表头**：ch07 源是 `Std Order | Run Order | Block | Factor A | …`，ch08 源是两行 `Factor`/`A:A`；
  译文统一压平为 `标准顺序 | 运行顺序 | 区组 | 因子 A | …`（`-1` 用 U+2212 `−`）。
- **表格编号**：补充材料用具内局部编号（表 1、表 2…），译成 `**表 1 …**` 加粗表题放在表格**上方**；
  ch02 源文里没有独立 `Table N.` 表题行，故不杜撰表题。
- **图片命名**（由 `place_supplement_images.py` 处理）：有编号图 → `figureS<章>.<图号>.jpg`；
  无编号图 → `figureS<章>.uf<序号>.jpg`（沿用图库 `cNNufNNM` 的 uf 词表）。
  ch01 三张（`figureS1.1/1.2/1.3.jpg`）、ch02 两张（`figureS2.uf001/002.jpg`）、
  ch06 三张（`figureS6.uf001/002/003.jpg`），ch03/04/05/07/08 无图。
  译文在 `translations_zh/` 下，引用一律写 `../images/…`。
- **ch01 原文事实**：第三张图的图题被 OCR 错印成 `Figure 2`，按正文
  “**Figure 3** [from Barton (1999)] shows a cause-and-effect diagram” 更正为**图 3**。
  ch01 的 `## Blank Guide Sheets…` 与 `## Interactions` 是 S1.1 的子部分，译文降级为 `###`。
- **ch08 原文事实（勿当 OCR 错误改）**：① 主效应别名链只列 3 项（省略 4/5/6 因子交互，PDF 亦然）；
  ② “a partial fold with only four runs over can be **constricted**” 是原文用词（应为 constructed），译文按词义译出并在答复里说明。
- **ch06 有几处原文印误，译文按原文保留并在答复中列出**（未加译者注）：S6.1 第 2 个正规方程首项印 `\\sum x_{i}`（应为 `x_{i1}`）；
  表 1「平方和」表头印 `(3)^{2}\\div n2^{k-1}` 而正文/表值为 `n2^{k}=16`；S6.4 加第 4 次重复处印 `\\frac{1}{12}` 而结果写 `\\sqrt{3/16}`；
  S6.7 的 `E(\\overline{y}_{F})` 第二行漏 `\\beta_{0}`、假设式印 `+\\dots\\beta_{kk}`（少一个 `+`）。
- **尚未翻译**：第 9–14 章补充材料、以及**第 15 章补充材料**（`ch15.md`，95 KB，含 S15.1–S15.4）。

## 并行翻译的教训（2026-09-22 实测）
- `runSubagent` 可以并行译各章（一次 4 个），**每个 prompt 里必须带**：源路径、输出路径、完整风格指南、
  以及「已完成的参考样例」路径（`Chapter 07/translations_zh/S7_补充材料_zh.md`），这样 8 篇风格才一致。
- **事后必须逐文件校验**，不能只看 subagent 的汇报：用 `verify_supplement_zh.py` 对比
  「源表数 → 译表行数」「源 `$$` 数 → 译 `$$` 数」「源图数 → 译图数」，本次 S3/S4/S5/S6 公式数**完全相等**才算过。
- **⚠️ subagent 会写仓库记忆文件，可能把它搞坏**：本次 ch06 的 subagent 把 `design-of-experiments-translation.md`
  **整段重复了一遍**（「## 命名与结构」出现两次）并丢了一处换行，还自己加了一条"本文件历史上被误编辑"的假告警。
  → **以后给 subagent 的 prompt 里要明确写「禁止修改 /memories/ 下的任何文件」**，并定期 `memory view` 抽查。
"""


def main() -> int:
    lines = MEM.read_text(encoding="utf-8").split("\n")
    print(f"原件：{len(lines)} 行")

    # 1..2 = 标题 + 空行；跳过 3..5（subagent 伪造的告警）；6..352 = 正文（到第 15 章 markdown 记录末尾）
    head = lines[0:2]
    body = lines[5:352]
    assert head[0].startswith("# 试验设计与分析翻译计划"), head[0]
    assert body[0].strip() == "", repr(body[0])
    assert body[1].startswith("## 命名与结构"), repr(body[1])
    assert "未处理（留给翻译阶段" in "\n".join(body[-6:]), body[-6:]

    # 716..end = 独有的尾部两节（第 15 章排版清理 / PowerShell 踩坑、第 6 章补充材料 S6）
    tail = lines[715:]
    assert tail[1].startswith("## 第 15 章排版清理"), repr(tail[1])
    assert any("第 6 章补充材料 S6" in x for x in tail), "尾部缺少 S6 小节"

    new_text = "\n".join(head + body + [SUPPLEMENT_SECTION] + tail)

    backup = MEM.with_suffix(".md.broken")
    shutil.copy2(MEM, backup)
    print(f"损坏件备份：{backup}")

    MEM.write_text(new_text, encoding="utf-8", newline="\n")
    out = MEM.read_text(encoding="utf-8")
    print(f"修复后：{len(out.splitlines())} 行，{len(out)} 字节")
    print(f"「## 命名与结构」出现次数：{out.count('## 命名与结构')}（应为 1）")
    print(f"伪告警残留：{'是' if '历史上被一次误编辑' in out else '否'}（应为「否」）")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
