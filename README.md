# 试验设计与分析翻译计划

《Design and Analysis of Experiments》（Douglas C. Montgomery）中文翻译。在线版：<https://souyerz1440.github.io/DAoE-translation/>

阶段一（拆分与格式转换）已完成：按书签把原书拆成 15 个章 PDF，用 MinerU 转成 Markdown，再按二级书签（Section）切分为逐节文件。
阶段二（逐节翻译）进行中：**第 1～15 章已译完并收录**（入口见 `book/index.md`），译文放在各章 `translations_zh/`。
暂不包含附录与 Index。

## 内容范围

新增“习题与答案”栏目：入口为 `Exercises/index.md`，顶部菜单和前 3 章章首页均可进入。按第 9 版编号完整收录前 3 章 143 题，题目、数据与答案解析均为中文；原答案题号错位及必要勘误在对应题旁注明。

习题按“章节 → 题号 → 独立页面”组织：`Exercises/translations_zh/chapter02.md` 是章目录，`Exercises/translations_zh/chapter02/exercise-2-31.md` 是单题页。每页包含完整中文背景、各小问、数据、直接展示的答案，以及上一题、下一题和关联题链接。

数据为 `Exercises/data/chapter01-03.json`，单题目录清单为 `Exercises/data/pages.json`，图表位于 `Exercises/images/`。维护时编辑 `_tools/exercise_templates/` 中的中文稿，再运行 `python _tools/build_exercises.py` 复算数据、生成单题页与章节索引，并同步 `toc.yml`（需要 numpy、scipy、matplotlib）。仅调整文字或页面组织时，可用 `python _tools/build_exercises.py --reuse-figures` 复用已有图表。常规电子书构建直接使用已生成文件，无需这些计算依赖。可运行 `python _tools/verify_exercises.py` 检查题号、页面导航、数据、引用和关键复算结果。

公式使用主题对应的 KaTeX 0.15.2 样式与字体，本地副本位于 `book/static/katex/`，保留原 MIT 许可证，避免 CDN 无法访问时公式重影或排版错乱。

| 项目 | 值 |
| --- | --- |
| 原书 PDF 总页数 | 682 页（扫描件，无文本层，163 ppi） |
| 对照用 PDF | `DesignandAnalysisofExperiments9thEdition.pdf`（第 9 版，**有文本层**，749 页） |
| 已转换 | 第 1～15 章 |
| 章数 | 15 |
| 二级 Section 数 | 100 |
| Markdown 文件数 | 115（100 个 Section + 15 个章首 intro） |
| 章节内图片数 | 517 |

根目录两本原书 PDF（第 10 版扫描件、第 9 版对照本）只在本地，不入库。

Section 文件名格式为 `{序号}_{编号}_{标题}.md`（与 `slp3翻译计划` 一致），`sections/` 内的图片引用自动改写为 `../images/…`。
另有一套一次性校验脚本在工作区根目录的 `tools/`（`d:\Code for VS\tools\`，不在本项目内，带 `_` 前缀）：`_contentcheck.py` 核对图片是否与高清插图逐字节一致，`_tblcheck*.py` 比对表格数值，`_tbl826_zh.py` / `_tbl823_check.py` 为第 8 章专项。

进度：

| 章节 | 英文文件数 | 已译 |
| --- | --- | --- |
| 第 1 章 | 7 | 已译完（已上线） |
| 第 2 章 | 7 | 已译完（已上线） |
| 第 3 章 | 12 | 已译完（已上线） |
| 第 4 章 | 5 | 已译完（已上线） |
| 第 5 章 | 7 | 已译完（已上线） |
| 第 6 章 | 10 | 已译完（已上线） |
| 第 7 章 | 9 | 已译完（已上线） |
| 第 8 章 | 10 | 已译完（已上线） |
| 第 9 章 | 7 | 已译完（已上线） |
| 第 10 章 | 9 | 已译完（已上线） |
| 第 11 章 | 8 | 已译完（已收录，含补充材料） |
| 第 12 章 | 6 | 已译完（已收录，含补充材料） |
| 第 13 章 | 7 | 已译完（已收录，含补充材料） |
| 第 14 章 | 6 | 已译完（已收录，含补充材料） |
| 第 15 章 | 5 | 已译完（已收录，含补充材料） |

“英文文件数”含章首 intro 文件。

新增整章译文上线时，需要同步更新**四处**：`toc.yml`（加章条目与 children）、
`myst.yml` 中 `project.exclude` 的对应排除项、`book/progress.md` 的状态，
以及 `book/index.md` 的“已上线章节”列表与导语。

## 已知的原文（OCR）问题

译文一般按原文译出；会影响理解的地方加“译者注”说明，扫描件中丢失的脚注按原书补出。
