# 试验设计与分析翻译计划

《Design and Analysis of Experiments》（Douglas C. Montgomery）中文翻译计划。

阶段一（结构拆分与格式转换）已完成：按书签把原书拆成 14 个章 PDF，用 MinerU 转成 Markdown，
再按二级书签（Section）切分为逐节文件。阶段二（逐节翻译）进行中：**第 1～8 章译文已完成并上线**
（入口见 `book/index.md`），译文放在各章 `translations_zh/`，命名与格式对照 `slp3翻译计划`，
详见下文“翻译”一节。暂不包含附录与 Index，不包含第 15 章及之后内容。

## 内容范围

| 项目 | 值 |
| --- | --- |
| 原书 PDF 总页数 | 682 页（扫描件，无文本层，163 ppi） |
| 对照用 PDF | `DesignandAnalysisofExperiments9thEdition.pdf`（第 9 版，**有文本层**，749 页） |
| 已转换 | 第 1～14 章（PDF 第 11～560 页） |
| 章数 | 14 |
| 二级 Section 数 | 96 |
| Markdown 文件数 | 110（96 个 Section + 14 个章首 intro） |
| 章节内图片数 | 339 |

## 目录结构

```
试验设计与分析翻译计划/
├─ Design and analysis of experiments -- Douglas C_ Montgomery.pdf   原书（第 10 版，扫描件；未入库）
├─ DesignandAnalysisofExperiments9thEdition.pdf   第 9 版（有文本层，用于交叉核对；未入库）
├─ .github/workflows/deploy.yml      GitHub Actions：构建并发布到 GitHub Pages
├─ README.md                         本文件
├─ Chapter 01/                        … Chapter 14/
│  ├─ 01_Introduction.md              该章完整 Markdown（MinerU 输出）
│  ├─ images/                         该章引用的图片
│  ├─ sections/                       按 Section 切分的结果
│  │  ├─ 00_Introduction.md           章首（书名标题 + 本章学习目标）
│  │  ├─ 01_1.1_Strategy_of_Experimentation.md
│  │  ├─ 02_1.2_Some_Typical_Applications_of_Experimental_Design.md
│  │  ├─ …
│  │  └─ sections_manifest.json       本节目录（Order/Kind/Label/Title/File/ImageCount）
│  └─ translations_zh/                逐节中文译文（文件名 = 同名英文文件 + _zh）
├─ MinerU-Skill/
│  └─ Design_and_analysis_of_experiments_--_Douglas_C_Montgomery_130e1e/
│     ├─ split_pdf/chapters/          拆分出的 15 个章 PDF（本地对照用，未入库）
│     ├─ markdown/                    MinerU 原始输出目录（图片被移动后仅剩未引用图片）
│     ├─ mineru_input.txt             供 MinerU 批处理的章 PDF 列表
│     └─ manifest.json                全部元数据：章页范围 + 每节页范围
├─ Image Gallery/                     出版商级高清插图（按章分目录，cNNfNNM = 图 N.MM）
├─ myst.yml                          Jupyter Book（MyST）项目配置
├─ toc.yml                           电子书目录
├─ requirements-book.txt             电子书依赖（jupyter-book）
├─ LXGWWenKai-Regular.ttf            电子书正文字体（稿子美楷体）
├─ start.bat                         一键本地预览（自动建 venv、装依赖并启动）
├─ .gitignore                        忽略 `_build/`、`.venv-book/`、**全部原书 PDF**、`start.bat` 与 `*.ps1` 等
├─ .venv-book/                       本地预览用的虚拟环境（start.bat 自动创建）
├─ book/                             电子书首页、翻译进度与关于页
│  ├─ index.md                       首页（已上线章节入口）
│  ├─ progress.md                    翻译进度
│  ├─ about.md                       关于本项目
│  └─ static/book.css                电子书样式（含 @font-face 字体声明）
├─ organize_markdown_by_chapter.ps1
├─ split_chapters_by_bookmarks.ps1
├─ split_markdown_by_sections.ps1
├─ _build/                            构建产物（已被 .gitignore 忽略）
├─ _tools/
│  ├─ dump_outlines.ps1               导出 PDF 书签树，便于核对
│  ├─ check_headings.ps1              核对 Markdown 标题 vs 书签
│  ├─ verify_output.ps1               端到端校验（结构 / 图片 / 无损拆分）
│  ├─ verify_split.ps1                校验拆分结果与清单一致性
│  ├─ clean_math.py                   清理 MinerU 在数学模式里留下的多余空格
│  └─ html_table_txt.py               把 MinerU 的单行 <table> blob 渲染成 Markdown 表
```

每个 Section 文件名格式为 `{序号}_{编号}_{标题}.md`，与 `slp3翻译计划` 保持一致；
`sections/` 内的图片引用会自动改写为 `../images/…`。

另有一套**一次性校验脚本**放在工作区根目录的 `tools/`（即 `d:\Code for VS\tools\`，不在本项目目录内）：
`_contentcheck.py` 逐字节核对图片与 `Image Gallery/` 是否一致，`_tblcheck*.py` 比对表格数值，
`_tbl826_zh.py` / `_tbl823_check.py` 等用于第 8 章的专项校验。这些脚本带 `_` 前缀，属辅助工具。

## 翻译

译文放在各章的 `translations_zh/` 中，命名与 `slp3翻译计划` 一致：与 `sections/` 中的英文文件同名，
再加 `_zh` 后缀（如 `01_1.1_Strategy_of_Experimentation_zh.md`）。格式约定：

- 每个文件首行用一级标题写该节中文标题（保留节号，如 `# 1.1 试验的策略`）；
  章首文件 `00_Introduction_zh.md` 在标题后附一行原书出处（版次、ISBN 与出版社链接）。
- 术语首次出现写作 `**中文**（English）`（如 `**部分因子试验**（fractional factorial experiment）`）。
- 图题在图片下方、表题在表格上方，均用粗体（`**图 1.5 …**`、`**表 1.1 …**`）；
  图片单独成段，引用沿用 `../images/…`。
- 公式统一改写为 `$$ … $$` 的 LaTeX（OCR 结果中的多余空格一并整理）。
- 脚注用 MyST 原生语法（`[^1]` / `[^1] …`）；扫描件中丢失的原书脚注改以“译者注”说明。
- 出现可疑数值时用**独立约束反算**校验（合计、平方和、$F$ 比、别名矩阵等），而不是只对照 OCR 文本。

每节译完后建议跑一遍 `python _tools/clean_math.py <文件…> --apply`：它会清理 MinerU 在数学模式里
留下的系统性空格（`9 0 + 1 0 0` → `90 + 100`、`\frac {1}` → `\frac{1}`、`S S _ {A}` → `SS_{A}` 等）。
不加 `--apply` 时只输出 diff。注意该脚本的替换是**大小写敏感**的，不要用 PowerShell `-replace` 代替
（它默认忽略大小写，会把处理组合 `a b` 误改成 `AB`）。

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
| 第 9～14 章 | 43 | 未开始 |

“英文文件数”含章首 intro 文件。

新增整章译文上线时，需要同步更新**四处**：`toc.yml`（加章条目与 children）、
`myst.yml` 中 `project.exclude` 的对应排除项、`book/progress.md` 的状态，
以及 `book/index.md` 的“已上线章节”列表与导语。

## 已知的原文（OCR）问题

翻译过程中发现的原文印刷或扫描问题记录如下。译文一般按原文译出，其中会影响理解的地方已加
“译者注”说明。

### 第 1～3 章

- **1.1 节**：原书脚注（淬火介质与炉次效应混杂）在扫描件中文字丢失，已按提供的原文补出。
- **2.5 节**：正文引用“表 2.5”，但硬度数据实际在表 2.6。
- **2.6 节**：正文写“例 2.2”，实际为**例 2.3**；同处公式印作 0.225，按上下文应为 0.255。
- **3.3 节**：两条脚注（指向第 3 章补充材料）文字丢失，已按提供的原文补出；
  例 3.3 的置信区间计算与处理效应估计被 OCR 交错，译文已重排。
- **3.4 节**：多处沿用较早版本的“拉伸强度”表述（该章例 3.1 实为刻蚀速率）；
  例 3.4 中原文写作“all five variances”（五个方差），但试验只有四个处理，正文已改为四个方差并加注。
- **3.5 节**：图 3.11 及正文中的 $MS_E$ 印作 330.70，应为 333.70（= 表 3.4 的误差均方）。
- **3.6 节**：Design-Expert 输出中 Pure Error 平方和印作 5338.20（译文按原文保留），
  但按表内一致性、误差均方 × 16 以及同节 Minitab/JMP 输出，应为 5339.20，已加注。
- **3.7 节**：均值列表中末一个印作 $\mu_{1}=675$，应为 $\mu_{4}=675$。
- **3.8 节**：图 3.18 的图题在扫描件中丢失，已按原书补出
  （含尺度因子 $\sqrt{MS_E/n}=\sqrt{0.094/6}=0.125$）。
- **3.9 节**：脚注（关于 $\{\tau_i\}$ 相互独立、故固定效应模型的 $\sum\tau_i=0$ 约束不适用）在扫描件中
  文字丢失，已按原书补出；此外 $\sigma_\tau^2$ 的备择假设条件印作 $\sigma_\tau^2=0$（应为 $>0$）、
  期望均方推导中的 $an\sigma_\tau^2$ 印作 $an^2$、$MS_E$ 印作 190（应为 1.90）、
  $\chi_{0.025,12}^2$ 印作 23,3367（应为 23.3367）。
- **3.10 节**：`N μ̂ = y..` 处原文括号不配对。

### 第 4～5 章

- **4.1 节**：表 4.5 之后一大段正文在扫描件中只剩末句，已按提供的原文补齐；该节脚注（重复测量设
  计，第 15 章）同样按原文补出；正文引用 `Section 3.9.2` 应为 **3.10.2**，缺失值算例中的
  `Equation 4.16` 应为 **4.20**。
- **4.2 节**：表 4.9 被 OCR 挪到 4.1 节末尾（并夹断了“处理与区组不正交”一句），已移回首次引用处；
  两处正文写 `Table 4.8`，应为**表 4.9**；表 4.13 的脚注标记印作 `4`，应为 `a`。
- **4.2 节**：原书两处都写作 `Example 4.3`（重复拉丁方、标准拉丁方），译文按原文保留并加说明注。
- **4.3 节**：表 4.19 中列平方和的总和记号印作 $y_{..}^2$，应为 $y_{..l.}$。
- **4.4 节**：式 4.41 应为 **4.40**；式 4.37 的系数 1/4 应为 $1/k$；$V(\tilde\tau_i)$ 的说明印作“区组内”，
  应为“区组间”；表 4.26 被 OCR 挪到 4.4.2 之后，已移回 4.4.1。含对照平方和
  $SS_c=k(\sum c_iQ_i)^2/(\lambda a\sum c_i^2)$ 的那一段已按提供的原文补齐。
- **5.3 节**：式 5.4 括号内印错；图 5.10a 第 9 个观测的残差印作 26.00，应为 **−6.00**；式 5.17 之后的
  拟合值印成 $\overline{y}_{j..}$；式 5.20 的方括号位置印错（可用例 5.2 的数值校验）。该节第 2 条脚注
  已按提供的原文译出。
- **5.5 节**：例 5.4 正文把温度写成 °C，应为 **°F**。
- **5.6 节**：表 5.24 第 3 天第 3 格印作 $G(\cdot)$，应为 $F(f_2g_3=95)$；表 5.23、5.24 被 OCR 插在
  正文之前。
- **表 5.16（刀具寿命）版式完全错乱**：用行/列合计、表 5.17–5.19 的 JMP 效应平方和
  （24.3333 / 25.3333 / 61.3333）、纯误差平方和 13.00000 等约束反解，复原出 18 个观测值，
  8 项独立校验全部通过。

### 第 6 章

- **6.3 节**：图 6.5 的分栏哈希切块图（原书整图被扫描切成几块）已删去；哈希图 `f9083d76…`
  实为**图 6.8**，已改用 `figure6.8.jpg` 并移到散度效应段之后（原书正文说的是“这些极差标在
  图 6.8 的立方体上”）。
- **6.3 节**：表 6.7 中约简模型的 `A-Gap` 95% 上限印作 28.10，按其标准误 10.42、df = 12、
  $t=2.179$ 应为 **−28.10**。
- **6.5 节**：表 6.14、6.15 的整块 JMP 输出与图 6.17 被 OCR 挪到正文之后，已按正文顺序重排；
  图 6.16 的图题被塞进表 6.14 的输出块且图片引用丢失，已按图片库补回；$\hat{y}=46.22$ 应为 **46.25**。
- **6.6 节**：正文 `Table 6.15` 应为**表 6.17**；`y ≠ ln y` 系 $y^{*}=\ln y$ 的 OCR 误读。
- **6.7 节**：模型式中 $\beta_{12}$ 被印成 $\beta_2$，已改正。
- **6.8 节**：正文 `Thus, in Table 6.22` 应为**表 6.24**；表 6.24 的图题写作 `Example 6.6`，
  应为**例 6.7**；图 6.37 与图 6.38 的图题在扫描件中完全相同。
- **6.9 节**：$t_0=-0.7935$ 后印作 `P = 0.76`，按 df = 4 复算约为 0.47，疑为扫描件笔误，
  译文按原文保留、未加注。

### 第 7 章

- **7.2 节**：源文件把 $SS_{Blocks}$ 的算式排在解释它的句子之前，译文按书序重排。
- **7.4 节**：图 7.1 的图题排在图片之前；**图 7.3 的图片引用丢失**，已改用 `figure7.3.jpg` 并补图题；
  正文一句重复文字已按提供的原文更正；JMP 输出中 `Parameter Estimates` 与 `Fixed Effect Tests`
  两个小标题在扫描件中互换，已按列标题改正。
- **7.6 节**：**图 7.6 的图片引用丢失**，该处只剩一段被严重破坏的四区组排列表，已按正文的
  $(L_1,L_2)$ 组合恢复并接入 `figure7.6.jpg`。
- **7.8 节**：例 7.3 的 $SS_{ABC}$、$SS_{AB}$ 算式在扫描件中数字被空格拆开，已清理。

### 第 8 章

- **8.3 节**：表 8.12（*Calculation of Dispersion Effects for Example 8.4*）在扫描件中被压缩成
  一个无表题、无列标题、无各次试验符号行的三行数组；已按原表恢复为完整表格（符号行由表 8.10
  的因子水平相乘得到，残差由回归式算得，三行汇总值用书中的原值）。
- **8.4 节**：表 8.14 被 OCR 插在例 8.5 的正文中间（原书此处跨页，残留“（续见第 300 页）”），
  已移回首次引用处；**图 8.18 的图片引用丢失**（只留下 `Block 1`、`Block 2` 两行文字），
  已接入 `figure8.18.jpg`。
- **8.6 节**：**表 8.26（别名矩阵，66 列）的扫描件 OCR 错误很多**（后半有 71 个空单元格、
  行宽被补成 34 而实际为 32）。译文按由表 8.25 的设计精确复算的
  $\mathbf{A}=(\mathbf{X}_1'\mathbf{X}_1)^{-1}\mathbf{X}_1'\mathbf{X}_2$ 列出，并与第 9 版 PDF 文本层逐格
  比对，**858 个单元格全部一致**；唯一差异是 $x_5$ 行、$x_3x_4$ 列原书印 0，复算为 0.2。
- **8.6 节**：表 8.23 的 $k=27,\ N=28$ 部分原书为 9 行 × 3 栏（每栏 9 个符号），扫描件丢了 1 行
  并压成单栏；$k=35,\ N=36$ 一行末尾少一个符号。均已按原书恢复。
- **8.6 节**：合并估计表中 $A$、$B$、$C$ 三行的第二列都印成 $BD+CE+FG=19.15$（别名链一并重复），
  按 $[B]'$、$[C]'$ 复算应为 0.33 与 1.53；**表 8.22 之后“折叠设计的定义关系”整段在扫描件中丢失**
  （只剩半句“and for the second fraction, they are…”），已按提供的原文补出。
- **8.6 节**：图 8.24 被扫描切成两块（整图 + 一个哈希切片），哈希图与整图下半的相关系数为 0.68，
  判为重复裁切，只保留整图。
- **8.9 节**：正文 `Table 8.28` 应为**表 8.30**。

## 本地预览

电子书用 Jupyter Book 2 和 MyST `book-theme` 构建，与 `slp3翻译计划` 相同。

直接运行仓库根目录的 `start.bat` 即可：首次运行会自动创建 `.venv-book` 虚拟环境、安装
`requirements-book.txt` 中的依赖，随后启动本地服务器并自动打开 <http://localhost:3000/>；
在控制台按 `Ctrl+C` 可停止服务。

手工等效命令：

```powershell
python -m venv .venv-book
.\.venv-book\Scripts\python.exe -m pip install -r requirements-book.txt
.\.venv-book\Scripts\jupyter.exe book start
```

只要构建静态站点（不启服务器）：

```powershell
.\.venv-book\Scripts\jupyter.exe book build --html
```

输出在 `_build/html`。新增整章译文后，需要同步更新 `toc.yml`、`myst.yml` 里
`project.exclude` 中未译章节的排除项、`book/progress.md` 与 `book/index.md`（四处）。

## 在线发布（GitHub Pages）

线上地址：<https://souyerz1440.github.io/DAoE-translation/>

`.github/workflows/deploy.yml` 会在每次向默认分支（`master` / `main`）推送后自动构建并发布，
也可以在 Actions 页面用 **Run workflow** 手动触发。首次启用只需做一次：仓库
Settings → Pages → Build and deployment 把 **Source** 选为 **GitHub Actions**。

工作流与本地 `start.bat` 用的是同一套工具链：装 Node 22 + Python 3.13 →
`pip install -r requirements-book.txt` → `jupyter-book build --html` →
把 `_build/html` 作为 Pages 制品上传并发布。

**关键点：站点部署在子路径上，构建时必须告诉 MyST 站点的位置**，否则所有 CSS/JS/图片都会 404。
工作流靠环境变量完成：

```yaml
env:
  BASE_URL: /${{ github.event.repository.name }}   # → /DAoE-translation
```

本地要复现线上路径效果，可先在 PowerShell 里 `$env:BASE_URL='/DAoE-translation'` 再构建；
日后绑定自定义域名并把站点放在域名根目录，则把该值改成 `''`。

另有两处与发布相关的约定：

- 各章首页不再链接原书 PDF，改为给出原书出处：Douglas C. Montgomery,
  *Design and Analysis of Experiments*, 10th ed., Wiley, 2019（ISBN 978-1-119-49244-3）
  与[出版社页面](https://www.wiley.com/en-us/Design+and+Analysis+of+Experiments%2C+10th+Edition-p-9781119492443)。
- **所有原书 PDF 都不入库**，仓库里只保留译文、图片、构建配置与工具脚本：
  - 根目录两本原书 PDF（第 10 版扫描件 56.8 MB、第 9 版对照本 108.9 MB）：后者超过 GitHub
    单文件 100 MB 硬上限，且两者都已被 `myst.yml` 的 `project.exclude` 排除出构建。
  - `MinerU-Skill/…/split_pdf/chapters/` 下拆分出的 15 个章 PDF（46.9 MB）：早期各章首页曾用它们
    做“原始 PDF”对照链接，现已改为链接出版社页面，这些 PDF 只留本地。
  - 入库内容合计约 88 MB，最大单文件是正文字体 `LXGWWenKai-Regular.ttf`（24.3 MB）。

  细节见 `.gitignore` 内注释。

## 依赖

- **qpdf**（读书签、拆页）：`D:\Program Files\qpdf 12.3.2\bin\qpdf.exe`
- **MinerU**：装在 conda 环境 `aitest` 中，须先 `conda activate aitest`；使用需 Token 的 `extract`
  命令与 `vlm` 模型（`flash-extract` 免费模式受 10 MB / 20 页限制，无法处理整章）
- **PowerShell 5.1**
- **Jupyter Book 2**：由 `start.bat` 在 `.venv-book` 虚拟环境中自动安装（`jupyter-book==2.1.6`）
- **LXGW WenKai（稿子美楷体）**：仓库根目录的 `LXGWWenKai-Regular.ttf`（24 MB），作为电子书正文字体。
  通过 `myst.yml` 的 `project.static_files` 发布到站点根目录，在 `book/static/book.css` 中用
  `@font-face` 引用；缺字体时会回退到微软雅黑等系统字体
- **pdftotext / pdftoppm**（TeX Live 自带，`D:\texlive\2026\bin\windows\`）：从对照 PDF 取文本层、
  渲染页面图片

## 已知问题

1. **原书为低分辨率扫描件**（163 ppi，无文本层），OCR 结果受扫描质量限制：正文偶有字符错认
   （如 `coefficients` → `Ceffi-cients`），个别公式会出现多余空格。
2. **MinerU 漏检 3 个二级标题**：`3.6 Sample Computer Output`、`3.7 Determining Sample Size`、
   `4.4 Balanced Incomplete Block Designs`。这三处已对照原书页面补回标题，使切分边界正确。
3. **第 4.4 节开头一句在 OCR 中丢失**：该节正文目前从 “…appear together an equal number of
   times.” 开始，前面缺失约 1～2 句。需要时可用更高分辨率重新渲染该页后再补。
4. **MinerU 原始输出目录保留了未被 Markdown 引用的图片**（约 1900 个，31 MB）。这些是与章节内容
   无关的中间产物，不影响 `Chapter NN/` 下的成果；如需瘦身可直接删除
   `MinerU-Skill/…/markdown/images/`。
5. **书签排序**：原 PDF 第 7 章的二级书签顺序错乱（7.2 排在 7.1 之前），脚本已按节号重排。
6. **表格比正文更容易出错**：宽表、只有符号的表、以及被排版压成一行的表，OCR 往往丢字符或错位
   （第 8 章的表 8.12 / 8.23 / 8.26 都是实例）。处理办法是用**独立约束反算**校验（合计、平方和、
   $F$ 比、别名矩阵、因子水平相乘…），而不是只对照 OCR 文本。`_tools/html_table_txt.py`
   把 MinerU 的单行 `<table>` blob 渲染成可读表格，是誊抄的第一步。
7. **对照 PDF 是另一版本**：根目录的 `DesignandAnalysisofExperiments9thEdition.pdf` 是第 9 版（有文本层），
   而本项目的扫描底本是第 10 版。两版在第 8 章的表格上完全一致（已逐格核对），但**其他章节的
   表号与内容可能不同**，引用前需先与 MinerU 正文的数值交叉确认。取文本用：
   `pdftotext -layout -f <页> -l <页> -enc UTF-8 <PDF> out.txt`。

## 后续工作

- 逐节翻译其余各章，译文同样放到各章 `translations_zh/`。
- 视需要补上附录、Index 与第 15 章。
- 新章上线时同步更新四处：`toc.yml`、`myst.yml` 的 `exclude`、`book/progress.md`、`book/index.md`。
- 推送后 GitHub Actions 会自动重新发布（见“在线发布”一节）；视情况再考虑自定义域名。
