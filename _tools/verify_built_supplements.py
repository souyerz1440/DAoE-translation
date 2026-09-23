"""核对构建产物里 8 篇补充材料页面是否生成、内容是否齐全（含 ch06 的译者注）。"""

from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
HTML = ROOT / "_build" / "html"

EXPECT = {
    1: ["图 1", "图 3", "试验前计划"],
    2: ["配对 t 检验"],
    3: ["可估函数"],
    4: ["尤登方"],
    5: ["模型层次"],
    6: ["译者注", "正规方程"],
    7: ["运行顺序很重要"],
    8: ["Yates"],
}


def strip_tags(text: str) -> str:
    out, skip = [], False
    for ch in text:
        if ch == "<":
            skip = True
        elif ch == ">":
            skip = False
        elif not skip:
            out.append(ch)
    return "".join(out)


def main() -> int:
    problems: list[str] = []
    print(f"{'章':>3}  {'页面':>10}  {'大小':>8}  内容检查")

    for n in range(1, 9):
        page = HTML / f"chapter-{n:02d}" / "translations-zh" / f"s{n}-zh" / "index.html"
        if not page.is_file():
            print(f"{n:>3}  {'缺失':>10}")
            problems.append(f"ch{n:02d}: 未生成 {page.relative_to(ROOT)}")
            continue

        raw = page.read_text(encoding="utf-8", errors="replace")
        text = strip_tags(raw)
        found = [kw for kw in EXPECT[n] if kw in text]
        missing = [kw for kw in EXPECT[n] if kw not in text]
        flag = "OK" if not missing else f"缺少 {missing}"
        if missing:
            problems.append(f"ch{n:02d}: 页面中未找到 {missing}")
        print(f"{n:>3}  {'s'+str(n)+'-zh':>10}  {page.stat().st_size:>8,}  {flag}  ({'、'.join(found)})")

    # 首页与进度页
    index = (HTML / "index.html").read_text(encoding="utf-8", errors="replace")
    progress = (HTML / "book" / "progress" / "index.html").read_text(encoding="utf-8", errors="replace")
    print()
    print(f"首页出现“补充材料”：{'是' if '补充材料' in index else '否'}")
    print(f"进度页出现“补充材料”：{'是' if '补充材料' in progress else '否'}")
    if "补充材料" not in index:
        problems.append("首页未提到补充材料")
    if "补充材料" not in progress:
        problems.append("进度页未提到补充材料")

    # 侧边栏目录：MyST 把站点 toc 写进 config.json（不是 index.json！）
    config_json = HTML / "config.json"
    if config_json.is_file():
        import json

        blob = json.dumps(json.loads(config_json.read_text(encoding="utf-8")), ensure_ascii=False)
        n_supp = blob.count("补充材料")
        n_pages = sum(1 for k in range(1, 9) if f"s{k}-zh" in blob)
        print(f"侧边栏 config.json：“补充材料” 出现 {n_supp} 次；sN-zh 页面条目 {n_pages}/8")
        if n_pages < 8:
            problems.append(f"config.json 只含 {n_pages}/8 个补充材料条目")
        # 每个页面的侧边栏都是内联的，逐页确认列全了 8 篇
        for k in range(1, 9):
            page = HTML / f"chapter-{k:02d}" / "translations-zh" / f"s{k}-zh" / "index.html"
            if not page.is_file():
                continue
            page_text = page.read_text(encoding="utf-8", errors="replace")
            if not all(f"s{j}-zh" in page_text for j in range(1, 9)):
                problems.append(f"ch{k:02d}: 侧边栏未列全 8 篇补充材料")
        print("逐页侧边栏检查完成。")
    else:
        problems.append("未找到 _build/html/config.json，无法核对侧边栏")

    # ch06 的 5 条译者注是否渲染出来
    s6 = strip_tags(
        (HTML / "chapter-06" / "translations-zh" / "s6-zh" / "index.html").read_text(
            encoding="utf-8", errors="replace"
        )
    )
    notes = ["第 2 个正规方程的首项", "平方和”一列的表头", "增加第四次重复后", "漏掉了", "省略号后均漏了一个加号"]
    got = [kw for kw in notes if kw in s6]
    print(f"ch06 译者注渲染：{len(got)}/{len(notes)} 条命中 {got}")
    if len(got) != len(notes):
        problems.append(f"ch06 译者注未全部渲染（缺 {[k for k in notes if k not in got]}）")

    print()
    if problems:
        print(f"发现 {len(problems)} 个问题：")
        for p in problems:
            print(f"  - {p}")
        return 1
    print("8 篇补充材料页面全部生成且内容齐全。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
