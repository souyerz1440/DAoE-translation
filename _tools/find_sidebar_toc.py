"""找出构建产物里承载「侧边栏目录」的文件，并核对其是否含 8 篇补充材料。"""

from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
HTML = ROOT / "_build" / "html"

SUPP = "补充材料"


def main() -> int:
    print("搜索同时含 s1-zh 与 s8-zh 的文件（即完整目录表）：\n")
    hits: list[tuple[Path, int]] = []
    for path in HTML.rglob("*"):
        if not path.is_file():
            continue
        if path.suffix.lower() not in {".html", ".json", ".js", ".txt"}:
            continue
        try:
            text = path.read_text(encoding="utf-8", errors="replace")
        except OSError:
            continue
        if "s1-zh" in text and "s8-zh" in text:
            hits.append((path, text.count(SUPP)))

    for path, n in hits:
        rel = path.relative_to(HTML)
        print(f"  {rel}   「{SUPP}」出现 {n} 次   大小 {path.stat().st_size:,}")

    if not hits:
        print("  （无）")
        return 1

    # 逐页检查侧边栏
    print("\n逐页检查（页面 HTML 中「补充材料」链接数）：")
    ok = 0
    for k in range(1, 9):
        page = HTML / f"chapter-{k:02d}" / "translations-zh" / f"s{k}-zh" / "index.html"
        if not page.is_file():
            print(f"  ch{k:02d}: 页面缺失")
            continue
        text = page.read_text(encoding="utf-8", errors="replace")
        n = text.count(SUPP)
        # 侧边栏应列出全部 8 篇
        has_all = all(f"s{j}-zh" in text for j in range(1, 9))
        print(f"  ch{k:02d}: 「{SUPP}」{n} 次；侧边栏含全部 8 篇 -> {'是' if has_all else '否'}")
        if has_all:
            ok += 1

    print(f"\n{ok}/8 个页面的侧边栏含全部 8 篇补充材料。")
    return 0 if ok == 8 else 1


if __name__ == "__main__":
    raise SystemExit(main())
