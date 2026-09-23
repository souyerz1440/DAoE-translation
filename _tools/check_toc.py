"""校验 toc.yml：解析 YAML、递归收集条目、检查文件是否存在、统计各章条目数。

用法：
    python _tools/check_toc.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parent.parent


def walk(entries: list, out: list[tuple[str, str]], depth: int = 0) -> None:
    for entry in entries:
        if not isinstance(entry, dict):
            continue
        if "file" in entry:
            out.append((entry["file"], entry.get("title", "")))
        if entry.get("children"):
            walk(entry["children"], out, depth + 1)


def main() -> int:
    toc_path = ROOT / "toc.yml"
    data = yaml.safe_load(toc_path.read_text(encoding="utf-8"))
    entries = data["project"]["toc"]

    flat: list[tuple[str, str]] = []
    walk(entries, flat)

    print(f"toc.yml 解析成功，共 {len(flat)} 个页面条目。\n")

    missing: list[str] = []
    per_chapter: dict[str, int] = {}
    for rel, title in flat:
        path = (ROOT / rel.replace("/", "\\")).resolve()
        label = rel.split("/")[0]
        per_chapter[label] = per_chapter.get(label, 0) + 1
        if not path.is_file():
            missing.append(rel)

    print("各目录下的页面数：")
    for key in sorted(per_chapter):
        print(f"  {key:>14}: {per_chapter[key]}")

    supplements = [rel for rel, _ in flat if "补充材料" in rel]
    print(f"\n补充材料条目（{len(supplements)}）：")
    for rel in supplements:
        print(f"  {rel}")

    print()
    if missing:
        print(f"有 {len(missing)} 个条目指向的文件不存在：")
        for rel in missing:
            print(f"  - {rel}")
        return 1
    print("全部条目指向的文件都存在。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
