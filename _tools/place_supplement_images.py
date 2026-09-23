"""把补充材料里被引用的图片拷进对应章的 images/ 目录，并改成规范文件名。

命名约定（沿用项目既有做法）：
    有编号图  -> figureS<章号>.<图号>.jpg      （S 表示 Supplemental）
    无编号图  -> figureS<章号>.uf<序号>.jpg    （uf = unnumbered figure）

用法：
    python _tools/place_supplement_images.py [--apply]
"""

from __future__ import annotations

import hashlib
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SUPP = ROOT / "Supplemental Materials" / "images"

# 章号 -> [(源图片文件名, 目标文件名), ...]
PLAN: dict[int, list[tuple[str, str]]] = {
    1: [
        ("9f9f0bc3b5d529dda2009c1916db624fef5f6b843a2145d3636e82dd877dee83.jpg", "figureS1.1.jpg"),
        ("a1231827762779f3c456d59abe78c69d50c3a4eb2f8605538b9a8212a3058f21.jpg", "figureS1.2.jpg"),
        ("4fc574b1648ad6553274b020fe440826c4bafb96d03273c125dd8302a84d306a.jpg", "figureS1.3.jpg"),
    ],
    2: [
        ("9972983185c91ba02e305cb9c81e4f4f69c097b2b5f454ae3e627f2764969541.jpg", "figureS2.uf001.jpg"),
        ("ec6c893dab4887b25995c7ddb51c27153effff05757951c9025c18d26318746c.jpg", "figureS2.uf002.jpg"),
    ],
    6: [
        ("1b958fb4f88e208816311248330c6b2dcd92611e437361e46df62956337b0aee.jpg", "figureS6.uf001.jpg"),
        ("2a30b30995b9a63eafb7ac184bc1acb7de866b0ffec2e6924440590e2875f15d.jpg", "figureS6.uf002.jpg"),
        ("7a52b65cfa79d5103781872225b58afd3879ca6cd38e2e27a2746cfacd6f0163.jpg", "figureS6.uf003.jpg"),
    ],
}


def md5(path: Path) -> str:
    return hashlib.md5(path.read_bytes()).hexdigest()


def main() -> int:
    apply = "--apply" in sys.argv
    problems = 0

    for chapter, entries in PLAN.items():
        images_dir = ROOT / f"Chapter {chapter:02d}" / "images"
        if not images_dir.is_dir():
            print(f"[错误] 目标目录不存在：{images_dir}")
            problems += 1
            continue

        print(f"--- Chapter {chapter:02d} -> {images_dir.relative_to(ROOT)}")
        for source_name, target_name in entries:
            source = SUPP / source_name
            target = images_dir / target_name
            if not source.is_file():
                print(f"    [缺失] 源图片不存在：{source_name}")
                problems += 1
                continue
            if target.exists():
                same = md5(source) == md5(target)
                print(f"    [已存在] {target_name}  内容{'一致' if same else '★不一致★'}")
                if not same:
                    problems += 1
                continue
            print(f"    {target_name}  <-  {source_name[:16]}…  ({source.stat().st_size} B)")
            if apply:
                shutil.copy2(source, target)

    if not apply:
        print("\n（预览模式，未拷贝；加 --apply 生效）")
    elif problems:
        print(f"\n有 {problems} 个问题，请检查。")
    else:
        print("\n完成。")
    return 1 if problems else 0


if __name__ == "__main__":
    raise SystemExit(main())
