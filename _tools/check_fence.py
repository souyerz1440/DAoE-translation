"""确认某个行号是否落在 ``` 围栏代码块内部（用于判断 clean_math 的告警是不是误报）。"""

from __future__ import annotations

import sys
from pathlib import Path

FENCE = "`" * 3


def fence_ranges(lines: list[str]) -> list[tuple[int, int]]:
    ranges: list[tuple[int, int]] = []
    start: int | None = None
    for i, line in enumerate(lines, 1):
        if line.strip().startswith(FENCE):
            if start is None:
                start = i
            else:
                ranges.append((start, i))
                start = None
    if start is not None:
        ranges.append((start, len(lines)))
    return ranges


def main() -> int:
    path = Path(sys.argv[1])
    targets = [int(x) for x in sys.argv[2:]]
    lines = path.read_text(encoding="utf-8").split("\n")
    ranges = fence_ranges(lines)

    print(f"文件：{path.name}")
    print(f"代码块区间：{ranges}")
    for t in targets:
        hit = next((r for r in ranges if r[0] <= t <= r[1]), None)
        if hit:
            print(f"  L{t:5d} 在代码块 {hit} 内  ->  clean_math 告警是误报，切勿 --apply")
        else:
            print(f"  L{t:5d} **不在**代码块内  ->  需人工确认")
            print(f"         {lines[t - 1][:90]!r}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
