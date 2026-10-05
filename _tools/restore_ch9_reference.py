"""Restore Chapter 9 tables/formulas from the local ninth-edition text layer.

Input: pdftotext -layout DesignandAnalysisofExperiments9thEdition.pdf
       _build/ch9_reference_9th.txt
Run without arguments to audit; pass --apply to update the translations.
"""
from pathlib import Path
from fractions import Fraction
from itertools import combinations
import re
import sys

ROOT = Path(__file__).resolve().parents[1]
PAGES = (ROOT / '_build/ch9_reference_9th.txt').read_text(encoding='utf-8').split('\f')
ZH = ROOT / 'Chapter 09/translations_zh'


def table(page, number, width, count, numbered=False):
    body = PAGES[page - 1].split(f'T A B L E 9 . {number}', 1)[1]
    rows = []
    for line in body.splitlines():
        tokens = line.split()
        # Printer marks can precede/follow the final row.
        tokens = [t for t in tokens if t != 'k']
        if len(tokens) == width + int(numbered):
            values = tokens[int(numbered):]
            if all(v in ('1', '−1', '+', '−') for v in values):
                if numbered:
                    assert int(tokens[0]) == len(rows) + 1
                rows.append(values)
                if len(rows) == count:
                    break
    assert len(rows) == count, (number, len(rows))
    return rows


TABLES = {
    15: table(443, 15, 15, 16, True),
    37: table(461, 37, 9, 18),
    38: table(463, 38, 6, 12),
}


def alias_rows(text):
    rows = {}
    current = None
    for line in text.splitlines():
        line = line.strip()
        if re.match(r'\[A[B-J]\] =', line):
            current, rhs = line.split(' = ', 1)
            current = current[1:-1]
            rows[current] = rhs
        elif current and re.match(r'^(?:[+−]|DE\b)', line):
            rows[current] += ' ' + line
        elif current:
            current = None
    return {k: re.sub(r'\s+', ' ', v).replace('−', '-') for k, v in rows.items()}


ALIASES = {
    37: alias_rows(PAGES[460] + '\n' + PAGES[461].split('We see that', 1)[0]),
    38: alias_rows(PAGES[461].split('The alias', 1)[1]),
}
assert len(ALIASES[37]) == 8
assert len(ALIASES[38]) == 5


def solve(matrix, rhs):
    """Exact Gauss-Jordan elimination; no rounding in coefficient verification."""
    n = len(matrix)
    a = [[Fraction(v) for v in row] + [Fraction(rhs[i])] for i, row in enumerate(matrix)]
    for col in range(n):
        pivot = next(i for i in range(col, n) if a[i][col])
        a[col], a[pivot] = a[pivot], a[col]
        scale = a[col][col]
        a[col] = [x / scale for x in a[col]]
        for i in range(n):
            if i != col:
                scale = a[i][col]
                a[i] = [x - scale * y for x, y in zip(a[i], a[col])]
    return [row[-1] for row in a]


def verify_aliases(number, letters):
    design = [[1 if v in ('1', '+') else -1 for v in row] for row in TABLES[number]]
    model = [[1] + row + [row[0] * x for x in row[1:]] for row in design]
    mismatches = []
    for i, j in combinations(range(1, len(letters)), 2):
        term = letters[i] + letters[j]
        coefficients = solve(model, [r[i] * r[j] for r in design])
        assert all(c == 0 for c in coefficients[1:len(letters) + 1])
        for index, (name, expression) in enumerate(ALIASES[number].items()):
            found = re.search(r'([+-])\s*(\d+(?:\.\d+)?)?\s*' + term + r'\b', expression)
            printed = float(found[2] or 1) * (1 if found[1] == '+' else -1) if found else 0
            exact = coefficients[len(letters) + 1 + index]
            if abs(printed - float(exact)) > 0.00051:
                mismatches.append((name, term, printed, str(exact)))
    return mismatches


def markdown_table(number, letters):
    numbered = number == 15
    header = (['试验号'] if numbered else []) + list(letters)
    rows = [([str(i)] if numbered else []) + row for i, row in enumerate(TABLES[number], 1)]
    return '\n'.join('| ' + ' | '.join(r) + ' |' for r in [header, ['---'] * len(header)] + rows)


def formulas(number):
    blocks = []
    for name, expression in ALIASES[number].items():
        if number == 37 and name == 'AC':
            expression = expression.replace('+ 10.143 FG', '+ 0.143 FG')
        terms = re.findall(r'(?:[+-]\s*)?(?:\d+(?:\.\d+)?\s*)?[A-J]{2}', expression)
        assert re.sub(r'\s', '', ''.join(terms)) == re.sub(r'\s', '', expression)
        lines = [' '.join(terms[i:i + 6]) for i in range(0, len(terms), 6)]
        lines[0] = f'[{name}] &= ' + lines[0]
        lines[1:] = ['&' + s for s in lines[1:]]
        blocks.append('$$\n\\begin{aligned}\n' + '\\\\\n'.join(lines) + '\n\\end{aligned}\n$$')
    return '\n\n'.join(blocks)


def replace_once(text, old, new):
    assert text.count(old) == 1, old
    return text.replace(old, new, 1)


def update():
    file5 = next(ZH.glob('05_*.md'))
    text = file5.read_text(encoding='utf-8')
    text = replace_once(text, '![](../images/table9.15_restored.png)', markdown_table(15, 'ABCDEFGHJKLMNPQ'))
    text = re.sub(r'^\[\^1\]:.*$', '[^1]: 译者注：英文分节稿漏掉了整张表 9.15；现按第 9 版对照 PDF 第 443 页（印刷页 428）恢复全部 16 次试验及 15 个因子的表格。第 6 次试验的 G 列在该版原表中印为 1，译表照录；后文表 9.20 使用的是 D、E、H、K、M、Q 列。', text, flags=re.M)
    # The original figure is vector content, rendered separately from page 449.
    text, count = re.subn(r'<table><tr><td>Lock</td>.*?</table>', '![表 9.24 无完全混杂设计的 JMP 逐步回归输出](../images/figure9.14_restored.png)', text, count=1, flags=re.S)
    assert count == 1
    text = re.sub(r'^\[\^2\]:.*$', '[^2]: 译者注：英文分节稿把图 9.14 识别为无表号 HTML 表，并误合并了 Lock、Entered 列的单元格；现从第 9 版对照 PDF 第 449 页（印刷页 434）补回完整图像。表 9.24 仍为试验安排与响应数据。', text, flags=re.M)
    file5.write_text(text, encoding='utf-8')

    file6 = next(ZH.glob('06_*.md'))
    text = file6.read_text(encoding='utf-8')
    for number, letters in [(37, 'ABCDEFGHJ'), (38, 'ABCDEF')]:
        text = replace_once(text, f'![](../images/table9.{number}_restored.png)', markdown_table(number, letters))
    text = replace_once(text, '**表 9.37 的部分别名关系（原书扫描图）**\n\n![](../images/figure9.37_aliases_part1.png)\n\n![](../images/figure9.37_aliases_part2.png)', '**表 9.37 的部分别名关系**\n\n' + formulas(37))
    text = replace_once(text, '**表 9.38 的部分别名关系（原书扫描图）**\n\n![](../images/figure9.38_aliases.png)', '**表 9.38 的部分别名关系**\n\n' + formulas(38))
    text = replace_once(text, '两因子交互作用的部分别名关系见下列原书扫描图。', '两因子交互作用的部分别名关系见下列公式。')
    text = replace_once(text, '两因子交互作用的部分别名关系见扫描图。', '两因子交互作用的部分别名关系见下列公式。')
    notes = {
        3: '译者注：两组别名式已按第 9 版对照 PDF 第 461—462 页（印刷页 446—447）恢复为公式。该版 [AC] 式也印有 $+10.143FG$；按表 9.37 的全部 18 次试验，以截距、九个主效应及 $AB,AC,\\ldots,AJ$ 为模型列复算，FG 在 [AC] 中的系数为 $1/7\\approx0.143$，译文据此改为 $+0.143FG$。表 9.38 的五条别名式也已逐项复算，与原书一致。',
        4: '译者注：英文分节稿中表 9.37 只提取出 17 行，后几行还重复误识；现按第 9 版对照 PDF 第 461 页（印刷页 446）恢复全部 18 行、九个因子的表格。',
        5: '译者注：英文分节稿中表 9.38 只提取出 11 行；现按第 9 版对照 PDF 第 463 页（印刷页 448）恢复全部 12 行、六个因子的表格。',
    }
    for n, note in notes.items():
        text, count = re.subn(r'^\[\^' + str(n) + r'\]:.*$', lambda m: f'[^{n}]: {note}', text, flags=re.M)
        assert count == 1
    file6.write_text(text, encoding='utf-8')


def verify_translations():
    """Compare the rendered source tables with every PDF cell, including signs."""
    total = 0
    for number, prefix in [(15, '05_'), (37, '06_'), (38, '06_')]:
        text = next(ZH.glob(prefix + '*.md')).read_text(encoding='utf-8')
        body = text.split(f'**表 9.{number}**', 1)[1].split('\n\n', 1)[1]
        lines = body.split('\n\n', 1)[0].splitlines()
        rows = [[v.strip() for v in line.strip('|').split('|')] for line in lines[2:]]
        if number == 15:
            assert [r[0] for r in rows] == [str(i) for i in range(1, 17)]
            rows = [r[1:] for r in rows]
        assert rows == TABLES[number], number
        total += sum(map(len, rows))
    for file in ZH.glob('*.md'):
        text = file.read_text(encoding='utf-8')
        for ref in re.findall(r'!\[[^\]]*\]\(([^)]+)\)', text):
            assert (file.parent / ref).is_file(), ref
    print(f'All {total} table cells match the PDF; all Chapter 9 image references exist.')


if __name__ == '__main__':
    for n, rows in TABLES.items():
        print(f'Table 9.{n}: {len(rows)} rows x {len(rows[0])} factors')
    issues37 = verify_aliases(37, 'ABCDEFGHJ')
    issues38 = verify_aliases(38, 'ABCDEF')
    print('Table 9.37 alias discrepancies:', issues37)
    print('Table 9.38 alias discrepancies:', issues38)
    assert issues37 == [('AC', 'FG', 10.143, '1/7')]
    assert issues38 == []
    # Check LaTeX reconstruction before making any edits.
    for number in (37, 38):
        formulas(number)
    if '--apply' in sys.argv:
        update()
        print('Updated sections 9.5 and 9.6.')
    if '--verify' in sys.argv:
        verify_translations()
