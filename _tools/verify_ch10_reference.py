"""Audit translated Chapter 10 against the local ninth edition and refit examples.

Requires pdftotext -layout output at _build/ch9_reference_9th.txt.
Uses exact rational arithmetic for all regression fits.
"""
from fractions import Fraction as F
from html.parser import HTMLParser
from pathlib import Path
import math
import re

ROOT = Path(__file__).resolve().parents[1]
ZH = ROOT / 'Chapter 10/translations_zh'
PAGES = (ROOT / '_build/ch9_reference_9th.txt').read_text(encoding='utf-8').split('\f')
TEXT = next(ZH.glob('03_*.md')).read_text(encoding='utf-8')


class TableParser(HTMLParser):
    def __init__(self):
        super().__init__()
        self.rows = []
        self.cell = None

    def handle_starttag(self, tag, attrs):
        if tag == 'tr':
            self.rows.append([])
        elif tag in ('td', 'th'):
            self.cell = ''

    def handle_data(self, data):
        if self.cell is not None:
            self.cell += data

    def handle_endtag(self, tag):
        if tag in ('td', 'th'):
            self.rows[-1].append(self.cell.strip())
            self.cell = None


def html_rows(block):
    parser = TableParser()
    parser.feed(block)
    return parser.rows


def table_after(caption):
    body = TEXT.split(caption, 1)[1]
    return html_rows(re.search(r'<table>.*?</table>', body, re.S)[0])


def reference_rows(page, marker, width, count):
    body = PAGES[page - 1].split(marker, 1)[1]
    rows = []
    for line in body.splitlines():
        tokens = line.replace('−', '-').split()
        if len(tokens) == width and tokens[0].isdigit() and all(re.fullmatch(r'<?-?\d+(?:\.\d+)?', t) for t in tokens):
            assert int(tokens[0]) == len(rows) + 1
            rows.append(tokens)
            if len(rows) == count:
                break
    assert len(rows) == count
    return rows


def inverse(matrix):
    n = len(matrix)
    a = [[F(x) for x in row] + [F(i == j) for j in range(n)] for i, row in enumerate(matrix)]
    for j in range(n):
        pivot = next(i for i in range(j, n) if a[i][j])
        a[j], a[pivot] = a[pivot], a[j]
        scale = a[j][j]
        a[j] = [v / scale for v in a[j]]
        for i in range(n):
            if i != j:
                scale = a[i][j]
                a[i] = [v - scale * w for v, w in zip(a[i], a[j])]
    return [row[n:] for row in a]


def fit(x, y):
    x = [[F(v) for v in row] for row in x]
    y = [F(v) for v in y]
    n, p = len(x), len(x[0])
    gram = [[sum(row[i] * row[j] for row in x) for j in range(p)] for i in range(p)]
    cross = [sum(row[j] * value for row, value in zip(x, y)) for j in range(p)]
    inv = inverse(gram)
    beta = [sum(a * b for a, b in zip(row, cross)) for row in inv]
    predicted = [sum(a * b for a, b in zip(row, beta)) for row in x]
    residuals = [v - q for v, q in zip(y, predicted)]
    sse = sum(e * e for e in residuals)
    mse = sse / (n - p)
    leverage = [sum(row[i] * inv[i][j] * row[j] for i in range(p) for j in range(p)) for row in x]
    press = sum((e / (1 - h)) ** 2 for e, h in zip(residuals, leverage))
    sst = sum(v * v for v in y) - sum(y) ** 2 / n
    return dict(gram=gram, cross=cross, inv=inv, beta=beta, predicted=predicted,
                residuals=residuals, mse=mse, sse=sse, leverage=leverage, press=press, sst=sst)


def compare_tables():
    total = 0
    for number, page, width, count, exceptions in [
        (2, 481, 4, 16, {}),
        (3, 481, 8, 16, {(13, 3): ('17.1', '7.1')}),
        (5, 486, 8, 12, {(5, 6): ('1.14', '1.4')}),
    ]:
        translated = [row for row in table_after(f'**表 10.{number}**') if row[0].isdigit()]
        original = reference_rows(page, f'T A B L E 10 . {number}', width, count)
        assert len(translated) == count
        for i, (r1, r2) in enumerate(zip(original, translated), 1):
            assert len(r1) == len(r2) == width
            for j, (a, b) in enumerate(zip(r1, r2)):
                if a != b:
                    assert exceptions.get((i, j)) == (a, b), (number, i, j, a, b)
            total += width
    print(f'{total} table cells checked against PDF; only two documented corrections differ.')


def numerical_audit():
    data = [r for r in table_after('**表 10.2**') if r[0].isdigit()]
    result = fit([[1, r[1], r[2]] for r in data], [r[3] for r in data])
    assert result['gram'] == [[16, 1458, 164], [1458, 133560, 14946], [164, 14946, 1726]]
    assert result['cross'] == [37577, 3429550, 385562]
    assert all(abs(float(a) - b) < 0.00001 for a, b in zip(result['beta'], [1566.07777, 7.62129, 8.58485]))
    assert abs(float(result['mse']) - 267.604) < 0.0005
    rows = [r for r in table_after('**表 10.3**') if r[0].isdigit()]
    for i, row in enumerate(rows):
        h, e = result['leverage'][i], result['residuals'][i]
        student = float(e) / math.sqrt(float(result['mse'] * (1 - h)))
        cook = student ** 2 / 3 * float(h / (1 - h))
        mse_deleted = (result['sse'] - e ** 2 / (1 - h)) / 12
        external = float(e) / math.sqrt(float(mse_deleted * (1 - h)))
        values = [result['predicted'][i], e, h, student, cook, external]
        tolerances = [0.051, 0.051, 0.00051, 0.0051, 0.00051, 0.0051]
        for j, (v, tolerance) in enumerate(zip(values, tolerances), 2):
            target = float(row[j].lstrip('<'))
            if row[j].startswith('<'):
                assert abs(float(v)) < target
            else:
                assert abs(float(v) - target) < tolerance, (i + 1, j, float(v), target)
    print('Example 10.1 refit and all 96 fitted/residual/diagnostic values verified.')
    print('Exact-data PRESS:', float(result['press']), 'predicted R^2:', float(1 - result['press'] / result['sst']))
    # Compare unrounded-data PRESS with PRESS based on the rounded table.
    rounded_press = sum((float(r[3]) / (1 - float(r[4]))) ** 2 for r in rows)
    print('Rounded-table PRESS:', rounded_press)
    se = math.sqrt(float(result['mse'] * result['inv'][1][1]))
    print('Example 10.7 standard error:', se)
    reduced = fit([[1, r[1]] for r in data], [r[3] for r in data])
    extra = reduced['sse'] - result['sse']
    print('Example 10.6 extra SS:', float(extra), 'F:', float(extra / result['mse']))
    example2 = html_rows(re.search(r'<table>.*?</table>', TEXT.split('## 例 10.2', 1)[1], re.S)[0])
    example2 = [r for r in example2 if r[0].isdigit()]
    for row in example2:
        assert (F(row[1]) - 140) / 20 == F(row[4])
        assert (F(row[2]) - 60) / 20 == F(row[5])
        assert (F(row[3]) - F('22.5')) / F('7.5') == F(row[6])
    result2 = fit([[1] + r[4:7] for r in example2], [r[7] for r in example2])
    assert result2['beta'] == [51, F('5.625'), F('10.625'), F('1.125')]
    kept = [r for r in example2 if r[0] != '8']
    result3 = fit([[1] + r[4:7] for r in kept], [r[7] for r in kept])
    assert result3['beta'] == [F(2655, 52), F(297, 52), F(557, 52), F(63, 52)]
    example3 = TEXT.split('## 例 10.3', 1)[1].split('## 例 10.4', 1)[0]
    vector = re.search(r'\\hat\{\\boldsymbol\\beta\}\\approx\s*\\begin\{bmatrix\}(.*?)\\end\{bmatrix\}', example3, re.S)[1]
    translated_beta = [F(v.strip()) for v in vector.split('\\\\')]
    assert all(abs(a - b) < F('0.0000051') for a, b in zip(translated_beta, result3['beta']))
    exact_inverse = [[F(v, 104) for v in row] for row in [[10, 2, 2, 2], [2, 16, 3, 3], [2, 3, 16, 3], [2, 3, 3, 16]]]
    assert result3['inv'] == exact_inverse
    print('Example 10.3 corrected coefficients:', [float(v) for v in result3['beta']])
    rows4 = [r for r in table_after('**表 10.5**') if r[0].isdigit()]
    result4 = fit([[1] + r[4:7] for r in rows4], [r[7] for r in rows4])
    assert all(abs(float(a) - b) < 0.00001 for a, b in zip(result4['beta'], [50.49391, 5.40996, 10.16316, 1.07245]))
    print('Examples 10.2–10.4 coding, normal equations, and coefficients verified.')
    print('Example 10.3 inverse:', [[float(v) for v in r] for r in result3['inv']])
    print('Example 10.4 inverse:', [[float(v) for v in r] for r in result4['inv']])
    # Independent check of the block/interaction dependency in Example 10.5.
    base = [[int(v) for v in row] for row in re.findall(r'^\| [1-8] \| (.*?) \|$', TEXT, re.M)
            for row in [row.replace('−', '-').split(' | ')]]
    assert len(base) == 8
    new = [1, -1, -1, -1, 1, 1, -1]
    for row, block in [(r, -1) for r in base] + [(new, 1)]:
        assert block == row[5] - row[6] - row[0]
    added = [[-1, -1, -1, 1], [1, -1, -1, -1], [-1, 1, 1, 1], [1, 1, 1, -1]]
    augmented = [[1] + r + [r[0] * r[1], r[2] * r[3]] for r in added]
    assert all(sum(r[j] for r in augmented) == 0 for j in range(1, 7))
    block_model = [r + [-1] for r in base] + [r + [2] for r in augmented]
    assert all(sum(r[j] * r[7] for r in block_model) == 0 for j in range(7))
    gram = [[sum(r[i] * r[j] for r in block_model) for j in range(8)] for i in range(8)]
    inverse(gram)  # Full rank: block and both interactions are estimable.
    print('Example 10.5 block dependency and four-run augmentation verified.')


def structural_audit():
    tags = []
    for file in ZH.glob('*.md'):
        text = file.read_text(encoding='utf-8')
        tags.extend(re.findall(r'\\tag\{(10\.[^}]+)\}', text))
        definitions = re.findall(r'^\[\^([^]]+)\]:', text, re.M)
        references = re.findall(r'\[\^([^]]+)\](?!:)', text)
        assert set(definitions) == set(references), file.name
        assert len(definitions) == len(set(definitions)), file.name
        for image in re.findall(r'!\[[^]]*\]\(([^)]+)\)', text):
            assert (file.parent / image).is_file(), image
    expected = [f'10.{i}' for i in range(1, 63) if i != 9] + ['10.9a', '10.9b']
    assert sorted(tags) == sorted(expected)
    print('All 63 numbered equations, image references, and footnotes checked.')


if __name__ == '__main__':
    compare_tables()
    numerical_audit()
    structural_audit()
