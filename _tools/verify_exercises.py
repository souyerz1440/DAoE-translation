"""Check exercise coverage, local links and independent numerical reference cases."""
from pathlib import Path
from urllib.parse import unquote
import importlib.util
import json
import math
import re

ROOT = Path(__file__).resolve().parents[1]
all_pages = [ROOT/'Exercises'/'index.md']
manifest=json.loads((ROOT/'Exercises/data/pages.json').read_text(encoding='utf8'))
expected_numbers=[f'{ch}.{n}' for ch,total in [(1,12),(2,60),(3,71)] for n in range(1,total+1)]
assert [item['number'] for item in manifest]==expected_numbers
assert len({item['file'] for item in manifest})==143
toc=(ROOT/'toc.yml').read_text(encoding='utf8')
for chapter,count in [(1,12),(2,60),(3,71)]:
    page=ROOT/'Exercises'/'translations_zh'/f'chapter{chapter:02}.md'
    text=page.read_text(encoding='utf8')
    numbers=[int(x) for x in re.findall(rf'^\| \[{chapter}\.(\d+)\]',text,re.M)]
    assert numbers==list(range(1,count+1)), (chapter,numbers)
    assert '## 答案与解析' not in text, 'Chapter pages must be indexes'
    assert '{{' not in text and '\ufffd' not in text
    all_pages.append(page)
    for n in range(1,count+1):
        rel=f'Exercises/translations_zh/chapter{chapter:02}/exercise-{chapter}-{n}.md'
        single=ROOT/rel; body=single.read_text(encoding='utf8')
        assert body.count('\n# 习题 ')==1
        assert body.count('\n## 题目\n')==1 and body.count('\n## 答案与解析\n')==1
        assert ':::{admonition}' not in body and '{{' not in body and '\ufffd' not in body
        assert '**答案来源**' in body and '[本章题目目录]' in body
        assert toc.count(f'- file: {rel}\n')==1,rel
        pos=expected_numbers.index(f'{chapter}.{n}')
        if pos: assert f'[上一题 {expected_numbers[pos-1]}]' in body
        if pos<142: assert f'[下一题 {expected_numbers[pos+1]}]' in body
        all_pages.append(single)
    print(f'Chapter {chapter}: {count} independent pages, directory entries and previous/next links')

for key,parts in [('2-31',7),('3-12',2),('3-62',6)]:
    chapter=int(key[0]); path=ROOT/f'Exercises/translations_zh/chapter{chapter:02}/exercise-{key}.md'
    statement=path.read_text(encoding='utf8').split('## 题目\n',1)[1].split('## 答案与解析',1)[0]
    assert len(re.findall(r'\*\*（[a-g]）\*\*',statement))==parts
    assert '第\n' not in statement, 'Inline references must not become new subquestions'
assert '溶剂' in (ROOT/'Exercises/translations_zh/chapter02/exercise-2-31.md').read_text(encoding='utf8')
assert '桡骨远端骨折' in (ROOT/'Exercises/translations_zh/chapter03/exercise-3-10.md').read_text(encoding='utf8')
print('Full backgrounds, source data for continuation questions and subquestion labels checked')

for page in all_pages:
    for url in re.findall(r'\]\(([^)]+)\)',page.read_text(encoding='utf8')):
        if '://' in url or url.startswith('#'): continue
        target=(page.parent/unquote(url.split('#')[0])).resolve()
        assert target.is_file(), (page,url)
print('All exercise links and figures exist')

css=ROOT/'book/static/katex/katex.min.css'
for url in re.findall(r'url\(([^)]+)\)',css.read_text(encoding='utf8')):
    assert (css.parent/url).is_file(),url
assert (css.parent/'LICENSE').is_file()
assert './katex/katex.min.css' in (ROOT/'book/static/book.css').read_text(encoding='utf8')
print('Local KaTeX stylesheet, font references and license verified')

data=json.loads((ROOT/'Exercises/data/chapter01-03.json').read_text(encoding='utf8'))
assert len(data)==43, len(data)
assert [len(g) for g in data['3.35']]==[5,4,6]
assert [len(g) for g in data['3.20']]==[5,4,5,4]
assert data['3.22'][3]==[129,127,132,129]
assert data['2.46'][1]==[16.62,16.75,17.37,17.12,16.98,16.87,17.34,17.02,17.08,17.27]
for groups in data.values():
    for group in groups: assert all(isinstance(x,(int,float)) and math.isfinite(x) for x in group)
print('Source data: all 43 data sets checked; unbalanced group sizes preserved')

if importlib.util.find_spec('scipy') and importlib.util.find_spec('matplotlib'):
    import runpy
    from scipy import stats
    import numpy as np
    m=runpy.run_path(str(ROOT/'_tools/build_exercises.py'))
    for key,expected in [('3.11',12.728107492613294),('3.26',16.082840236686387),('3.58',66.79707321945857)]:
        # scipy's independent implementation versus the rendering helper.
        actual=stats.f_oneway(*data[key]).statistic
        assert np.isclose(actual,expected) and np.isclose(actual,m['anova'](data[key])['F'])
    assert np.isclose(stats.f.sf(2.91,2,12),.09324837624538171)
    assert np.isclose(stats.kruskal(*data['3.28']).statistic,14.931024096385547)
    for variance,required in [(25,5),(36,7),(49,8)]:
        assert m['powera'](required,[50,60,50,60],variance)>=.9
        assert m['powera'](required-1,[50,60,50,60],variance)<.9
    for key,effect,required,one in [('2.28',1,1825,False),('2.31',1.5,28,True)]:
        sd=np.sqrt(m['anova'](data[key])['mse'])
        assert m['power2'](required,effect,sd,one=one)>=.9
        assert m['power2'](required-1,effect,sd,one=one)<.9
    print('Independent ANOVA, corrected P values, ranks and minimum sample sizes verified')
else:
    print('Numerical recomputation skipped: requires scipy and matplotlib')

print('PASS: 143 exercises')
