"""Check exercise coverage, local links and independent numerical reference cases."""
from pathlib import Path
from urllib.parse import unquote
import importlib.util
import json
import math
import re

ROOT = Path(__file__).resolve().parents[1]
import os
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'_out/exercises/matplotlib'))
all_pages = [ROOT/'Exercises'/'index.md']
manifest=json.loads((ROOT/'Exercises/data/pages.json').read_text(encoding='utf8'))
expected_numbers=[f'{ch}.{n}' for ch,total in [(1,12),(2,60),(3,71),(4,59)] for n in range(1,total+1)]
assert [item['number'] for item in manifest]==expected_numbers
assert len({item['file'] for item in manifest})==202
toc=(ROOT/'toc.yml').read_text(encoding='utf8')
for chapter,count in [(1,12),(2,60),(3,71),(4,59)]:
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
        assert not re.search(r'[\x00-\x09\x0b\x0c\x0e-\x1f]',body),(single,'Unexpected control character in equations')
        assert '**答案来源**' in body and '[本章题目目录]' in body
        question=body.split('\n## 题目\n',1)[1].split('\n## 答案与解析\n',1)[0]
        answer=body.split('\n## 答案与解析\n',1)[1]
        parts=list(dict.fromkeys(re.findall(r'\*\*（([a-g])）\*\*',question)))
        if f'{chapter}.{n}' not in ['2.12','3.5','3.71']:
            answers=re.findall(r'^### （([a-g])）',answer,re.M)
            assert parts==answers,(single,parts,answers)
        for lineno,line in enumerate(body.splitlines(),1):
            assert len(re.findall(r'(?<!\\)\$',line))%2==0,(single,lineno,'Unbalanced math delimiters')
        assert toc.count(f'- file: {rel}\n')==1,rel
        pos=expected_numbers.index(f'{chapter}.{n}')
        if pos: assert f'[上一题 {expected_numbers[pos-1]}]' in body
        if pos<201: assert f'[下一题 {expected_numbers[pos+1]}]' in body
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

anova_keys=set()
for template in (ROOT/'_tools/exercise_templates').glob('chapter*.md'):
    for block in re.finditer(r'^## ([123]\.\d+) [^\n]+\n([\s\S]*?)(?=^## [123]\.\d+ |\Z)',template.read_text(encoding='utf8'),re.M):
        if re.search(r'\{\{(?:anova|summary|fsummary|planning|levene)\|',block[2]): anova_keys.add(block[1])
for key in anova_keys:
    c,n=map(int,key.split('.')); text=(ROOT/f'Exercises/translations_zh/chapter{c:02}/exercise-{c}-{n}.md').read_text(encoding='utf8').split('## 答案与解析',1)[1]
    for required in ['模型','假设','H_0','H_1','SS_A','SS_E','MS_A','MS_E','| 来源 |','自由度']:
        assert required in text,(key,required)
diagnostics=json.loads((ROOT/'Exercises/data/diagnostics.json').read_text(encoding='utf8'))
for item in diagnostics.values():
    assert item['actual_run_order_available'] is False
    assert item['sequence_basis']==('题干明确的炉内连续加热次序' if item['exercise']=='3.35' else '原表组内记录顺序')
    assert item['n']>1
assert diagnostics['2.35|paired']['n']==10
assert diagnostics['3.35|raw']['n']==15
assert diagnostics['3.35|raw']['within_group_order_available'] is True
print(f'Requirements: subquestions match; {len(anova_keys)} ANOVA solutions contain all five elements; {len(diagnostics)} diagnostic sets document order availability')

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
    for hours,required in [(10,111431),(8,174110)]:
        d=hours/168
        assert m['powera'](required,[0,d/2,d],15.6)>=.9
        assert m['powera'](required-1,[0,d/2,d],15.6)<.9
    random_method=m['anova_method']('3.34')
    assert '不是随机效应的 BLUP' in random_method
    assert r'\hat\mu_i=\bar y_i' not in random_method
    print('Independent ANOVA, corrected P values, ranks and minimum sample sizes verified')
else:
    print('Numerical recomputation skipped: requires scipy and matplotlib')

if importlib.util.find_spec('numpy'):
    import numpy as np
    source=json.loads((ROOT/'Exercises/data/chapter04-source.json').read_text(encoding='utf8'))
    results=json.loads((ROOT/'Exercises/data/chapter04-results.json').read_text(encoding='utf8'))
    # Closed-form balanced RCBD residuals independently check the regression code.
    for key in ['8','9','12','13','14','15','16','18','19','26']:
        y=np.array(source['rectangular'][key],float)
        e=y-y.mean(axis=0)-y.mean(axis=1)[:,None]+y.mean()
        assert np.isclose(np.sum(e**2),results[key+'{}']['sse'])
        assert results[key+'{}']['df']==(y.shape[0]-1)*(y.shape[1]-1)
    # Independent BIBD information-matrix and adjusted-total calculations.
    for key in ['45','47','48']:
        arr=np.array(source['rectangular'][key],float);inc=np.isfinite(arr).astype(float)
        a,b=arr.shape;r=inc.sum(axis=1)[0];k=inc.sum(axis=0)[0]
        assert np.all(inc.sum(axis=1)==r) and np.all(inc.sum(axis=0)==k)
        pairs=inc@inc.T;lam=pairs[0,1]
        assert np.all((pairs-np.diag(np.diag(pairs)))[~np.eye(a,dtype=bool)]==lam)
        Q=np.nansum(arr,axis=1)-inc@np.nansum(arr,axis=0)/k
        treatment=k*(Q@Q)/(lam*a)
        total=np.nansum(arr**2)-np.nansum(arr)**2/np.isfinite(arr).sum()
        blocks=np.sum(np.nansum(arr,axis=0)**2)/k-np.nansum(arr)**2/np.isfinite(arr).sum()
        assert np.isclose(total-blocks-treatment,results[key+'{}']['sse'])
    assert np.isclose(results['48{}']['sse'],3.25)
    assert results["8{'missing': [(1, 2), (3, 3)]}"]['df']==10
    assert results["27{'dropcol': 4}"]['df']==8
    # Verify all pairs in the published binary affine construction.
    blocks=[[x for x in range(8) if (u&x).bit_count()%2==v] for u in range(1,8) for v in [0,1]]
    assert all(sum(i in bb and j in bb for bb in blocks)==3 for i in range(8) for j in range(i+1,8))
    text=(ROOT/'Exercises/translations_zh/chapter04/exercise-4-58.md').read_text(encoding='utf8')
    assert '| 误差 | 15 | 25.75 | 1.71667 |' in text
    assert '18.1457%' in text
    print('Chapter 4: independent RCBD and BIBD formulas, missing-value degrees of freedom and construction verified')
print('PASS: 202 exercises')
