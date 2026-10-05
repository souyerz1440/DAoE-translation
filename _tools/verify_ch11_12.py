"""Validate translated chapter assets, equation numbering and key numerical examples."""
from pathlib import Path
import re,math
from urllib.parse import unquote

ROOT=Path(__file__).resolve().parents[1]
files=[p for c in (11,12) for p in (ROOT/f'Chapter {c}/translations_zh').glob('*.md')]
assert len(files)==16
tags=[];images=0;tables=0;cells=set()
for p in files:
    s=p.read_text(encoding='utf-8')
    assert not [c for c in s if ord(c)<32 and c not in '\n\t'],p
    assert s.count('$$')%2==0,p
    tags+=re.findall(r'\\tag\{([^}]+)\}',s)
    tables+=s.count('<table>')
    for a,b in re.findall(r'(?:!\[[^]]*\]\(([^)]+)\)|<img[^>]*src="([^"]+)")',s):
        assert (p.parent/unquote(a or b)).is_file(),(p,a or b)
        images+=1
    for foot in re.findall(r'\[\^([^]]+)\](?!:)',s):
        assert f'[^{foot}]:' in s,(p,foot)
    cells.update(t for t in re.findall(r'<td[^>]*>(.*?)</td>',s) if re.search(r'[A-Za-z]{4}',t) and '<img' not in t and '$' not in t)
assert sorted(tags)==sorted([f'11.{i}' for i in range(1,31)]+[f'12.{i}' for i in range(1,7)])

# Check the conference matrix and its defining orthogonality.
C=[[0,1,1,1,1,1],[1,0,1,-1,-1,1],[1,1,0,1,-1,-1],[1,-1,1,0,1,-1],[1,-1,-1,1,0,1],[1,1,-1,-1,1,0]]
assert all(sum(C[k][i]*C[k][j] for k in range(6))==(5 if i==j else 0) for i in range(6) for j in range(6))

# Validate the corrected Gaussian process kernel against all printed observations.
a=[-1943.3447961328,3941.78888206788,3488.57543918861,-2040.39522592773,-742.642897583584,519.91871208163,-3082.85411601115,958.926988711818,80.468182554262,-1180.44117607546]
x=[.0560573769818389,.0947,.0765974898313444,.0947,.0898402482375096,.0717377150616494,.0644873310121405,.0499,.0499,.0790747191607881]
r=[.0618,.0126487944665913,.0618,.0608005210868486,.0367246615426894,.0377241897055609,.0148210408248663,0,.0347687447931648,0]
expected=[338.07,1613.04,335.91,327.82,449.23,440.58,1173.82,1140.36,453.83,1261.39]
predicted=[734.545842514493+sum(b*math.exp(-65.4025404276544*(t-u)**2-3603.24827558717*(v-z)**2) for b,u,z in zip(a,x,r)) for t,v in zip(x,r)]
assert max(abs(p-y) for p,y in zip(predicted,expected))<.005
assert abs(sum(a*a for a in [40.3,40.5,40.7,40.2,40.6])-202.3**2/5-.172)<1e-8

# Check corrected supplemental means directly against the 72 observations.
s=(ROOT/'Chapter 12/translations_zh/S12_补充材料_zh.md').read_text(encoding='utf-8')
for line in s.splitlines():
    values=[v.strip() for v in line.strip('|').split('|')]
    if len(values)==11 and values[0].isdigit():
        y=list(map(float,values[1:9]));mean=float(values[9]);sn=float(values[10])
        assert abs(sum(y)/8-mean)<=.000500001,(values[0],mean)
        assert abs(-10*math.log10(sum(1/t**2 for t in y)/8)-sn)<.001,(values[0],sn)

print(f'PASS: {len(files)} files, {len(tags)} numbered equations, {images} image references, {tables} HTML tables.')
print(f'PASS: conference matrix, pure error, supplemental means and signal-to-noise ratios; GP maximum error {max(abs(p-y) for p,y in zip(predicted,expected)):.6f}.')
print('Remaining software labels:',sorted(cells))
