"""Check translated assets, numbering and calculations against data in the files."""
from pathlib import Path
from urllib.parse import unquote
from collections import Counter
import re, math
ROOT=Path(__file__).resolve().parents[1]
files=[p for c in (13,14,15) for p in (ROOT/f'Chapter {c}/translations_zh').glob('*.md')]
assert len(files)==21
tags=[];images=0;tables=0;labels=set()
for p in files:
 s=p.read_text(encoding='utf-8')
 assert s.count('$$')%2==0,p
 for m in re.findall(r'\$\$(.*?)\$\$',s,re.S):
  m=re.sub(r'\\[{}]', '', m)
  assert m.count('{')==m.count('}'),p
 assert not re.search(r'■\s*(图|表|TABLE|FIGURE)',s),p
 assert not any(ord(c)<32 and c not in '\n\t' for c in s),p
 refs=list(dict.fromkeys(re.findall(r'\[\^(\d+)\](?!:)',s)))
 defs=re.findall(r'^\[\^(\d+)\]:',s,re.M)
 assert refs==[str(i+1) for i in range(len(refs))],(p,refs)
 assert defs==refs,(p,refs,defs)
 tags+=re.findall(r'\\tag\{([^}]+)\}',s)
 tables+=s.count('<table>')
 for a,b in re.findall(r'(?:!\[[^]]*\]\(([^)]+)\)|<img[^>]*src="([^"]+)")',s):
  assert (p.parent/unquote(a or b)).is_file(),(p,a or b)
  images+=1
 for t in re.findall(r'<td[^>]*>(.*?)</td>',s,re.S):
  if re.search(r'[a-zA-Z]{4}',t) and '$' not in t and '<img' not in t:labels.add(t)
 source=ROOT/f'Chapter {p.parent.parent.name.split()[-1]}/sections'/p.name.replace('_zh.md','.md')
 if source.is_file():
  src=source.read_text(encoding='utf-8')
  for path in re.findall(r'!\[[^]]*\]\(([^)]+)\)',src):assert path in s,(p,path)
expected=[f'13.{i}' for i in range(1,28)]+[f'14.{i}' for i in range(1,21)]
for i in range(1,49):
 expected.extend([f'15.{i}{a}' for a in ('abc' if i==37 else 'ab')] if i in (37,38,41) else [f'15.{i}'])
assert Counter(tags)==Counter(expected),(Counter(expected)-Counter(tags),Counter(tags)-Counter(expected))
def content(ch,prefix):return next((ROOT/f'Chapter {ch}/translations_zh').glob(prefix+'*.md')).read_text(encoding='utf-8')
def htmlrows(table):return [re.findall(r'<td[^>]*>(.*?)</td>',r,re.S) for r in re.findall(r'<tr>(.*?)</tr>',table,re.S)]
def approx(a,b,tol=.005):assert abs(a-b)<tol,(a,b)
# ANCOVA: recover the 15 pairs from table 15.10 rather than rounded sums.
s=content(15,'03_');rows=htmlrows(re.findall(r'<table>.*?</table>',s,re.S)[0]);groups=[[],[],[]]
for row in rows[2:7]:
 for i in range(3):groups[i].append((float(row[2*i+1]),float(row[2*i])))
allpairs=sum(groups,[]);xs=sum(x for x,y in allpairs);ys=sum(y for x,y in allpairs)
Sxx=sum(x*x for x,y in allpairs)-xs*xs/15;Syy=sum(y*y for x,y in allpairs)-ys*ys/15;Sxy=sum(x*y for x,y in allpairs)-xs*ys/15
Txx=sum(sum(x for x,y in g)**2/5 for g in groups)-xs*xs/15
Tyy=sum(sum(y for x,y in g)**2/5 for g in groups)-ys*ys/15
Txy=sum(sum(x for x,y in g)*sum(y for x,y in g)/5 for g in groups)-xs*ys/15
for a,b in zip((Sxx,Syy,Sxy,Txx,Tyy,Txy),(261.733333333,346.4,282.6,66.133333333,140.4,96)):approx(a,b,1e-6)
Exx,Eyy,Exy=Sxx-Txx,Syy-Tyy,Sxy-Txy
sse=Eyy-Exy**2/Exx;f=((Syy-Sxy**2/Sxx)-sse)/2/(sse/11)
approx(sse,27.9858895706,1e-8);approx(f,2.609,.002)
# Staggered nesting: error is the sum of the paired-sample deviations.
s=content(14,'S14');pairs=[]
for a,b in re.findall(r'^\| \d+ \| ([\d.]+)，([\d.]+) \|',s,re.M):pairs.append((float(a),float(b)))
assert len(pairs)==10;approx(sum((a-b)**2/2 for a,b in pairs),5.62,1e-8)
# Supplemental unbalanced battery example: fit all nine cell means via normal equations.
s=content(15,'S15');rows=[]
for line in s.splitlines():
 v=[t.strip() for t in line.strip('|').split('|')]
 if len(v)==9 and all(re.fullmatch(r'\d+',t) for t in v):rows.append(list(map(float,v)))
assert len(rows)==31
for y,x1,x2,x3,x4,x5,x6,x7,x8 in rows:assert (x5,x6,x7,x8)==(x1*x3,x1*x4,x2*x3,x2*x4)
X=[[1]+r[1:] for r in rows];Y=[r[0] for r in rows]
A=[[sum(x[i]*x[j] for x in X) for j in range(9)]+[sum(x[i]*y for x,y in zip(X,Y))] for i in range(9)]
for i in range(9):
 pivot=max(range(i,9),key=lambda k:abs(A[k][i]));A[i],A[pivot]=A[pivot],A[i]
 scale=A[i][i];A[i]=[a/scale for a in A[i]]
 for j in range(9):
  if j!=i:
   scale=A[j][i];A[j]=[a-scale*b for a,b in zip(A[j],A[i])]
coef=[r[-1] for r in A]
approx(sum((y-sum(a*b for a,b in zip(x,coef)))**2 for x,y in zip(X,Y)),9553.833333333,1e-6)
approx(coef[0],155,1e-8)
print(f'PASS: {len(files)} files; {len(tags)} numbered equations; {images} image references; {tables} HTML tables; matching consecutive footnotes.')
print('PASS: source image preservation, ANCOVA sums and treatment test, staggered error, 31-observation regression fit and interaction encoding.')
print('Software labels retained:',sorted(labels))
