"""Recompute exercise tables/figures and expand the Chinese editorial templates.

Run from the repository root: python _tools/build_exercises.py
Requires numpy, scipy and matplotlib. Book builds consume the generated Markdown,
so these packages are not required for readers or GitHub Pages.
"""
from pathlib import Path
import itertools
import json
import re
import argparse
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'Exercises' / 'translations_zh'
OUT.mkdir(parents=True, exist_ok=True)
IMAGES = ROOT / 'Exercises' / 'images'
IMAGES.mkdir(exist_ok=True)
FONT = FontProperties(fname=str(ROOT / 'LXGWWenKai-Regular.ttf'))
plt.rcParams.update({'axes.unicode_minus': False, 'font.size': 10})
FIGURES_DONE=set()
DIAGNOSTICS={}

CHAPTERS = {1: ('引言', 12), 2: ('简单比较试验', 60), 3: ('单因子试验与方差分析', 71)}
TEXTBOOK_PAGES = {1: '21～22', 2: '55～63', 3: '125～134'}
BODY_LINKS = {
    1: 'Chapter%2001/translations_zh/00_Introduction_zh.md',
    2: 'Chapter%2002/translations_zh/00_Simple_Comparative_Experiments_zh.md',
    3: 'Chapter%2003/translations_zh/00_Experiments_with_a_Single_Factor_-_The_Analysis_of_Variance_zh.md',
}
RELATED = {
    '2.23':['2.22'], '2.25':['2.24'], '2.38':['2.37'], '2.41':['2.40'],
    '2.43':['2.42'], '2.46':[], '2.48':['2.28'], '2.49':['2.26'],
    '2.50':['2.26'], '2.52':['2.51'], '2.54':['2.53'],
    '3.12':['3.11'], '3.13':['3.11'], '3.15':['3.14'], '3.21':['3.20'],
    '3.23':['3.22'], '3.37':['3.36'], '3.40':['3.38'], '3.41':['3.39'],
    '3.43':['3.40'], '3.47':['3.30'], '3.48':['3.30'], '3.49':['3.26'],
    '3.50':['3.30'], '3.51':['3.30'], '3.53':['3.52'], '3.55':['3.54'],
    '3.57':['3.26'], '3.59':['3.27'], '3.60':['3.28'], '3.61':['3.58'],
    '3.71':['3.70'],
}

def exercise_file(ch, n):
    return f'chapter{ch:02}/exercise-{ch}-{n}.md'

def exercise_link(key, current_ch):
    ch, n = map(int, key.split('.'))
    return f'./exercise-{ch}-{n}.md' if ch == current_ch else f'../chapter{ch:02}/exercise-{ch}-{n}.md'

def link_references(text, current_ch):
    # Link exercise references only; leave numerical values, equations and section numbers alone.
    pattern = r'(?P<prefix>习题\s*|沿用\s*|参见\s*|重做\s*|见\s*|对\s*|与\s*)(?P<key>[123]\.[1-9]\d*)(?![\d.])'
    def replace(m):
        ch, n = map(int, m['key'].split('.'))
        if n > CHAPTERS[ch][1]: return m.group(0)
        return m['prefix'] + f'[{m["key"]}]({exercise_link(m["key"],current_ch)})'
    # Do not touch math, code spans or existing Markdown links.
    return ''.join(part if i % 2 else re.sub(pattern, replace, part)
                   for i, part in enumerate(re.split(r'(\$\$[\s\S]*?\$\$|\$[^$\n]*\$|`[^`]*`|!?\[[^\]]*\]\([^)]*\))', text)))

def write_pages(ch, text):
    chapter_title, count = CHAPTERS[ch]
    blocks = list(re.finditer(r'^## (\d+)\.(\d+) (.+)\n([\s\S]*?)(?=^## \d+\.\d+ |\Z)', text, re.M))
    assert len(blocks) == count
    target_dir = OUT / f'chapter{ch:02}'
    target_dir.mkdir(exist_ok=True)
    index_rows = []
    manifest = []
    all_keys = [f'{c}.{n}' for c,(_,total) in CHAPTERS.items() for n in range(1,total+1)]
    for block in blocks:
        n = int(block[2]); key = f'{ch}.{n}'; title = block[3]
        question, answer = block[4].split(':::{admonition} 答案与解析\n:class: dropdown\n',1)
        question = question.strip().removeprefix('**题目**：').strip()
        answer = answer.strip()
        assert answer.endswith(':::')
        answer = answer[:-3].strip()
        question = link_references(question,ch)
        answer = link_references(answer,ch)
        # A numbered subquestion gets its own paragraph, preserving the book's labels.
        question = re.sub(r'(?<!第)（([a-g])）',r'\n\n**（\1）** ',question).strip()
        # Present the source data after the background and before the questions.
        trailing_table = re.search(r'\n\n(\|[^\n]*\n(?:\|[^\n]*(?:\n|$))+)$',question)
        first_part = question.find('**（a）**')
        if trailing_table and 0 <= first_part < trailing_table.start():
            source_table = trailing_table[1].strip()
            question = question[:trailing_table.start()].strip()
            question = question[:first_part].rstrip()+'\n\n'+source_table+'\n\n'+question[first_part:]
        # Copies of the original source data keep continuation exercises readable on their own.
        for related in RELATED.get(key,[]):
            if related in DATA and not re.search(r'^\|',question,re.M):
                supplement = f'**所用数据（来自习题 [{related}]({exercise_link(related,ch)})）**\n\n' + raw_table(related)
                pos = question.find('**（a）**')
                if pos >= 0: question = question[:pos].rstrip()+'\n\n'+supplement+'\n\n'+question[pos:]
                else: question += '\n\n'+supplement
                break
        related_links = [f'[{ref}]({exercise_link(ref,ch)})' for ref in RELATED.get(key,[])]
        meta = f'第 9 版 · 第 {ch} 章 · 原书章末习题见印刷页 {TEXTBOOK_PAGES[ch]}。'
        position = all_keys.index(key); navigation = []
        if position: navigation.append(f'← [上一题 {all_keys[position-1]}]({exercise_link(all_keys[position-1],ch)})')
        navigation.append(f'[本章题目目录](../chapter{ch:02}.md)')
        if position+1 < len(all_keys): navigation.append(f'[下一题 {all_keys[position+1]}]({exercise_link(all_keys[position+1],ch)}) →')
        content = f'''---
title: 习题 {key}：{title}
short_title: {key} {title}
description: 第9版习题{key}的完整中文题目、数据与答案解析
---

(exercise-{ch}-{n})=
# 习题 {key}：{title}

[习题与答案](../../index.md) / [第 {ch} 章　{chapter_title}](../chapter{ch:02}.md)

{meta}

## 题目

{question}

## 答案与解析

{answer}

---

**答案来源**：{source_note(ch,n)}
'''
        if related_links: content += '\n**关联题目**：' + ' · '.join(related_links) + '\n'
        content += '\n' + ' · '.join(navigation) + '\n'
        # Individual pages sit one directory deeper than the former chapter-long pages.
        content = content.replace('(../images/', '(../../images/')
        (target_dir/f'exercise-{ch}-{n}.md').write_text(content,encoding='utf8')
        index_rows.append([f'[{key}](./{exercise_file(ch,n)})',title])
        manifest.append({'number':key,'title':title,'file':f'Exercises/translations_zh/{exercise_file(ch,n)}'})
    intro = f'''---
title: 第 {ch} 章习题：{chapter_title}
short_title: 第 {ch} 章习题
description: 第9版第{ch}章的{count}道习题目录，每道题提供独立的完整中文页面
---

# 第 {ch} 章习题：{chapter_title}

[全部章节](../index.md) · [本章正文](../../{BODY_LINKS[ch]})

本章收录第 **9 版 {ch}.1～{ch}.{count}，共 {count} 题**。点击题号进入独立页面，阅读完整题目背景、各小问、原始数据和中文答案解析。

## 题目目录

'''
    (OUT/f'chapter{ch:02}.md').write_text(intro+table(['题号','题目'],index_rows)+'\n',encoding='utf8')
    return manifest

def fmt(x):
    return f'{x:.6g}'

def table(headers, rows):
    return '\n'.join(['| ' + ' | '.join(headers) + ' |', '| ' + ' | '.join(['---']*len(headers)) + ' |'] + ['| ' + ' | '.join(map(str,r)) + ' |' for r in rows])

DATA = {
 '2.8': [[9.37,13.04,11.69,8.21,11.18,10.41,13.15,11.51,13.21,7.75]],
 '2.22': [[108,124,124,106,115,138,163,159,134,139]],
 '2.24': [[159,224,222,149,280,379,362,260,101,179,168,485,212,264,250,170]],
 '2.26': [[16.03,16.04,16.05,16.05,16.02,16.01,15.96,15.98,16.02,15.99],[16.02,15.97,15.96,16.01,15.99,16.03,16.04,16.02,16.01,16]],
 '2.28': [[65,81,57,66,82,82,67,59,75,70],[64,71,83,59,65,56,69,74,82,79]],
 '2.29': [[2.7,4.6,2.6,3,3.2,3.8],[4.6,3.4,2.9,3.5,4.1,5.1]],
 '2.31': [[11.176,7.089,8.097,11.739,11.291,10.759,6.467,8.315],[5.263,6.748,7.461,7.015,8.133,7.418,3.772,8.963]],
 '2.32': [[1,2,1,3,5,1,5,2,3,5,3,6,5,3,2,1,6,8,2,3],[7,8,5,9,5,8,6,4,6,7,6,9,5,7,4,6,8,5,8,7]],
 '2.33': [[5.34,6.65,4.76,5.98,7.25,6,7.55,5.54,5.62,6.21,5.97,7.35,5.44,4.39,4.98,5.25,6.35,4.61,6,5.32]],
 '2.34': [[.265,.265,.266,.267,.267,.265,.267,.267,.265,.268,.268,.265],[.264,.265,.264,.266,.267,.268,.264,.265,.265,.267,.268,.269]],
 '2.35': [[6.08,6.22,7.99,7.44,6.48,7.99,6.32,7.60,6.03,7.52],[5.73,5.80,8.42,6.84,6.43,8.76,6.32,7.62,6.59,7.67]],
 '2.36': [[1.186,1.151,1.322,1.339,1.2,1.402,1.365,1.537,1.559],[1.061,.992,1.063,1.062,1.065,1.178,1.037,1.086,1.052]],
 '2.37': [[206,193,192,188,207,210,205,185,194,187,189,178],[177,176,198,197,185,188,206,200,189,201,197,203]],
 '2.39': [[9.9,9.4,10,10.3,10.6,10.3,9.3,9.8],[10.2,10,10.7,10.5,10.6,10.2,10.4,10.3]],
 '2.46': [[16.85,16.40,17.21,16.35,16.52,17.04,16.96,17.15,16.59,16.57],[16.62,16.75,17.37,17.12,16.98,16.87,17.34,17.02,17.08,17.27]],
 '3.10': [[4.51,7.95,4.97,3,7.97,2.23,3.95,5.64,9.35,6.52,4.96,6.1,7.19,4.03,2.72,9.19,5.17,5.7,5.85,6.45],[5.32,6,5.12,7.08,5.48,6.52,4.09,6.28,7.77,5.68,8.47,4.58,4.11,5.72,5.91,6.89,6.99,4.98,9.94,6.38],[4.73,5.81,5.69,3.86,4.06,6.56,8.34,3.01,6.71,6.51,1.7,5.89,6.55,5.34,5.88,7.5,3.28,5.38,7.3,5.46],[7.03,4.65,6.65,5.49,6.98,4.85,7.26,5.92,5.58,7.91,4.9,4.54,8.18,5.42,6.03,7.04,5.17,7.6,7.9,7.91]],
 '3.11': [[3129,3000,2865,2890],[3200,3300,2975,3150],[2800,2900,2985,3050],[2600,2700,2600,2765]],
 '3.14': [[7,7,15,11,9],[12,17,12,18,18],[14,19,19,18,18],[19,25,22,19,23],[7,10,11,15,11]],
 '3.16': [[24,28,37,30],[37,44,31,35],[42,47,52,38]],
 '3.17': [[3,5,3,7,6,5,3,2,1,6],[1,3,4,7,5,6,3,2,1,7],[4,1,3,5,7,1,2,4,2,7],[3,5,7,5,10,3,4,7,2,7]],
 '3.18': [[83,85,85,87,90,88,88,84,91,90],[91,87,84,87,85,86,83],[94,91,87,85,87,91,92,86]],
 '3.19': [[1000,1500,1200,1800,1600,1100,1000,1250],[1500,1800,2000,1200,2000,1700,1800,1900],[900,1000,1200,1500,1200,1550,1000,1100]],
 '3.20': [[21.8,21.9,21.7,21.6,21.7],[21.7,21.4,21.5,21.4],[21.9,21.8,21.8,21.6,21.5],[21.9,21.7,21.8,21.4]],
 '3.22': [[143,141,150,146],[152,149,137,143],[134,136,132,127],[129,127,132,129]],
 '3.24': [[1530,1530,1440],[1610,1650,1500],[1560,1730,1530],[1500,1490,1510]],
 '3.25': [[80,83,83,85],[75,75,79,79],[74,73,76,77],[67,72,74,74],[62,62,67,69],[60,61,64,66]],
 '3.26': [[9,12,10,8,15],[20,21,23,17,30],[6,5,8,16,7]],
 '3.27': [[17.6,18.9,16.3,17.4,20.1,21.6],[16.9,15.3,18.6,17.1,19.5,20.3],[21.4,23.6,19.4,18.5,20.5,22.3],[19.3,21.1,16.9,17.5,18.3,19.8]],
 '3.28': [[19,20,19,30,8],[80,61,73,56,80],[47,26,25,35,50],[95,46,83,78,97]],
 '3.29': [[84.99,84.04,84.38],[85.15,85.13,84.88],[84.72,84.48,85.16],[84.20,84.10,84.55]],
 '3.30': [[100,96,92,96,92],[76,80,75,84,82],[108,100,96,98,100]],
 '3.31': [[58.2,57.2,58.4,55.8,54.9],[56.3,54.5,57,55.3],[50.1,54.2,55.4],[52.9,49.9,50,51.7]],
 '3.32': [[110,157,194,178],[1,2,4,18],[880,1256,5276,4355],[495,7040,5307,10050],[7,5,29,2]],
 '3.33': [[31,10,21,4,1],[62,40,24,30,35],[53,27,120,97,68]],
 '3.34': [[23.46,23.48,23.56,23.39,23.40],[23.59,23.46,23.42,23.49,23.50],[23.51,23.64,23.46,23.52,23.49],[23.28,23.40,23.37,23.46,23.39],[23.29,23.46,23.37,23.32,23.38]],
 '3.35': [[491.50,498.30,498.10,493.50,493.60],[488.50,484.65,479.90,477.35],[490.10,484.80,488.25,473,471.85,478.65]],
 '3.36': [[2.76,5.67,4.49],[1.43,1.70,2.19],[2.34,1.97,1.47],[.94,1.36,1.65]],
 '3.42': [[77.199,74.466,92.746,76.208,82.876],[80.522,79.306,81.914,80.346,73.385],[79.417,78.017,91.596,80.802,80.626],[78.001,78.358,77.544,77.364,77.386]],
 '3.54': [[4.93,4.86,4.75,4.95,4.79,4.88],[4.85,4.91,4.79,4.85,4.75,4.85],[4.83,4.88,4.9,4.75,4.82,4.9],[4.89,4.77,4.94,4.86,4.79,4.76]],
 '3.55': [[.05,.04,.05,.06,.03,.05],[.04,.02,.03,.05,.03,.02],[.09,.13,.11,.15,.08,.12],[.03,.04,.05,.05,.03,.02]],
 '3.56': [[.09,.10,.13,.08,.07],[.06,.09,.12,.07,.12],[.11,.08,.08,.05,.06],[.19,.13,.15,.20,.11]],
 '3.58': [[575,542,530,539,570],[565,593,590,579,610],[600,651,610,637,629],[725,700,715,685,710]],
 '3.62': [[14,14.1,14.2,14,14.1],[13.9,13.8,13.9,14,14],[14.1,14.2,14.1,14,13.9],[13.6,13.8,14,13.9,13.7],[13.8,13.6,13.9,13.8,14]],
}

def anova(groups):
    gs = [np.array(g,dtype=float) for g in groups]
    ns = np.array([len(g) for g in gs]); means=np.array([g.mean() for g in gs])
    grand=np.concatenate(gs).mean(); a=len(gs); N=sum(ns)
    ssa=sum(ns*(means-grand)**2); sse=sum(sum((g-g.mean())**2) for g in gs)
    mse=sse/(N-a); msa=ssa/(a-1); F=msa/mse
    return dict(ns=ns,means=means,grand=grand,a=a,N=N,ssa=ssa,sse=sse,mse=mse,msa=msa,F=F,p=stats.f.sf(F,a-1,N-a))

def analysis(key, transform='raw'):
    groups=[np.array(g,dtype=float) for g in DATA[key]]
    if transform=='log': groups=[np.log(g) for g in groups]
    if transform=='sqrt': groups=[np.sqrt(g) for g in groups]
    if transform=='neglog': groups=[-np.log(g) for g in groups]
    r=anova(groups); a=r['a']; N=r['N']; df=N-a
    txt=anova_method(key,transform)+'\n\n'+table(['来源','自由度','平方和','均方','$F$','$P$ 值'],[
        ['处理',a-1,fmt(r['ssa']),fmt(r['msa']),fmt(r['F']),fmt(r['p'])],
        ['误差',df,fmt(r['sse']),fmt(r['mse']),'—','—'],
        ['总计',N-1,fmt(r['ssa']+r['sse']),'—','—','—']])
    txt+='\n\n'+table(['组别','样本量','均值','样本标准差'],[[i+1,len(g),fmt(g.mean()),fmt(g.std(ddof=1))] for i,g in enumerate(groups)])
    return txt

def anova_method(key='', transform='raw', random=None):
    if random is None: random=key in ['3.34','3.35','3.36','3.37','3.40','3.41','3.42','3.43','3.62']
    response={'raw':'原响应','log':r'$z_{ij}=\ln(y_{ij})$','sqrt':r'$z_{ij}=\sqrt{y_{ij}}$','neglog':r'$z_{ij}=-\ln(s_{ij})$'}.get(transform,'原响应')
    if random:
        model=r'''**模型与模型假设（随机效应）**：$y_{ij}=\mu+u_i+\epsilon_{ij}$，$u_i\overset{\mathrm{iid}}\sim N(0,\sigma_\tau^2)$，$\epsilon_{ij}\overset{\mathrm{iid}}\sim N(0,\sigma^2)$；随机效应与误差相互独立，抽取的水平代表目标水平总体。不同组独立，同组观测的边际协方差为 $\sigma_\tau^2$，不能把同组观测的边际独立性当作假设。

**原假设和备择假设**：$H_0:\sigma_\tau^2=0$；$H_1:\sigma_\tau^2>0$。'''
    else:
        model=r'''**模型与模型假设（固定效应）**：$y_{ij}=\mu_i+\epsilon_{ij}=\mu+\tau_i+\epsilon_{ij}$，$i=1,\ldots,a$，$j=1,\ldots,n_i$；$\epsilon_{ij}\overset{\mathrm{iid}}\sim N(0,\sigma^2)$。误差独立、正态且同方差，处理分配及抽样符合题意；不平衡设计使用实际 $n_i$。取 $\sum_i n_i\tau_i=0$ 识别参数；平衡时等价于 $\sum_i\tau_i=0$。

**原假设和备择假设**：$H_0:\mu_1=\cdots=\mu_a$；$H_1$：至少一个总体均值不同。'''
    formulas=r'''

**平方和公式与检验统计量**：令 $N=\sum_i n_i$，$\bar y_i=\sum_jy_{ij}/n_i$，$\bar y_{\cdot\cdot}=\sum_i\sum_jy_{ij}/N$。

$$SS_A=\sum_i n_i(\bar y_i-\bar y_{\cdot\cdot})^2,\qquad SS_E=\sum_i\sum_j(y_{ij}-\bar y_i)^2,$$
$$SS_T=\sum_i\sum_j(y_{ij}-\bar y_{\cdot\cdot})^2=SS_A+SS_E,$$
$$MS_A=\frac{SS_A}{a-1},\quad MS_E=\frac{SS_E}{N-a},\quad F_0=\frac{MS_A}{MS_E}\overset{H_0}{\sim}F_{a-1,N-a}.$$

当 $F_0>F_{1-\alpha,a-1,N-a}$ 时拒绝 $H_0$，$P=P(F_{a-1,N-a}\ge F_0)$。本文分位数均为**下侧累积概率分位数**。

**点估计量**：$\hat\mu_i=\bar y_i$，$\hat\sigma^2=MS_E$；各组样本方差为 $s_i^2=\sum_j(y_{ij}-\bar y_i)^2/(n_i-1)$。以下为 ANOVA Table 与计算结果。'''
    if transform!='raw': model=f'本次分析的响应为 {response}；以下模型与平方和中的 $y$ 均代表该变换后的响应。\n\n'+model
    if random:
        formulas=formulas.replace(r'**点估计量**：$\hat\mu_i=\bar y_i$，$\hat\sigma^2=MS_E$；各组样本方差为 $s_i^2=\sum_j(y_{ij}-\bar y_i)^2/(n_i-1)$。',r'**方差分量估计量**：$\hat\sigma^2=MS_E$，$\hat\sigma_\tau^2=(MS_A-MS_E)/n_0$，其中 $n_0=[N-\sum_i n_i^2/N]/(a-1)$，平衡时 $n_0=n$。表中组样本均值仅作描述，不是随机效应的 BLUP；各组样本方差为 $s_i^2=\sum_j(y_{ij}-\bar y_i)^2/(n_i-1)$。')
    return model+formulas

def summary_anova(a,n,F=None,total=None,treatment=None,error=None,error_df=None,total_df=None,random=False):
    N=a*n if n else total_df+1
    dfa=a-1; dfe=N-a
    if total is not None and treatment is not None:
        ssa=treatment; sse=total-treatment; mse=sse/dfe; msa=ssa/dfa; f=msa/mse
        rows=[['处理',dfa,fmt(ssa),fmt(msa),fmt(f),fmt(stats.f.sf(f,dfa,dfe))],['误差',dfe,fmt(sse),fmt(mse),'—','—'],['总计',N-1,fmt(total),'—','—','—']]
        note=''
    else:
        rows=[['处理',dfa,f'{fmt(F*dfa)}m',f'{fmt(F)}m',fmt(F),fmt(stats.f.sf(F,dfa,dfe))],['误差',dfe,f'{dfe}m','$m$','—','—'],['总计',N-1,f'{fmt(F*dfa+dfe)}m','—','—','—']]
        note='原题未给原始观测或平方和，无法唯一恢复平方和的绝对数值。令 $m=MS_E>0$，下表保留其尺度，不虚构平方和。\n\n'
    return anova_method(random=random)+'\n\n'+note+table(['来源','自由度','平方和','均方','$F$','$P$ 值'],rows)

def planning_anova(a):
    return anova_method()+'\n\n样本量规划尚未产生观测值，实际平方和未知，以下 ANOVA Table 仅列结构；不能把非中心参数当作已观测平方和。\n\n'+table(['来源','自由度','平方和','均方'],[['处理',str(a-1),'$SS_A$',f'$SS_A/{a-1}$'],['误差',f'${a}(n-1)$','$SS_E$',f'$SS_E/[{a}(n-1)]$'],['总计',f'${a}n-1$','$SS_T$','—']])

def plot(key,groups,tag='raw'):
    fn=f'ch{key.replace(".","-")}-{tag}.png'
    caption=f'![习题 {key}：观测值、残差 Q–Q 图、残差记录顺序图、残差 vs 拟合值图（{tag}）](../images/{fn})'
    if fn in FIGURES_DONE: return caption
    fig,axes=plt.subplots(2,2,figsize=(11,7),layout='constrained'); axs=axes.ravel()
    res=np.concatenate([g-g.mean() for g in groups]); fitted=np.concatenate([np.full(len(g),g.mean()) for g in groups])
    rng=np.random.default_rng(20261102)
    for i,g in enumerate(groups):
        axs[0].scatter(np.full(len(g),i+1)+rng.uniform(-.09,.09,len(g)),g,s=17,alpha=.75)
        axs[0].plot([i+.75,i+1.25],[g.mean()]*2,color='#a33c27',lw=2)
    axs[0].set_xticks(range(1,len(groups)+1)); axs[0].set_xlabel('组别',fontproperties=FONT); axs[0].set_ylabel('响应值',fontproperties=FONT)
    axs[0].set_title('观测值与组均值',fontproperties=FONT)
    stats.probplot(res,plot=axs[1]); qqtitle='组内残差 Q–Q 图'; qqylabel='有序残差'
    axs[1].set_title(qqtitle,fontproperties=FONT); axs[1].set_xlabel('理论正态分位数',fontproperties=FONT); axs[1].set_ylabel(qqylabel,fontproperties=FONT)
    if key=='3.35':
        for i,g in enumerate(groups): axs[2].plot(range(1,len(g)+1),g-g.mean(),'o-',ms=4,lw=1,label=f'炉 {i+1}')
        axs[2].legend(prop=FONT); axs[2].set_xlabel('炉内连续加热序号（各炉分别）',fontproperties=FONT); axs[2].set_title('各炉残差时序图（非统一时间轴）',fontproperties=FONT)
    else:
        axs[2].plot(np.arange(1,len(res)+1),res,'o-',ms=3,lw=.8)
        for boundary in np.cumsum([len(g) for g in groups])[:-1]: axs[2].axvline(boundary+.5,color='gray',ls=':',lw=.7)
        axs[2].set_xlabel('记录序号（非实际运行时间顺序）',fontproperties=FONT); axs[2].set_title('残差记录顺序图；虚线为组界',fontproperties=FONT)
    axs[2].axhline(0,color='#a33c27',lw=1); axs[2].set_ylabel('组内残差',fontproperties=FONT)
    axs[3].scatter(fitted,res,s=20,alpha=.7); axs[3].axhline(0,color='#a33c27',lw=1); axs[3].set_xlabel('拟合值',fontproperties=FONT); axs[3].set_ylabel('残差',fontproperties=FONT); axs[3].set_title('残差 vs 拟合值',fontproperties=FONT)
    for ax in axs: ax.grid(alpha=.15)
    fig.savefig(IMAGES/fn,dpi=140); plt.close(fig)
    FIGURES_DONE.add(fn)
    return caption

def mean_plot(key):
    gs=[np.array(g) for g in DATA[key]]; r=anova(gs); fn=f'ch{key.replace(".","-")}-means.png'
    if fn not in FIGURES_DONE:
        fig,axes=plt.subplots(1,2,figsize=(10,3.8),layout='constrained'); rng=np.random.default_rng(20261103)
        for i,g in enumerate(gs): axes[0].scatter(np.full(len(g),i+1)+rng.uniform(-.07,.07,len(g)),g,s=18)
        half=stats.t.ppf(.975,r['N']-r['a'])*np.sqrt(r['mse']/r['ns']); axes[1].errorbar(range(1,len(gs)+1),r['means'],yerr=half,fmt='o',capsize=5)
        for ax,title in zip(axes,['观测值散点图','组均值及逐组 95% 区间']):
            ax.set_xlabel('组别',fontproperties=FONT); ax.set_ylabel('响应值',fontproperties=FONT); ax.set_title(title,fontproperties=FONT); ax.set_xticks(range(1,len(gs)+1)); ax.grid(alpha=.15)
        fig.savefig(IMAGES/fn,dpi=140); plt.close(fig); FIGURES_DONE.add(fn)
    return r'图中均值估计量为 $\bar y_i$，误差条为 $\bar y_i\pm t_{0.975,N-a}\sqrt{MS_E/n_i}$；均值区间的重叠情况不能直接代替均值差检验。'+f'\n\n![习题 {key}：散点图和处理均值区间图](../images/{fn})'

def diagnostic_view(key,view):
    gs=[np.array(g) for g in DATA[key]]; residual=np.concatenate([g-g.mean() for g in gs]); fitted=np.concatenate([np.full(len(g),g.mean()) for g in gs]); fn=f'ch{key.replace(".","-")}-{view}.png'
    titles={'qqres':'残差 Q–Q 图','sequence':'残差记录顺序图（非实际运行时序）','resfit':'残差 vs 拟合值图'}
    DIAGNOSTICS.setdefault(key+'|raw',{'exercise':key,'scale':'raw','actual_run_order_available':False,'sequence_basis':'原表组内记录顺序','residual':'y_ij - group_mean','n':len(residual)})
    if fn not in FIGURES_DONE:
        fig,ax=plt.subplots(figsize=(7,4),layout='constrained')
        if view=='qqres': stats.probplot(residual,plot=ax); ax.set_xlabel('理论正态分位数',fontproperties=FONT); ax.set_ylabel('有序残差',fontproperties=FONT)
        elif view=='sequence': ax.plot(range(1,len(residual)+1),residual,'o-',ms=4); ax.axhline(0,color='gray',lw=1); ax.set_xlabel('原表组内记录序号',fontproperties=FONT); ax.set_ylabel('残差',fontproperties=FONT)
        else: ax.scatter(fitted,residual); ax.axhline(0,color='gray',lw=1); ax.set_xlabel('拟合值',fontproperties=FONT); ax.set_ylabel('残差',fontproperties=FONT)
        ax.set_title(titles[view],fontproperties=FONT); fig.savefig(IMAGES/fn,dpi=140); plt.close(fig); FIGURES_DONE.add(fn)
    text=f'![习题 {key}：{titles[view]}](../images/{fn})'
    if view=='sequence': text+='\n\n原资料未给实际运行时间和随机运行次序。本图仅按原表记录顺序绘制，**不能当作实际残差时序图**，也不能用来证明时间独立性；取得真实运行顺序后须重画。'
    return text

def diagnostics(key,transform='raw'):
    gs=[np.array(g,dtype=float) for g in DATA[key]]
    if transform=='log': gs=[np.log(g) for g in gs]
    elif transform=='sqrt': gs=[np.sqrt(g) for g in gs]
    elif transform=='neglog': gs=[-np.log(g) for g in gs]
    elif transform=='paired': gs=[gs[0]-gs[1]]
    DIAGNOSTICS[key+'|'+transform]={'exercise':key,'scale':transform,'actual_run_order_available':False,'within_group_order_available':key=='3.35','sequence_basis':'题干明确的炉内连续加热次序' if key=='3.35' else '原表组内记录顺序','residual':'paired_difference - mean_difference' if transform=='paired' else 'transformed_y_ij - transformed_group_mean','n':sum(len(g) for g in gs)}
    residual=r'$\hat y_{ij}=\bar y_i$，$e_{ij}=y_{ij}-\bar y_i$'
    if transform=='paired': residual=r'对配对差 $d_i=y_{i1}-y_{i2}$ 拟合常数均值：$\hat d_i=\bar d$，$e_i=d_i-\bar d$'
    qualifier='随机效应题这里采用组内中心化残差检查组内误差；组均值是诊断基准，不冒充 REML 的 BLUP 拟合值。' if key in ['3.34','3.35','3.36','3.42','3.62'] else ''
    figures=plot(key,gs,transform)
    if key=='3.35':
        return f'**残差定义**：{residual}。{qualifier}\n\n'+figures+'\n\n**三项诊断**：右上为组内残差 Q–Q 图，左下按各炉连续加热序号分别画残差时序，右下为残差 vs 拟合值图。题干明确了各炉内部的连续加热先后，因此可检查炉内随顺序的变化；没有跨炉统一运行时间，不能把三炉拼成真实的共同时间轴，也不能据此验证跨炉时间关系。'
    return f'**残差定义**：{residual}。变换分析时使用相应变换尺度的响应。{qualifier}\n\n'+figures+r'''

**三项诊断**：Q–Q 图用于检查尾部、偏态及异常残差；残差记录顺序图用于展示原表记录的相邻变化；残差 vs 拟合值图用于观察方差是否随处理均值变化。

**时序资料的限制**：原题只提供分组观测，没有实际运行时间和随机运行次序。左下图按原表组内记录顺序排列，虚线表示组界，**不是实际残差时序图**；组界附近的变化不能解释为时间效应。若要完成真正的时序诊断，须取得原始运行顺序后按该顺序重画。当前不能仅凭此图断言误差在时间上独立。单样本或配对差的拟合值为常数时，残差 vs 拟合值图也不能检验复杂均值结构。'''

def qq_groups(key):
    gs=[np.array(g) for g in DATA[key]]; fn=f'ch{key.replace(".","-")}-groups-qq.png'
    fig,axes=plt.subplots(1,len(gs),figsize=(5*len(gs),3.8),squeeze=False,layout='constrained')
    for i,(g,ax) in enumerate(zip(gs,axes.ravel())):
        stats.probplot(g,plot=ax); ax.set_title(f'第 {i+1} 组正态 Q–Q 图',fontproperties=FONT); ax.set_xlabel('理论正态分位数',fontproperties=FONT); ax.set_ylabel('有序观测值',fontproperties=FONT)
    fig.savefig(IMAGES/fn,dpi=140); plt.close(fig)
    return f'![习题 {key}：各样本分别绘制的正态 Q–Q 图](../images/{fn})'

DATA_LABELS = {
    '2.26':['机器 1','机器 2'], '2.28':['配方 1','配方 2'],
    '2.29':['125 SCCM','200 SCCM'], '2.31':['95°C','100°C'],
    '2.32':['10 秒','20 秒'], '2.37':['配方 1','配方 2'],
    '2.39':['溶液 1','溶液 2'], '2.46':['改性砂浆','未改性砂浆'],
    '3.10':['假处理','PEMF 每日 1 小时','PEMF 每日 2 小时','PEMF 每日 4 小时'],
    '3.14':['15%','20%','25%','30%','35%'], '3.16':['20 g','30 g','40 g'],
    '3.17':['微型车','小型车','中型车','全尺寸车'], '3.18':['夏季','过渡季','冬季'],
    '3.20':['100','125','150','175'], '3.24':['插捣 10 次','插捣 15 次','插捣 20 次','插捣 25 次'],
    '3.25':['0.37','0.51','0.71','1.02','1.40','1.99'],
    '3.56':['10 in/min','12 in/min','14 in/min','16 in/min'],
    '3.58':['160 W','180 W','200 W','220 W'],
}

def raw_table(key):
    paired = {'2.34': ('检验员','卡尺 1','卡尺 2'),
              '2.35': ('双胞胎对','先出生','后出生'),
              '2.36': ('梁的顺序','Karlsruhe 法','Lehigh 法')}
    if key in paired:
        return table(list(paired[key]),[[i+1,fmt(x),fmt(y)] for i,(x,y) in enumerate(zip(*DATA[key]))])
    labels=DATA_LABELS.get(key,list(range(1,len(DATA[key])+1)))
    return table(['组别／处理条件','观测值（同一行属于同一组）'], [[labels[i],'，'.join(fmt(x) for x in g)] for i,g in enumerate(DATA[key])])

def pairs(key,alpha=.05,method='tukey'):
    r=anova(DATA[key]); df=r['N']-r['a']; rows=[]
    for i,j in itertools.combinations(range(r['a']),2):
        se=np.sqrt(r['mse']*(1/r['ns'][i]+1/r['ns'][j])); diff=r['means'][i]-r['means'][j]
        half=(stats.studentized_range.ppf(1-alpha,r['a'],df)/np.sqrt(2) if method=='tukey' else stats.t.ppf(1-alpha/2,df))*se
        rows.append([f'{i+1} − {j+1}',fmt(diff),f'({fmt(diff-half)}, {fmt(diff+half)})','是' if abs(diff)>half else '否'])
    formula=r'**估计量与区间公式**：$\widehat{\mu_i-\mu_j}=\bar y_i-\bar y_j$，$SE_{ij}=\sqrt{MS_E(1/n_i+1/n_j)}$。'
    if method=='tukey': formula+=r' Tukey–Kramer 同时区间为 $(\bar y_i-\bar y_j)\pm q_{1-\alpha,a,N-a}SE_{ij}/\sqrt2$；区间不含零则相应均值对显著，控制全部两两比较的家族错误率。'
    else: formula+=r' Fisher LSD 的逐对区间为 $(\bar y_i-\bar y_j)\pm t_{1-\alpha/2,N-a}SE_{ij}$；它不是全部比较的同时置信区间。'
    return formula+'\n\n'+table(['均值对','均值差',f'{100*(1-alpha):g}% 区间','显著'],rows)

def mean_ci(key,level=.95):
    r=anova(DATA[key]); q=stats.t.ppf((1+level)/2,r['N']-r['a'])
    return r'**估计量与区间公式**：$\hat\mu_i=\bar y_i$，$\widehat{SE}(\bar y_i)=\sqrt{MS_E/n_i}$，区间为 $\bar y_i\pm t_{1-\alpha/2,N-a}\sqrt{MS_E/n_i}$。以下各区间均为逐个区间。'+'\n\n'+table(['组别','均值点估计','估计标准误',f'{level*100:g}% 均值置信区间（逐个区间）'],[[i+1,fmt(m),fmt(np.sqrt(r['mse']/n)),f'({fmt(m-q*np.sqrt(r["mse"]/n))}, {fmt(m+q*np.sqrt(r["mse"]/n))})'] for i,(m,n) in enumerate(zip(r['means'],r['ns']))])

def two(key,paired=False,known=None,delta=0,alternative='two-sided',level=.95,welch=False,mode='all'):
    g1,g2=map(np.array,DATA[key]); n1=len(g1); n2=len(g2); diff=g1.mean()-g2.mean()
    if paired:
        d=g1-g2; se=d.std(ddof=1)/np.sqrt(n1); df=n1-1
    elif known:
        se=np.sqrt(known[0]**2/n1+known[1]**2/n2); df=np.inf
    elif welch:
        v1=g1.var(ddof=1)/n1; v2=g2.var(ddof=1)/n2; se=np.sqrt(v1+v2); df=(v1+v2)**2/(v1**2/(n1-1)+v2**2/(n2-1))
    else:
        sp2=((n1-1)*g1.var(ddof=1)+(n2-1)*g2.var(ddof=1))/(n1+n2-2); se=np.sqrt(sp2*(1/n1+1/n2)); df=n1+n2-2
    dist=stats.norm if np.isinf(df) else stats.t(df)
    t=(diff-delta)/se; p=2*dist.sf(abs(t)) if alternative=='two-sided' else (dist.sf(t) if alternative=='greater' else dist.cdf(t))
    q=dist.ppf((1+level)/2); lo=diff-q*se; hi=diff+q*se
    formula=r'$\hat\Delta=\bar y_1-\bar y_2$'
    if paired: formula=r'$d_i=y_{i1}-y_{i2}$，$\hat\Delta=\bar d$，$s_d^2=\sum_i(d_i-\bar d)^2/(n-1)$，$SE=s_d/\sqrt n$，$\nu=n-1$'
    elif known: formula+=r'，$SE=\sqrt{\sigma_1^2/n_1+\sigma_2^2/n_2}$，使用标准正态分布'
    elif welch: formula+=r'，$SE=\sqrt{s_1^2/n_1+s_2^2/n_2}$，$\nu=\dfrac{(s_1^2/n_1+s_2^2/n_2)^2}{(s_1^2/n_1)^2/(n_1-1)+(s_2^2/n_2)^2/(n_2-1)}$'
    else: formula+=r'，$s_p^2=\dfrac{(n_1-1)s_1^2+(n_2-1)s_2^2}{n_1+n_2-2}$，$SE=s_p\sqrt{1/n_1+1/n_2}$，$\nu=n_1+n_2-2$'
    estimate=f'**估计量公式**：{formula}。\n\n**计算结果**：两组均值为 {fmt(g1.mean())}、{fmt(g2.mean())}；差估计为 {fmt(diff)}，标准误为 {fmt(se)}。'
    distribution='标准正态分布' if known else f'自由度 {fmt(df)} 的 t 分布'
    test=f'检验统计量 ${"Z" if known else "t"}_0=(\\hat\\Delta-\\Delta_0)/SE={fmt(t)}$，按{distribution}取'+{'two-sided':'双侧','greater':'右侧','less':'左侧'}[alternative]+f'尾概率，$P={fmt(p)}$。'
    ci=r'**区间公式**：$\hat\Delta\pm '+(r'z_{1-\alpha/2}' if known else r't_{1-\alpha/2,\nu}')+r'SE$。'+f'\n\n**结果**：均值差的 {level*100:g}% 双侧置信区间为 $({fmt(lo)}, {fmt(hi)})$。'
    if mode=='p': return f'$P={fmt(p)}$。'
    return estimate+'\n\n'+(test if mode=='test' else ci if mode=='ci' else test+'\n\n'+ci)

def two_summary(n,m1,m2,s1,s2,welch=False):
    se=np.sqrt((s1*s1+s2*s2)/n); df=(s1*s1+s2*s2)**2/((s1**4+s2**4)/(n-1)) if welch else 2*n-2
    diff=m1-m2; t=diff/se; p=2*stats.t.sf(abs(t),df); h=stats.t.ppf(.975,df)*se
    model='两总体独立、正态，方差不假定相等，使用 Welch–Satterthwaite 近似。' if welch else '两总体独立、正态、同方差，使用合并方差 t 方法。'
    formula=r'$\hat\Delta=\bar y_1-\bar y_2$，$SE=\sqrt{s_1^2/n+s_2^2/n}$。'
    formula+=(r'$\nu=(s_1^2/n+s_2^2/n)^2/[(s_1^2/n)^2/(n-1)+(s_2^2/n)^2/(n-1)]$。' if welch else r'$s_p^2=[(n-1)s_1^2+(n-1)s_2^2]/(2n-2)$，$\nu=2n-2$。')
    return model+'\n\n**估计量与检验公式**：'+formula+r'$t_0=\hat\Delta/SE$ 检验 $H_0:\Delta=0$ 对 $H_1:\Delta\ne0$，区间公式为 $\hat\Delta\pm t_{0.975,\nu}SE$。'+f'\n\n**结果**：差估计 {fmt(diff)}，标准误 {fmt(se)}，自由度 {fmt(df)}，$t_0={fmt(t)}$，$P={fmt(p)}$；95% 区间为 $({fmt(diff-h)}, {fmt(diff+h)})$。'

def random_components(key):
    r=anova(DATA[key]); n0=(r['N']-sum(r['ns']**2)/r['N'])/(r['a']-1)
    v=(r['msa']-r['mse'])/n0; rho=v/(v+r['mse'])
    txt=r'**估计量公式**：$n_0=[N-\sum_i n_i^2/N]/(a-1)$（平衡时为 $n$），$\hat\sigma^2=MS_E$，$\hat\sigma_\tau^2=(MS_A-MS_E)/n_0$；$\widehat{\sigma_T^2}=\hat\sigma^2+\hat\sigma_\tau^2$，$\hat\rho=\hat\sigma_\tau^2/\widehat{\sigma_T^2}$。'+'\n\n'+f'**结果**：$n_0={fmt(n0)}$，矩估计为 $\\hat\\sigma^2={fmt(r["mse"])}$，$\\hat\\sigma_\\tau^2={fmt(v)}$。'
    if v>=0: txt+=f'总方差为 {fmt(v+r["mse"])}，组间方差占比为 {rho:.2%}。'
    else: txt+='组间分量的无约束矩估计为负，不能解释为有效的方差或方差占比。'
    if len(set(r['ns']))==1:
        L=(r['F']/stats.f.ppf(.975,r['a']-1,r['N']-r['a'])-1)/n0
        U=(r['F']/stats.f.ppf(.025,r['a']-1,r['N']-r['a'])-1)/n0
        txt+=r' 组内相关区间由 $r=MS_A/MS_E$，$L=(r/F_{0.975,a-1,N-a}-1)/n$，$U=(r/F_{0.025,a-1,N-a}-1)/n$，再将 $[L/(1+L),U/(1+U)]$ 与 $[0,1]$ 相交求得。'+f'结果为 $({fmt(max(0,L/(1+L)))}, {fmt(min(1,U/(1+U)))})$。'
    lo=(r['N']-r['a'])*r['mse']/stats.chi2.ppf(.975,r['N']-r['a']); hi=(r['N']-r['a'])*r['mse']/stats.chi2.ppf(.025,r['N']-r['a'])
    txt+=r' 误差方差区间公式为 $[(N-a)MS_E/\chi^2_{0.975,N-a},(N-a)MS_E/\chi^2_{0.025,N-a}]$，'+f'精确 95% 卡方区间结果为 $({fmt(lo)}, {fmt(hi)})$。'
    return txt

def power2(n,d,s,alpha=.05,one=False,known=False):
    nc=d/(s*np.sqrt(2/n))
    if known:
        q=stats.norm.ppf(1-alpha/(1 if one else 2)); return stats.norm.sf(q-nc)+(0 if one else stats.norm.cdf(-q-nc))
    df=2*n-2; q=stats.t.ppf(1-alpha/(1 if one else 2),df)
    return stats.nct.sf(q,df,nc)+(0 if one else stats.nct.cdf(-q,df,nc))

def powera(n,means,s2,alpha=.05):
    a=len(means); nc=n*sum((np.array(means)-np.mean(means))**2)/s2
    return stats.ncf.sf(stats.f.ppf(1-alpha,a-1,a*(n-1)),a-1,a*(n-1),nc)

def component(key,which):
    r=anova(DATA[key]); n0=(r['N']-sum(r['ns']**2)/r['N'])/(r['a']-1); v=(r['msa']-r['mse'])/n0
    if which=='tau': return r'$\hat\sigma_\tau^2=(MS_A-MS_E)/n_0$，其中 $n_0=[N-\sum_i n_i^2/N]/(a-1)$；'+f'$n_0={fmt(n0)}$，结果为 $\\hat\\sigma_\\tau^2={fmt(v)}$。'+('无约束估计为负，应限制在非负参数空间中重新估计，不能解释为负的实际方差。' if v<0 else '')
    if which=='error': return r'$\hat\sigma^2=MS_E=SS_E/(N-a)$，'+f'结果为 $\\hat\\sigma^2={fmt(r["mse"])}$。'
    if which=='total': return r'$\widehat{\sigma_T^2}=\hat\sigma_\tau^2+\hat\sigma^2$，'+f'结果为 {fmt(v+r["mse"])}。'
    if which=='rho': return r'$\hat\rho=\hat\sigma_\tau^2/(\hat\sigma_\tau^2+\hat\sigma^2)$，'+f'结果为 {v/(v+r["mse"]):.4%}。'
    if which=='rhoci':
        assert len(set(r['ns']))==1
        lo=(r['F']/stats.f.ppf(.975,r['a']-1,r['N']-r['a'])-1)/n0; hi=(r['F']/stats.f.ppf(.025,r['a']-1,r['N']-r['a'])-1)/n0
        return r'令 $r=MS_A/MS_E$、$\theta=\sigma_\tau^2/\sigma^2$，由 $r/(1+n\theta)\sim F_{a-1,N-a}$，得 $L=(r/F_{0.975,a-1,N-a}-1)/n$、$U=(r/F_{0.025,a-1,N-a}-1)/n$。对 $\rho=\theta/(1+\theta)$ 作单调变换，并截断到参数空间：'+f'\n\n$$CI_{{0.95}}(\\rho)=({fmt(max(0,lo/(1+lo)))}, {fmt(min(1,hi/(1+hi)))})$$'
    if which=='errorci':
        df=r['N']-r['a']; low=df*r['mse']/stats.chi2.ppf(.975,df); high=df*r['mse']/stats.chi2.ppf(.025,df)
        return r'$CI_{0.95}(\sigma^2)=[(N-a)MS_E/\chi^2_{0.975,N-a},(N-a)MS_E/\chi^2_{0.025,N-a}]$，'+f'结果为 $({fmt(low)},{fmt(high)})$。'
    raise ValueError(which)

def source_note(ch,n):
    if ch==1: return f'参考：`ch01.pdf`，'+(f'保留题 1.{n-10}R。' if n>10 else f'原题 1.{n}'+('S' if n in [1,9,10] else '')+'。')
    if ch==2:
        if n in [5,15,20,30,33,56,59]: return '参考答案文件无直接对应题；解析由编者补充。'
        if n<=4: src=n
        elif n<=14: src=n-1
        elif n<=19: src=n-2
        elif n<=29: src=n-3
        elif n<=32: src=n-4
        elif n<=55: src=n-5
        elif n<=58: src=n-6
        else: src=n-7
        suffix='S' if src in [2,4,9,12,19,26,28,31,34,43] else ''
    else:
        if n in [7,10,24,32,34,42,43,62,70]: return '参考答案文件无直接对应题；解析由编者补充。'
        if n<=6: src=n
        elif n<=9: src=n-1
        elif n<=23: src=n-2
        elif n<=31: src=n-3
        elif n==33: src=29
        elif n<=41: src=n-5
        elif n<=58: src=n-7
        elif n in [59,60]: src=52
        elif n==61: src=53
        elif n<=69: src=n-9
        else: src=61
        suffix='S' if src in [3,6,11,14,16,21,22,24,27,30,38,48] or n==59 else ''
    return f'参考：`ch0{ch}.pdf`，原题 {ch}.{src}{suffix}；本页按第 9 版编号。'

def expand(match):
    op,key,*args=match.group(1).split('|')
    if op=='data': return raw_table(key)
    if op=='anova': return analysis(key,args[0] if args else 'raw')
    if op=='pairs': return pairs(key,float(args[0]) if args else .05,args[1] if len(args)>1 else 'tukey')
    if op=='ci': return mean_ci(key,float(args[0]) if args else .95)
    if op=='random': return random_components(key)
    if op in ['plot','diagnostics']: return diagnostics(key,args[0] if args else 'raw')
    if op=='qq': return qq_groups(key)
    if op=='meansplot': return mean_plot(key)
    if op in ['qqres','sequence','resfit']: return diagnostic_view(key,op)
    if op=='model': return anova_method(key,args[0] if args else 'raw')
    if op=='component': return component(key,args[0])
    if op=='twosummary': return two_summary(int(key),*map(float,args[:4]),welch=len(args)>4 and args[4]=='welch')
    if op=='pvalue':
        r=anova(DATA[key]); return f'$P=P(F_{{{r["a"]-1},{r["N"]-r["a"]}}}\\ge{fmt(r["F"])})={fmt(r["p"])}$。'
    if op=='levene':
        groups=[np.abs(np.array(g)-np.median(g)) for g in DATA[key]]; r=anova(groups)
        text=anova_method().replace('**模型与模型假设（固定效应）**','**绝对偏差的工作模型与近似假设**')
        text=text.replace('以下为 ANOVA Table 与计算结果。','以下对 $z_{ij}=|y_{ij}-\\tilde y_i|$ 计算 ANOVA Table。F 校准为近似；绝对偏差未必严格正态，同组中位数估计还引入约束。')
        return '以下模型和平方和中的 $y$ 均指绝对偏差 $z$，不是原始电池寿命。\n\n'+text+'\n\n'+table(['来源','自由度','平方和','均方','$F$','$P$'],[['组别',r['a']-1,fmt(r['ssa']),fmt(r['msa']),fmt(r['F']),fmt(r['p'])],['误差',r['N']-r['a'],fmt(r['sse']),fmt(r['mse']),'—','—'],['总计',r['N']-1,fmt(r['ssa']+r['sse']),'—','—','—']])
    if op=='corrupted':
        gs=[g[:] for g in DATA[key]]; i,j=max(((i,j) for i,g in enumerate(gs) for j in range(len(g))),key=lambda ij:gs[ij[0]][ij[1]]); gs[i][j]=250; r=anova(gs)
        return table(['来源','自由度','平方和','均方','$F$','$P$'],[['处理',r['a']-1,fmt(r['ssa']),fmt(r['msa']),fmt(r['F']),fmt(r['p'])],['误差',r['N']-r['a'],fmt(r['sse']),fmt(r['mse']),'—','—'],['总计',r['N']-1,fmt(r['ssa']+r['sse']),'—','—','—']])
    if op=='fsummary': return summary_anova(int(key),int(args[0]),F=float(args[1]))
    if op=='planning': return planning_anova(int(key))
    if op=='summary': return summary_anova(int(key),int(args[0]),total=float(args[1]),treatment=float(args[2]),random=len(args)>3 and args[3]=='random')
    if op=='diffci':
        i,j=int(args[0])-1,int(args[1])-1; level=float(args[2]) if len(args)>2 else .95; r=anova(DATA[key]); diff=r['means'][i]-r['means'][j]; se=np.sqrt(r['mse']*(1/r['ns'][i]+1/r['ns'][j])); h=stats.t.ppf((1+level)/2,r['N']-r['a'])*se
        return r'$\widehat{\mu_i-\mu_j}=\bar y_i-\bar y_j$，$SE=\sqrt{MS_E(1/n_i+1/n_j)}$，区间为 $(\bar y_i-\bar y_j)\pm t_{1-\alpha/2,N-a}SE$。'+f'\n\n对第 {i+1} 组减第 {j+1} 组，点估计 {fmt(diff)}，标准误 {fmt(se)}；{100*level:g}% 区间为 $({fmt(diff-h)}, {fmt(diff+h)})$。'
    if op=='box':
        fn=f'ch{key.replace(".","-")}-box.png'; fig,ax=plt.subplots(figsize=(5.5,3.4),layout='constrained'); ax.boxplot(DATA[key],tick_labels=['125 SCCM','200 SCCM']); ax.set_ylabel('刻蚀均匀性',fontproperties=FONT); ax.set_title('两种流量的箱线图',fontproperties=FONT); fig.savefig(IMAGES/fn,dpi=140); plt.close(fig)
        return f'![125 与 200 SCCM 流量的刻蚀均匀性箱线图](../images/{fn})'
    if op=='wald':
        r=anova(DATA[key]); se=r['mse']*np.sqrt(2/(r['N']-r['a'])); return f'平衡设计且方差分量为正时，REML 点估计与上述 ANOVA 矩估计相同。作为明确的近似比较，误差方差的正态 Wald 区间为 $\\hat\\sigma^2\\pm1.959964\\hat\\sigma^2\\sqrt{{2/(N-a)}}$，得 $({fmt(max(0,r["mse"]-1.959964*se))}, {fmt(r["mse"]+1.959964*se)})$。它与精确卡方区间不同，且小样本覆盖率可能较差；软件若使用轮廓似然区间，也会得到不同端点。'
    if op=='normal':
        groups=[np.array(g) for g in DATA[key]]
        if args and args[0]=='paired': groups=[groups[0]-groups[1]]
        return table(['组别（配对时指差值）','Shapiro–Wilk 正态性检验 $P$ 值'],[[i+1,fmt(stats.shapiro(g).pvalue)] for i,g in enumerate(groups)])+'\n\n正态性检验不显著仅表示未发现明显偏离；须同时检查概率图、异常值及试验的独立性。'
    if op=='ftest':
        g1,g2=map(np.array,DATA[key]); f=g1.var(ddof=1)/g2.var(ddof=1); p=min(1,2*min(stats.f.cdf(f,len(g1)-1,len(g2)-1),stats.f.sf(f,len(g1)-1,len(g2)-1)))
        return f'方差比 $F={fmt(f)}$，自由度为 $({len(g1)-1},{len(g2)-1})$，双侧 $P={fmt(p)}$。'
    if op=='one':
        g=np.array(DATA[key][0]); mu=float(args[0]); level=float(args[1]) if len(args)>1 else .95; alt=args[2] if len(args)>2 else 'greater'; n=len(g); se=g.std(ddof=1)/np.sqrt(n); t=(g.mean()-mu)/se; p=2*stats.t.sf(abs(t),n-1) if alt=='two-sided' else stats.t.sf(t,n-1); half=stats.t.ppf((1+level)/2,n-1)*se
        mode=args[3] if len(args)>3 else 'all'
        if mode=='p': return f'$P={fmt(p)}$。'
        point=r'$\hat\mu=\bar y=\sum_i y_i/n$，$s^2=\sum_i(y_i-\bar y)^2/(n-1)$，$SE=s/\sqrt n$。'+f'$n={n}$，$\\bar y={fmt(g.mean())}$，$s={fmt(g.std(ddof=1))}$。'
        test=r'正态独立样本下，$t_0=(\bar y-\mu_0)/(s/\sqrt n)$，在 $H_0$ 下服从 $t_{n-1}$。'+f'统计量 {fmt(t)}，自由度 {n-1}，'+('双侧' if alt=='two-sided' else '右侧')+f'$P={fmt(p)}$。'
        ci=r'均值区间公式为 $\bar y\pm t_{1-\alpha/2,n-1}s/\sqrt n$。'+f'{100*level:g}% 双侧置信区间为 $({fmt(g.mean()-half)}, {fmt(g.mean()+half)})$。'
        return point+'\n\n'+(ci if mode=='ci' else test if mode=='test' else test+'\n\n'+ci)
    if op=='variance':
        g=np.array(DATA[key][0]); df=len(g)-1; s2=g.var(ddof=1); chi=df*s2; p=min(1,2*min(stats.chi2.cdf(chi,df),stats.chi2.sf(chi,df)))
        mode=args[0] if args else 'all'
        point=r'$\hat\sigma^2=s^2=\sum_i(y_i-\bar y)^2/(n-1)$，'+f'$s^2={fmt(s2)}$。'
        ci=r'正态总体下，$CI_{0.95}(\sigma^2)=[(n-1)s^2/\chi^2_{0.975,n-1},(n-1)s^2/\chi^2_{0.025,n-1}]$，'+f'结果为 $({fmt(chi/stats.chi2.ppf(.975,df))}, {fmt(chi/stats.chi2.ppf(.025,df))})$。'
        test=r'$H_0:\sigma^2=1$ 对 $H_1:\sigma^2\ne1$；$\chi_0^2=(n-1)s^2/\sigma_0^2$，'+f'结果为 {fmt(chi)}，自由度 {df}，双侧 $P={fmt(p)}$。'
        return point+'\n\n'+(ci if mode=='ci' else test if mode=='test' else ci+'\n\n'+test)
    if op=='two': return two(key,**json.loads('|'.join(args)) if args else {})
    raise ValueError(match.group(0))

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--reuse-figures',action='store_true',help='Reuse figures when changing only page structure or wording')
    args=parser.parse_args()
    if args.reuse_figures:
        original_plot=plot
        def plot(key,groups,tag='raw'):
            fn=f'ch{key.replace(".","-")}-{tag}.png'
            assert (IMAGES/fn).is_file(),fn
            return f'![习题 {key}：观测值、残差 Q–Q 图、残差记录顺序图、残差 vs 拟合值图（{tag}）](../images/{fn})'
    (ROOT/'Exercises'/'data').mkdir(exist_ok=True)
    (ROOT/'Exercises'/'data'/'chapter01-03.json').write_text(json.dumps(DATA,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
    manifest=[]
    for p in sorted((ROOT/'_tools'/'exercise_templates').glob('chapter*.md')):
        text=re.sub(r'\{\{([^{}]+(?:\{[^{}]*\}[^{}]*)?)\}\}',expand,p.read_text(encoding='utf8'))
        assert '{{' not in text, p
        manifest.extend(write_pages(int(p.stem[-2:]),text))
        print(p.name)
    (ROOT/'Exercises/data/pages.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
    (ROOT/'Exercises/data/diagnostics.json').write_text(json.dumps(DIAGNOSTICS,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
    toc_path=ROOT/'toc.yml'
    toc=toc_path.read_text(encoding='utf8')
    start=toc.index('    - file: Exercises/index.md')
    end=toc.index('    - file: Chapter 01/',start)
    lines=['    - file: Exercises/index.md','      title: 习题与答案','      children:']
    for ch,(title,_) in CHAPTERS.items():
        lines.extend([f'        - file: Exercises/translations_zh/chapter{ch:02}.md',f'          title: 第 {ch} 章　{title}习题','          children:'])
        for item in manifest:
            if item['number'].startswith(f'{ch}.'):
                lines.extend([f'            - file: {item["file"]}',f'              title: "{item["number"]} {item["title"]}"'])
    # Preserve the existing repository's CRLF convention.
    toc_path.write_bytes((toc[:start]+'\n'.join(lines)+'\n'+toc[end:]).replace('\r\n','\n').replace('\n','\r\n').encode('utf8'))
    # Chapter 4 has its own blocked-design computations and publishing step.
    import runpy
    runpy.run_path(str(ROOT/'_tools/build_chapter04_exercises.py'),run_name='__main__')
