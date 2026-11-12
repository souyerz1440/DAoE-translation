"""Recompute and publish ninth-edition Chapter 4 exercises. Run with .venv-book Python."""
from pathlib import Path
import json, re, itertools
import os
os.environ.setdefault('MPLCONFIGDIR',str(Path(__file__).resolve().parents[1]/'_out/exercises/matplotlib'))
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'Exercises/translations_zh/chapter04'
OUT.mkdir(exist_ok=True)
IMG=ROOT/'Exercises/images'
FONT=FontProperties(fname=str(ROOT/'LXGWWenKai-Regular.ttf'))
plt.rcParams.update({'font.family':FONT.get_name(),'axes.unicode_minus':False})
from matplotlib import font_manager
font_manager.fontManager.addfont(str(ROOT/'LXGWWenKai-Regular.ttf'))
def f(x): return f'{float(x):.6g}'
def tab(head,rows): return '\n'.join(['| '+' | '.join(head)+' |','| '+' | '.join(['---']*len(head))+' |']+['| '+' | '.join(map(str,r))+' |' for r in rows])

DATA={
8:[[73,68,74,71,67],[73,67,75,72,70],[75,68,78,73,68],[73,71,75,75,69]],
9:[[13,22,18,39],[16,24,17,44],[5,4,1,22]],
12:[[9.3,9.4,9.6,10],[9.4,9.3,9.8,9.9],[9.2,9.4,9.5,9.7],[9.7,9.6,10,10.2]],
13:[[250,350,219,375],[400,525,390,580],[275,340,200,310]],
14:[[.500,.634,.487,.329,.512],[.535,.675,.520,.435,.540],[.513,.595,.488,.400,.510]],
15:[[.78,.80,.81,.75,.77,.78],[.85,.85,.92,.86,.81,.83],[.93,.92,.95,.89,.89,.83],[1.14,.97,.98,.88,.86,.83],[.97,.86,.78,.76,.76,.75]],
16:[[1244,21,82,2221,905,839],[281,129,396,1306,336,910],[220,84,458,543,300,794],[225,83,425,552,291,826],[19,11,-34,121,15,103],[-20,35,-53,170,104,199]],
17:[[334.5,31.6,701,41.2,61.2,69.6,67.5,66.6,120.7,881.9],[919.4,404.2,1024.8,54.1,62.8,671.6,882.1,354.2,321.9,91.1],[108.4,26.1,240.8,191.1,69.7,242.8,62.7,396.9,23.6,290.4]],
18:[[4.93,4.86,4.75,4.95,4.79,4.88],[4.85,4.91,4.79,4.85,4.75,4.85],[4.83,4.88,4.90,4.75,4.82,4.90],[4.89,4.77,4.94,4.86,4.79,4.76]],
19:[[8,4,5,6],[14,5,6,9],[14,6,9,2],[17,9,3,6]],
26:[[10,6,6,6,6],[7,6,6,1,6],[5,3,3,2,5],[6,4,4,2,3]],
45:[[None,17,14,13,12],[14,14,None,13,10],[12,None,13,12,9],[13,11,11,12,None],[11,12,10,None,8]],
47:[[114,None,None,None,120,None,117],[126,120,None,None,None,119,None],[None,137,117,None,None,None,134],[141,None,129,149,None,None,None],[None,145,None,150,143,None,None],[None,None,120,None,118,123,None],[None,None,None,136,None,130,127]],
48:[[73,74,None,71],[None,75,67,72],[73,75,68,None],[75,None,72,75]],
}
NOISE=np.array([[.05,.04,.05,.06,.03,.05],[.04,.02,.03,.05,.03,.02],[.09,.13,.11,.15,.08,.12],[.03,.04,.05,.05,.03,.02]])
LATIN={27:('ABDCE CEADB BACED DCEBA EDBAC'.split(),[[8,7,1,7,3],[11,2,7,3,8],[4,9,10,1,5],[6,8,6,6,10],[4,2,3,8,8]]),28:('CDAB BCDA ABCD DABC'.split(),[[10,14,7,8],[7,18,11,8],[5,10,11,9],[10,10,12,14]]),40:('ABCDE BCDEA CDEAB DEABC EABCD'.split(),[[26,16,19,16,13],[18,21,18,11,21],[20,12,16,25,13],[15,15,22,14,17],[10,24,17,17,14]]),41:('CBDA BCAD ADBC DACB'.split(),[[11,10,14,8],[8,12,10,12],[9,11,7,15],[9,8,18,6]])}
GREEK={40:np.array([[0,1,2,3,4],[2,3,4,0,1],[4,0,1,2,3],[1,2,3,4,0],[3,4,0,1,2]]),41:np.array([[1,2,3,0],[0,3,2,1],[3,0,1,2],[2,1,0,3]])}
LABELS={8:['化学剂1','化学剂2','化学剂3','化学剂4'],9:['洗涤液1','洗涤液2','洗涤液3'],12:['压头1','压头2','压头3','压头4'],13:['设计1','设计2','设计3'],14:['油1','油2','油3'],15:['喷嘴1','喷嘴2','喷嘴3','喷嘴4','喷嘴5'],16:['SLIM','COCOMO-A','COCOMO-R','COCOMO-C','FUNCTION POINTS','ESTIMALS'],17:['MP','MP＋HDMTX','MP＋LDMTX'],19:['5 rpm','10 rpm','15 rpm','20 rpm'],26:['4 ft','6 ft','8 ft','10 ft']}
def raw(k):
    if k in LATIN:
        letters,y=LATIN[k]; greek=GREEK.get(k)
        return tab(['行／列']+[str(j+1) for j in range(len(y))],[[i+1]+[f'{letters[i][j]}'+('αβγδε'[greek[i,j]] if greek is not None else '')+' = '+f(v) for j,v in enumerate(row)] for i,row in enumerate(y)])
    y=DATA[k]
    return tab(['处理／区组（4.17 为组内受试者序号）']+[str(j+1) for j in range(len(y[0]))],[[LABELS.get(k,[str(i+1) for i in range(len(y))])[i]]+['—' if v is None else f(v) for v in row] for i,row in enumerate(y)])

def records(k,missing=(),dropcol=None,greek=True):
    if k in LATIN:
        ls,y=LATIN[k]; p=len(y)
        obs=[(float(y[i][j]),ord(ls[i][j])-65,i,j) for i in range(p) for j in range(p) if (i,j) not in missing and j!=dropcol]
        factors=[np.array([o[t] for o in obs]) for t in [1,2,3]]
        names=['处理','行','列']
        if k in GREEK and greek:
            factors.append(np.array([GREEK[k][o[2],o[3]] for o in obs]));names.append('希腊字母因子')
        return np.array([o[0] for o in obs]),factors,names
    a=np.array(DATA[k],dtype=float)
    obs=[(a[i,j],i,j) for i in range(a.shape[0]) for j in range(a.shape[1]) if np.isfinite(a[i,j]) and (i,j) not in missing]
    return np.array([o[0] for o in obs]),[np.array([o[t] for o in obs]) for t in [1,2]],['处理','区组']

def matrix(factors,n):
    return np.column_stack([np.ones(n)]+[(v==level).astype(float) for v in factors for level in np.unique(v)[1:]])
def fit(y,factors):
    X=matrix(factors,len(y)); coef=np.linalg.lstsq(X,y,rcond=None)[0];yh=X@coef;e=y-yh
    return float(e@e),len(y)-np.linalg.matrix_rank(X),yh,e,X
RESULT={}
def compute(k,missing=(),dropcol=None,greek=True,oneway=False,transform=None):
    y,fac,names=records(k,missing,dropcol,greek)
    if oneway:fac=fac[:1];names=names[:1]
    if transform:y=transform(y)
    sse,df,yh,e,X=fit(y,fac);mse=sse/df
    rows=[]
    # Sequential sums add to total; adjusted treatment SS tests the requested effect.
    prev=float(np.sum((y-y.mean())**2))
    order=list(range(1,len(fac)))+[0]
    included=[]
    for z in order:
        included.append(fac[z]);cur,d,_,_,_=fit(y,included)
        dfn=len(np.unique(fac[z]))-1;ss=prev-cur;prev=cur
        F=(ss/dfn)/mse
        nonorthogonal=bool(missing) or dropcol is not None or k in [45,47,48]
        rows.append([names[z]+('（调整后）' if z==0 and len(fac)>1 else ''),dfn,ss,ss/dfn,None if nonorthogonal and z!=0 else F,None if nonorthogonal and z!=0 else stats.f.sf(F,dfn,df)])
    rows += [['误差',df,sse,mse,None,None],['总计',len(y)-1,float(np.sum((y-y.mean())**2)),None,None,None]]
    return dict(y=y,fac=fac,names=names,sse=sse,df=df,mse=mse,yh=yh,e=e,X=X,rows=rows)

RCBD=r'''**模型与假设**：$y_{ij}=\mu+\tau_i+\beta_j+\epsilon_{ij}$，$\sum_i\tau_i=\sum_j\beta_j=0$；误差独立且服从 $N(0,\sigma^2)$，处理与区组效应可加，没有处理×区组交互作用。处理在各区组内随机分配。

**假设与检验**：$H_0:\tau_1=\cdots=\tau_a=0$；$H_1$：至少一个处理效应非零。令 $C=G^2/(ab)$，$T_i$、$B_j$ 分别为处理及区组总和，

$$SS_T=\sum_{ij}y_{ij}^2-C,\quad SS_A=\frac{\sum_iT_i^2}{b}-C,\quad SS_B=\frac{\sum_jB_j^2}{a}-C,$$
$$SS_E=SS_T-SS_A-SS_B,\quad F_0=\frac{SS_A/(a-1)}{SS_E/[(a-1)(b-1)]}\sim F_{a-1,(a-1)(b-1)}\quad(H_0).$$

除题目另有说明，取 $\alpha=0.05$，以 $P<0.05$ 拒绝原假设。本文 $t_q,F_q$ 均指下侧累积概率为 $q$ 的分位数。区组作为控制因素，不将其显著性作为保留区组的前提。
'''
LS=r'''**模型与假设**：$y_{ij}=\mu+\tau_{L(i,j)}+\rho_i+\kappa_j+\epsilon_{ij}$，各组效应和为零；各误差独立同分布为 $N(0,\sigma^2)$，所有因子效应可加，没有交互作用。$L(i,j)$ 是该单元的处理标签。行、列及处理均为固定效应。

**原假设与备择假设**：$H_0:\tau_1=\cdots=\tau_p=0$；$H_1$：至少一个处理效应不为零。设 $C=G^2/p^2$，各因子水平总和为 $U_l$，

$$SS_T=\sum y^2-C,\qquad SS_{\text{因子}}=\frac1p\sum_lU_l^2-C,$$
$$SS_E=SS_T-SS_A-SS_R-SS_C,\qquad F_0=\frac{SS_A/(p-1)}{SS_E/[(p-1)(p-2)]}.$$

显著性水平取 $0.05$；处理、行、列自由度各为 $p-1$，总自由度 $p^2-1$。
'''
GENERAL=r'''**缺失／非正交资料的精确检验**：在保留下来的真实观测上拟合上述加性模型。令全模型和删去待检验因子的简化模型设计矩阵分别为 $X_F,X_R$，$H=X X^+$ 为正交投影，

$$SS_E=y'(I-H_F)y,\quad SS_{A\mid\text{其他}}=y'(H_F-H_R)y,$$
$$F_0=\frac{[SSE_R-SSE_F]/[\operatorname{rank}X_F-\operatorname{rank}X_R]}{SSE_F/[N-\operatorname{rank}X_F]}.$$

参数用最小二乘估计，检验 $H_0$：所有可估计的处理差为零；备择为至少一个不为零。误差仍须独立、正态、同方差。表中先列区组因子的顺序平方和，最后列调整后的处理平方和，使各项相加等于总平方和；它们不能全部当作偏平方和。
'''
def anova(r):
    return tab(['来源','自由度','平方和','均方','$F$','$P$ 值'],[[row[0],int(row[1])]+['—' if v is None else f(v) for v in row[2:]] for row in r['rows']])
def report(k,**kwargs):
    r=compute(k,**kwargs); RESULT[str(k)+str(kwargs)]=dict(sse=r['sse'],df=int(r['df']),mse=r['mse'])
    prefix=LS if k in LATIN else RCBD
    if k in GREEK and kwargs.get('greek',True):
        prefix=LS+'\n再加入希腊字母因子 $\\gamma_{G(i,j)}$（效应和为零），从总平方和再减去 $SS_G$；误差自由度改为 $(p-1)(p-3)$，其余假设与检验方式相同。\n'
    if kwargs.get('oneway'):
        prefix=r'''**模型与假设**：$y_{ij}=\mu+\tau_i+\epsilon_{ij}$，误差独立且服从 $N(0,\sigma^2)$。$H_0$：各处理均值相等；备择为至少一个不同。

$$SS_A=\sum_i n_i(\bar y_i-\bar y)^2,\quad SS_E=\sum_{ij}(y_{ij}-\bar y_i)^2,\quad SS_T=SS_A+SS_E,\quad F_0=\frac{SS_A/(a-1)}{SS_E/(N-a)}.$$
'''
    if kwargs.get('missing') or kwargs.get('dropcol') is not None or k in [45,47,48]:prefix+='\n'+GENERAL
    if k in [45,47,48]:
        prefix=r'''**模型与假设**：对实际安排的处理—区组单元拟合 $y_{ij}=\mu+\tau_i+\beta_j+\epsilon_{ij}$，效应和为零，误差独立同分布为 $N(0,\sigma^2)$，处理与区组效应可加。$H_0$：各处理总体均值相等；备择为至少一个不同。取 $\alpha=0.05$。

设总运行数为 $N_0=ar=bk$，处理数为 $a$，区组数为 $b$，每个区组运行 $k$ 个处理，各处理重复 $r$ 次，每对不同处理同区组出现 $\lambda$ 次。令 $C_0=G^2/N_0$，

$$Q_i=T_i-\frac1k\sum_jn_{ij}B_j,\quad SS_T=\sum y_{ij}^2-C_0,\quad SS_B=\frac1k\sum_jB_j^2-C_0,$$
$$SS_{A\mid B}=\frac{k\sum_iQ_i^2}{\lambda a},\quad SS_E=SS_T-SS_B-SS_{A\mid B},\quad F_0=\frac{SS_{A\mid B}/(a-1)}{SS_E/(N_0-a-b+1)}.$$

表中区组平方和未调整处理，用于加法分解，不据此检验区组效应。处理平方和已调整区组。
'''+GENERAL
    treatment=next(row for row in r['rows'] if row[0].startswith('处理'))
    conclusion=f"$F_0={f(treatment[4])}$，$P={f(treatment[5])}$，"+('拒绝处理均值相等的原假设，有处理差异的证据。' if treatment[5]<.05 else '不拒绝处理均值相等的原假设；尚无充分证据认定处理均值不同。')
    return prefix+'\n\n'+anova(r)+'\n\n'+conclusion,r

def diagnostics(k,r,suffix=''):
    name=f'ch4-{k}{suffix}-diagnostics.png'
    fig,ax=plt.subplots(1,3,figsize=(14,4));stats.probplot(r['e'],plot=ax[0]);ax[0].set(title='残差正态 Q–Q 图',xlabel='理论正态分位数',ylabel='残差')
    ax[1].plot(np.arange(1,len(r['e'])+1),r['e'],'o-',lw=.7);ax[1].set(title='残差记录顺序图',xlabel='原表逐行记录序号',ylabel='残差')
    ax[2].scatter(r['yh'],r['e']);ax[2].set(title='残差 vs 拟合值',xlabel='拟合值',ylabel='残差')
    for a in ax[1:]:a.axhline(0,color='gray',lw=.7)
    fig.tight_layout();fig.savefig(IMG/name,dpi=160);plt.close(fig)
    sw=stats.shapiro(r['e'])
    assessment=('正态性检验提示残差偏离正态，应结合 Q–Q 图中的尾部和异常值谨慎解释 F 检验。' if sw.pvalue<.05 else '未检出明显的正态性偏离；小样本下检验功效有限，仍应结合 Q–Q 图判断。')
    return f'![习题 4.{k}：残差 Q–Q 图、记录顺序图与残差 vs 拟合值图](../../images/{name})\n\n'+r'''残差为 $e=y-X\hat\theta$。左图检查正态性及极端残差，右图检查方差和均值结构。中图按表中逐行记录顺序绘制；题目未给出实际随机运行顺序，**它不是真实时序图，不能据此确认时间独立性**。需获得运行日志后再按实际时间重画。未见明显图形异常也不能证明假设成立。'''+f'\n\n残差 Shapiro–Wilk 统计量 $W={f(sw.statistic)}$、$P={f(sw.pvalue)}$（作为辅助描述，拟合残差相互相关，不能当成独立原始样本的精确检验）。'+assessment
def means(k,r):
    levels=np.unique(r['fac'][0])
    labels=LABELS.get(k,[str(j+1) for j in range(len(levels))])
    return tab(['处理','样本均值'],[[labels[i],f(r['y'][r['fac'][0]==i].mean())] for i in levels])
def lsd(k,r):
    a=len(np.unique(r['fac'][0]));b=len(r['y'])/a;L=stats.t.ppf(.975,r['df'])*np.sqrt(2*r['mse']/b)
    mu=np.array([r['y'][r['fac'][0]==i].mean() for i in range(a)])
    return r'估计量为 $\hat\mu_i=\bar y_{i\cdot}$，$\widehat{\mu_i-\mu_j}=\bar y_i-\bar y_j$；'+f'$LSD=t_{{0.975,{int(r["df"])}}}\\sqrt{{2MS_E/b}}={f(L)}$。差值的逐对 95% 区间为估计差 $\\pm LSD$。\n\n'+means(k,r)+'\n\n'+tab(['处理对','均值差','95% 区间','LSD 判定'],[[f'{i+1} − {j+1}',f(mu[i]-mu[j]),f'[{f(mu[i]-mu[j]-L)}, {f(mu[i]-mu[j]+L)}]','显著' if abs(mu[i]-mu[j])>L else '不显著'] for i in range(a) for j in range(i+1,a)])+'\n\nFisher LSD 为逐对比较，不能将其解释为所有区间联合具有 95% 覆盖率。整体检验不显著时，不以保护性 LSD 宣布处理差异。'
def scaled(k,r):
    mu=np.array([r['y'][r['fac'][0]==i].mean() for i in np.unique(r['fac'][0])]);b=len(r['y'])/len(mu);s=np.sqrt(r['mse']/b)
    fig,ax=plt.subplots(figsize=(8,4));grid=np.linspace(stats.t.ppf(.001,r['df']),stats.t.ppf(.999,r['df']),400)
    for i,m in enumerate(mu):ax.plot(m+s*grid,stats.t.pdf(grid,r['df'])/s,label=f'处理 {i+1}')
    ax.scatter(mu,np.zeros_like(mu),color='black');ax.set(xlabel='处理均值尺度',ylabel='缩放 t 密度',title='各均值处的缩放 t 分布');ax.legend();fig.tight_layout();fig.savefig(IMG/f'ch4-{k}-scaled-t.png',dpi=160);plt.close(fig)
    return f'![处理均值与缩放 t 分布](../../images/ch4-{k}-scaled-t.png)\n\n'+f'绘制 $m_i+s_{{\\bar y}} t_{{{int(r["df"])}}}$ 的密度曲线，其中 $s_{{\\bar y}}=\\sqrt{{MS_E/b}}={f(s)}$，$m_i=\\bar y_i$。曲线重叠仅作直观参考，正式逐对判断使用 $SE(\\bar y_i-\\bar y_j)=\\sqrt{{2MS_E/b}}$，不能简单按两条均值区间是否重叠作显著性判断。\n\n'+lsd(k,r)

P={}
def page(n,title,q,a,source=None):P[n]=(title,q,a,source)
def part(letter,title,text):return f'### （{letter}）{title}\n\n{text}'

page(1,'忽略区组后的误差方差','单因子有 4 个水平、6 次重复，采用区组设计。误差平方和为 500，区组平方和为 250。若相同数据按完全随机设计分析，误差方差估计为多少？选项：（a）25.0；（b）25.5；（c）35.0；（d）37.5；（e）以上皆非。',r'完全随机模型为 $y_{ij}=\mu+\tau_i+\epsilon_{ij}$，假设误差独立、正态且同方差。忽略区组后，$SS_E^*=500+250=750$，$df_E^*=4(6-1)=20$，故 $\hat\sigma^2=MS_E^*=750/20=37.5$，选 **（d）**。这里是对相同数据重新分析；不意味着实际取消区组后必然得到相同数据。')
page(2,'区组化损失的误差自由度','单因子有 5 个水平、3 次重复，原采用完全随机设计。改为随机完全区组设计后，纯误差自由度减少多少？选项：（a）3；（b）5；（c）2；（d）4；（e）以上皆非。',r'完全随机设计误差自由度为 $5(3-1)=10$；随机完全区组设计为 $(5-1)(3-1)=8$，减少 $3-1=2$，选 **（c）**。两个设计都采用可加、独立同方差误差模型；单次观测的区组设计无法单独估计交互作用。')
page(3,'区组控制干扰变异','判断正误：区组化可以控制试验中不可控干扰因子所传递的变异。','**正确。** 只要能够识别干扰来源并据此形成内部相对均匀的区组，就可从误差中分离相应变异；干扰因子未必能调到任意指定值。每个区组内仍应随机化。')
page(4,'区组数与处理数','判断正误：随机完全区组设计中，区组数必须始终等于处理数或因子水平数。','**错误。** 有 $a$ 个处理、$b$ 个区组时，每区组含全部 $a$ 个处理，共 $ab$ 次运行；$a$ 与 $b$ 无须相等。例如 4.8 有 4 个处理、5 个区组。')
page(5,'区组化与随机化原则','“能区组化就区组化，不能区组化就随机化”的关键含义是什么？（a）区组内通常最好不随机化；（b）区组化违反等方差假设；（c）将干扰因子各水平作为区组，并在区组内随机化；（d）随机化运行优于随机化区组。','选 **（c）**。能够识别并控制的干扰来源进入区组结构；其余影响通过随机化分散。区组化与区组内随机化配合使用。')

q6='随机完全区组试验的 ANOVA 表如下。\n\n'+tab(['来源','自由度','平方和','均方','$F$'],[['处理',4,1010.56,'?',29.84],['区组','?','?',64.765,'?'],['误差',20,169.33,'?','—'],['总计',29,1503.71,'—','—']])+'\n\n**（a）** 填写空白，可给出 P 值界限。\n\n**（b）** 使用了多少个区组？\n\n**（c）** 可以得出什么结论？'
ssb=1503.71-1010.56-169.33;mse=169.33/20
r6={'rows':[['处理',4,1010.56,1010.56/4,(1010.56/4)/mse,stats.f.sf((1010.56/4)/mse,4,20)],['区组',5,ssb,ssb/5,(ssb/5)/mse,stats.f.sf((ssb/5)/mse,5,20)],['误差',20,169.33,mse,None,None],['总计',29,1503.71,None,None,None]]}
page(6,'补全区组试验方差分析表',q6,part('a','补全表格',RCBD+'\n'+anova(r6))+ '\n\n'+part('b','区组数',r'$N=30$，$a=5$，所以 $b=30/5=6$，区组自由度为 5。')+'\n\n'+part('c','结论','处理检验显著，至少两个处理总体均值不同；仅凭汇总表不能确认具体处理对，也不能进行残差诊断。'), '无直接对应题；按第 9 版汇总表重新计算。')
page(7,'将完全随机设计改按区组分析','对习题 [3.8](../chapter03/exercise-3-8.md) 的单因子完全随机试验，假定实际采用随机完全区组设计，区组平方和为 80.00。修改 ANOVA 表。\n\n原分析：处理自由度 4、平方和 987.71；误差自由度 25、平方和 186.53；总自由度 29、平方和 1174.24。',RCBD+'\n\n'+anova({'rows':[['处理',4,987.71,987.71/4,(987.71/4)/(106.53/20),stats.f.sf((987.71/4)/(106.53/20),4,20)],['区组',5,80,16,None,None],['误差',20,106.53,106.53/20,None,None],['总计',29,1174.24,None,None,None]]})+'\n\n区组数为 6，误差平方和减少 80、自由度减少 5。处理效应显著。原数据未指定各观测的区组归属，不能据此恢复残差图。')

BACKGROUND={8:'化学家研究 4 种化学剂对某种布料抗拉强度的影响。为控制布卷差异，选取 5 卷布作为区组，每卷均以随机顺序使用全部 4 种化学剂。',9:'比较 3 种洗涤液抑制 5 加仑牛奶容器中细菌生长的效果。每天只能进行 3 次试验，以 4 个试验日为区组，每天测试全部洗涤液。',12:'硬度试验使用 4 种压头和 4 块试片，以试片为区组。数据为 Rockwell C 标度读数减去 40。',13:'公司比较 3 种直邮宣传册。对东北、西北、东南、西南 4 个地区，每种设计各寄出 5000 份，以地区为区组，记录回应数。',14:'比较 3 种润滑油对柴油卡车燃油经济性的影响，响应为发动机运行 15 分钟后的制动比油耗，以 5 台发动机为区组。',15:'《Fire Safety Journal》（第 4 卷，1981 年 8 月）研究喷嘴设计对湍流水射流稳定性和性能的影响。在出口速度 11.73、14.37、16.59、20.43、23.46、28.74 m/s 下测得形状因子，以速度为干扰因子（区组），比较 5 种喷嘴。',16:'《Communications of the ACM》（30 卷 5 期，1987 年）比较 6 种软件开发成本估计算法。对 6 个项目应用全部算法，响应为成本估计的百分误差，以项目为区组。教材第 4 行印作 COCONO-C，这里按通常名称标为 COCOMO-C。',17:'《Nature Genetics》（34 卷，2003 年，85～90 页）研究白血病治疗与基因表达。3 个组各含 10 名受试者：仅 MP、MP＋高剂量甲氨蝶呤（HDMTX）、MP＋低剂量甲氨蝶呤（LDMTX）。表中列出某基因的响应。',18:'第 3.8 节的比率控制算法试验实际按随机区组设计进行，以 6 个时间段为区组，每个时间段测试全部 4 种算法。记录平均槽电压和电压标准差（pot noise）。',19:'铝母合金生产商使用 4 台具有不同运行特性的熔炉，以熔炉为区组，研究 5、10、15、20 rpm 搅拌速度对晶粒尺寸的影响。',26:'工业工程师研究物体距眼睛的距离对聚焦时间的影响。比较 4、6、8、10 ft，以 5 名受试者为区组，每人接受全部距离条件。'}
TITLES={8:'化学剂与布料抗拉强度',9:'洗涤液与细菌生长',12:'压头与硬度读数',13:'直邮宣传册回应数',14:'润滑油与比油耗',15:'喷嘴设计与形状因子',16:'软件成本估计算法',17:'白血病治疗与基因表达',18:'比率控制算法与槽电压',19:'搅拌速度与晶粒尺寸',26:'距离与眼睛聚焦时间'}
for n in [8,9,26,27,28,40,41,45,47]:
    if n in LATIN:
        bg={27:'研究 5 种配料 A～E 对化学反应时间的影响。每批原料仅够运行 5 次，每次约需 1.5 小时，每天最多运行 5 次，以原料批次为行、日期为列，采用拉丁方控制这两个干扰来源。',28:'工程师比较 4 种电视部件装配方法 A～D。选 4 名操作员，并用装配次序控制疲劳影响，以次序为行、操作员为列安排拉丁方，响应为装配时间。',40:'化学过程产率试验包含 5 批原料（行）、5 种酸浓度（列）、5 种静置时间 A～E 和 5 种催化剂浓度 α～ε，采用下列希腊—拉丁方。',41:'在习题 4.28 的研究中，工程师进一步以工作场所 α～δ 控制潜在变异。重新试验得到下列希腊—拉丁方：行为装配次序，列为操作员，A～D 为方法，α～δ 为工作场所。'}[n]
        title={27:'配料反应时间的拉丁方',28:'装配方法的拉丁方',40:'化学产率的希腊—拉丁方',41:'增加工作场所的装配试验'}[n]
    elif n in [45,47]:
        bg={45:'工程师比较 5 种汽油添加剂对行驶里程的影响，以汽车为区组，因时间限制每车只能测试 4 种添加剂，得到下列平衡不完全区组设计。破折号表示未安排该组合。',47:'研究 7 种硬木浓度（依次为 2%、4%、6%、8%、10%、12%、14%）对纸张强度的影响。试验厂每天只能运行 3 次，以日期为区组，采用下列平衡不完全区组设计。破折号表示未安排该组合。'}[n]
        title={45:'汽油添加剂的不完全区组试验',47:'硬木浓度与纸张强度'}[n]
    else:bg=BACKGROUND[n];title=TITLES[n]
    ans,r=report(n)
    if n in [45,47]:
        ans+='\n\n'+r'平衡不完全区组设计中 $Q_i=T_i-k^{-1}\sum_j n_{ij}B_j$，$\hat\tau_i=kQ_i/(\lambda a)$，$SS_{A\mid B}=k\sum_iQ_i^2/(\lambda a)$，$SS_B=\sum_j B_j^2/k-G^2/N$。误差自由度为 $N-a-b+1$。处理平方和必须调整区组，不能使用完全随机设计的原始处理总和公式。'
    ans+='\n\n'+diagnostics(n,r)
    page(n,title,bg+'\n\n'+raw(n)+'\n\n以 $\\alpha=0.05$ 分析数据并给出结论。',ans)

for n,parent in [(10,8),(11,9)]:
    r=compute(parent)
    page(n,'处理均值与缩放 t 分布',f'绘制习题 [4.{parent}](./exercise-4-{parent}.md) 中各处理的样本均值，与适当缩放的 t 分布比较，并解释结论。\n\n'+raw(parent),scaled(n,r))

for n in [12,13,14,15,16]:
    ans,r=report(n)
    if n in [12,13,14]:
        q='**（a）** 分析试验数据。\n\n**（b）** 用 Fisher LSD 法比较各处理，确定哪些总体均值不同。\n\n**（c）** 分析残差。'
        answer=part('a','方差分析',ans)+'\n\n'+part('b','Fisher LSD 比较',lsd(n,r))+'\n\n'+part('c','残差分析',diagnostics(n,r))
    elif n==15:
        q='**（a）** 用散点图和方差分析判断喷嘴设计是否影响形状因子，取 $\\alpha=0.05$。\n\n**（b）** 分析残差。\n\n**（c）** 绘制均值与缩放 t 分布，判断哪些喷嘴有差异，并与 Duncan 多重极差检验比较。'
        fig,ax=plt.subplots(figsize=(7,4))
        for i,row in enumerate(DATA[n]):ax.scatter([i+1]*len(row),row)
        ax.set(xlabel='喷嘴设计',ylabel='形状因子');fig.tight_layout();fig.savefig(IMG/'ch4-15-scatter.png',dpi=160);plt.close(fig)
        mu=np.array(DATA[15]).mean(axis=1);idx=np.argsort(mu);sr=np.sqrt(r['mse']/6);dr=[]
        for lo in range(4):
            for hi in range(lo+1,5):
                p=hi-lo+1;crit=stats.studentized_range.ppf((1-.05)**(p-1),p,r['df'])*sr;d=mu[idx[hi]]-mu[idx[lo]]
                dr.append([f'{idx[lo]+1} 与 {idx[hi]+1}',p,f(d),f(crit),'显著' if d>crit else '不显著'])
        duncan=r'Duncan 法对排序均值中跨度为 $p$ 的一对，采用 $\alpha_p=1-(1-0.05)^{p-1}$，$R_p=q_{1-\alpha_p;p,\nu}\sqrt{MS_E/b}$，按逐步极差程序比较。'+'\n\n'+tab(['喷嘴对','跨度 p','均值差','临界极差','比较'],dr)+'\n\n本题 Duncan 与保护性 LSD 的逐对判定一致：1–2、1–3、1–4、2–4、3–5、4–5 显著，其余处理对不显著。一般情况下两法可能不同；Duncan 不控制全部逐对比较的家族错误率。图中的 t 曲线只是直观展示，应按各方法自身临界值判断。'
        answer=part('a','散点图与方差分析','![喷嘴观测值散点图](../../images/ch4-15-scatter.png)\n\n'+ans)+'\n\n'+part('b','残差分析',diagnostics(n,r)+'\n\n第 4、5 种喷嘴的响应随速度变化较明显，提示喷嘴×速度交互作用或非恒定变异；单次观测的可加模型不能将交互作用与误差分开，显著性结论须结合这一限制解释。')+'\n\n'+part('c','处理比较',scaled(n,r)+'\n\n'+duncan)
    else:
        q='**（a）** 算法的平均成本估计准确性是否不同？\n\n**（b）** 分析残差。\n\n**（c）** 实际中推荐哪种算法？'
        answer=part('a','方差分析',ans)+'\n\n'+part('b','残差分析',diagnostics(n,r)+'\n\nSLIM 在个别项目上的巨大误差导致重尾及不均匀变异。原始尺度的正态同方差模型值得怀疑；存在负数，不能直接对全部误差取对数。')+'\n\n'+part('c','推荐算法',means(n,r)+'\n\n准确性应结合接近零的误差和误差离散度评估，不能仅追求最小的有符号均值。FUNCTION POINTS 的绝对百分误差均值为 '+f(np.mean(np.abs(DATA[16][4])))+'，ESTIMALS 为 '+f(np.mean(np.abs(DATA[16][5])))+'；本组项目上优先考虑 FUNCTION POINTS，再以新项目验证。')
    page(n,TITLES[n],BACKGROUND[n]+'\n\n'+raw(n)+'\n\n'+q,answer)

# Independent subjects in 4.17: the column index is not a pairing identifier.
ar,rr=report(17,oneway=True);al,rl=report(17,oneway=True,transform=np.log)
page(17,TITLES[17],BACKGROUND[17]+'\n\n'+raw(17)+'\n\n**（a）** 有证据认为治疗均值不同吗？\n\n**（b）** 检查正态假设，能否认为来自正态总体？\n\n**（c）** 对原数据取对数，重新检验治疗均值。\n\n**（d）** 分析变换后残差并评价模型。',part('a','原始尺度分析','原题各组为不同受试者，未说明跨组配对，因此采用独立组单因子模型，不把表中同列序号当作区组。\n\n'+ar)+'\n\n'+part('b','原始尺度正态性',diagnostics(17,rr)+'\n\n原始数据明显右偏，包含较大极值，正态近似值得怀疑；不能由正态性检验未显著推断总体必然正态。')+'\n\n'+part('c','自然对数尺度分析','所有观测为正，令 $z=\\ln y$。\n\n'+al)+'\n\n'+part('d','变换后诊断',diagnostics(17,rl,'-log')+'\n\n对数尺度减轻了右偏与方差不齐。变换后检验针对平均对数响应，而非原始尺度算术均值。'), '对应参考答案 4.15；参考答案将各列当作区组，原题未给出配对依据，这里按独立组重算。')

a18,r18=report(18)
DATA[180]=NOISE.tolist();an,rn=report(180)
page(18,TITLES[18],BACKGROUND[18]+'\n\n**平均电压**\n\n'+raw(18)+'\n\n**电压标准差（与上表逐项对应）**\n\n'+raw(180)+'\n\n**（a）** 平均电压是否受算法影响？取 $\\alpha=0.05$。\n\n**（b）** 分析电压标准差，算法是否影响 pot noise？\n\n**（c）** 进行适当残差分析。\n\n**（d）** 若要同时降低平均电压与噪声，选择哪种算法？',part('a','平均电压',a18)+'\n\n'+part('b','槽噪声',an)+'\n\n'+part('c','两个响应的诊断',diagnostics(18,r18)+'\n\n'+diagnostics(18,rn,'-noise')+'\n\n噪声响应在算法 3 下较大且离散更强，方差结构须谨慎评价，可另考察对数标准差模型。')+'\n\n'+part('d','算法选择',means(18,r18)+'\n\n噪声均值依次为 '+ '、'.join(f(v) for v in NOISE.mean(axis=1))+'。算法 2 的平均电压较低且噪声最低，是兼顾两响应的合理选择；平均电压检验不显著，不能声称它在总体上电压显著最小。'))
a19,r19=report(19)
fig,axes=plt.subplots(1,2,figsize=(10,4))
for ax,fac,label in zip(axes,r19['fac'][::-1],['熔炉','搅拌速度水平']):ax.scatter(fac+1,r19['e']);ax.axhline(0,color='gray');ax.set(xlabel=label,ylabel='残差')
fig.tight_layout();fig.savefig(IMG/'ch4-19-factor-residuals.png',dpi=160);plt.close(fig)
page(19,TITLES[19],BACKGROUND[19]+'\n\n'+raw(19)+'\n\n**（a）** 有证据认为搅拌速度影响晶粒尺寸吗？\n\n**（b）** 绘制残差正态概率图并解释。\n\n**（c）** 绘制残差对熔炉及速度图，解释其信息。\n\n**（d）** 若希望晶粒较小，应推荐什么速度和熔炉？',part('a','方差分析',a19)+'\n\n'+part('b','残差正态图',diagnostics(19,r19))+'\n\n'+part('c','残差对因子','![残差对熔炉和搅拌速度](../../images/ch4-19-factor-residuals.png)\n\n部分熔炉的响应对速度变化不同，提示可加性可能不足。每个组合只有一次观测，无法独立估计交互作用；需要重复试验验证。')+'\n\n'+part('d','工艺建议',means(19,r19)+'\n\n熔炉均值依次为 '+ '、'.join(f(v) for v in np.array(DATA[19]).mean(axis=0))+'。样本上 5 rpm 和熔炉 3 的边际均值较低，可作为后续验证起点；速度主效应未显著，不能据此宣布该速度最优，也不能把单个最小单元等同于最优工艺。'))

a20,r20=report(9)
page(20,'用一般回归显著性检验分析洗涤液','用一般回归显著性检验重新分析习题 [4.9](./exercise-4-9.md)。\n\n'+raw(9),GENERAL+'\n\n全模型含截距、2 个处理指示变量、3 个日期指示变量，秩为 6；简化模型仅含截距及日期，秩为 4。\n\n'+a20+'\n\n'+f'$SSE_R={f(fit(r20["y"],r20["fac"][1:])[0])}$，$SSE_F={f(r20["sse"])}$。分子自由度为 2，分母自由度为 6，与随机完全区组 ANOVA 一致。')
y8=np.array(DATA[8]);mu=y8.mean()
page(21,'固定效应参数估计','在习题 [4.8](./exercise-4-8.md) 中，假定化学剂和布卷均为固定效应，估计 $\\tau_i$ 与 $\\beta_j$。\n\n'+raw(8),r'采用 $y_{ij}=\mu+\tau_i+\beta_j+\epsilon_{ij}$，约束 $\sum_i\tau_i=\sum_j\beta_j=0$。最小二乘估计为'+ '\n\n'+r'$$\hat\mu=\bar y_{\cdot\cdot},\quad \hat\tau_i=\bar y_{i\cdot}-\bar y_{\cdot\cdot},\quad \hat\beta_j=\bar y_{\cdot j}-\bar y_{\cdot\cdot}.$$'+'\n\n'+f'$\\hat\\mu={f(mu)}$；$\\hat\\tau_i$ 依次为 '+ '、'.join(f(v) for v in y8.mean(axis=1)-mu)+'；$\\hat\\beta_j$ 依次为 '+ '、'.join(f(v) for v in y8.mean(axis=0)-mu)+'。两组估计的和均为零。')
r9=compute(9);dgrid=np.linspace(0,3,151);beta=stats.ncf.cdf(stats.f.ppf(.95,2,6),2,6,4*dgrid**2)
fig,ax=plt.subplots(figsize=(7,4));ax.plot(dgrid,beta);ax.set(xlabel='标准化效应 Δ = √(Στ²) / σ',ylabel='第二类错误概率 β',title='4.9 设计的 OC 曲线（α = 0.05）');fig.tight_layout();fig.savefig(IMG/'ch4-22-oc.png',dpi=160);plt.close(fig)
page(22,'洗涤液试验的操作特性曲线','为习题 [4.9](./exercise-4-9.md) 的设计绘制 OC 曲线。该检验对处理效应的小差异敏感吗？',r'固定效应 RCBD 中 $a=3,b=4$，非中心参数 $\Lambda=b\sum_i\tau_i^2/\sigma^2$，自由度为 $(2,6)$。定义 $\Delta=\sqrt{\sum_i\tau_i^2}/\sigma$，则 $\Lambda=4\Delta^2$，'+ '\n\n'+r'$$\beta(\Delta)=P\{F_{2,6}(4\Delta^2)\le F_{0.95;2,6}\},\qquad \mathrm{Power}=1-\beta.$$'+'\n\n![OC 曲线](../../images/ch4-22-oc.png)\n\n'+tab(['Δ','功效'],[[f(d),f(stats.ncf.sf(stats.f.ppf(.95,2,6),2,6,4*d*d))] for d in [.5,1,1.5,2]])+'\n\n当效应相对误差较小时，功效较低；只有 4 个区组，不能认为它能灵敏检出任意小差异。')

def missing_analysis(k,missing):
    y,fac,names=records(k,missing);sse,df,yh,e,X=fit(y,fac)
    full_y,full_fac,_=records(k)
    # Build full-level reference coding for prediction (all levels survive these deletions).
    coef=np.linalg.lstsq(X,y,rcond=None)[0];pred=matrix(full_fac,len(full_y))@coef
    shape=np.array(LATIN[k][1] if k in LATIN else DATA[k]).shape
    vals=[pred[i*shape[1]+j] for i,j in missing]
    a,r=report(k,missing=missing)
    filled=np.array(LATIN[k][1] if k in LATIN else DATA[k],dtype=float)
    for (i,j),v in zip(missing,vals):filled[i,j]=v
    if k in LATIN:
        ls,_=LATIN[k];p=len(ls);tot=filled.sum();ssA=sum(filled[np.array([[c==chr(65+t) for c in row] for row in ls])].sum()**2/p for t in range(p))-tot**2/p**2
    else:ssA=shape[1]*sum((filled.mean(axis=1)-filled.mean())**2)
    apF=ssA/(shape[0]-1)/r['mse']
    extra=f'缺失值的最小二乘预测为 '+ '、'.join(f(v) for v in vals)+f'。补值后误差自由度必须减少 {len(missing)}。补值传统法给出处理平方和 {f(ssA)}、$F={f(apF)}$；精确法则在真实观测上使用调整平方和，见下表。两法的误差平方和相同，但处理平方和一般不同，正式结论以精确检验为准。'
    return extra+'\n\n'+a,vals,r
for n,k,miss in [(23,8,[(1,2)]),(24,12,[(1,2)]),(34,27,[(2,3)]),(44,12,[(1,2)])]:
    text,vals,r=missing_analysis(k,miss)
    q={23:'假设化学剂 2、布卷 3 的观测缺失。先估计缺失值并作分析，再作精确分析并比较。',24:'假设压头 2、试片 3 的观测缺失。估计缺失值并分析。',34:'假设批次 3、日期 4 的观测缺失。估计缺失值并用其分析。',44:'对习题 4.24 的缺失值问题，使用第 4.1.4 节的一般回归方法作精确分析，并与近似分析比较。'}[n]
    formula=r'RCBD 的单个缺失值估计为 $\hat x=(aT_i^\prime+bB_j^\prime-G^\prime)/[(a-1)(b-1)]$，各带撇总和均不含缺失值。' if k not in LATIN else r'拉丁方单个缺失值估计为 $\hat x=[p(R^\prime+C^\prime+T^\prime)-2G^\prime]/[(p-1)(p-2)]$，各带撇总和均不含缺失值。'
    page(n,'缺失观测的估计与分析',f'沿用习题 [4.{k}](./exercise-4-{k}.md) 的数据。'+q+'\n\n下表为原完整数据；分析时删除指定观测。\n\n'+raw(k),formula+'\n\n'+text)

text,vals,r=missing_analysis(8,[(1,2),(3,3)])
iterfilled=np.array(DATA[8],float);iterfilled[1,2]=iterfilled[3,3]=np.mean([v for i,row in enumerate(DATA[8]) for j,v in enumerate(row) if (i,j) not in [(1,2),(3,3)]])
for iteration in range(1000):
    previous=iterfilled.copy()
    for i,j in [(1,2),(3,3)]:
        iterfilled[i,j]=0
        iterfilled[i,j]=(4*iterfilled[i,:].sum()+5*iterfilled[:,j].sum()-iterfilled.sum())/12
    if np.max(np.abs(iterfilled-previous))<1e-10:break
assert np.allclose([iterfilled[1,2],iterfilled[3,3]],vals)
page(25,'随机区组中两个缺失值','习题 [4.8](./exercise-4-8.md) 中，化学剂 2／布卷 3、化学剂 4／布卷 4 两个观测缺失。\n\n'+raw(8)+'\n\n**（a）** 迭代估计两缺失值并分析。\n\n**（b）** 对 SSE 分别求导，联立求解并分析。\n\n**（c）** 推导两缺失值位于不同区组时的一般公式。\n\n**（d）** 推导两缺失值位于同一区组时的一般公式。',part('a','迭代补值',r'用剩余观测的均值初始化两缺失值，交替采用单缺失值公式更新，直到两次变化均小于 $10^{-10}$。收敛值为 '+ '、'.join(f(v) for v in vals)+'。\n\n'+text)+'\n\n'+part('b','联立方程',r'令 $D=(a-1)(b-1)$，$h_{uv}=aI(i_u=i_v)+bI(j_u=j_v)-1$，$A_u=aT_{i_u}^\prime+bB_{j_u}^\prime-G^\prime$。对补值后的 $SSE=\sum y^2-\sum_iT_i^2/b-\sum_jB_j^2/a+G^2/(ab)$ 求导，得到'+ '\n\n'+r'$$Dx_1-h_{12}x_2=A_1,\qquad Dx_2-h_{12}x_1=A_2.$$'+'\n\n'+f'本题 $a=4,b=5,D=12,h=-1$，解为 $x_1={f(vals[0])},x_2={f(vals[1])}$，与迭代法一致。使用实际 18 个观测，误差自由度为 10，分析见（a）。')+'\n\n'+part('c','不同区组',r'若处理也不同，则 $h=-1$；若处理相同，则 $h=a-1$。统一解为'+ '\n\n'+r'$$\hat x_1=\frac{DA_1+hA_2}{D^2-h^2},\qquad \hat x_2=\frac{hA_1+DA_2}{D^2-h^2}.$$'+'\n\n带撇总和排除两个缺失值；仅在剩余设计保持连通且分母非零时适用。')+'\n\n'+part('d','同一区组',r'同一区组中的两个不同单元必属不同处理，因此 $h=b-1$。代入（c）的统一解即可。缺失值是拟合预测，不是新增独立观测，误差自由度减少 2。'), '参考答案无对应题；由 SSE 求导补充。')

for n,k in [(29,9),(30,12),(31,14),(32,16)]:
    a,r=report(k);arr=np.array(DATA[k]);A,B=arr.shape;msb=A*np.sum((arr.mean(axis=0)-arr.mean())**2)/(B-1);v=(msb-r['mse'])/A
    page(n,'随机区组方差分量估计',f'对习题 [4.{k}](./exercise-4-{k}.md)，假定区组是随机抽取的，估计区组方差分量。\n\n'+raw(k),r'**混合模型**：$y_{ij}=\mu+\tau_i+B_j+\epsilon_{ij}$，$\sum_i\tau_i=0$；$B_j\overset{iid}\sim N(0,\sigma_B^2)$，$\epsilon_{ij}\overset{iid}\sim N(0,\sigma^2)$，二者独立，无处理×区组交互作用。'+ '\n\n'+r'$$E(MS_E)=\sigma^2,\quad E(MS_B)=\sigma^2+a\sigma_B^2,\quad \hat\sigma^2=MS_E,\quad \hat\sigma_B^2=\frac{MS_B-MS_E}{a}.$$'+'\n\n'+f'$MS_B={f(msb)}$，$MS_E={f(r["mse"])}$，$a={A}$，故 $\\hat\\sigma_B^2={f(v)}$。'+ ('若要求非负方差分量，可将负的矩估计截为零；这会改变估计量性质。' if v<0 else '')+'\n\n下表用于展示平方和及计算基础；处理检验为 $H_0:\\tau_i=0$，备择为至少一个非零，采用 RCBD 处理 F 检验。\n\n'+RCBD+'\n\n'+anova(r))
page(33,'受试者方差分量的可识别性','对习题 [4.17](./exercise-4-17.md)，假定受试者是随机的，估计区组方差分量。\n\n'+raw(17),r'题目存在设计信息缺口：4.17 说明三个治疗组各有 10 名受试者，并未说明同一受试者接受三种治疗，也未说明 10 个匹配组。单纯“受试者随机”不能构成跨治疗区组。若每名受试者仅测一次，其随机效应与测量误差合并，$\sigma_S^2$ 与 $\sigma^2$ 无法分别识别。'+'\n\n'+r'**只有确认每列是一组匹配区组后**，才可用 $y_{ij}=\mu+\tau_i+B_j+\epsilon_{ij}$，假定 $B_j$ 与误差独立、正态，且方差分别为 $\sigma_B^2,\sigma^2$，计算 $\hat\sigma_B^2=(MS_B-MS_E)/3$。'+'\n\n'+f'作为明确附条件的演算，若强行把同列当作区组，$MS_E={f(compute(17)["mse"])}$，矩估计为 '+f((3*np.sum((np.array(DATA[17]).mean(axis=0)-np.array(DATA[17]).mean())**2)/9-compute(17)['mse'])/3)+'。这不是由现有题意可支持的区组方差结论。', '对应参考答案 4.30；补充说明其区组解释所需的配对前提。')

page(35,'拉丁方固定效应的最小二乘估计',r'对 $p\times p$ 拉丁方，行效应 $\alpha_i$、列效应 $\beta_k$ 和处理效应 $\tau_j$ 均固定，求这些模型参数的最小二乘估计。',LS+'\n\n'+r'在 $\sum_i\alpha_i=\sum_k\beta_k=\sum_j\tau_j=0$ 下，因每个处理在各行各列各出现一次，正规方程为'+ '\n\n'+r'$$p^2\hat\mu=G,\quad p\hat\mu+p\hat\alpha_i=R_i,\quad p\hat\mu+p\hat\beta_k=C_k,\quad p\hat\mu+p\hat\tau_j=T_j.$$'+'\n\n'+r'故 $\hat\mu=G/p^2$，$\hat\alpha_i=R_i/p-\hat\mu$，$\hat\beta_k=C_k/p-\hat\mu$，$\hat\tau_j=T_j/p-\hat\mu$。这是点估计公式；题目无数值数据，不虚构 ANOVA 数值。')
page(36,'拉丁方缺失值公式的推导','推导拉丁方单个缺失值的估计公式（正文式 4.28）。',r'设缺失值为 $x$；所在行、列、处理的已知总和为 $R^\prime,C^\prime,T^\prime$，所有已知观测总和为 $G^\prime$。令其余观测的平方和为 $K$，补值后的 SSE 中与 $x$ 有关的项为'+ '\n\n'+r'$$SSE(x)=K+x^2-\frac{(R^\prime+x)^2+(C^\prime+x)^2+(T^\prime+x)^2}{p}+\frac{2(G^\prime+x)^2}{p^2}+\text{与 }x\text{ 无关的项}.$$'+'\n\n'+r'求导并令零：$x-(R^\prime+C^\prime+T^\prime+3x)/p+2(G^\prime+x)/p^2=0$。所以'+ '\n\n'+r'$$\hat x=\frac{p(R^\prime+C^\prime+T^\prime)-2G^\prime}{(p-1)(p-2)}.$$'+'\n\n需 $p>2$ 且剩余设计参数可估计。补值仅用于拟合计算，误差自由度比完整拉丁方减少 1。')

q37=r'''为了增加处理重复，运行 $n$ 个 $p\times p$ 拉丁方。模型为

$$y_{ijkh}=\mu+\rho_h+\alpha_{i(h)}+\tau_j+\beta_{k(h)}+(\tau\rho)_{jh}+\epsilon_{ijkh},$$

其中 $\rho_h$ 为方效应，行列嵌套于方。采用 $\sum_h\rho_h=0$，每方行、列效应分别和为零，$\sum_j\tau_j=0$，交互作用按处理及方分别求和均为零。

**（a）** 写出正规方程并求参数估计。

**（b）** 写出方差分析表。
'''
a37=r'''误差独立同分布为 $N(0,\sigma^2)$；除显式加入的处理×方交互作用外，不含其他交互作用。所有效应暂按固定效应解释。令 $G_h,R_{ih},C_{kh},T_j,T_{jh}$ 分别为方、方内行、方内列、处理、处理×方总和。正规方程可写为

$$np^2\hat\mu=G,\quad p^2(\hat\mu+\hat\rho_h)=G_h,$$
$$p(\hat\mu+\hat\rho_h+\hat\alpha_{i(h)})=R_{ih},\quad p(\hat\mu+\hat\rho_h+\hat\beta_{k(h)})=C_{kh},$$
$$np(\hat\mu+\hat\tau_j)=T_j,\quad p[\hat\mu+\hat\rho_h+\hat\tau_j+\widehat{(\tau\rho)}_{jh}]=T_{jh}.$$

因此 $\hat\mu=\bar y$，$\hat\rho_h=\bar y_h-\bar y$，$\hat\tau_j=\bar y_j-\bar y$，$\hat\alpha_{i(h)}=\bar y_{ih}-\bar y_h$，$\hat\beta_{k(h)}=\bar y_{kh}-\bar y_h$，$\widehat{(\tau\rho)}_{jh}=\bar y_{jh}-\bar y_j-\bar y_h+\bar y$。
'''
b37=r'''设 $C=G^2/(np^2)$，$SS_T=\sum y^2-C$，各平方和及自由度如下。均方等于平方和除以自由度。
'''+ '\n\n'+tab(['来源','平方和','自由度','固定效应 F'],[['处理',r'$\sum_jT_j^2/(np)-C$',r'$p-1$',r'$MS_A/MS_E$'],['方',r'$\sum_hG_h^2/p^2-C$',r'$n-1$',r'$MS_S/MS_E$'],['行（方内）',r'$\sum_{ih}R_{ih}^2/p-\sum_hG_h^2/p^2$',r'$n(p-1)$',r'$MS_R/MS_E$'],['列（方内）',r'$\sum_{kh}C_{kh}^2/p-\sum_hG_h^2/p^2$',r'$n(p-1)$',r'$MS_C/MS_E$'],['处理×方',r'$\sum_{jh}T_{jh}^2/p-C-SS_A-SS_S$',r'$(p-1)(n-1)$',r'$MS_{AS}/MS_E$'],['误差','由总平方和减去以上各项',r'$n(p-1)(p-2)$','—'],['总计',r'$\sum y^2-C$',r'$np^2-1$','—']])+'\n\n'+r'处理检验 $H_0:\tau_j=0$ 对至少一个非零；交互作用检验 $H_0:(\tau\rho)_{jh}=0$ 对至少一个非零。若方为随机抽取并把处理×方也作为随机项，处理检验分母应改为 $MS_{AS}$，自由度为 $(p-1,(p-1)(n-1))$。有交互作用时还应解释各方中的处理表现。'
page(37,'多个拉丁方的联合分析',q37,part('a','模型、正规方程与估计',a37)+'\n\n'+part('b','符号 ANOVA 表',b37))
page(38,'拉丁方的样本量与功效','讨论如何确定拉丁方设计所需样本量。',r'固定处理效应模型下，$p$ 阶单个拉丁方有 $N=p^2$，处理自由度 $p-1$、误差自由度 $(p-1)(p-2)$，非中心参数为'+ '\n\n'+r'$$\Lambda=\frac{p\sum_j\tau_j^2}{\sigma^2},\quad \mathrm{Power}=P\{F_{p-1,(p-1)(p-2)}(\Lambda)>F_{1-\alpha;p-1,(p-1)(p-2)}\}.$$'+'\n\n预先指定有实际意义的处理差配置、误差方差估计、显著性水平和目标功效，再比较可行方案。处理数既定时不能随意增大 $p$ 来增加重复；应增加独立拉丁方数 $n$。若可假定共同处理效应且无处理×方交互作用，可按相应模型计算功效；若保留随机方交互作用，则必须用该方差分量及交互作用自由度规划，不能仍用单方的误差分母。')
a39,r39=report(27,dropcol=4)
page(39,'丢失整个日期的拉丁方','习题 [4.27](./exercise-4-27.md) 中，日期 5 的数据分析有误，必须全部舍弃。为余下数据建立适当分析。\n\n'+raw(27),a39+'\n\n删除第 5 列后剩 20 个观测，设计矩阵秩为 12，误差自由度为 8。行、列、处理不再完全正交，不能沿用完整拉丁方平方和，不能把错误的 5 个值重新用于估计。')

# Hyper-Graeco-Latin construction over GF(5).
hyper=np.array([[(i+j)%5,(i+2*j)%5,(i+3*j)%5] for i in range(5) for j in range(5)])
page(42,'五因子五水平超方的构造','构造用于研究 5 个因子的 $5\\times5$ 超方，并列出 ANOVA 表。',r'取行因子 $R=i$、列因子 $C=j$，$i,j\in\{0,1,2,3,4\}$；三个符号因子为 $A=i+j$、$B=i+2j$、$D=i+3j$，运算模 5。下表每格依次列 A、B、D（用 1～5 标记水平）。'+'\n\n'+tab(['行／列']+[str(i) for i in range(1,6)],[[i+1]+['/'.join(str((i+m*j)%5+1) for m in [1,2,3]) for j in range(5)] for i in range(5)])+'\n\n'+r'任意两因子的 25 个有序水平对都恰出现一次，因而两两正交。模型为 $y_{ij}=\mu+R_i+C_j+A_{i+j}+B_{i+2j}+D_{i+3j}+\epsilon_{ij}$；各效应和为零，误差独立正态、同方差，无交互作用。检验各因子 $H_0$：所有水平效应为零，对至少一个非零。令 $C_0=G^2/25$，每个因子的平方和均为 $\sum_lT_l^2/5-C_0$。'+ '\n\n'+tab(['来源','自由度','平方和','F'],[[v,4,r'$SS_'+v+'$',r'$(SS_'+v+r'/4)/(SS_E/4)$'] for v in ['R','C','A','B','D']]+[['误差',4,r'$\sum y^2-C_0-SS_R-SS_C-SS_A-SS_B-SS_D$','—'],['总计',24,r'$\sum y^2-C_0$','—']])+'\n\n每项检验分子、分母自由度均为 4。无响应数据，故此处给符号表。')

# Two actual squares; suppress Greek letters in the second.
def combined_squares():
    y=[];fac=[[],[],[],[],[]]
    for h,k in enumerate([28,41]):
        ls,values=LATIN[k]
        for i in range(4):
            for j in range(4):
                y.append(values[i][j]);t=ord(ls[i][j])-65
                for dest,v in zip(fac,[h,h*4+i,h*4+j,t,h*4+t]):dest.append(v)
    return np.array(y,float),[np.array(v) for v in fac]
yc,fc=combined_squares();sse,df,*_=fit(yc,fc);mse=sse/df
grand=yc.mean();ssS=16*sum((yc[fc[0]==h].mean()-grand)**2 for h in range(2));ssA=8*sum((yc[fc[3]==t].mean()-grand)**2 for t in range(4));ssR=4*sum((yc[fc[1]==v].mean()-yc[fc[0]==v//4].mean())**2 for v in range(8));ssC=4*sum((yc[fc[2]==v].mean()-yc[fc[0]==v//4].mean())**2 for v in range(8));ssAS=4*sum((yc[fc[4]==h*4+t].mean()-yc[fc[0]==h].mean()-yc[fc[3]==t].mean()+grand)**2 for h in range(2) for t in range(4))
rows=[]
for label,ss,d in [('方',ssS,1),('行（方内）',ssR,6),('列（方内）',ssC,6),('处理',ssA,3),('处理×方',ssAS,3)]:rows.append([label,d,ss,ss/d,(ss/d)/mse,stats.f.sf((ss/d)/mse,d,df)])
rows += [['误差',df,sse,mse,None,None],['总计',31,float(np.sum((yc-grand)**2)),None,None,None]]
page(43,'两个装配拉丁方的联合分析','将习题 [4.28](./exercise-4-28.md) 与 [4.41](./exercise-4-41.md) 联合分析。忽略 4.41 中希腊字母，使用习题 4.37 的方法。\n\n**第一个方**\n\n'+raw(28)+'\n\n**第二个方**\n\n'+raw(41),a37+'\n\n'+b37+'\n\n**代入本题数据（两方暂按固定效应）**\n\n'+anova({'rows':rows})+'\n\n'+f'处理×方检验 $P={f(stats.f.sf((ssAS/3)/mse,3,df))}$；处理主效应的固定效应检验 $P={f(stats.f.sf((ssA/3)/mse,3,df))}$。若把两方看作随机样本，处理的 $F=MS_A/MS_{{AS}}={f(ssA/ssAS)}$，自由度 $(3,3)$，$P={f(stats.f.sf(ssA/ssAS,3,3))}$。解释取决于方的抽样方式；第二方工作场所效应在忽略后可能进入残差，应谨慎比较。')

def bibd(k):
    arr=np.array(DATA[k],float);N=np.isfinite(arr).astype(float);a,b=arr.shape;rr=int(N.sum(axis=1)[0]);kk=int(N.sum(axis=0)[0]);lam=int((N@N.T)[0,1]);T=np.nansum(arr,axis=1);B=np.nansum(arr,axis=0);Q=T-N@B/kk;tau=kk*Q/(lam*a)
    return arr,N,a,b,rr,kk,lam,T,B,Q,tau
arr,N,a,b,rep,k,lam,T,B,Q,tau=bibd(45);contrasts=[[-1,-1,0,1,1],[1,-1,0,0,0],[0,0,0,-1,1],[-1,-1,4,-1,-1]];r45=compute(45)
cr=[]
for c in np.array(contrasts):
    L=c@tau;ss=L**2/((k/(lam*a))*(c@c));F=ss/r45['mse'];cr.append([','.join(map(str,c)),f(L),f(ss),f(F),f(stats.f.sf(F,1,r45['df']))])
page(46,'不完全区组的正交对比','为习题 [4.45](./exercise-4-45.md) 构造一组正交对比，计算各对比平方和。\n\n'+raw(45),r'本设计 $a=b=5,r=k=4,\lambda=3$。采用调整后效应 $\hat\tau_i=kQ_i/(\lambda a)$。对比系数满足 $\sum_i c_i=0$，估计量 $\hat L=c^T\hat\tau$，'+ '\n\n'+r'$$\operatorname{Var}(\hat L)=\frac{k\sigma^2}{\lambda a}\sum_i c_i^2,\quad SS_L=\frac{\hat L^2}{[k/(\lambda a)]\sum_i c_i^2},\quad F_L=SS_L/MS_E.$$'+'\n\n'+tab(['系数 c（处理1～5）','对比估计','平方和','F','P'],cr)+'\n\n'+f'各系数向量两两内积为零，平方和之和为 {f(sum(float(row[2]) for row in cr))}，等于调整后处理平方和（末位受显示精度影响）。检验 $H_0:L=0$ 对 $H_1:L\\ne0$，误差自由度为 11。未作多重检验校正时，对比 1 在 1% 水平显著；对比 2 在 5% 水平显著。')
a48,r48=report(48)
page(48,'催化剂例题的一般回归检验','用一般回归显著性检验分析第 9 版例 4.4 的催化剂不完全区组试验。\n\n'+raw(48),a48+'\n\n'+f'全模型秩为 7，$SSE_F={f(r48["sse"])}$；仅区组的简化模型秩为 4，$SSE_R={f(fit(r48["y"],r48["fac"][1:])[0])}$。相减得到处理调整平方和 22.75。参考答案 4.45S 在相减一行混入其他题的 138.78、99.25、39.53；本题按实际数据重算。')
page(49,'BIBD 调整处理平方和的证明',r'证明 BIBD 的调整处理平方和为 $k\sum_iQ_i^2/(\lambda a)$。',r'设处理—区组关联矩阵为 $N$，$T$、$B$ 是处理、区组总和向量。消去区组后的正规方程为 $C\hat\tau=Q$，其中'+ '\n\n'+r'$$C=rI-\frac1kNN^T,\quad Q=T-\frac1kNB.$$'+'\n\n'+r'BIBD 中 $NN^T=(r-\lambda)I+\lambda J$，并有 $r(k-1)=\lambda(a-1)$，故 $C=(\lambda a/k)(I-J/a)$。又 $\mathbf1^TQ=0$，在 $\sum_i\tau_i=0$ 下 $\hat\tau=kQ/(\lambda a)$。因此'+ '\n\n'+r'$$SS_{A\mid B}=\hat\tau^TQ=\frac{k}{\lambda a}\sum_iQ_i^2.$$'+'\n\n该量等于加入处理后减少的残差平方和，自由度为 $a-1$。')
blocks50=list(itertools.combinations(range(1,5),2))
page(50,'四处理二次运行区组的 BIBD','比较 4 个处理，每区组只能运行 2 次。构造一个含 6 个区组的 BIBD。',tab(['区组','处理'],[[j+1,','.join(map(str,v))] for j,v in enumerate(blocks50)])+'\n\n每个无序处理对恰好安排一次，所以 $a=4,b=6,k=2,r=3,\\lambda=1$。检验 $ar=bk=12$、$r(k-1)=\\lambda(a-1)=3$ 均成立。区组内运行顺序随机化。')
blocks51=[]
for u in range(1,8):
    for v in [0,1]:blocks51.append([x+1 for x in range(8) if (u&x).bit_count()%2==v])
page(51,'八处理十四区组的 BIBD','比较 8 个处理，每区组有 4 次运行。构造 $b=14,\\lambda=3$ 的 BIBD。',r'把 8 个处理编号对应于三维二元向量 $x\in GF(2)^3$。对每个非零 $u$，取 $u\cdot x=0$ 与 $u\cdot x=1$ 两个区组。'+ '\n\n'+tab(['区组','处理'],[[i+1,','.join(map(str,bb))] for i,bb in enumerate(blocks51)])+'\n\n每个区组有 4 个处理，每个处理出现 7 次；任意两个不同处理同处一个区组的次数为 3。因此 $a=8,b=14,k=4,r=7,\\lambda=3$，$ar=bk=56$。')

def interblock(n,parent):
    arr,N,a,b,rep,k,lam,T,B,Q,tau=bibd(parent);r=compute(parent)
    # Block variance via adjusted block SS, and GLS recovery from block totals.
    sse_t=fit(r['y'],r['fac'][:1])[0];ssbadj=sse_t-r['sse'];msbadj=ssbadj/(b-1)
    v=(msbadj-r['mse'])*(b-1)/(a*(rep-1));v=max(0,v)
    M=N@N.T;bi=np.linalg.pinv(M)@N@(B-k*B.mean()/k);bi-=bi.mean()
    # Symmetric BIBD: interblock information eigenvalue r-lambda.
    vintra=k*r['mse']/(lam*a);vinter=(k*r['mse']+k*k*v)/(rep-lam)
    comb=(tau/vintra+bi/vinter)/(1/vintra+1/vinter)
    return r'''**模型**：$y_{ij}=\mu+\tau_i+B_j+\epsilon_{ij}$，$\sum_i\tau_i=0$；随机区组 $B_j\sim N(0,\sigma_B^2)$，独立误差 $\epsilon_{ij}\sim N(0,\sigma^2)$，二者独立。区组总和向量满足

$$B=k\mu\mathbf1+N^T\tau+u,\quad \operatorname{Var}(u)=(k\sigma^2+k^2\sigma_B^2)I.$$

先由组内分析估计 $\hat\sigma^2=MS_E$，由调整区组平方和估计

$$\hat\sigma_B^2=\frac{(b-1)(MS_{B\mid A}-MS_E)}{a(r-1)}.$$

本题为对称 BIBD，$NN^T=(r-\lambda)I+\lambda J$，所以组间估计为 $\hat\tau^{\rm inter}=(r-\lambda)^{-1}N(B-\bar B\mathbf1)$。在对比子空间中，组内与组间估计的协方差分别为 $v_I(I-J/a)$、$v_B(I-J/a)$，其中

$$v_I=\frac{k\sigma^2}{\lambda a},\quad v_B=\frac{k\sigma^2+k^2\sigma_B^2}{r-\lambda},\quad \hat\tau^{\rm combined}=\frac{\hat\tau^{\rm intra}/\hat v_I+\hat\tau^{\rm inter}/\hat v_B}{1/\hat v_I+1/\hat v_B}.$$
'''+'\n\n'+f'$MS_{{B\\mid A}}={f(msbadj)}$，$\\hat\\sigma^2={f(r["mse"])}$，$\\hat\\sigma_B^2={f(v)}$。\n\n'+tab(['处理','组内效应','组间效应','合并效应'],[[i+1,f(tau[i]),f(bi[i]),f(comb[i])] for i in range(a)])+'\n\n组间信息可能很弱；合并权重由方差决定，不能简单平均两套估计。方差代入后属于经验估计，不能直接声称保留已知方差时的精确分布。'
for n,parent in [(52,45),(53,47)]:page(n,'不完全区组的组间信息恢复',f'对习题 [4.{parent}](./exercise-4-{parent}.md) 的设计进行组间分析。\n\n'+raw(parent),interblock(n,parent))
page(54,'BIBD 不存在的必要条件',r'验证 $a=8,r=8,k=4,b=16$ 的 BIBD 不存在。',r'虽然 $ar=64=bk$，另一必要条件要求 $\lambda=r(k-1)/(a-1)=24/7$。同一区组出现的处理对次数必须是非负整数，因此不可能存在这些参数的 BIBD。')
page(55,'BIBD 组内效应估计量的方差',r'证明组内处理效应估计量的方差为 $\operatorname{Var}(\hat\tau_i)=k(a-1)\sigma^2/(\lambda a^2)$。',r'由区组消元的正规方程，$C=(\lambda a/k)(I-J/a)$，$\hat\tau=C^+Q$，且 $\operatorname{Var}(Q)=\sigma^2C$，所以'+ '\n\n'+r'$$\operatorname{Cov}(\hat\tau)=\sigma^2C^+=\frac{k\sigma^2}{\lambda a}(I-J/a).$$'+'\n\n'+r'对角元为 $k(a-1)\sigma^2/(\lambda a^2)$；非对角元为 $-k\sigma^2/(\lambda a^2)$。因此 $\operatorname{Var}(\hat\tau_i-\hat\tau_j)=2k\sigma^2/(\lambda a)$，各处理效应估计并非相互独立。')
page(56,'扩展不完全区组设计',r'当 $a<k<2a$ 时，每区组先安排所有处理各一次，再增加一个区组大小 $k^*=k-a$ 的不完全区组设计。平衡时有 $r^*=r-b$ 和 $\lambda^*$。写出其统计分析，提示 $\lambda=2r-b+\lambda^*$。',r'''设附加设计关联矩阵为 $N^*$，完整设计的出现次数矩阵 $N=J+N^*$。每行总和为 $r=b+r^*$，每列总和为 $k=a+k^*$。注意单元出现次数可能为 2，应按每次真实运行记录响应。

**模型与假设**：$y_{ij\ell}=\mu+\tau_i+\beta_j+\epsilon_{ij\ell}$，误差独立同分布 $N(0,\sigma^2)$，处理、区组效应可加，效应和为零。$H_0:\tau_i=0$ 对至少一个非零。

因为 $NN^T=(r^*-\lambda^*)I+(b+2r^*+\lambda^*)J$，定义 $\lambda=b+2r^*+\lambda^*=2r-b+\lambda^*$。则

$$C=rI-NN^T/k=\frac{\lambda a}{k}(I-J/a),\quad Q=T-NB/k,\quad \hat\tau_i=\frac{kQ_i}{\lambda a}.$$

设总运行数 $N_0=ar=bk$，$G$ 为总响应，$SS_T=\sum y^2-G^2/N_0$，$SS_B=\sum_jB_j^2/k-G^2/N_0$，$SS_{A\mid B}=k\sum_iQ_i^2/(\lambda a)$，$SS_E=SS_T-SS_B-SS_{A\mid B}$。
'''+ '\n\n'+tab(['来源','自由度','平方和','F'],[['区组',r'$b-1$',r'$SS_B$','—'],['调整处理',r'$a-1$',r'$SS_{A\mid B}$',r'$[SS_{A\mid B}/(a-1)]/[SS_E/(N_0-a-b+1)]$'],['误差',r'$N_0-a-b+1$',r'$SS_E$','—'],['总计',r'$N_0-1$',r'$SS_T$','—']])+'\n\n这里的 λ 是出现次数矩阵两个不同处理行的内积，不能理解为每区组只计一次的处理对出现次数。重复单元若保留真实重复，还可进一步分离纯误差与不适合度。', '参考答案无直接对应题；按广义出现次数矩阵推导。')
page(57,'区组化的自由度代价','单因子有 5 个水平、3 次重复，原采用完全随机设计；改按区组设计后，纯误差自由度减少多少？选项：（a）3；（b）5；（c）2；（d）4；（e）以上皆非。',r'本题与 [4.2](./exercise-4-2.md) 相同。完全随机设计 $df_E=15-5=10$，随机完全区组设计 $df_E=(5-1)(3-1)=8$，减少 2，选 **（c）**。')

ms0=37.75/18;F0=(70.88/5)/ms0;p0=stats.f.sf(F0,5,18);ms1=25.75/15;F1=(70.88/5)/ms1
q58='物理研究生 Laura 完成一项单因子完全随机试验，ANOVA 汇总为：处理均方 14.18、误差平方和 37.75、总平方和 108.63、总自由度 23。\n\n'
for letter,question in zip('abcdefghijk',['求处理平方和。','求处理自由度。','求误差自由度。','求误差均方。','求检验统计量。','在 0.05 水平，“不拒绝原假设”是否正确？','给出 P 值上界。','给出 P 值下界。','有多少个处理水平？','每个水平重复多少次？','假设实际采用随机完全区组设计，区组平方和为 12，重建 ANOVA 表。区组化使误差估计减少多少？']):q58+=f'**（{letter}）** {question}\n\n'
answers=[r'$SS_A=108.63-37.75=70.88$。',r'$df_A=SS_A/MS_A\approx70.88/14.18\approx5$，差异来自均方四舍五入。',r'$df_E=23-5=18$。',f'$MS_E=37.75/18={f(ms0)}$。',RCBD.replace('RCBD','RCBD')+'\n\n本部分采用完全随机模型 $y_{ij}=\\mu+\\tau_i+\\epsilon_{ij}$，误差独立正态、同方差。$H_0$ 为各均值相等，备择为至少一个不同；$F=(SS_A/5)/(SS_E/18)$。\n\n'+anova({'rows':[['处理',5,70.88,70.88/5,F0,p0],['误差',18,37.75,ms0,None,None],['总计',23,108.63,None,None,None]]}),f'**不正确**。$P={f(p0)}<0.05$，拒绝原假设。',f'$P={f(p0)}<0.005$。',f'$P={f(p0)}>0.001$。不能将 0.001 当作上界。',r'$a=df_A+1=6$。',r'$N=24$，平衡重复数 $b=24/6=4$。',RCBD+'\n\n'+anova({'rows':[['处理',5,70.88,70.88/5,F1,stats.f.sf(F1,5,15)],['区组',3,12,4,None,None],['误差',15,25.75,ms1,None,None],['总计',23,108.63,None,None,None]]})+'\n\n'+f'区组化后误差自由度为 $18-3=15$，误差均方从 {f(ms0)} 降至 {f(ms1)}，减少 {f(ms0-ms1)}，相对减少 {f(100*(ms0-ms1)/ms0)}%。参考答案误保留 18 个误差自由度，从而高估了误差减少比例。']
# Remove irrelevant RCBD block from the completely randomized part.
answers[4]=answers[4].split('本部分采用',1)[1];answers[4]='本部分采用'+answers[4]
page(58,'Laura 的 ANOVA 表与区组化',q58,'\n\n'.join(part(c,'计算与解释',a) for c,a in zip('abcdefghijk',answers)), '对应参考答案 4.58；更正 P 值界限、区组化后的误差自由度及方差减少比例。')
a59,r59=report(13,oneway=True);rb=compute(13)
page(59,'直邮试验忽略地区区组的影响','假定习题 [4.13](./exercise-4-13.md) 以完全随机设计进行，忽略地区差异但得到完全相同的数据。重新分析，比较忽略区组对结论的影响。\n\n'+raw(13),a59+'\n\n'+f'区组模型 $MS_E={f(rb["mse"])}$，完全随机模型 $MS_E={f(r59["mse"])}$；后者把地区差异并入误差。处理平方和保持不变，但误差自由度由 6 增为 9。完全随机模型的处理检验仍显著；区组化显著提高了比较精度。\n\n'+diagnostics(59,r59), '无直接对应题；按第 9 版题意重新计算。')

# Match by content, not by numbering.
SOURCES={1:'4.1',2:'4.2S',3:'4.3',4:'4.4',5:'4.5',7:'4.6',8:'4.7S',9:'4.8',10:'4.9',11:'4.10',12:'4.11',13:'4.12S',14:'4.13',16:'4.14',17:'4.15',18:'4.16',19:'4.17S',20:'4.18',21:'4.19',22:'4.20',23:'4.21',24:'4.22',26:'4.23',27:'4.24',28:'4.25S',29:'4.26',30:'4.27S',31:'4.28',32:'4.29',33:'4.30',34:'4.31',35:'4.32',36:'4.33',37:'4.34',38:'4.35',39:'4.36',40:'4.37',41:'4.38',42:'4.39',43:'4.40',44:'4.41',45:'4.42S',46:'4.43S',47:'4.44S',48:'4.45S',49:'4.46S',50:'4.47S',51:'4.48',52:'4.49',53:'4.50',54:'4.51S',55:'4.52S',57:'4.53',58:'4.58'}
assert sorted(P)==list(range(1,60))
manifest=[]
for n,(title,q,a,source) in sorted(P.items()):
    if source is None:source=f'按题意对应 `Exercises/ch04.pdf` 的 {SOURCES[n]}；公式与数值由编者复核，补充完整推导与诊断。' if n in SOURCES else '参考答案无直接对应题；由编者补充推导、计算及诊断。'
    nav=[]
    nav.append(f'← [上一题 {"3.71" if n==1 else "4."+str(n-1)}]({"../chapter03/exercise-3-71.md" if n==1 else "./exercise-4-"+str(n-1)+".md"})')
    nav.append('[本章题目目录](../chapter04.md)')
    if n<59:nav.append(f'[下一题 4.{n+1}](./exercise-4-{n+1}.md) →')
    content=f'''---
title: 习题 4.{n}：{title}
short_title: 4.{n} {title}
description: 第9版习题4.{n}的中文题目、数据与答案解析
---

(exercise-4-{n})=
# 习题 4.{n}：{title}

[习题与答案](../../index.md) / [第 4 章　随机区组、拉丁方与相关设计](../chapter04.md)

第 9 版 · 第 4 章 · 原书章末习题见印刷页 171～178。题目依据第 9 版，参考答案按题意匹配。

## 题目

{q}

## 答案与解析

{a}

---

**答案来源**：{source}

{' · '.join(nav)}
'''
    (OUT/f'exercise-4-{n}.md').write_text(content,encoding='utf8')
    manifest.append(dict(number=f'4.{n}',title=title,file=f'Exercises/translations_zh/chapter04/exercise-4-{n}.md'))
(OUT.parent/'chapter04.md').write_text('''---
title: 第 4 章习题：随机区组、拉丁方与相关设计
short_title: 第 4 章习题
description: 第9版第四章59道习题及解答
---

# 第 4 章习题：随机区组、拉丁方与相关设计

[全部章节](../index.md) · [本章正文](../../Chapter%2004/translations_zh/00_Randomized_Blocks_Latin_Squares_and_Related_Designs_zh.md)

本章收录第 **9 版 4.1～4.59，共 59 题**。各题给出题目、数据、公式与解析；有原始数据且涉及模型评价时提供残差图。参考答案文件按题意匹配，不直接沿用其错位编号。

4.17 未说明跨治疗组配对，因此按独立组分析；4.33 说明随机受试者方差分量的可识别性限制。4.58 更正参考答案的误差自由度与 P 值界限。缺失观测同时给出补值法与真实观测上的精确检验。

## 题目目录

'''+tab(['题号','题目'],[[f'[4.{n}](./chapter04/exercise-4-{n}.md)',P[n][0]] for n in range(1,60)])+'\n',encoding='utf8')
(ROOT/'Exercises/data/chapter04-pages.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
(ROOT/'Exercises/data/chapter04-results.json').write_text(json.dumps(RESULT,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
(ROOT/'Exercises/data/chapter04-source.json').write_text(json.dumps({'rectangular':DATA,'latin_squares':LATIN,'greek_symbols':{k:v.tolist() for k,v in GREEK.items()},'source_mapping':SOURCES},ensure_ascii=False,indent=2)+'\n',encoding='utf8')
# Integrate without rewriting earlier chapter content.
manifest_path=ROOT/'Exercises/data/pages.json'
existing=json.loads(manifest_path.read_text(encoding='utf8'))
manifest_path.write_text(json.dumps([v for v in existing if not v['number'].startswith('4.')]+manifest,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
tocpath=ROOT/'toc.yml';toc=tocpath.read_text(encoding='utf8')
start=toc.find('        - file: Exercises/translations_zh/chapter04.md')
end=toc.index('    - file: Chapter 01/')
if start>=0:toc=toc[:start]+toc[end:]
end=toc.index('    - file: Chapter 01/')
lines=['        - file: Exercises/translations_zh/chapter04.md','          title: 第 4 章　随机区组、拉丁方与相关设计习题','          children:']
for v in manifest:lines.extend([f'            - file: {v["file"]}',f'              title: "{v["number"]} {v["title"]}"'])
tocpath.write_text(toc[:end]+'\n'.join(lines)+'\n'+toc[end:],encoding='utf8')
last=ROOT/'Exercises/translations_zh/chapter03/exercise-3-71.md';s=last.read_text(encoding='utf8')
if '[下一题 4.1]' not in s:last.write_text(s.rstrip()+' · [下一题 4.1](../chapter04/exercise-4-1.md) →\n',encoding='utf8')
print('Generated 59 exercise pages, chapter index, numerical results and figures.')
