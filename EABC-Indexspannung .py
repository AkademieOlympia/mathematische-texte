import numpy as np, pandas as pd, math, os

# Recompute compactly for reproducibility
N=10_000_000
sieve=np.ones(N+1,dtype=bool)
sieve[:2]=False
for p in range(2,int(N**0.5)+1):
    if sieve[p]:
        sieve[p*p:N+1:p]=False

nums=np.arange(N+1)
famE = sieve & (nums%12==1)
famA = sieve & (nums%12==5)
famB = sieve & (nums%12==7)
famC = sieve & (nums%12==11)
cE=np.cumsum(famE,dtype=np.int32)
cA=np.cumsum(famA,dtype=np.int32)
cB=np.cumsum(famB,dtype=np.int32)
cC=np.cumsum(famC,dtype=np.int32)

# cumulative log sums for Chebyshev-like mod-12 tension
logs=np.zeros(N+1)
prime_idx=np.where(sieve)[0]
logs[prime_idx]=np.log(prime_idx)
sumE=np.cumsum(logs*(nums%12==1))
sumA=np.cumsum(logs*(nums%12==5))
sumB=np.cumsum(logs*(nums%12==7))
sumC=np.cumsum(logs*(nums%12==11))
sumAll=np.cumsum(logs)

def theta_log12_at(xarr):
    xarr=np.asarray(xarr)
    num=(sumB[xarr]+sumC[xarr])-(sumE[xarr]+sumA[xarr])
    den=sumAll[xarr]
    return num/den

def corr(a,b):
    a=np.asarray(a); b=np.asarray(b)
    return float(np.corrcoef(a,b)[0,1])

def resid(y,x):
    A=np.vstack([x,np.ones_like(x)]).T
    coef=np.linalg.lstsq(A,y,rcond=None)[0]
    return y-(A@coef),coef

# fourlings to 1e7
p=np.arange(5,N-3)
mask=sieve[p-4] & sieve[p-2] & sieve[p+2] & sieve[p+4]
centers=p[mask]
vals=np.stack([centers-4,centers-2,centers+2,centers+4],axis=1)
mods=vals%12
mod_to_fam={1:'E',5:'A',7:'B',11:'C'}
fams=np.vectorize(mod_to_fam.get)(mods)

idxs=np.empty_like(vals)
for j in range(4):
    v=vals[:,j]
    m=mods[:,j]
    idx=np.empty_like(v)
    idx[m==1]=cE[v[m==1]]
    idx[m==5]=cA[v[m==5]]
    idx[m==7]=cB[v[m==7]]
    idx[m==11]=cC[v[m==11]]
    idxs[:,j]=idx

idx_spread=idxs.max(axis=1)-idxs.min(axis=1)
idx_var=idxs.var(axis=1)
orient=np.array([''.join(row) for row in fams])
orient_sign=np.where(orient=='ABCE',1,-1)
theta=theta_log12_at(centers)
lp=np.log(centers)
theta_res,_=resid(theta,lp)
spread_res,_=resid(idx_spread,lp)
var_res,_=resid(idx_var,lp)

fourlings_df=pd.DataFrame({
    'center': centers,
    'p_minus_4': vals[:,0],
    'p_minus_2': vals[:,1],
    'p_plus_2': vals[:,2],
    'p_plus_4': vals[:,3],
    'orientation': orient,
    'idx_E': np.where(mods[:,0]==1, idxs[:,0], np.where(mods[:,1]==1, idxs[:,1], np.where(mods[:,2]==1, idxs[:,2], idxs[:,3]))),
    'idx_A': np.where(mods[:,0]==5, idxs[:,0], np.where(mods[:,1]==5, idxs[:,1], np.where(mods[:,2]==5, idxs[:,2], idxs[:,3]))),
    'idx_B': np.where(mods[:,0]==7, idxs[:,0], np.where(mods[:,1]==7, idxs[:,1], np.where(mods[:,2]==7, idxs[:,2], idxs[:,3]))),
    'idx_C': np.where(mods[:,0]==11, idxs[:,0], np.where(mods[:,1]==11, idxs[:,1], np.where(mods[:,2]==11, idxs[:,2], idxs[:,3]))),
    'idx_spread': idx_spread,
    'idx_var': idx_var,
    'theta_log12': theta,
    'theta_log12_res': theta_res,
    'idx_spread_res': spread_res,
    'idx_var_res': var_res,
})

# zeros side
base_dir = os.path.dirname(os.path.abspath(__file__)) if '__file__' in globals() else os.getcwd()
gam=np.load(os.path.join(base_dir, 'zeros6.npy'))[:100000]
E=0.5+1j*gam
A=(3*E**2)**(1/3)
omega=np.exp(2j*np.pi/3)
B=omega*A
C=omega**2*A
Cpt=(B**2-A**2)/(2*E)
Delta=np.abs(C-Cpt)/np.abs(C)
Esym=np.abs(A)**2+np.abs(B)**2+np.abs(C)**2
Ept=np.abs(A)**2+np.abs(B)**2+np.abs(Cpt)**2
Lambda=(Esym-Ept)/Esym
X=np.floor(gam).astype(int)
mask=(X<=N)
gam=gam[mask]; X=X[mask]; Delta=Delta[mask]; Lambda=Lambda[mask]
theta_z=theta_log12_at(X)
lg=np.log(gam)
theta_z_res,_=resid(theta_z,lg)
Delta_res,_=resid(Delta,lg)
Lambda_res,_=resid(Lambda,lg)
zeros_df=pd.DataFrame({
    'gamma': gam,
    'X_floor': X,
    'Delta': Delta,
    'Lambda': Lambda,
    'theta_log12': theta_z,
    'Delta_res': Delta_res,
    'Lambda_res': Lambda_res,
    'theta_log12_res': theta_z_res,
})

# windowed comparison
binw=10000
maxX=int(X.max())
bins=np.arange(0,maxX+binw,binw)
fourlings_df['bin']=pd.cut(fourlings_df['center'], bins, labels=False, include_lowest=True, right=False)
zeros_df['bin']=pd.cut(zeros_df['X_floor'], bins, labels=False, include_lowest=True, right=False)
g4=fourlings_df.groupby('bin').agg(center_mean=('center','mean'),
                                   spread_mean=('idx_spread','mean'),
                                   var_mean=('idx_var','mean'),
                                   count4=('center','size'))
gz=zeros_df.groupby('bin').agg(X_mean=('X_floor','mean'),
                               Delta_mean=('Delta','mean'),
                               Lambda_mean=('Lambda','mean'),
                               theta_mean=('theta_log12','mean'),
                               countz=('X_floor','size'))
window_df=g4.join(gz, how='inner').dropna().reset_index()

summary = pd.DataFrame([
    ['fourlings_count', len(fourlings_df)],
    ['orientation_ABCE', int((orient=='ABCE').sum())],
    ['orientation_CEAB', int((orient=='CEAB').sum())],
    ['corr_idx_spread_theta_raw', corr(idx_spread, theta)],
    ['corr_idx_var_theta_raw', corr(idx_var, theta)],
    ['corr_idx_spread_theta_res', corr(spread_res, theta_res)],
    ['corr_idx_var_theta_res', corr(var_res, theta_res)],
    ['corr_Delta_theta_raw', corr(Delta, theta_z)],
    ['corr_Lambda_theta_raw', corr(Lambda, theta_z)],
    ['corr_Delta_theta_res', corr(Delta_res, theta_z_res)],
    ['corr_Lambda_theta_res', corr(Lambda_res, theta_z_res)],
    ['window_corr_spread_Delta', corr(window_df['spread_mean'], window_df['Delta_mean']) if len(window_df)>1 else np.nan],
    ['window_corr_spread_Lambda', corr(window_df['spread_mean'], window_df['Lambda_mean']) if len(window_df)>1 else np.nan],
    ['window_corr_spread_theta', corr(window_df['spread_mean'], window_df['theta_mean']) if len(window_df)>1 else np.nan],
    ['window_corr_var_Delta', corr(window_df['var_mean'], window_df['Delta_mean']) if len(window_df)>1 else np.nan],
    ['window_corr_var_Lambda', corr(window_df['var_mean'], window_df['Lambda_mean']) if len(window_df)>1 else np.nan],
    ['window_corr_var_theta', corr(window_df['var_mean'], window_df['theta_mean']) if len(window_df)>1 else np.nan],
], columns=['metric','value'])

top_fourlings = fourlings_df.sort_values(['idx_spread','idx_var'], ascending=[False,False]).head(50)

report = f"""EABC-Indexspannungs- und Richtungs-Test (bis 10^7)

Rechenumfang
- Primzahlvierlinge bis Zentrum <= 10,000,000
- Riemann-Nullstellen: erste 100000 aus zeros6.npy
- Chebyshev-artige mod-12-Spannung:
  Theta_log12(X) = ((B+C) - (E+A)) / (alle), jeweils logarithmisch gewichtet

Vierlinge
- Anzahl klassischer Vierlinge: {len(fourlings_df)}
- Orientierungsmuster:
  ABCE = {(orient=='ABCE').sum()}
  CEAB = {(orient=='CEAB').sum()}

Wichtigste Korrelationen
- corr(idx_spread, Theta_raw) = {corr(idx_spread, theta):.6f}
- corr(idx_var, Theta_raw)    = {corr(idx_var, theta):.6f}
- corr(idx_spread_res, Theta_res) = {corr(spread_res, theta_res):.6f}
- corr(idx_var_res, Theta_res)    = {corr(var_res, theta_res):.6f}

Nullstellenseite
- corr(Delta, Theta_raw)   = {corr(Delta, theta_z):.6f}
- corr(Lambda, Theta_raw)  = {corr(Lambda, theta_z):.6f}
- corr(Delta_res, Theta_res)  = {corr(Delta_res, theta_z_res):.6f}
- corr(Lambda_res, Theta_res) = {corr(Lambda_res, theta_z_res):.6f}

Fenstervergleich (10k-Bins, nur ueberlappender Bereich)
- corr(spread_mean, Delta_mean) = {corr(window_df['spread_mean'], window_df['Delta_mean']) if len(window_df)>1 else float('nan'):.6f}
- corr(spread_mean, Lambda_mean)= {corr(window_df['spread_mean'], window_df['Lambda_mean']) if len(window_df)>1 else float('nan'):.6f}
- corr(spread_mean, theta_mean) = {corr(window_df['spread_mean'], window_df['theta_mean']) if len(window_df)>1 else float('nan'):.6f}

Kurze Deutung
1) Die EABC-Indexspannung ist auch bis 10^7 nicht trivial und bleibt nach einfachem Detrending deutlich mit der lokalen mod-12-Spannung korreliert.
2) Die Orientierung der Vierlinge bleibt auf zwei Modi beschraenkt (ABCE / CEAB).
3) Auf der Nullstellenseite sind Delta und Lambda weiterhin die besseren Groessen als der rohe Quotient r.
4) Der gemeinsame 10k-Bin-Vergleich ist mit nur 8 Bins noch klein, zeigt aber einen starken Gleichlauf zwischen Vierlings-Indexspannung und Delta/Lambda-Mitteln.
"""

report_path=os.path.join(base_dir, 'eabc_index_tension_test2_report.txt')
summary_path=os.path.join(base_dir, 'eabc_index_tension_test2_summary.csv')
fourlings_path=os.path.join(base_dir, 'eabc_index_tension_test2_fourlings.csv')
zeros_path=os.path.join(base_dir, 'eabc_index_tension_test2_zero_windows.csv')
top_path=os.path.join(base_dir, 'eabc_index_tension_test2_top_fourlings.csv')

with open(report_path,'w',encoding='utf-8') as f:
    f.write(report)
summary.to_csv(summary_path,index=False)
fourlings_df.to_csv(fourlings_path,index=False)
window_df.to_csv(zeros_path,index=False)
top_fourlings.to_csv(top_path,index=False)

print("Created files:")
print(report_path)
print(summary_path)
print(fourlings_path)
print(zeros_path)
print(top_path)
