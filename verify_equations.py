"""Independent numerical check in Python; this is not a MATLAB execution.
Requires numpy/scipy/matplotlib only for developer verification, NOT MATLAB use.
"""
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent
# Same manually digitized pixel positions and axes as the MATLAB file.
q=np.array([[210,519],[259,262],[271,267],[284,219],[295,170],[307,189],[318,170],[330,171],[341,197],[354,179]])
C=np.c_[(q[:,0]-210)*200/98,(558-q[:,1])*5000/54]
Q=[[[172,351],[193,334],[213,310],[233,263],[253,217],[273,185],[293,176],[311,194],[332,197],[352,183]],[[172,733],[193,732],[213,725],[233,716],[253,699],[273,602],[293,568],[304,550],[321,549],[343,571]],[[172,1123],[193,1127],[213,1129],[233,1132],[253,1131],[273,1123],[293,1092],[304,1050],[312,1015],[371,936],[383,957],[392,928]]]
cs=[]
for q,d,y0 in zip(Q,[156.25,312.5,625],[381,771,1160]):
 q=np.array(q);cs.append(dict(dose=d,t=(q[:,0]-172)*200/82.5,y=(y0-q[:,1])*10000/82.5))
kEF=3.5e-6;R=250/.004629

def logistic(g,t,n):return g[1]/(1+(g[1]/n-1)*np.exp(-g[0]*t))
g=np.exp(least_squares(lambda z:(logistic(np.exp(z),C[:,0],C[0,1])-C[:,1])/max(C[:,1]),np.log([.03,35000])).x)
def simulate(p,t,c,rtol=1e-8):
 f,fb,b,beta=p;t=np.asarray(t);n0=c['y'][0]
 def rhs(ce):
  def fun(t,y):
   return [kEF*R*ce-(f+fb)*y[0]+b*y[1],fb*y[0]-b*y[1],g[0]*(1-np.exp(y[2])/g[1])-beta*y[1]]
  return fun
 s1=solve_ivp(rhs(c['dose']),[0,6],[0,0,np.log(n0)],method='LSODA',rtol=rtol,atol=rtol/100,dense_output=True)
 assert s1.success
 y=np.empty((len(t),3));mask=t<=6;y[mask]=s1.sol(t[mask]).T
 if np.any(~mask):
  s2=solve_ivp(rhs(0),[6,max(t)],s1.y[:,-1],method='LSODA',rtol=rtol,atol=rtol/100,dense_output=True)
  assert s2.success;y[~mask]=s2.sol(t[~mask]).T
 y[:,2]=np.exp(y[:,2]);return y

def residual(p):return np.concatenate([(simulate(p,c['t'],c)[:,2]-c['y'])/max(c['y']) for c in cs])
lo=np.array([1e-4,1e-4,1e-7,1e-9]);hi=np.array([10,10,10,.1])
def fit(mapper,guess,lower,upper):
 fit=least_squares(lambda z:residual(mapper(np.exp(z))),np.log(guess),bounds=(np.log(lower),np.log(upper)),diff_step=1e-4,max_nfev=250,ftol=1e-9,xtol=1e-9,gtol=1e-8)
 return mapper(np.exp(fit.x)),sum(fit.fun**2)
pA,a=fit(lambda x:np.r_[1.25,.21,0,x[0]],[.0005],lo[3:],hi[3:]);print('A',pA,a,flush=True)
pB,b=fit(lambda x:np.array([x[0],x[1],0,x[2]]),[1.25,.21,pA[3]],lo[[0,1,3]],hi[[0,1,3]]);print('B',pB,b,flush=True)
rev=[]
for guess in [[1.25,.21,.005,.0005],[.03,.03,.005,.0001],[.1,.1,.02,.0005]]:
 pp,ss=fit(lambda x:x,guess,lo,hi);rev.append((ss,pp));print('C',pp,ss,flush=True)
c,pC=min(rev,key=lambda x:x[0]);P=np.array([pA,pB,pC])
# Exact equivalent family from invariant a,b,c.
aSum=sum(pC[:3]);bProd=pC[0]*pC[2];gain=pC[3]*pC[1];fast=(aSum+np.sqrt(aSum*aSum-4*bProd))/2;slow=bProd/fast
family=[pC]
for bf in np.geomspace(slow*1.01,fast*.99,30):
 fe=bProd/bf;fb=aSum-bf-fe;be=gain/fb;pp=np.array([fe,fb,bf,be])
 if np.all(pp>=lo) and np.all(pp<=hi):family.append(pp)
family=np.array(family);t=np.unique(np.r_[np.linspace(0,24,100),np.linspace(24,600,200)]);ref=simulate(pC,t,cs[2],1e-10)
maxdiff=0
for pp in family:
 assert np.allclose([sum(pp[:3]),pp[0]*pp[2],pp[3]*pp[1]],[aSum,bProd,gain],rtol=1e-9)
 Y=simulate(pp,t,cs[2],1e-10);maxdiff=max(maxdiff,np.max(np.abs(Y[:,2]-ref[:,2]))/max(ref[:,2]))
assert maxdiff<1e-5
# Washout states continuous by construction. Drug flux identity holds pointwise.
y=np.array([12.0,25.0,1000]);fe,fb,bf,be=pC
cf_dot=kEF*R*625-(fe+fb)*y[0]+bf*y[1];cb_dot=fb*y[0]-bf*y[1]
assert np.isclose(cf_dot+cb_dot,kEF*R*625-fe*y[0])
assert np.all(ref[:,:2]>=-1e-7) and np.all(ref[:,2]>0)
# With beta=0, drug dynamics must not affect logistic growth.
null=pC.copy();null[3]=0
nullY=simulate(null,t,cs[2]);assert np.max(np.abs(nullY[:,2]-logistic(g,t,cs[2]['y'][0])))/g[1]<1e-5
# Tighter tolerances should not materially change the objective.
precise=np.concatenate([(simulate(pC,cc['t'],cc,1e-10)[:,2]-cc['y'])/max(cc['y']) for cc in cs])
assert abs(sum(precise**2)-c)<1e-4
fig,axs=plt.subplots(2,2,figsize=(11,7));ax=axs.flat[0];tt=np.linspace(0,300,200);ax.plot(C[:,0],C[:,1],'ko');ax.plot(tt,logistic(g,tt,C[0,1]),'k-');ax.set_title('SUM-149PT control')
for ax,cc in zip(list(axs.flat)[1:],cs):
 tt=np.linspace(0,max(cc['t']),300);ax.plot(cc['t'],cc['y'],'ko',label='Digitized subset')
 for pp,style,label in zip(P,[':', '--','-'],['Old PK / direct PD','Joint irreversible','Joint reversible']):ax.plot(tt,simulate(pp,tt,cc)[:,2],style,label=label)
 ax.set_title(f"{cc['dose']:g} nM, 6 h");ax.legend(fontsize=8)
for ax in axs.flat:ax.set_xlabel('Hours from treatment start');ax.set_ylabel('Cell count');ax.grid(alpha=.2)
fig.suptitle('Independent Python check (not MATLAB output)');fig.tight_layout();fig.savefig(ROOT/'verification_preview.png',dpi=140)
report=f'''Independent Python numerical verification, NOT MATLAB execution.
Growth kp={g[0]:.8g}, theta={g[1]:.8g}
Scaled SSE A/B/C: {a:.8g}, {b:.8g}, {c:.8g}
One reversible representative [kFE,kFB,kBF,beta]: {pC}
(a,b,c) combinations: {aSum}, {bProd}, {gain*kEF*R}
Equivalent sets checked: {len(family)}
Maximum relative N discrepancy between equivalent sets: {maxdiff:.5g}
PASS: mass-balance cancellation; positive trajectories; beta=0 logistic limit;
tolerance sensitivity; exact combination invariance; equal N predictions.
Individual PK rates are not uniquely estimated. Python optimizer differs from
MATLAB fminsearch, so representative values may differ even with similar fits.
'''
(ROOT/'verification_report.txt').write_text(report);print(report)
