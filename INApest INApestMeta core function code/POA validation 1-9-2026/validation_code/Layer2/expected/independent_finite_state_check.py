"""Independent finite-state check of the Layer-2 trajectory fixtures.
This does not call INApestPoA.R. It verifies that the synthetic trajectory sets
encode the intended reinvasion/extinction/management-feedback probabilities.
"""
import math, csv, os

def impose_prior(absent, prior):
    na=sum(absent); np=len(absent)-na
    if na==0 or np==0: raise ValueError('both classes required')
    return [prior/na if a else (1-prior)/np for a in absent]

def posterior_zero(weights, absent, q):
    L=q
    z=sum(w*l for w,l in zip(weights,L))
    post=[w*l/z for w,l in zip(weights,L)]
    return post, sum(w for w,a in zip(post,absent) if a)

def rescale_for_propagation(priorw, match, absent, postw):
    ta=sum(w for w,a in zip(postw,absent) if a)
    tp=1-ta
    out=[0.0]*len(priorw)
    for want_abs,target in [(True,ta),(False,tp)]:
        idx=[i for i,(m,a) in enumerate(zip(match,absent)) if m and a==want_abs]
        den=sum(priorw[i] for i in idx)
        if target>0 and den<=0: raise RuntimeError('unsupported class')
        for i in idx: out[i]=priorw[i]/den*target
    s=sum(out); return [x/s for x in out]

def run_two_round(present, events, prior, pdet):
    n=len(present)
    absent1=[not bool(present[i][0]) for i in range(n)]
    w=impose_prior(absent1,prior)
    q1=[1.0 if absent1[i] else 1-pdet for i in range(n)]
    post1,poa1=posterior_zero(w,absent1,q1)
    match=[events[i][0]==0 for i in range(n)]
    w2=rescale_for_propagation(w,match,absent1,post1)
    absent2=[not bool(present[i][1]) for i in range(n)]
    prior2=sum(x for x,a in zip(w2,absent2) if a)
    q2=[1.0 if absent2[i] else 1-pdet for i in range(n)]
    post2,poa2=posterior_zero(w2,absent2,q2)
    return poa1,prior2,poa2

rows=[]
def chk(id,got,expected,label,tol=1e-12):
    err=abs(got-expected); ok=err<=tol
    rows.append(dict(test_id=id,quantity=label,observed=got,expected=expected,abs_error=err,status='PASS' if ok else 'FAIL'))
    if not ok: raise AssertionError((id,label,got,expected))

# L2-03 reinvasion 10%
present=[[0,0] for _ in range(9)]+[[0,1]]+[[1,1] for _ in range(10)]
events=[[0,0] for _ in range(20)]
a,b,c=run_two_round(present,events,.8,.8)
e1=.8/(.8+.2*.2); e2=e1*.9; e3=e2/(e2+(1-e2)*.2)
for got,exp,lbl in [(a,e1,'posterior_round_1'),(b,e2,'prior_round_2'),(c,e3,'posterior_round_2')]: chk('L2-03',got,exp,lbl)

# L2-04 extinction 30%, reinvasion 20%
present=[[0,0] for _ in range(8)]+[[0,1] for _ in range(2)]+[[1,0] for _ in range(3)]+[[1,1] for _ in range(7)]
events=[[0,0] for _ in range(20)]
a,b,c=run_two_round(present,events,.7,.6)
e1=.7/(.7+.3*.4); e2=e1*.8+(1-e1)*.3; e3=e2/(e2+(1-e2)*.4)
for got,exp,lbl in [(a,e1,'posterior_round_1'),(b,e2,'prior_round_2'),(c,e3,'posterior_round_2')]: chk('L2-04',got,exp,lbl)

# L2-06 detected present trajectory becomes absent, but observed zero excludes it
present=[[0,0],[0,0],[1,0],[1,1]]
events=[[0,0],[0,0],[1,0],[0,0]]
a,b,c=run_two_round(present,events,.5,.5)
for got,exp,lbl in [(a,2/3,'posterior_round_1'),(b,2/3,'prior_round_2'),(c,.8,'posterior_round_2')]: chk('L2-06',got,exp,lbl)

# Unsupported present class check.
present=[[0,0],[0,0],[1,1],[1,1]]; events=[[0,0],[0,0],[1,0],[1,0]]
try:
    run_two_round(present,events,.5,.5)
    raise AssertionError('L2-08 failed to trigger support error')
except RuntimeError:
    rows.append(dict(test_id='L2-08',quantity='unsupported_present_class',observed=1,expected=1,abs_error=0,status='PASS'))

path=os.path.join(os.path.dirname(__file__),'independent_finite_state_check.csv')
with open(path,'w',newline='',encoding='utf-8') as f:
    w=csv.DictWriter(f,fieldnames=rows[0].keys()); w.writeheader(); w.writerows(rows)
print('PASS: independent finite-state fixtures')
for r in rows: print(r)
