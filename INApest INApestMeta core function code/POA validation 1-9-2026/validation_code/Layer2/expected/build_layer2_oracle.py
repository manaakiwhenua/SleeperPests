import csv, math, os
out=[]
def add(test_id, quantity, value, formula, note=''):
    out.append(dict(test_id=test_id, quantity=quantity, exact_value=float(value), formula=formula, note=note))

# L2-01
prior=.8; q=.2
post=prior/(prior+(1-prior)*q)
add('L2-01','posterior_poa',post,'0.8 / (0.8 + 0.2*0.2)','single zero-detection round')
add('L2-01','background_sse',.8,'1 - 0.2')

# L2-02 repeated rounds
for k in (1,2,3):
    po=prior/(prior+(1-prior)*(q**k))
    add('L2-02',f'posterior_round_{k}',po,f'0.8 / (0.8 + 0.2*0.2^{k})')

# L2-03 reinvasion
post1=prior/(prior+(1-prior)*q)
r=.1
prior2=post1*(1-r)
post2=prior2/(prior2+(1-prior2)*q)
add('L2-03','posterior_round_1',post1,'0.8 / (0.8 + 0.2*0.2)')
add('L2-03','prior_round_2',prior2,'posterior_round_1 * (1 - 0.1)')
add('L2-03','posterior_round_2',post2,'prior_round_2 / (prior_round_2 + (1-prior_round_2)*0.2)')

# L2-04 extinction + reinvasion
prior=.7; q=.4; r=.2; e=.3
post1=prior/(prior+(1-prior)*q)
prior2=post1*(1-r)+(1-post1)*e
post2=prior2/(prior2+(1-prior2)*q)
add('L2-04','posterior_round_1',post1,'0.7 / (0.7 + 0.3*0.4)')
add('L2-04','prior_round_2',prior2,'post1*(1-0.2) + (1-post1)*0.3')
add('L2-04','posterior_round_2',post2,'prior2 / (prior2 + (1-prior2)*0.4)')

# L2-05 two surveillance pathways
prior=.5
q_bg=.4
q_info_mean=.75
q_comb=.5*(.4*1)+.5*(.4*.5)
post=prior/(prior+(1-prior)*q_comb)
add('L2-05','background_sse',1-q_bg,'1 - 0.4')
add('L2-05','info_triggered_sse',1-q_info_mean,'1 - mean(1,0.5)')
add('L2-05','combined_sse',1-q_comb,'1 - mean(0.4,0.2)')
add('L2-05','posterior_poa',post,'0.5 / (0.5 + 0.5*0.3)')

# L2-06 management/history benchmark
prior=.5; q=.5
post1=prior/(prior+(1-prior)*q)
prior2=post1
post2=prior2/(prior2+(1-prior2)*q)
naive_prior2=post1+(1-post1)*.5
naive_post2=naive_prior2/(naive_prior2+(1-naive_prior2)*q)
add('L2-06','posterior_round_1',post1,'0.5 / (0.5 + 0.5*0.5)')
add('L2-06','prior_round_2_compatible_history',prior2,'posterior_round_1; detected-present histories are incompatible')
add('L2-06','posterior_round_2',post2,'(2/3) / ((2/3) + (1/3)*0.5)')
add('L2-06','naive_prior_round_2',naive_prior2,'post1 + (1-post1)*0.5','wrong comparator: propagates detected histories')
add('L2-06','naive_posterior_round_2',naive_post2,'(5/6) / ((5/6)+(1/6)*0.5)','wrong comparator: should NOT match PoA core')

# L2-07 high PoA exact target and discretisation sequence
alpha=math.sqrt(2)/2
qbar=alpha*.8+(1-alpha)*(.8**20)
prior=.999
post=prior/(prior+(1-prior)*qbar)
add('L2-07','alpha_n1',alpha,'sqrt(2)/2')
add('L2-07','mean_q_present',qbar,'alpha*0.8 + (1-alpha)*0.8^20')
add('L2-07','exact_posterior_poa',post,'0.999 / (0.999 + 0.001*mean_q_present)')
for n in (1000,10000,50000,200000):
    m=max(2,int(math.sqrt(n)))
    k=round(m*alpha)
    ah=k/m
    qa=ah*.8+(1-ah)*(.8**20)
    pa=prior/(prior+(1-prior)*qa)
    add('L2-07',f'approx_posterior_N{n}',pa,f'alpha_hat={k}/{m}; prior=0.999',f'absolute error={abs(pa-post):.12g}')

# L2-08 support safeguard has expected error, no numeric oracle
# L2-09 simulation-mode oracle
prior=.8; q=.2
post=prior/(prior+(1-prior)*q)
add('L2-09','exact_posterior_poa',post,'0.8 / (0.8 + 0.2*0.2)','simulation-mode target')

# L2-10 positive binary detection
add('L2-10','posterior_poa',0.0,'P(detection|absent)=0, therefore posterior absence=0')

path=os.path.join(os.path.dirname(__file__),'layer2_exact_oracle.csv')
with open(path,'w',newline='',encoding='utf-8') as f:
    w=csv.DictWriter(f,fieldnames=['test_id','quantity','exact_value','formula','note'])
    w.writeheader(); w.writerows(out)
print(path)
for row in out:
    print(f"{row['test_id']} {row['quantity']}: {row['exact_value']:.12g}")
