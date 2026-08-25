from pathlib import Path
import re, hashlib
root=Path('/mnt/data/INApest_definitive_2026-08-25/src')
expected={
'INApest.R':['INApest <- function','Pathogen'],
'INApestParallel.R':['INApestParallel <- function'],
'INApestPathogen.R':['INApestPathogen <- function'],
'INApestMeta.r':['INApestMeta','Pathogen'],
'INApestMetaParallel.r':['INApestMetaParallel','Pathogen'],
'INApestMetaMultipleLandUse.r':['INApestMetaMultipleLandUse','Pathogen'],
'INApestMetaParallelMultipleLandUse.r':['INApestMetaParallelMultipleLandUse','Pathogen'],
'INApestMetaTransitionMatrix.r':['INApestMetaTransitionMatrix','Pathogen','InitialPathogenState'],
'INApestMetaTransitionMatrixParallel.r':['INApestMetaTransitionMatrixParallel','Pathogen','InitialPathogenState'],
'INApestPathogenTransitionMatrix.R':['local.dynamics.transition.matrix.pathogen'],
'INApestPointPathogen.R':['INApestPointPathogenInteraction'],
'INApestMetaPoint.R':['INApestMetaPoint <- function','Pathogen','PathogenEvents'],
'INApestMetaPointParallel.R':['INApestMetaPointParallel <- function','PathogenEvents'],
'INApestPointTransitionMatrix.R':['INApestPointTransitionMatrix <- function','Pathogen','PathogenEvents'],
'INApestPointTransitionMatrixParallel.R':['INApestPointTransitionMatrixParallel <- function','PathogenEvents'],
}

def scan(text):
    stack=[]; pairs={')':'(',']':'[','}':'{'}; opens=set(pairs.values())
    i=0; line=1; state='code'; quote=None
    while i<len(text):
        c=text[i]; n=text[i+1] if i+1<len(text) else ''
        if c=='\n': line+=1
        if state=='comment':
            if c=='\n': state='code'
        elif state=='string':
            if c=='\\': i+=1
            elif c==quote: state='code'; quote=None
        else:
            if c=='#': state='comment'
            elif c in ('"',"'",'`'): state='string'; quote=c
            elif c in opens: stack.append((c,line))
            elif c in pairs:
                if not stack or stack[-1][0]!=pairs[c]: return False,f'mismatch {c} line {line}, stack={stack[-3:]}'
                stack.pop()
        i+=1
    if state=='string': return False,'unterminated string'
    return (not stack, 'ok' if not stack else f'unclosed {stack[-5:]}')

rows=[]
for fn,needles in expected.items():
    p=root/fn
    if not p.exists(): rows.append((fn,'MISSING','','')); continue
    s=p.read_text(errors='replace')
    ok,msg=scan(s)
    missing=[x for x in needles if x not in s]
    sha=hashlib.sha256(p.read_bytes()).hexdigest()
    rows.append((fn,'PASS' if ok and not missing else 'FAIL', msg, ', '.join(missing),len(s.splitlines()),sha))
for r in rows: print('\t'.join(map(str,r)))
if any(r[1]!='PASS' for r in rows): raise SystemExit(1)
