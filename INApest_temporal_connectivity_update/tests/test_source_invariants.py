#!/usr/bin/env python3
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1] / 'modified'
FILES = [
    'INApest.R',
    'INApestParallel.R',
    'INApestParallelInAScene.R',
    'INApestMetaParallel.r',
    'INApestMetaMultipleLandUse.r',
    'INApestMetaParallelMultipleLandUse.r',
    'INApestMetaTransitionMatrix.r',
    'INApestMetaTransitionMatrixParallel.r',
]


def delimiter_check(path: Path):
    text = path.read_text()
    stack = []
    pairs = {')': '(', ']': '[', '}': '{'}
    openers = set(pairs.values())
    quote = None
    escaped = False
    comment = False
    line = 1
    for ch in text:
        if ch == '\n':
            line += 1
            comment = False
            if quote is None:
                escaped = False
            continue
        if comment:
            continue
        if quote is not None:
            if escaped:
                escaped = False
                continue
            if ch == '\\' and quote in ('\"', "'"):
                escaped = True
                continue
            if ch == quote:
                quote = None
            continue
        if ch == '#':
            comment = True
            continue
        if ch in ('\"', "'", '`'):
            quote = ch
            continue
        if ch in openers:
            stack.append((ch, line))
        elif ch in pairs:
            if not stack or stack[-1][0] != pairs[ch]:
                raise AssertionError(f'{path.name}:{line}: unmatched {ch}')
            stack.pop()
    if quote is not None:
        raise AssertionError(f'{path.name}: unterminated quote {quote}')
    if stack:
        ch, ln = stack[-1]
        raise AssertionError(f'{path.name}:{ln}: unclosed {ch}')


for name in FILES:
    path = ROOT / name
    s = path.read_text()
    delimiter_check(path)
    assert '3D array (nodes x nodes x timesteps)' in s, name
    assert 'SDDprob 3D array must have dimensions nodes x nodes x Ntimesteps' in s, name
    assert 'LDDprob 3D array must have dimensions nodes x nodes x Ntimesteps' in s, name
    assert re.search(r'SDDprob\[\s*,\s*,\s*timestep\s*\]', s), f'{name}: missing SDD timestep slice'
    assert re.search(r'LDDprob\[\s*,\s*,\s*timestep\s*\]', s), f'{name}: missing LDD timestep slice'

# Occupancy functions must recompute DispProb from selected matrices.
for name in ['INApest.R', 'INApestParallel.R', 'INApestParallelInAScene.R']:
    s = (ROOT/name).read_text()
    assert 'DispProb = NodeSDDprob' in s
    assert '1-(1-NodeSDDprob)*(1-NodeLDDprob)' in s or '1-(1-NodeSDDprob)*(1-NodeLDDprob)' in s.replace(' ', '')
    assert 'BPAM = sweep(DispProb,2,EnvEstabProb' in s

# Meta functions must pass selected matrices into existing local dynamics or
# use them directly in the established propagule-dispersal calculations.
for name in ['INApestMetaMultipleLandUse.r', 'INApestMetaTransitionMatrix.r']:
    s = (ROOT/name).read_text()
    assert re.search(r'sddprob\s*=\s*NodeSDDprob', s), name
    assert re.search(r'lddprob\s*=\s*NodeLDDprob', s), name

for name in ['INApestMetaParallel.r', 'INApestMetaParallelMultipleLandUse.r']:
    s = (ROOT/name).read_text()
    assert 'rowSums(NodeSDDprob)' in s, name
    assert '%*% NodeSDDprob' in s, name
    assert 'rowSums(NodeLDDprob)' in s, name
    assert '%*% NodeLDDprob' in s, name

s = (ROOT/'INApestMetaParallel.r').read_text()
assert 'sweep(NodeSDDprob,1,Pout' in s
assert 'sweep(NodeLDDprob,1,Qout' in s

s = (ROOT/'INApestMetaTransitionMatrixParallel.r').read_text()
assert re.search(r'sddprob\s*=\s*NodeSDDprob', s)
assert re.search(r'lddprob\s*=\s*NodeLDDprob', s)

# Guarded initialization migration: all older sibling families should now have
# the no-input case explicitly represented instead of attempting sample(..., NA).
for name in ['INApest.R','INApestParallel.R','INApestParallelInAScene.R','INApestMetaParallel.r','INApestMetaMultipleLandUse.r','INApestMetaParallelMultipleLandUse.r']:
    s = (ROOT/name).read_text()
    assert 'Infested = integer(0)' in s, name
    assert 'risk = NULL' in s, name

# The serial transition mortality timestep must match the already-correct parallel implementation.
s = (ROOT/'INApestMetaTransitionMatrix.r').read_text()
assert 'NodeMortalityProb[, s,timestep] * Managing' in s
assert 'NodeMortalityProb[, s,t] * Managing' not in s

print(f'PASS: source invariants and delimiter balance for {len(FILES)} modified R files')
