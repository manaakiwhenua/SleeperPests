from pathlib import Path
import re
root=Path('/mnt/data/INApest_info_persistence_work/modified')
files=[
'INApest.R','INApestParallel.R','INApestParallelInAScene.R','INApestMetaParallel.r',
'INApestMetaMultipleLandUse.r','INApestMetaParallelMultipleLandUse.r',
'INApestMetaTransitionMatrix.r','INApestMetaTransitionMatrixParallel.r']
problems=[]
for name in files:
    s=(root/name).read_text()
    for needle in [
        'InfoPersistenceSteps = NA',
        'InfoPersistenceSteps matrix must have dimensions nodes x Ntimesteps',
        'InfoPersistenceSteps values must be non-negative whole numbers or NA',
        'Programmed stopping takes priority where InfoPersistenceSteps is not NA',
        'LastKnownPresence',
        'ProgrammedInfoNodes',
        'is.na(NodeInfoPersistenceSteps) & NodeInfoRetentionProb < 1',
    ]:
        if needle not in s: problems.append(f'{name}: missing {needle}')

for name in ['INApest.R','INApestParallel.R','INApestParallelInAScene.R']:
    s=(root/name).read_text()
    if 'HaveInfo == 1 & Invaded == 1' not in s:
        problems.append(f'{name}: occupancy known-presence reset missing')

for name in ['INApestMetaParallel.r','INApestMetaMultipleLandUse.r','INApestMetaParallelMultipleLandUse.r']:
    s=(root/name).read_text()
    if 'ManagementKillProb' not in s or 'ConditionalManagementMortality' not in s:
        problems.append(f'{name}: management-kill attribution missing')
    if 'InitialKnownPresence = which(InitDetection == 1)' not in s:
        problems.append(f'{name}: initial direct-detection clock missing')

for name in ['INApestMetaTransitionMatrix.r','INApestMetaTransitionMatrixParallel.r']:
    s=(root/name).read_text()
    if not re.search(r'rowSums\(N\s*-\s*N0\)\s*>\s*0',s):
        problems.append(f'{name}: direct management-kill reset missing')
    if 'InitialKnownPresence' not in s or 'InitDetection == 1' not in s:
        problems.append(f'{name}: initial direct-detection clock missing')

serial=(root/'INApestMetaTransitionMatrix.r').read_text()
if re.search(r'DetectionProbPerStage\s*<-.*NodeDetectionProb\[,\s*,\s*t\]',serial):
    problems.append('INApestMetaTransitionMatrix.r: stale timestep-loop detection indexing remains')
if not re.search(r'NodeDetectionProb\[,\s*,\s*timestep\]',serial):
    problems.append('INApestMetaTransitionMatrix.r: current-timestep detection indexing missing')

mlu=(root/'INApestMetaMultipleLandUse.r').read_text()
if 'SEAM*apply(Detected,1,max)' not in mlu:
    problems.append('INApestMetaMultipleLandUse.r: node-level SEAM detection collapse missing')

ina=(root/'INApestParallelInAScene.R').read_text()
for needle in ['UseExternalInfoRefresh', 'seam = if(UseExternalInfoRefresh == T) 0 else SEAM*Detected', 'SEAM*Detected']:
    if needle not in ina: problems.append(f'INApestParallelInAScene.R: missing {needle}')

if problems:
    print('SOURCE INVARIANT FAIL')
    print('\n'.join(problems))
    raise SystemExit(1)
print('SOURCE INVARIANT PASS')
for name in files: print('PASS',name)
