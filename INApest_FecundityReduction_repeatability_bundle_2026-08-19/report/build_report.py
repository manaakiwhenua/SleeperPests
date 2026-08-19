#!/usr/bin/env python3
from pathlib import Path
import base64, csv, html, hashlib

ROOT = Path(__file__).resolve().parents[1]
ARCH = ROOT / 'results' / 'archived'
LATEST = ROOT / 'results' / 'latest_rerun'
FIG = ROOT / 'figures'
OUT = ROOT / 'report' / 'FecundityReduction_repeatability_report.html'

def read_csv(path):
    with path.open(newline='', encoding='utf-8-sig') as f:
        return list(csv.DictReader(f))

def esc(x): return html.escape(str(x))

def fmt(x, dp=3):
    if x is None or x == '': return ''
    try:
        v=float(x)
        if v != v: return 'NA'
        s=f'{v:.{dp}f}'.rstrip('0').rstrip('.')
        return s
    except Exception:
        return esc(x)

def table(rows, cols, labels=None, classes=''):
    labels = labels or cols
    h='<div class="table-wrap '+classes+'"><table><thead><tr>' + ''.join(f'<th>{esc(l)}</th>' for l in labels) + '</tr></thead><tbody>'
    for r in rows:
        h += '<tr>' + ''.join(f'<td>{r.get(c, "")}</td>' for c in cols) + '</tr>'
    h += '</tbody></table></div>'
    return h

def image_data(path):
    return 'data:image/png;base64,' + base64.b64encode(path.read_bytes()).decode('ascii')

def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()

# Evidence tables
equiv_raw = read_csv(ARCH/'equivalence_results.csv')
equiv=[]
for r in equiv_raw:
    equiv.append({'model':esc(r['model']), 'output':esc(r['output'].replace('LargeOut.rds','')), 'equal':'PASS' if r['all_equal'].upper()=='TRUE' else 'FAIL', 'diff':fmt(r['max_abs_difference'])})
smoke_raw=read_csv(ARCH/'parameterisation_smoke_results.csv')
smoke=[{'model':esc(r['model'].replace('INApestMetaMultipleLandUse','Meta MLU').replace('INApestMetaTransitionMatrix','Meta transition matrix').replace('INApestMeta','Meta')), 'shape':esc(r['parameterisation'].replace('_',' ')), 'status':esc(r['status'])} for r in smoke_raw]
shapeeq_raw=read_csv(ARCH/'parameterisation_parallel_equivalence.csv')
shapeeq=[{'model':esc(r['model'].replace('INApestMetaMultipleLandUse','Meta MLU').replace('INApestMetaTransitionMatrix','Meta transition matrix')), 'shape':esc(r['parameterisation'].replace('_',' ')), 'equal':'PASS' if r['all_equal'].upper()=='TRUE' else 'FAIL', 'diff':fmt(r['max_abs_difference'])} for r in shapeeq_raw]
meta_current_raw=read_csv(LATEST/'current_meta_serial_parallel.csv')
meta_current=[{'shape':esc(r['parameterisation'].replace('_',' ')), 'equal':'PASS' if r['all_equal'].upper()=='TRUE' else 'FAIL', 'diff':fmt(r['max_abs_difference'])} for r in meta_current_raw]
contract_raw=read_csv(LATEST/'current_local_dynamics_contract.csv')
contract=[{'check':esc(r['check'].replace('_',' ')), 'pass':'PASS' if r['pass'].upper()=='TRUE' else 'FAIL', 'message':esc(r['message'])} for r in contract_raw]
compare_raw=read_csv(LATEST/'archived_vs_latest_comparison.csv')
compare=[{'file':esc(r['file']), 'same':'YES' if r['byte_identical'].lower()=='true' else 'NO'} for r in compare_raw]

mort_raw=read_csv(ARCH/'mortality_by_fecundity_results.csv')
mort=[]
for r in mort_raw:
    mort.append({'mort':fmt(r['mortality'],2),'fr':fmt(r['fecundity_reduction'],2),'mean':fmt(r['final_population_mean'],2),'sd':fmt(r['final_population_sd'],2),'q':f"{fmt(r['q05'],2)}–{fmt(r['q95'],2)}",'ratio':fmt(float(r['ratio_to_no_fecundity_reduction'])*100,1)+'%'})
lu_raw=read_csv(ARCH/'landuse_mortality_fecundity_results.csv')
lu03=[]
for r in lu_raw:
    if abs(float(r['mortality'])-.3)<1e-9:
        lu03.append({'pattern':esc(r['pattern'].replace('_',' ')),'fr':f"{fmt(r['FR_LU1'],1)} / {fmt(r['FR_LU2'],1)} / {fmt(r['FR_LU3'],1)}",'total':fmt(r['final_total_mean'],1),'lu1':fmt(r['node1_LU1_mean'],1),'lu2':fmt(r['node2_LU2_mean'],1),'lu3':fmt(r['node3_LU3_mean'],1),'mixed':fmt(r['node4_mixed_mean'],1)})

manifest = read_csv(ROOT/'docs'/'source_snapshot_manifest.csv')
latest_meta_files=[r for r in manifest if r['snapshot']=='latest' and r['file'].startswith('INApestMeta')]
manifest_rows=[{'file':esc(r['file']),'bytes':esc(r['bytes']),'sha':esc(r['sha256'][:16]+'…')} for r in latest_meta_files]

# Work example
abund=[40,30,30]; fr=[.1,.5,.9]; contrib=[n*(1-f) for n,f in zip(abund,fr)]; tot=sum(contrib); shares=[100*c/tot for c in contrib]

img1=image_data(FIG/'mortality_fecundity_effect.png')
img2=image_data(FIG/'landuse_fecundity_gradient_effect.png')

CSS=r'''
:root{--ink:#20323f;--muted:#60717c;--line:#d6e2e5;--page:#eef4f5;--teal:#0b7285;--teal2:#075766;--blue:#2878c8;--purple:#6d5ac7;--orange:#d97706;--green:#2b8a5a;--coral:#c85f4a;--soft:#f7fafb;--amber:#fff7e6;--red:#b63b32}
*{box-sizing:border-box}html{background:var(--page);scroll-behavior:smooth}body{font-family:system-ui,-apple-system,'Segoe UI',Arial,sans-serif;color:var(--ink);line-height:1.57;max-width:1500px;margin:auto;padding:26px 42px 82px;background:#fff;box-shadow:0 0 28px rgba(30,70,80,.08)}
header{background:linear-gradient(120deg,#075766,#0b7285 55%,#2f87a2);color:#fff;padding:32px 36px;border-radius:17px}header h1{margin:.05rem 0 .45rem;font-size:2.32rem;line-height:1.1;color:#fff;border:0;padding:0}header p{margin:.3rem 0;color:#e3f6f8;max-width:118ch}header code{background:rgba(255,255,255,.14);color:#fff;padding:.05rem .22rem}
h1,h2,h3,h4{line-height:1.22}h1{color:var(--teal2);border-bottom:2px solid #d9e8eb;padding-bottom:.32rem;margin-top:2.55rem}h2{color:#315e75;margin-top:1.9rem}h3{color:#3f5f6b}.lede{font-size:1.08rem;color:#324c57;max-width:118ch}
.toc{background:#f7fbfc;border:1px solid var(--line);padding:14px 18px;border-radius:12px;margin:15px 0}.toc a{color:#136f91;text-decoration:none;margin-right:.72rem;display:inline-block}.toc a:hover{text-decoration:underline}
.callout{border-left:5px solid var(--teal);background:#edf7f8;padding:13px 16px;border-radius:0 9px 9px 0;margin:15px 0}.callout.good{border-left-color:var(--green);background:#edf8f1}.callout.warn{border-left-color:var(--coral);background:#fff0ec}.callout.amber{border-left-color:var(--orange);background:#fff7e6}
.score-grid{display:grid;grid-template-columns:repeat(5,minmax(160px,1fr));gap:10px;margin:16px 0}.score{border:1px solid var(--line);border-radius:12px;padding:14px;background:#fff}.score .n{font-size:1.48rem;font-weight:850;color:var(--teal2)}.score .k{font-size:.81rem;color:var(--muted)}
.cards{display:grid;grid-template-columns:repeat(3,1fr);gap:12px;margin:15px 0}.card{border:1px solid var(--line);border-radius:12px;padding:14px 16px;background:#fff;min-width:0}.card .head{font-weight:850;color:#315e75;margin-bottom:.35rem}.card p{margin:.35rem 0}.card.mort{background:#fff4ee;border-top:5px solid var(--coral)}.card.fec{background:#eef7f8;border-top:5px solid var(--teal)}.card.spread{background:#f3f0ff;border-top:5px solid var(--purple)}
.badge{display:inline-block;padding:.13rem .5rem;border-radius:999px;background:#eef5fc;border:1px solid #bdd6e9;color:#226399;font-size:.76rem;font-weight:780;margin:.08rem .2rem .08rem 0}
.flow{display:grid;grid-template-columns:repeat(7,minmax(105px,1fr));align-items:stretch;gap:7px;margin:16px 0}.flow .box{border:1px solid var(--line);border-radius:11px;background:#fff;padding:11px 9px;text-align:center;min-height:108px;display:flex;flex-direction:column;justify-content:center}.flow .box strong{color:#315e75}.flow .arrow{display:flex;align-items:center;justify-content:center;color:#77929b;font-size:1.55rem;font-weight:900;position:absolute}.process{display:grid;grid-template-columns:1fr 38px 1fr 38px 1fr 38px 1fr;gap:6px;align-items:center;margin:14px 0}.process .arrow{text-align:center;font-size:1.5rem;color:#6f8790}.process .box{border:1px solid var(--line);border-radius:12px;padding:13px;text-align:center;background:#fff;min-height:105px;display:flex;flex-direction:column;justify-content:center}
.formula{font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;background:#f3f7f8;border:1px solid var(--line);border-radius:10px;padding:12px 14px;overflow:auto;margin:10px 0;font-size:.91rem}.notation{font-family:Georgia,'Times New Roman',serif;font-size:1.03rem}
.mechanism{border:1px solid var(--line);border-radius:15px;overflow:hidden;margin:18px 0;box-shadow:0 3px 12px rgba(30,70,80,.05)}.mechanism-head{display:grid;grid-template-columns:1fr 1fr 1fr}.panel{min-width:0}.phead{color:#fff;font-weight:800;padding:.52rem .72rem}.phead.blue{background:var(--blue)}.phead.orange{background:var(--orange)}.phead.green{background:var(--green)}.pbody{padding:14px 16px;background:#fff;min-height:158px}.pbody p{margin:.42rem 0}
.landuse-node{display:grid;grid-template-columns:repeat(3,1fr);gap:10px;margin:14px 0}.lu{border-radius:11px;padding:12px;border:1px solid var(--line);background:#fff}.lu:nth-child(1){border-top:5px solid var(--blue)}.lu:nth-child(2){border-top:5px solid var(--orange)}.lu:nth-child(3){border-top:5px solid var(--green)}.lu .title{font-weight:850;color:#315e75}.big{font-size:1.48rem;font-weight:850;color:var(--teal2)}
.two{display:grid;grid-template-columns:1fr 1fr;gap:13px}.implementation-card{border:1px solid var(--line);border-radius:12px;padding:14px 16px;background:#fbfdfe}.implementation-card h3{margin:.05rem 0 .5rem}
table{border-collapse:separate;border-spacing:0;width:100%;font-size:.84rem;border:1px solid var(--line);border-radius:10px;overflow:hidden}th,td{padding:.52rem .62rem;border-right:1px solid var(--line);border-bottom:1px solid var(--line);vertical-align:top}th:last-child,td:last-child{border-right:0}tr:last-child td{border-bottom:0}th{background:#176f81;color:white;text-align:left}tbody tr:nth-child(even) td{background:#f8fbfc}.table-wrap{overflow-x:auto;margin:12px 0}.caption{font-size:.88rem;color:#405963;margin:.4rem 0 .8rem;padding:.65rem .78rem;background:#f7fafb;border-left:4px solid #9bc5cc}
.figure{border:1px solid var(--line);border-radius:12px;padding:12px;background:#fff;margin:16px 0}.figure img{display:block;width:100%;max-width:1120px;height:auto;margin:auto}.figure .caption{border-left-color:var(--teal);margin-top:10px}.small{font-size:.82rem;color:var(--muted)}code{font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;background:#f1f5f6;padding:.08rem .28rem;border-radius:4px}
details{border:1px solid var(--line);border-radius:10px;margin:11px 0;background:#fbfdfe}summary{cursor:pointer;padding:.7rem .85rem;font-weight:760;color:#315e75}details>div{padding:0 .85rem .85rem}.check{font-weight:850;color:var(--green)}.pending{font-weight:850;color:var(--orange)}
.popoutable{cursor:zoom-in}.overlay{position:fixed;inset:0;background:rgba(15,27,33,.78);z-index:9999;display:flex;align-items:center;justify-content:center;padding:22px}.overlaybox{background:#fff;max-width:95vw;max-height:94vh;overflow:auto;border-radius:14px;padding:18px;box-shadow:0 18px 60px rgba(0,0,0,.35);position:relative}.overlaybox img{width:min(1250px,88vw);height:auto}.close{position:sticky;top:0;float:right;border:0;background:#20323f;color:#fff;border-radius:999px;width:34px;height:34px;font-size:18px;cursor:pointer;z-index:2}
.filetree{font-family:ui-monospace,SFMono-Regular,Menlo,monospace;background:#f6f8f9;border:1px solid var(--line);border-radius:11px;padding:14px;white-space:pre-wrap;font-size:.84rem}.footer{margin-top:50px;border-top:1px solid var(--line);padding-top:14px;color:var(--muted);font-size:.84rem}
@media(max-width:1150px){.score-grid{grid-template-columns:repeat(3,1fr)}.process{grid-template-columns:1fr}.process .arrow{transform:rotate(90deg)}}@media(max-width:850px){body{padding:14px}.cards,.two,.mechanism-head,.landuse-node{grid-template-columns:1fr}.score-grid{grid-template-columns:1fr 1fr}}@media print{html{background:white}body{max-width:none;box-shadow:none;padding:10mm}header{break-inside:avoid}.popoutable{cursor:default}}
'''

html_doc=f'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>INApest FecundityReduction — implementation, validation and repeatability</title><style>{CSS}</style></head><body>
<header><h1>INApest <code>FecundityReduction</code></h1><p><strong>Model implementation, validation, demonstration and full repeatability record</strong></p><p>Audience: expert population, invasion and spatial modellers with little or no prior knowledge of INApest · 19 August 2026</p></header>
<div class="toc"><strong>Contents</strong><br><a href="#summary">Executive summary</a><a href="#primer">INApest primer</a><a href="#purpose">Purpose</a><a href="#ordinary">Ordinary Meta</a><a href="#mlu">Multiple land use</a><a href="#tm">Transition matrix</a><a href="#contract">Parameter contract</a><a href="#engineering">Engineering</a><a href="#repeatability">Repeatability</a><a href="#validation">Validation</a><a href="#trials">Demonstration simulations</a><a href="#reproduce">How to reproduce</a><a href="#limits">Limits</a></div>

<h1 id="summary">Executive summary</h1>
<p class="lede">The INApestMeta model family already represented management that removes individuals and management that reduces onward dispersal. It lacked a separate mechanism for management that leaves individuals alive but lowers their reproductive output. The <code>FecundityReduction</code> upgrade adds that mechanism while preserving the previous behaviour at the default value zero.</p>
<div class="score-grid"><div class="score"><div class="n">3</div><div class="k">Meta model structures support fecundity management</div></div><div class="score"><div class="n">14 / 14</div><div class="k">parameterisation and edge-case smoke checks passed</div></div><div class="score"><div class="n">12 / 12</div><div class="k">saved serial/parallel outputs matched exactly on the executable one-core path</div></div><div class="score"><div class="n">9 / 9</div><div class="k">MLU and transition-matrix input shapes matched across serial/parallel-function paths</div></div><div class="score"><div class="n">6 / 6</div><div class="k">archived evidence CSVs reproduced byte-for-byte with the latest source defaults</div></div></div>
<div class="callout good"><strong>Mechanistic definition.</strong> Management mortality is realised first. <code>FecundityReduction</code> then reduces the per-capita reproductive output of the survivors, and only in units that are currently managed. It is therefore not a second mortality term and is not a global rescaling of baseline fecundity.</div>
<div class="callout"><strong>Reproducibility result.</strong> The full original fecundity validation suite was rerun against the current source files—which now also contain restored custom <code>LocalDynamics</code> hooks and optional transition movement. With those later features at their defaults, every archived CSV evidence file was reproduced byte-for-byte.</div>
<div class="callout amber"><strong>Parallel caveat.</strong> Exact serial/parallel equality was demonstrated on the parallel functions' one-core fallback, where the RNG stream is shared. A genuine multi-worker PSOCK run intentionally uses independent worker RNG streams, so exact realisation-by-realisation equality is neither expected nor the correct criterion. A separate native-R Monte Carlo equivalence test is included and must be run on a machine where PSOCK workers can actually be created.</div>

<h1 id="primer">1. INApest model family: the minimum needed to read this report</h1>
<p>INApest is a family of discrete-time stochastic invasion-response models defined on a network of spatial nodes. A node can represent a farm, management unit, grid cell or other spatial unit. Replicate simulations (<code>Nperm</code>) propagate biological state and response state through <code>Ntimesteps</code>. Connectivity matrices describe short-distance dispersal (SDD), optional long-distance dispersal (LDD), and—separately—information exchange where supplied.</p>
<div class="cards"><div class="card"><div class="head">Ordinary Meta</div><p>One abundance value per node. Local survival, management mortality, reproduction, dispersal, recruitment and carrying capacity update node abundance.</p><span class="badge">state: node</span></div><div class="card"><div class="head">Multiple-land-use Meta</div><p>Each node contains an explicit abundance vector over land uses. Management effects can differ among those land uses before the node's propagules are pooled.</p><span class="badge">state: node × land use</span></div><div class="card"><div class="head">Transition-matrix Meta</div><p>Each node contains a stage vector. A stage transition matrix controls fecundity, stasis/progression and survival, with spatial dispersal and carrying capacities layered around that structure.</p><span class="badge">state: node × stage</span></div></div>
<h2>Where the new mechanism sits in a timestep</h2>
<div class="process"><div class="box"><strong>Existing population</strong><span>node / land-use / stage state</span></div><div class="arrow">→</div><div class="box"><strong>Detection &amp; management</strong><span>information controls whether management is adopted</span></div><div class="arrow">→</div><div class="box"><strong>Survival &amp; mortality</strong><span>natural survival and management mortality create surviving state <em>N</em><sup>0</sup></span></div><div class="arrow">→</div><div class="box"><strong>Fecundity</strong><span><code>FecundityReduction</code> modifies reproductive output of managed survivors</span></div></div>
<div class="process"><div class="box"><strong>Propagule draw</strong><span>Poisson production from effective reproductive contribution</span></div><div class="arrow">→</div><div class="box"><strong>Dispersal</strong><span>SDD and optional LDD; spread reduction acts on the management-sensitive LDD pathway</span></div><div class="arrow">→</div><div class="box"><strong>Establishment &amp; capacity</strong><span>recruitment probability and available capacity determine additions</span></div><div class="arrow">→</div><div class="box"><strong>Next state</strong><span>saved and passed to the next timestep</span></div></div>
<p class="small">The transition-matrix model contains additional stage-transition steps. The diagram isolates the location of fecundity management in the update sequence rather than reproducing every line of each model.</p>

<h1 id="purpose">2. Why add <code>FecundityReduction</code>?</h1>
<div class="cards"><div class="card mort"><div class="head">Mortality</div><p>Changes the number of individuals that survive management. In ordinary Meta, the survivor draw combines environmental survival and management mortality.</p><div class="formula">N⁰ ~ Binomial(N, Survival × (1 − MortalityProb × Managing))</div></div><div class="card fec"><div class="head">Fecundity reduction</div><p>Changes the reproductive contribution of individuals that survived. It is applied only when the relevant node or land-use unit is managed.</p><div class="formula">effective reproduction = N⁰ × (1 − FecundityReduction × Managing)</div></div><div class="card spread"><div class="head">Spread reduction</div><p>Acts later, on movement of propagules. In the current Meta defaults it reduces the management-sensitive LDD output, rather than changing how many propagules were produced.</p><div class="formula">LDD output ∝ propagules × LDDrate × (1 − SpreadReduction × Managing)</div></div></div>
<p>Keeping these pathways separate matters whenever management is sublethal, sterilising, growth-regulating, host-modifying, mating-disrupting or otherwise capable of lowering reproduction without killing the treated individuals. Collapsing such an effect into mortality changes abundance immediately; collapsing it into baseline <code>PropaguleProduction</code> changes reproduction even where management is absent. <code>FecundityReduction</code> instead represents a conditional management effect on per-capita reproduction.</p>

<h1 id="ordinary">3. Ordinary INApestMeta implementation</h1>
<p>For node <span class="notation">i</span> at timestep <span class="notation">t</span>, let <span class="notation">N⁰<sub>it</sub></span> be abundance after natural survival and management mortality, <span class="notation">M<sub>it</sub></span> the binary management state, <span class="notation">f<sub>it</sub></span> the fecundity reduction, and <span class="notation">b<sub>it</sub></span> baseline propagule production per survivor.</p>
<div class="formula">Nrep_it = N⁰_it × (1 − f_it × M_it)</div><div class="formula">Propagules_it ~ Poisson( b_it × Nrep_it )</div>
<p>This is implemented in the default <code>local.dynamics</code> function as <code>EffectiveReproductivePopulation</code>, followed by the existing Poisson propagule draw. The transition from survivors to reproduction is therefore explicit: if management is absent, the multiplier is one; if <code>f = 1</code> and management is active, survivors remain in the population but contribute zero new propagules during that timestep.</p>
<div class="callout"><strong>Default compatibility.</strong> <code>FecundityReduction = 0</code> makes the multiplier exactly one. The original implementation check found zero numerical difference from the pre-upgrade MLU and transition-matrix serial models in the tested scenarios; the current-source scientific tables also reproduce the archived results exactly.</div>

<h1 id="mlu">4. Multiple-land-use implementation: keep reproductive refuges explicit</h1>
<p>The multiple-land-use model is the most important place not to collapse the new parameter too early. A node contains separate abundance components for each land use. Management status and fecundity suppression are therefore applied to each component <strong>before</strong> those contributions are summed to obtain node-level propagule production.</p>
<div class="formula">C_iℓt = N⁰_iℓt × (1 − f_iℓt × M_iℓt)</div><div class="formula">Propagules_it ~ Poisson( b_it × Σℓ C_iℓt )</div>
<p>Here <span class="notation">C<sub>iℓt</sub></span> is the surviving reproductive contribution of land use <span class="notation">ℓ</span>. This is analogous to the refuge mechanism used elsewhere in the modelling work: the weakly affected component is kept separate long enough for the relevant process to act, rather than averaging its protection away at node level.</p>
<h2>Worked reproductive-refuge example</h2>
<div class="landuse-node"><div class="lu"><div class="title">Land use 1</div><p>Post-mortality abundance: <strong>{abund[0]}</strong><br>Fecundity reduction: <strong>{fr[0]:.1f}</strong></p><div class="big">{contrib[0]:.0f}</div><div class="small">effective reproductive contribution ({shares[0]:.0f}% of node total)</div></div><div class="lu"><div class="title">Land use 2</div><p>Post-mortality abundance: <strong>{abund[1]}</strong><br>Fecundity reduction: <strong>{fr[1]:.1f}</strong></p><div class="big">{contrib[1]:.0f}</div><div class="small">effective reproductive contribution ({shares[1]:.0f}% of node total)</div></div><div class="lu"><div class="title">Land use 3</div><p>Post-mortality abundance: <strong>{abund[2]}</strong><br>Fecundity reduction: <strong>{fr[2]:.1f}</strong></p><div class="big">{contrib[2]:.0f}</div><div class="small">effective reproductive contribution ({shares[2]:.0f}% of node total)</div></div></div>
<p>The node contains 100 survivors, but effective reproductive abundance is only <strong>{tot:.0f}</strong>. Land use 1 holds 40% of survivors yet supplies about <strong>{shares[0]:.0f}%</strong> of effective reproduction because its fecundity is weakly suppressed. At baseline propagule production 2, the node's Poisson mean is 108 rather than 200. This is the sense in which a land use can function as a <strong>reproductive refuge</strong>.</p>
<div class="two"><div class="implementation-card"><h3>SDD</h3><p>The land-use contributions are summed to generate one node-level propagule pool. The natural/SDD component is then dispersed using the node-level SDD matrix. Land-use identity is not carried as a separate SDD network.</p></div><div class="implementation-card"><h3>LDD</h3><p>Before the LDD component is reduced by land-use-specific <code>SpreadReduction</code>, the model weights the node's propagules by each land use's share of <code>FecundityContrib</code>. Thus land-use identity is retained where it is needed to combine fecundity and spread management.</p></div></div>
<p>Recruitment is subsequently allocated among land uses according to unoccupied land-use capacity. The crucial design choice is that <code>FecundityReduction</code> affects the land-use reproductive contribution before node pooling, not after it.</p>

<h1 id="tm">5. Transition-matrix implementation</h1>
<p>In the stage-structured model, fecundity is encoded in the first row of the transition matrix. For stages <span class="notation">s = 2,…,S</span>, the first-row coefficient <span class="notation">a<sub>1s</sub></span> is treated as per-capita fecundity of stage <span class="notation">s</span>. The supplied transition matrix is not modified. Instead, its fecundity coefficients are multiplied by the management-specific stage multipliers when the expected propagule count is calculated.</p>
<div class="formula">a′_1s,it = a_1s,it × (1 − f_ist × M_it)</div><div class="formula">λ_it = Σ(s=2…S) a′_1s,it × N⁰_ist&nbsp;&nbsp;&nbsp;&nbsp;; Propagules_it ~ Poisson(λ_it)</div>
<p>This permits a treatment to suppress reproduction strongly in one reproductive stage but weakly in another. A node vector applies the same reduction to all stages within each node; a stage vector applies the same stage profile across nodes; the full array can vary by node, stage and timestep.</p>
<div class="callout amber"><strong>Stage 1 convention.</strong> A length-<code>Nstages</code> vector is accepted so indexing aligns with biological stage numbering, but element 1 is unused in the fecundity calculation because the active first-row fecundity coefficients correspond to stages <code>2:Nstages</code>. A dedicated test confirmed that changing only stage 1 reduction changes the output by exactly zero.</div>
<div class="callout"><strong>Ambiguity protection.</strong> If the number of nodes equals the number of stages, an unnamed vector of that shared length could mean either nodes or stages. The function rejects that case rather than guessing; users can supply a dimensioned node × timestep matrix or node × stage × timestep array.</div>
<p>The current transition-matrix source also supports optional movement during stage progression. That feature was added after the fecundity work and is not used in these fecundity demonstrations; its default-off setting was present in the latest-source rerun that reproduced the archived fecundity evidence exactly.</p>

<h1 id="contract">6. User-facing parameter contract</h1>
{table([
 {'model':'Ordinary Meta','forms':'scalar; node vector','time':'node × timestep matrix','resolution':'current node value at timestep','note':'no stage/land-use dimension'},
 {'model':'Multiple-land-use Meta','forms':'scalar; land-use vector; node × land-use matrix','time':'node × land-use × timestep array','resolution':'current node × land-use slice','note':'land-use contribution reduced before node pooling'},
 {'model':'Transition-matrix Meta','forms':'scalar; node vector; stage vector','time':'node × timestep matrix; node × stage × timestep array','resolution':'current node × stage values','note':'stage 1 unused; ambiguous node/stage vector rejected'}
], ['model','forms','time','resolution','note'], ['Model','Static forms','Time-varying form','Resolved within timestep','Special rule'])}
<p>All values must be finite and within [0,1]. There is deliberately no <code>FecundityReductionSD</code>: the parameter expresses the management effect itself rather than adding another stochastic layer to the management-efficacy input.</p>

<h1 id="engineering">7. Source architecture, parallel processing and custom local dynamics</h1>
<h2>Default dynamics remain inspectable and replaceable</h2>
<p>All current Meta serial and parallel source files expose a <code>LocalDynamics</code> argument. The current default dynamics function is defined in the same source file: <code>local.dynamics</code>, <code>local.dynamicsLU</code> or <code>local.dynamics.transition.matrix</code>. Users can therefore replace local biology without editing the main simulation driver.</p>
<p>For backward compatibility, an older custom function that lacks <code>nodefecundityreduction</code> still runs when <code>FecundityReduction = 0</code>. If fecundity management is active, the driver stops with an informative error unless the custom function accepts <code>nodefecundityreduction</code> or <code>...</code>. This prevents the management effect from being silently ignored.</p>
{table(contract,['check','pass','message'],['Current custom-LocalDynamics check','Status','Observed message'])}
<h2>Parallel implementation</h2>
<p>The active parallel source family was standardised on base-R PSOCK clusters and <code>parallel::parLapply()</code>. The selected <code>LocalDynamics</code> is forced before the worker closure is created, so the default or user-supplied function is part of the worker dependency. The transition-matrix parallel function no longer relies on a long manually maintained <code>clusterExport()</code> list that could omit new parameters.</p>
<p>The exact equivalence evidence below comes from the one-core fallback used by WebR; that path still executes the parallel function's worker body and result-combination code. The bundle's native PSOCK script is the test for the genuinely multi-process branch.</p>

<h1 id="repeatability">8. Repeatability design and provenance</h1>
<p>The bundle deliberately retains two source states:</p>
<div class="two"><div class="implementation-card"><h3>Fecundity validation snapshot</h3><p>The exact source snapshot used to generate the original archived fecundity evidence. It records the implementation state before the later transition-movement and full custom-<code>LocalDynamics</code> refinements.</p><p><code>source/fecundity_validation_snapshot/</code></p></div><div class="implementation-card"><h3>Current source snapshot</h3><p>The latest source set at bundle assembly. It includes the fecundity upgrade, harmonised parallel processing, restored custom local dynamics, and the later transition-movement extension.</p><p><code>source/latest/</code></p></div></div>
<div class="callout good"><strong>Strong regression result.</strong> Running the original fecundity validation workflow against <em>current</em> sources produced the same six CSV evidence files byte-for-byte. This is stronger than rounded numerical agreement: it shows the later default-off/default-equivalent refinements did not alter these recorded fecundity tests.</div>
{table(compare,['file','same'],['Evidence file','Byte-identical: archived vs latest rerun'])}
<h2>Current Meta source checksums</h2>
{table(manifest_rows,['file','bytes','sha'],['File','Bytes','SHA-256 (prefix)'])}
<p class="small">Full SHA-256 values for both snapshots are in <code>docs/source_snapshot_manifest.csv</code>; the complete bundle also has a root checksum file.</p>

<h1 id="validation">9. Validation evidence</h1>
<h2>Serial/parallel-function agreement in the executable one-core environment</h2>
{table(equiv,['model','output','equal','diff'],['Model','Saved output','Status','Max absolute difference'])}
<p>All 12 saved outputs compared exactly. This covered binary INApest as a reference for the common parallel architecture, MLU Meta, and transition-matrix Meta.</p>
<h2>Input-shape smoke tests</h2>
{table(smoke,['model','shape','status'],['Model','FecundityReduction parameterisation / edge case','Status'])}
<p>All agreed parameter forms executed. The two transition-matrix edge checks—unused stage-1 value and ambiguous node/stage vector—also behaved as designed.</p>
<details><summary>Serial/parallel equivalence for every MLU and transition-matrix input shape</summary><div>{table(shapeeq,['model','shape','equal','diff'],['Model','Shape','Status','Max absolute difference'])}</div></details>
<h2>Current ordinary Meta serial counterpart</h2>
<p>The serial ordinary <code>INApestMeta</code> function was restored/constructed after the initial fecundity snapshot. A supplemental current-source test now checks its three supported fecundity shapes against <code>INApestMetaParallel</code> on the same one-core RNG path.</p>
{table(meta_current,['shape','equal','diff'],['Shape','Status','Max absolute difference'])}
<h2>Backward compatibility and true multi-worker scope</h2>
<div class="two"><div class="implementation-card"><h3 class="check">Fecundity default</h3><p>The original executed regression check reported <code>max_abs_difference = 0</code> against the pre-change serial MLU and transition-matrix models when <code>FecundityReduction = 0</code>. The current source rerun further reproduces all archived evidence tables exactly.</p></div><div class="implementation-card"><h3 class="pending">Native PSOCK</h3><p>WebR cannot create native R subprocesses, so actual multi-worker PSOCK execution has not been run in this environment. <code>scripts/05_native_psock_validation_current.R</code> forces more than one worker and compares serial and parallel Monte Carlo summaries using a four-combined-SE criterion.</p></div></div>

<h1 id="trials">10. Demonstration simulations</h1>
<div class="callout"><strong>These are implementation demonstrations, not calibrated ecological predictions.</strong> Connectivity is deliberately simplified (identity SDD, no LDD) so the population response can be attributed cleanly to mortality and fecundity management rather than spatial network effects.</div>
<h2>Mortality × fecundity reduction</h2>
<p>The ordinary Meta demonstration used four nodes, six timesteps and 60 stochastic replicates for each of 15 combinations: management mortality 0, 0.3 or 0.6 crossed with fecundity reductions 0, 0.25, 0.5, 0.75 or 1. Baseline propagule production was 0.6 per survivor, establishment was 0.0015 and carrying capacity was 500 per node.</p>
<div class="figure popoutable" data-title="Mortality × fecundity reduction"><img src="{img1}" alt="Final population under mortality and fecundity-reduction combinations"><div class="caption"><strong>Figure 1. Mortality and fecundity suppression act as separate, compounding mechanisms.</strong> Lines/points show mean final total population; intervals show the archived 5th–95th percentile across 60 stochastic replicates. At zero management mortality, complete fecundity suppression leaves the initial 100 survivors but prevents new recruitment; when mortality is also active, both survivor number and subsequent reproduction decline.</div></div>
{table(mort,['mort','fr','mean','sd','q','ratio'],['Mortality','Fecundity reduction','Mean final population','SD','5th–95th percentile','Mean relative to FR=0'])}
<p>The response is monotonic within each mortality level but is not expected to be linear at the final-population scale: suppression acts upstream of stochastic propagule production and recruitment, so effects compound through six population updates. At mortality 0.3, increasing fecundity reduction from 0 to 0.5 reduces the mean final population from 101.1 to 37.3; at complete suppression it falls to 11.4.</p>
<h2>Land-use-specific fecundity reduction</h2>
<p>The MLU experiment used three land uses. The first three nodes each contained one land use; the fourth was mixed. This design makes the land-use-specific mapping visible while also retaining a node in which the land-use contributions interact. The gradient treatment set fecundity reductions to 0, 0.5 and 0.9 across land uses.</p>
<div class="figure popoutable" data-title="Land-use fecundity gradient"><img src="{img2}" alt="Land-use-specific final populations under fecundity reduction gradient"><div class="caption"><strong>Figure 2. The land-use vector is applied to separate reproductive contributions.</strong> At management mortality 0.3, otherwise comparable single-land-use nodes finished at mean populations 20.5, 10.5 and 3.9 for fecundity reductions 0, 0.5 and 0.9. The ordering confirms that the vector is not collapsed to one node-level fecundity-reduction value before reproduction is calculated.</div></div>
{table(lu03,['pattern','fr','total','lu1','lu2','lu3','mixed'],['Pattern','FR: LU1 / LU2 / LU3','Mean final total','LU1 node','LU2 node','LU3 node','Mixed node'])}
<p>The gradient produces a larger final total than uniform 0.5 suppression because the unsuppressed land use remains a strong source. That is the operational meaning of a reproductive refuge: residual population in a weakly suppressed land use can contribute disproportionately to future propagules even when other parts of the same node are strongly affected.</p>

<h1 id="reproduce">11. How to reproduce the work</h1>
<h2>Bundle structure</h2>
<div class="filetree">INApest_FecundityReduction_repeatability_bundle_2026-08-19/
├── source/
│   ├── fecundity_validation_snapshot/   # exact source state behind archived evidence
│   └── latest/                          # current source state
├── inputs/                              # parameter contract and scenario manifest
├── scripts/
│   ├── 00_run_all.R
│   ├── 01_validate_fecundity_reduction_current.R
│   ├── 02_current_meta_serial_parallel_smoke.R
│   ├── 03_compare_archived_and_rerun.R
│   ├── 04_plot_results.R
│   ├── 05_native_psock_validation_current.R
│   └── 06_reproduce_scientific_results_serial.R
├── results/
│   ├── archived/                        # frozen original evidence
│   ├── latest_rerun/                    # current-source rerun in WebR
│   ├── native_serial_reproduction/      # serial current-source reproduction
│   └── webr_exact_rerun/                # exact rerun made from bundled WebR runtime
├── figures/
├── environment/
│   ├── webr_runtime/                    # working WebR 0.6.0 browser runtime
│   ├── run_webr_exact_reproduction.py
│   └── DEPENDENCIES.md
├── docs/                                # original implementation report, patch, manifests
└── report/                              # this HTML and its build script</div>
<h2>Route A — ordinary native R: reproduce the scientific results</h2>
<div class="formula">Rscript scripts/00_run_all.R</div>
<p>On a normal multi-core machine, the entry point uses serial current-source functions to reproduce the scientific smoke-test and demonstration tables without parallel RNG differences. The archived parameterisation, mortality × fecundity and land-use tables were independently reproduced exactly using this serial route under the recorded R runtime.</p>
<h2>Route B — validate actual multi-worker PSOCK</h2>
<div class="formula">Rscript scripts/05_native_psock_validation_current.R</div>
<p>This requires at least three detected cores and intentionally forces the parallel functions to use more than one PSOCK worker. It tests distributional/Monte Carlo agreement with serial simulations rather than exact RNG identity.</p>
<h2>Route C — reproduce the archived one-core evidence exactly</h2>
<div class="formula">python environment/run_webr_exact_reproduction.py</div>
<p>The WebR runtime used for the archived validation is bundled. The Python harness requires Playwright and a usable Chromium/Chrome installation. In this assembled bundle it reproduced all six archived evidence CSVs byte-for-byte. This route is primarily for provenance; most modelling work should use native R.</p>
<h2>Rebuild figures and report</h2>
<p><code>scripts/04_plot_results.R</code> regenerates the two demonstration figures from CSV evidence using base R graphics. <code>report/build_report.py</code> rebuilds this self-contained HTML from the saved result tables and figures.</p>

<h1 id="limits">12. Interpretation, scope and limitations</h1>
<ul><li><strong>The demonstration effect sizes are not species estimates.</strong> Their purpose is to show that the implemented mechanism responds in the expected direction and that different parameter shapes are mapped correctly.</li><li><strong>Management is conditional on the model's management state.</strong> A high <code>FecundityReduction</code> has no effect in a node or land-use component that is not currently managed.</li><li><strong>The default ordinary/MLU biology uses one node-level Poisson propagule draw.</strong> In the MLU model, land-use contributions determine its mean; the upgrade deliberately preserved this stochastic structure rather than introducing independent Poisson draws per land use.</li><li><strong>Spread reduction is a different estimand.</strong> Fecundity reduction changes the number of propagules produced; spread reduction changes the LDD fraction successfully leaving managed source components in the default local dynamics.</li><li><strong>Current custom local dynamics can replace the defaults.</strong> Once users supply their own dynamics, they are responsible for implementing any active fecundity-management mechanism they want; the driver prevents silent omission when the new parameter is non-zero.</li><li><strong>True PSOCK remains an external runtime check.</strong> The code path and test are supplied, but this report does not claim a multi-process execution was performed inside WebR.</li></ul>
<div class="callout good"><strong>Overall assessment.</strong> The upgrade has a clear mechanistic interpretation, preserves the zero-effect baseline, supports the requested spatial/stage/land-use/time dimensions, is exposed to custom local dynamics, and is accompanied by archived evidence plus a current-source rerun that reproduces that evidence exactly.</div>

<h1>Appendix A. Key evidence files</h1>
<ul><li><code>results/archived/mortality_by_fecundity_results.csv</code> — 15 mortality × fecundity combinations.</li><li><code>results/archived/landuse_mortality_fecundity_results.csv</code> — 12 mortality × land-use-pattern combinations.</li><li><code>results/archived/parameterisation_smoke_results.csv</code> — all supported forms and transition-matrix edge cases.</li><li><code>results/archived/equivalence_results.csv</code> — saved serial/parallel outputs on the one-core executable path.</li><li><code>results/latest_rerun/archived_vs_latest_comparison.csv</code> — byte-level regression result against current sources.</li><li><code>inputs/scenario_manifest.csv</code> — compact reconstruction of all main validation/demonstration scenarios.</li><li><code>docs/source_snapshot_manifest.csv</code> — source-file SHA-256 hashes for both frozen source states.</li></ul>
<div class="footer">Prepared as the repeatability record for the INApest <code>FecundityReduction</code> upgrade. Final HTML is self-contained; source files, scripts, archived outputs and runtime materials are packaged alongside it.</div>
<script>
function pop(el){{const ov=document.createElement('div');ov.className='overlay';const box=document.createElement('div');box.className='overlaybox';const btn=document.createElement('button');btn.className='close';btn.textContent='×';btn.onclick=()=>ov.remove();box.appendChild(btn);const clone=el.cloneNode(true);clone.classList.remove('popoutable');box.appendChild(clone);ov.appendChild(box);ov.onclick=e=>{{if(e.target===ov)ov.remove();}};document.body.appendChild(ov);}}
document.querySelectorAll('.popoutable').forEach(el=>el.addEventListener('dblclick',()=>pop(el)));
</script></body></html>'''
OUT.write_text(html_doc, encoding='utf-8')
print(OUT)
