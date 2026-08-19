#!/usr/bin/env python3
"""Re-run the archived one-core FecundityReduction validation in bundled WebR.

Requirements outside the bundle: Python 3, playwright, and a Chromium/Chrome
binary usable by Playwright. The WebR 0.6.0 runtime itself is bundled.
"""
from pathlib import Path
from playwright.sync_api import sync_playwright
import hashlib, json, shutil, sys

bundle = Path(__file__).resolve().parents[1]
runtime_dir = bundle / 'environment' / 'webr_runtime'
if not (runtime_dir / 'webr.mjs').exists():
    raise SystemExit('Bundled WebR runtime is missing')

srcdir = bundle / 'source' / 'latest'
script = bundle / 'scripts' / '01_validate_fecundity_reduction_current.R'
outdir = bundle / 'results' / 'webr_exact_rerun'
outdir.mkdir(parents=True, exist_ok=True)
output_names = [
    'equivalence_results.csv',
    'parameterisation_smoke_results.csv',
    'parameterisation_parallel_equivalence.csv',
    'mortality_by_fecundity_results.csv',
    'landuse_mortality_fecundity_results.csv',
    'validation_environment.csv',
]
source_files = [p.name for p in srcdir.iterdir() if p.is_file() and p.suffix.lower() == '.r']

def content_type(path):
    return {
        '.mjs':'text/javascript','.js':'text/javascript','.wasm':'application/wasm',
        '.so':'application/octet-stream','.r':'text/plain','.css':'text/css',
        '.gz':'application/gzip','.metadata':'application/octet-stream'
    }.get(path.suffix.lower(),'application/octet-stream')

browser_exe = shutil.which('chromium') or shutil.which('chromium-browser') or shutil.which('google-chrome')
with sync_playwright() as p:
    launch = {'headless': True, 'args':['--no-sandbox','--disable-gpu','--disable-web-security']}
    if browser_exe:
        launch['executable_path'] = browser_exe
    browser = p.chromium.launch(**launch)
    page = browser.new_page()
    def route_handler(route, request):
        url = request.url
        if url.startswith('https://webr.local/'):
            path = runtime_dir / url[len('https://webr.local/'):]
        elif url.startswith('https://src.local/'):
            path = srcdir / url[len('https://src.local/'):]
        elif url == 'https://script.local/validation.R':
            path = script
        else:
            route.fulfill(status=404, body='not found'); return
        if path.is_file():
            route.fulfill(status=200, body=path.read_bytes(), content_type=content_type(path), headers={'Access-Control-Allow-Origin':'*'})
        else:
            route.fulfill(status=404, body='not found')
    page.route('https://webr.local/**', route_handler)
    page.route('https://src.local/**', route_handler)
    page.route('https://script.local/**', route_handler)
    html = f'''<!doctype html><html><body><pre id="out">START</pre><script type="module">
import {{ WebR }} from 'https://webr.local/webr.mjs';
const out=document.getElementById('out'); const log=x=>out.textContent+='\\n'+x;
try {{
 const w=new WebR({{baseUrl:'https://webr.local/'}}); await w.init(); log('INIT_OK');
 await w.evalR("dir.create('/home/web_user/bundle/source/latest',recursive=TRUE,showWarnings=FALSE);dir.create('/home/web_user/bundle/results/out',recursive=TRUE,showWarnings=FALSE)");
 for (const f of {json.dumps(source_files)}) {{
   const txt=await (await fetch('https://src.local/'+f)).text();
   await w.FS.writeFile('/home/web_user/bundle/source/latest/'+f,new TextEncoder().encode(txt));
 }}
 const st=await (await fetch('https://script.local/validation.R')).text();
 await w.FS.writeFile('/home/web_user/bundle/validation.R',new TextEncoder().encode(st));
 await w.evalR("setwd('/home/web_user/bundle');Sys.setenv(INAPEST_BUNDLE_ROOT='/home/web_user/bundle',INAPEST_TEST_OUT='/home/web_user/bundle/results/out');source('validation.R')");
 window.__files={{}};
 for (const f of {json.dumps(output_names)}) {{
   const by=await w.FS.readFile('/home/web_user/bundle/results/out/'+f);
   window.__files[f]=new TextDecoder().decode(by);
 }}
 log('DONE');
}} catch(e) {{ log('ERROR '+e.stack); }}
</script></body></html>'''
    page.set_content(html, wait_until='domcontentloaded')
    page.wait_for_function("document.querySelector('#out').textContent.includes('DONE') || document.querySelector('#out').textContent.includes('ERROR')", timeout=240000)
    status = page.locator('#out').inner_text()
    print(status)
    if 'DONE' not in status:
        raise SystemExit('WebR validation failed')
    files = page.evaluate('window.__files')
    browser.close()

rows = []
for name, text in files.items():
    p = outdir / name
    p.write_text(text)
    a = (bundle / 'results' / 'archived' / name).read_bytes()
    b = p.read_bytes()
    rows.append((name, a == b, hashlib.sha256(a).hexdigest(), hashlib.sha256(b).hexdigest()))

cmp = outdir / 'archived_vs_webr_exact.csv'
with cmp.open('w', encoding='utf-8', newline='') as f:
    f.write('file,byte_identical,archived_sha256,rerun_sha256\n')
    for r in rows:
        f.write(','.join([r[0], str(r[1]).upper(), r[2], r[3]]) + '\n')
if not all(r[1] for r in rows):
    raise SystemExit('At least one evidence CSV differs from the archived reference')
print('PASS: all six archived evidence CSVs reproduced byte-for-byte.')
print(cmp)
