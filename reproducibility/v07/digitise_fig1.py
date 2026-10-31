"""Digitise bar heights in de Zelicourt et al. 2018 PLoS Genet Fig 1 (alfalfa biomass, t/ha).
Calibration: y-axis tick marks (0..14 t/ha, step 2) detected at x=156; bars detected by fill colour
(mock grey, SA187 green); top = first fill row from above (excluding error-bar pixels, which are darker grey)."""
import hashlib, json, sys
from pathlib import Path
import numpy as np
from PIL import Image

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / 'data/external/de_zelicourt_2018/fig1.png'
OUT = ROOT / 'data/external/de_zelicourt_2018/fig1_digitised.json'
SHA = '2ec2c32c63f7cdafea3c43916138593af5d571989ee13269cd2e0227a35eb839'
# percentages printed above the bars in Fig 1 (SA187 vs mock), used as an internal check
PRINTED = {'low_3.12': {'fresh_s1': .16, 'fresh_s2': .12, 'dry_s1': .14, 'dry_s2': .17},
           'high_7.81': {'fresh_s1': .15, 'fresh_s2': .12, 'dry_s1': .13, 'dry_s2': .17}}
if hashlib.sha256(SRC.read_bytes()).hexdigest() != SHA:
    sys.exit('fig1.png does not match the pinned sha256')
im = np.asarray(Image.open(SRC).convert('RGB')).astype(int)
dark = im.sum(2) < 300

def runs(ys):
    g = []
    for y in ys:
        if g and y - g[-1][-1] <= 2: g[-1].append(y)
        else: g.append([y])
    return [float(np.mean(x)) for x in g]

ticks = runs(np.where(dark[:, 156])[0])
panels = [ticks[:8], ticks[8:16]]          # 8 ticks (0..14) per panel, top to bottom
out = {}
grey = (abs(im[..., 0] - 192) < 12) & (abs(im[..., 1] - 192) < 12) & (abs(im[..., 2] - 192) < 12)
green = (im[..., 1] > 180) & (im[..., 0] > 120) & (im[..., 0] < 170) & (im[..., 2] < 110)
labels = ['fresh_s1', 'fresh_s2', 'dry_s1', 'dry_s2']
for pname, tk in zip(['low_3.12', 'high_7.81'], panels):
    tk = sorted(tk)                         # y pixel of 14,12,...,0
    vals = np.arange(14, -1, -2)
    slope, icpt = np.polyfit(tk, vals, 1)
    resid = float(np.max(np.abs(np.polyval([slope, icpt], tk) - vals)))
    y0, y1 = int(tk[0]) - 5, int(tk[-1])
    for name, mask in (('mock', grey), ('sa187', green)):
        colsum = mask[y0:y1].sum(0)
        xs = np.where(colsum > 20)[0]
        bars = []
        for x in xs:
            if bars and x - bars[-1][-1] <= 1: bars[-1].append(x)
            else: bars.append([x])
        bars = [b for b in bars if len(b) > 30]
        for lab, b in zip(labels, bars):
            # sample a 6-px strip near the left edge, away from the central error bar
            xc = slice(b[0] + 6, b[0] + 12)
            col = mask[y0:y1, xc].all(1)
            # first row that starts a run of >= 15 fill rows (skips legend swatches and anti-aliasing)
            run = np.array([col[i:i + 15].all() for i in range(len(col) - 15)])
            first_fill = y0 + int(np.argmax(run))
            # bar value = middle of the black top border just above the first fill row
            border = [y for y in range(first_fill - 8, first_fill) if dark[y, b[0] + 9]]
            top = float(np.mean(border)) if border else first_fill - 1.5
            out[f'{pname}/{lab}/{name}'] = round(float(np.polyval([slope, icpt], top)), 3)
    out[f'{pname}/calibration_max_residual_t_ha'] = round(resid, 4)
checks = {}
for pname, d in PRINTED.items():
    for lab, pct in d.items():
        got = out[f'{pname}/{lab}/sa187'] / out[f'{pname}/{lab}/mock'] - 1
        checks[f'{pname}/{lab}'] = {'printed': pct, 'digitised': round(got, 4), 'abs_diff': round(abs(got - pct), 4)}
worst = max(c['abs_diff'] for c in checks.values())
if worst > .015:
    sys.exit(f'digitised effects disagree with printed percentages by {worst:.3f}')
ratios = {t: round(out[f'high_7.81/{t}_s1/mock'] / out[f'low_3.12/{t}_s1/mock'], 4) for t in ('fresh', 'dry')}
result = {'source': 'de Zelicourt A et al. (2018) PLoS Genet 14:e1007273, doi:10.1371/journal.pgen.1007273, Fig 1 (CC BY 4.0)',
          'image_sha256': SHA, 'unit': 't/ha', 'method': __doc__.strip(),
          'values': out, 'check_against_printed_percentages': checks, 'max_abs_diff': worst,
          'season1_mock_ratio_high_over_low': ratios,
          'note': 'Bar means read from pixels; per-plot data are not published. Season 2 is confounded by rainfall (authors).'}
OUT.write_text(json.dumps(result, indent=1, sort_keys=True) + '\n')
print(json.dumps({'max_abs_diff': worst, 'ratios': ratios}))
