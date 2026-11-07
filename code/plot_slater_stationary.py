"""Presentation-only render from sealed numerical data; no science rerun."""
from pathlib import Path
import json, hashlib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter, NullLocator
ROOT=Path(__file__).resolve().parent.parent
DIR=ROOT/'results/susskind-slater-stationary'
SOURCE=DIR/'stationary-checks.json'
FIG=DIR/'finite-time-states-layout.png'
RECORD=DIR/'layout-evidence.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
assert sha(SOURCE)=='ec0c798d12dbb09f6f8020d939636f105637f98007b402284f9e5f6c4317f9f9'
if FIG.exists() or RECORD.exists():raise FileExistsError('Preserve successful rendering')
data=json.loads(SOURCE.read_text())
fig,axes=plt.subplots(1,2,figsize=(10,3.5),layout='constrained')
colors=['#3465a4','#c06c21','#2a8b66']
for row,color in zip(data['rows'],colors):
    vals=row['quadratures'][1]['time_rows'];ts=[d['T'] for d in vals]
    for ax,key in zip(axes,('T_times_projection','T_times_purity')):
        ax.plot(ts,[d[key] for d in vals],marker='o',ms=4,color=color,label=f"q={row['q']}")
        ax.axhline(row['kappa'],ls=':',alpha=.7,color=color)
axes[0].set(ylabel='T times vacuum probability')
axes[1].set(ylabel='T times purity')
for ax in axes:
    ax.set_xscale('log',base=2)
    ax.set_xticks([2,4,8,16,32,64])
    ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.xaxis.set_minor_locator(NullLocator())
    ax.set_xlabel('Averaging half-width T')
    ax.grid(alpha=.2);ax.legend(frameon=False,fontsize=9)
fig.savefig(FIG,dpi=180);plt.close(fig)
record=dict(source_sha256=sha(SOURCE),script_sha256=sha(Path(__file__)),figure_sha256=sha(FIG),
            matplotlib=matplotlib.__version__,scope='Layout only: readable base-two ticks; values from sealed JSON unchanged.')
with RECORD.open('x') as f:json.dump(record,f,indent=2);f.write('\n')
print(json.dumps(record))
