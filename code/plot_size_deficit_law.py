"""Plot the saved normalized limiting law; no finite-N sampling or recomputation."""
from pathlib import Path
import hashlib
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / 'results/susskind-memory-size-bridge'


if __name__ == '__main__':
    target = OUT / 'deficit-law.png'
    metadata = OUT / 'deficit-law-plot.json'
    if target.exists() or metadata.exists():
        raise FileExistsError('Plot output already exists')
    source = OUT / 'rank-law-checks.json'
    data = json.loads(source.read_text())
    rows = data['probabilities']['rows']
    selected = [next(r for r in rows if r['j'] == j and r['lambda_value'] == '1/4')
                for j in [1, 2, 4, 6]]
    fig, axes = plt.subplots(2, 2, figsize=(9, 5.4), constrained_layout=True)
    plot_records = []
    for ax, row, color in zip(axes.flat, selected, ['#246b8e', '#16735b', '#7c519c', '#ad532d']):
        values = [float(v) for v in row['normalization']['pmf']]
        shown = 28
        ax.bar(range(min(shown+1, len(values))), values[:shown+1], width=.77, color=color)
        ax.set_xlim(-.8, shown+.8)
        ax.set_xlabel('Deficit d')
        ax.set_ylabel('Probability')
        ax.set_title(f"j = {row['j']}   |   mean = {float(row['mean']):.3f}")
        ax.grid(axis='y', alpha=.2)
        ax.set_axisbelow(True)
        for side in ['top', 'right']:
            ax.spines[side].set_visible(False)
        plot_records.append(dict(j=row['j'], lambda_value='1/4', display_cutoff=shown,
                                 saved_summation_cutoff=row['normalization']['cutoff'],
                                 saved_probability_beyond_display=sum(values[shown+1:])))
    fig.suptitle('Limiting rank-deficit law of the highest Gaussian chaos', fontsize=13)
    fig.savefig(target, dpi=200)
    plt.close(fig)
    record = dict(scope='Visualization of the verified limiting probability law, not finite-N simulation.',
                  rows=plot_records, matplotlib_version=matplotlib.__version__,
                  input_sha256={source.relative_to(ROOT).as_posix():hashlib.sha256(source.read_bytes()).hexdigest(),
                                Path(__file__).relative_to(ROOT).as_posix():hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},
                  output_sha256=hashlib.sha256(target.read_bytes()).hexdigest())
    with metadata.open('x', encoding='utf-8') as handle:
        handle.write(json.dumps(record, indent=2)+'\n')
    print(json.dumps(dict(status='saved', plot=str(target), rows=plot_records)))
