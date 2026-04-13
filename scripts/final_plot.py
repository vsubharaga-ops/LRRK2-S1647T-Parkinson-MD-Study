import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def read_xvg(filepath):
    x, y = [], []
    with open(filepath) as f:
        for line in f:
            if line.startswith(('@','#')):
                continue
            cols = line.split()
            if len(cols) >= 2:
                try:
                    x.append(float(cols[0]))
                    y.append(float(cols[1]))
                except ValueError:
                    continue
    return np.array(x), np.array(y)

BASE = '/home/raga/pd_project/analysis'
WT  = '#2196F3'
MUT = '#E53935'

fig, axes = plt.subplots(3, 2, figsize=(14, 15))
fig.suptitle('WT vs S1647T LRRK2 — MD Comparative Analysis\n(AMBER99SB-ILDN, TIP3P, 300K, 1 bar)',
             fontsize=13, fontweight='bold')

# RMSD
ax = axes[0, 0]
t, r = read_xvg(f'{BASE}/rmsd_wt.xvg')
ax.plot(t, r, color=WT,  label=f'WT      mean={r.mean():.3f}±{r.std():.3f} nm', lw=1.5)
t, r = read_xvg(f'{BASE}/rmsd_mutant.xvg')
ax.plot(t, r, color=MUT, label=f'S1647T  mean={r.mean():.3f}±{r.std():.3f} nm', lw=1.5)
ax.set_xlabel('Time (ns)'); ax.set_ylabel('RMSD (nm)')
ax.set_title('Backbone RMSD'); ax.legend(fontsize=9); ax.grid(alpha=0.3)

# RMSF
ax = axes[0, 1]
res, f = read_xvg(f'{BASE}/rmsf_wt.xvg')
ax.plot(res, f, color=WT,  label='WT', lw=1.2)
res, f = read_xvg(f'{BASE}/rmsf_mutant.xvg')
ax.plot(res, f, color=MUT, label='S1647T', lw=1.2)
ax.axvline(x=1647, color='black', linestyle='--', lw=1.2, label='Res 1647')
ax.set_xlabel('Residue number'); ax.set_ylabel('RMSF (nm)')
ax.set_title('Per-residue RMSF (Cα)'); ax.legend(fontsize=9); ax.grid(alpha=0.3)

# Radius of gyration
ax = axes[1, 0]
t, g = read_xvg(f'{BASE}/gyrate_wt.xvg')
ax.plot(t, g, color=WT,  label=f'WT      mean={g.mean():.3f} nm', lw=1.5)
t, g = read_xvg(f'{BASE}/gyrate_mutant.xvg')
ax.plot(t, g, color=MUT, label=f'S1647T  mean={g.mean():.3f} nm', lw=1.5)
ax.set_xlabel('Time (ps)'); ax.set_ylabel('Rg (nm)')
ax.set_title('Radius of gyration'); ax.legend(fontsize=9); ax.grid(alpha=0.3)

# H-bonds
ax = axes[1, 1]
t, h = read_xvg(f'{BASE}/hbond_wt.xvg')
ax.plot(t, h, color=WT,  label=f'WT      mean={h.mean():.0f}', lw=1.5)
t, h = read_xvg(f'{BASE}/hbond_mutant.xvg')
ax.plot(t, h, color=MUT, label=f'S1647T  mean={h.mean():.0f}', lw=1.5)
ax.set_xlabel('Time (ps)'); ax.set_ylabel('Number of H-bonds')
ax.set_title('Intramolecular H-bonds'); ax.legend(fontsize=9); ax.grid(alpha=0.3)

# SASA
ax = axes[2, 0]
t, s = read_xvg(f'{BASE}/sasa_wt.xvg')
ax.plot(t, s, color=WT,  label=f'WT      mean={s.mean():.2f} nm²', lw=1.5)
t, s = read_xvg(f'{BASE}/sasa_mutant.xvg')
ax.plot(t, s, color=MUT, label=f'S1647T  mean={s.mean():.2f} nm²', lw=1.5)
ax.set_xlabel('Time (ps)'); ax.set_ylabel('SASA (nm²)')
ax.set_title('Solvent accessible surface area'); ax.legend(fontsize=9); ax.grid(alpha=0.3)

# Summary statistics table
ax = axes[2, 1]
ax.axis('off')
metrics = ['RMSD (nm)', 'Rg (nm)', 'H-bonds', 'SASA (nm²)']
files   = [('rmsd_wt.xvg','rmsd_mutant.xvg'),
           ('gyrate_wt.xvg','gyrate_mutant.xvg'),
           ('hbond_wt.xvg','hbond_mutant.xvg'),
           ('sasa_wt.xvg','sasa_mutant.xvg')]
rows = []
for metric, (fw, fm) in zip(metrics, files):
    _, yw = read_xvg(f'{BASE}/{fw}')
    _, ym = read_xvg(f'{BASE}/{fm}')
    diff = ym.mean() - yw.mean()
    rows.append([metric,
                 f'{yw.mean():.3f}±{yw.std():.3f}',
                 f'{ym.mean():.3f}±{ym.std():.3f}',
                 f'{diff:+.3f}'])
table = ax.table(
    cellText=rows,
    colLabels=['Metric', 'WT', 'S1647T', 'Δ (mut−WT)'],
    cellLoc='center', loc='center',
    bbox=[0, 0.2, 1, 0.7]
)
table.auto_set_font_size(False)
table.set_fontsize(10)
for (r, c), cell in table.get_celld().items():
    if r == 0:
        cell.set_facecolor('#DDDDDD')
        cell.set_text_props(fontweight='bold')
ax.set_title('Summary statistics', fontweight='bold')

plt.tight_layout()
outpath = f'{BASE}/WT_vs_S1647T_FINAL.png'
plt.savefig(outpath, dpi=300, bbox_inches='tight')
print(f"Saved: {outpath}")

# Print to terminal
print("\n=== Summary Statistics ===")
for row in rows:
    print(f"{row[0]:15s}  WT: {row[1]}   S1647T: {row[2]}   Δ: {row[3]}")
