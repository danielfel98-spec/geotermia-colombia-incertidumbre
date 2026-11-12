"""08 · Figuras 1 y 2 del articulo.

Figura 1 · diagrama ternario Na-K-Mg de Giggenbach con los 328 manantiales
Figura 2 · distribucion del potencial nacional bajo las dos estructuras de
           agregacion, con la banda del intervalo oficial

Salida: figuras/fig1_giggenbach.png y figuras/fig2_intervalo.png
"""
import pathlib, sys
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np, pandas as pd
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
cadena = __import__('05_cadena_oficial')
prop   = __import__('06_propagacion')

RAIZ = pathlib.Path(__file__).resolve().parents[1]
DER, FIG = RAIZ/'datos'/'derivados', RAIZ/'figuras'
S1, S2, S3 = '#2a78d6', '#eb6834', '#1baf7a'
INK, INK2, GRID = '#0b0b0b', '#52514e', '#d6d5d0'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.edgecolor':INK2,
    'axes.labelcolor':INK,'text.color':INK,'xtick.color':INK2,'ytick.color':INK2,
    'axes.linewidth':0.8,'figure.facecolor':'white','axes.facecolor':'white',
    'axes.spines.top':False,'axes.spines.right':False})

def figura1():
    d = pd.read_csv(DER/'manantiales_geoquimica_sgc.csv')
    d = d[d[['na','k','mg']].notna().all(axis=1)]
    a, b, c = d.na/1000.0, d.k/100.0, np.sqrt(d.mg.clip(lower=0))
    s = a + b + c; a, b, c = a/s, b/s, c/s
    x = 0.5*(2*b + c)/(a + b + c); y = (np.sqrt(3)/2)*c/(a + b + c)
    col = np.select([d.clase.values=='inmadura', d.clase.values=='parcial'], [S2, S1], default=S3)
    fig, ax = plt.subplots(figsize=(5.2, 4.6))
    for p, q in ((( 0,0),(1,0)), ((1,0),(0.5,np.sqrt(3)/2)), ((0.5,np.sqrt(3)/2),(0,0))):
        ax.plot([p[0],q[0]], [p[1],q[1]], color=INK2, lw=0.9)
    ax.scatter(x, y, s=11, c=col, alpha=0.75, edgecolor='none')
    for frac, et in ((0.5,'inmadura'), (0.2,'equilibrio parcial')):
        ax.plot([0.5*frac, 1-0.5*frac], [np.sqrt(3)/2*frac]*2, color=GRID, lw=0.8, ls='--')
        ax.text(1-0.5*frac+0.01, np.sqrt(3)/2*frac, f'√Mg = {frac}', fontsize=7, color=INK2, va='center')
    ax.text(-0.03,-0.03,'Na/1000',ha='right',fontsize=8); ax.text(1.03,-0.03,'K/100',fontsize=8)
    ax.text(0.5,np.sqrt(3)/2+0.03,'√Mg',ha='center',fontsize=8)
    from matplotlib.lines import Line2D
    n = d.clase.value_counts()
    ax.legend(handles=[Line2D([],[],marker='o',ls='',color=c,markersize=5,
              label=f'{e} ({n.get(e,0)})') for e,c in
              (('inmadura',S2),('parcial',S1),('equilibrio',S3))],
              frameon=False, fontsize=8, loc='upper left', labelcolor=INK)
    ax.set_axis_off(); ax.set_aspect('equal')
    fig.tight_layout(); fig.savefig(FIG/'fig1_giggenbach.png', dpi=200, bbox_inches='tight')
    print(f"figura 1: {len(d)} manantiales con Na, K y Mg")

def figura2():
    obs = pd.read_csv(DER/'observaciones_21_areas.csv'); N = prop.N
    T = prop.escenario_C(N, obs)
    tc = cadena.correr(T, N, correlacionado=True)
    ti = cadena.correr(T, N, correlacionado=False)
    fig, ax = plt.subplots(figsize=(6.2, 3.5))
    ax.hist(tc, bins=80, range=(0,4200), color=S1, alpha=0.45, edgecolor='white',
            linewidth=0.3, label='Agregación correlacionada')
    ax.hist(ti, bins=80, range=(0,4200), color=S3, alpha=0.45, edgecolor='white',
            linewidth=0.3, label='Agregación independiente')
    ymax = ax.get_ylim()[1]
    ax.axvspan(1138.81, 1201.58, color=S2, alpha=0.85, zorder=5)
    ax.annotate('intervalo oficial\n1.138,8 a 1.201,6\n(63 MWe de ancho)', xy=(1170, ymax*0.62),
                xytext=(1900, ymax*0.80), fontsize=8, color=S2,
                arrowprops=dict(arrowstyle='->', color=S2, lw=1.2))
    ax.set_xlabel('Potencial geotérmico nacional (MWe)'); ax.set_ylabel('Realizaciones')
    ax.set_xlim(0, 4200); ax.grid(axis='y', color=GRID, lw=0.6); ax.set_axisbelow(True)
    ax.legend(frameon=False, fontsize=8, loc='upper right', labelcolor=INK)
    fig.tight_layout(); fig.savefig(FIG/'fig2_intervalo.png', dpi=200, bbox_inches='tight')
    print(f"figura 2: P5-P95 correlacionado {np.percentile(tc,5):.0f} a {np.percentile(tc,95):.0f} MWe")

if __name__ == '__main__':
    FIG.mkdir(exist_ok=True); figura1(); figura2()
