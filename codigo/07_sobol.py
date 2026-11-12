"""07 · Indices de Sobol agrupados, de primer orden y totales (Saltelli).

Grupos de entrada: factor de recuperacion, calor almacenado y nivel termico.
Se agrupan porque bajo independencia el factor y el calor aportan 21 variables
cada uno, de modo que el diseno lleva 21+21+1 = 43 columnas.

Escenario C (con geoquimica propia), que es el unico en que las tres entradas
aportan varianza.

Salida: datos/derivados/sobol_indices.csv
"""
import pathlib, sys
import numpy as np, pandas as pd
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
cadena = __import__('05_cadena_oficial')

RAIZ = pathlib.Path(__file__).resolve().parents[1]
DER  = RAIZ/'datos'/'derivados'
N, SEMILLAS = 32_768, (11, 29, 47, 83)

obs = pd.read_csv(DER/'observaciones_21_areas.csv').set_index('id_area')
P   = np.load(DER/'posterior_mu_tau_21obs.npy')
MUg = np.load(DER/'rejilla_mu.npy'); TAUg = np.load(DER/'rejilla_tau.npy')
cdf = np.cumsum(P.ravel()/P.sum())
NA  = len(cadena.AREAS)
Y   = np.array([obs.loc[a,'T_observada_C'] for a in cadena.AREAS])
S   = np.array([obs.loc[a,'sigma_C']       for a in cadena.AREAS])
QE  = np.array([cadena.Q[a][0] for a in cadena.AREAS])
QL  = np.array([cadena.Q[a][1] for a in cadena.AREAS])
QH  = np.array([cadena.Q[a][2] for a in cadena.AREAS])

def temperatura(u):
    """u -> (mu,tau) por inversa de la cdf -> media condicional por area."""
    i = np.searchsorted(cdf, np.clip(u, 1e-12, 1-1e-12))
    MU  = MUg[i//P.shape[1]]; TAU = np.clip(TAUg[i % P.shape[1]], 1e-6, None)
    return (Y[None,:]/S[None,:]**2 + MU[:,None]/TAU[:,None]**2) \
         / (1/S[None,:]**2 + 1/TAU[:,None]**2)

def modelo(U, correlacionado):
    U = np.atleast_2d(U)
    urg, uq, T = U[:,:NA], U[:,NA:2*NA], temperatura(U[:,2*NA])
    if correlacionado:                       # una sola decision para las 21
        urg = np.repeat(urg[:,:1], NA, axis=1); uq = np.repeat(uq[:,:1], NA, axis=1)
    rg = cadena._tri_ppf(urg, *cadena.RG_TRI)
    q  = cadena._tri_ppf(uq, QL[None,:], QE[None,:], QH[None,:])
    return (q*1e18*rg*cadena.EC(T)/(cadena.FPS*cadena.VIDA)/1e6).sum(axis=1)

GRUPOS = {'Factor de recuperación': np.arange(0, NA),
          'Calor almacenado':       np.arange(NA, 2*NA),
          'Nivel térmico':          np.array([2*NA])}

def saltelli(correlacionado, semilla):
    g = np.random.default_rng(semilla); k = 2*NA + 1
    A, B = g.random((N, k)), g.random((N, k))
    yA, yB = modelo(A, correlacionado), modelo(B, correlacionado)
    V = np.concatenate([yA, yB]).var()
    out = {}
    for nom, cols in GRUPOS.items():
        AB = A.copy(); AB[:, cols] = B[:, cols]
        yAB = modelo(AB, correlacionado)
        out[nom] = (np.mean(yB*(yAB - yA))/V, np.mean((yA - yAB)**2)/(2*V))
    return out

def spearman(correlacionado, n=100_000, semilla=7):
    """Correlacion de rango entre cada entrada y el total, para comparar."""
    from scipy.stats import spearmanr
    g = np.random.default_rng(semilla)
    u_rg = g.random(n); u_q = g.random(n); u_T = g.random(n)
    T = temperatura(u_T); tot = np.zeros(n)
    RG = np.empty((n, NA)); QQ = np.empty((n, NA))
    for k, a in enumerate(cadena.AREAS):
        RG[:, k] = cadena._tri_ppf(u_rg if correlacionado else g.random(n), *cadena.RG_TRI)
        QQ[:, k] = cadena._tri_ppf(u_q if correlacionado else g.random(n), QL[k], QE[k], QH[k])
        tot += cadena.potencia(QQ[:, k], RG[:, k], T[:, k])
    # bajo independencia cada entrada es un grupo de 21 variables: la correlacion
    # de rango se calcula contra la media del grupo, que es su resumen escalar
    return {'Factor de recuperación': spearmanr(RG.mean(axis=1), tot).statistic,
            'Calor almacenado':       spearmanr(QQ.mean(axis=1), tot).statistic,
            'Nivel térmico':          spearmanr(T.mean(axis=1), tot).statistic}


def main():
    fil = []
    for corr, et in ((True, 'correlacionada'), (False, 'independiente')):
        acc = {}; sp = spearman(corr)
        for s in SEMILLAS:
            for n, v in saltelli(corr, s).items(): acc.setdefault(n, []).append(v)
        print(f"\n{et}  (N={N:,}, {len(SEMILLAS)} semillas)".replace(',','.'))
        for n, v in acc.items():
            a = np.array(v)
            print(f"  {n:24} Spearman {sp[n]:+.3f}   S1 {a[:,0].mean():.3f} ±{a[:,0].std():.3f}"
                  f"   ST {a[:,1].mean():.3f} ±{a[:,1].std():.3f}")
            fil.append([et, n, round(sp[n],3), round(a[:,0].mean(),4), round(a[:,0].std(),4),
                        round(a[:,1].mean(),4), round(a[:,1].std(),4)])
        print(f"  suma ST = {sum(np.array(v)[:,1].mean() for v in acc.values()):.3f}")
    pd.DataFrame(fil, columns=['dependencia','entrada','spearman','S1','S1_ee','ST','ST_ee']) \
      .to_csv(DER/'sobol_indices.csv', index=False)

if __name__ == '__main__':
    main()
