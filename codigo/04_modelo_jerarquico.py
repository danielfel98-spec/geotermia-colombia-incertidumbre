"""04 · Modelo jerarquico bayesiano de temperatura de reservorio por area.

    T_i   ~ Normal(mu, tau^2)            nivel poblacional
    T_obs ~ Normal(T_i, sigma_i^2)       observacion con su incertidumbre
    mu    ~ Normal(180, 60)
    tau   ~ HalfNormal(60)

T_i se marginaliza de forma analitica, de modo que la verosimilitud de cada
observacion es Normal(mu, tau^2 + sigma_i^2). El posterior de (mu, tau) se
obtiene por cuadratura sobre una rejilla, sin MCMC, lo que lo hace exactamente
reproducible. Los posteriores por area salen de la formula conjugada.

Entrada : datos/derivados/observaciones_21_areas.csv
Salida  : datos/derivados/posterior_21_areas.csv  y los parametros poblacionales
"""
import pathlib
import numpy as np, pandas as pd
from scipy import stats

RAIZ = pathlib.Path(__file__).resolve().parents[1]
DER  = RAIZ/'datos'/'derivados'
MU_RANGO, TAU_RANGO = (80, 330, 901), (1e-4, 160, 1101)

def ajustar(y, s, prior_tau=60.0):
    mu  = np.linspace(*MU_RANGO); tau = np.linspace(*TAU_RANGO)
    MU, TAU = np.meshgrid(mu, tau, indexing='ij')
    lp = stats.norm.logpdf(MU, 180, 60) + stats.halfnorm.logpdf(TAU, 0, prior_tau)
    for yi, si in zip(y, s):
        lp += stats.norm.logpdf(yi, MU, np.sqrt(TAU**2 + si**2))
    lp -= lp.max(); P = np.exp(lp); P /= P.sum()
    return mu, tau, P

def resumen(x, p, masa=0.94):
    m = (x*p).sum(); i = np.argsort(-p); c = np.cumsum(p[i]); k = i[c <= masa]
    return m, x[k].min(), x[k].max()

def ajuste_y_guardado(o, etiqueta):
    """Ajusta sobre el subconjunto o y guarda la rejilla con su etiqueta."""
    y = o.T_observada_C.values.astype(float); s = o.sigma_C.values.astype(float)
    mu, tau, P = ajustar(y, s)
    m_mu, lo_mu, hi_mu = resumen(mu, P.sum(1)); m_tau, lo_tau, hi_tau = resumen(tau, P.sum(0))
    print(f"[{etiqueta}] n={len(o):2d}  mu = {m_mu:6.1f} HDI94 [{lo_mu:.1f}, {hi_mu:.1f}]"
          f"   tau = {m_tau:5.1f} HDI94 [{lo_tau:.1f}, {hi_tau:.1f}]")
    np.save(DER/f'posterior_mu_tau_{etiqueta}.npy', P)
    return mu, tau, P


def main():
    o = pd.read_csv(DER/'observaciones_21_areas.csv')
    # ajuste con las 12 areas que tenian temperatura en la literatura: es el
    # contrafactual sin geoquimica primaria (escenario B de la propagacion)
    ajuste_y_guardado(o[o.nivel_evidencia != 'geoquimica propia'].reset_index(drop=True), '12obs')
    # ajuste con las 21: escenario C
    mu, tau, P = ajuste_y_guardado(o, '21obs')
    # posterior por area, integrando sobre (mu, tau)
    pm, pt = P.sum(1), P.sum(0)
    MU, TAU = np.meshgrid(mu, tau, indexing='ij')
    post = []
    for _, r in o.iterrows():
        si = r.sigma_C
        m = (r.T_observada_C/si**2 + MU/TAU**2)/(1/si**2 + 1/TAU**2)
        post.append((m*P).sum())
    o['T_posterior_C'] = np.round(post, 1)
    o['desplazamiento_C'] = np.round(o.T_posterior_C - o.T_observada_C, 1)
    o.to_csv(DER/'posterior_21_areas.csv', index=False)
    np.save(DER/'rejilla_mu.npy', mu); np.save(DER/'rejilla_tau.npy', tau)
    print(f"\n{len(o)} areas escritas en datos/derivados/posterior_21_areas.csv")
    print(o[['area','T_observada_C','sigma_C','T_posterior_C','desplazamiento_C']].to_string(index=False))

if __name__ == '__main__':
    main()
