"""09 · Robustez numerica y al supuesto distribucional.

(a) Error de Monte Carlo de los percentiles, por bootstrap sobre las
    realizaciones del escenario A.
(b) Sensibilidad del resultado a la familia distribucional del factor de
    recuperacion, que es el parametro dominante: triangular declarada,
    lognormal con los mismos momentos y uniforme sobre el rango declarado.

Salida: datos/derivados/robustez.csv
"""
import pathlib, sys
import numpy as np, pandas as pd
from scipy.stats import norm
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
cadena = __import__('05_cadena_oficial'); prop = __import__('06_propagacion')

DER = pathlib.Path(__file__).resolve().parents[1]/'datos'/'derivados'
ANCHO_OFIC, BOOT = 62.77, 2000
LO, MO, HI = cadena.RG_TRI
MEDIA_T = (LO + MO + HI)/3
VAR_T   = (LO**2 + HI**2 + MO**2 - LO*HI - LO*MO - HI*MO)/18
S2 = np.log(1 + VAR_T/MEDIA_T**2); MU_LN = np.log(MEDIA_T) - S2/2

def rg(fam, u):
    if fam == 'triangular': return cadena._tri_ppf(u, LO, MO, HI)
    if fam == 'uniforme':   return LO + u*(HI - LO)
    return np.exp(MU_LN + np.sqrt(S2)*norm.ppf(np.clip(u, 1e-9, 1-1e-9)))

def correr(fam, correlacionado, T, n, semilla=7):
    g = np.random.default_rng(semilla); tot = np.zeros(n)
    if correlacionado: u1, u2 = g.random(n), g.random(n)
    for a in cadena.AREAS:
        if not correlacionado: u1, u2 = g.random(n), g.random(n)
        qe, ql, qh = cadena.Q[a]
        tot += cadena.potencia(cadena._tri_ppf(u2, ql, qe, qh), rg(fam, u1), T[a])
    return tot

def main():
    fil = []
    print("(a) error de Monte Carlo de los percentiles, escenario A")
    print(f"    {'n':>8}{'dependencia':>16}{'ancho':>9}{'ee(ancho)':>11}{'ee(P5)':>9}{'ee(P95)':>9}")
    for n in (10_000, 100_000):
        T = prop.escenario_A(n)
        for corr, et in ((False,'independiente'), (True,'correlacionada')):
            x = cadena.correr(T, n, correlacionado=corr); g = np.random.default_rng(3)
            B = np.array([np.percentile(x[g.integers(0,n,n)], [5,50,95]) for _ in range(BOOT)])
            p5, _, p95 = np.percentile(x, [5,50,95])
            print(f"    {n:>8,}{et:>16}{p95-p5:>9.1f}{(B[:,2]-B[:,0]).std():>11.1f}"
                  f"{B[:,0].std():>9.1f}{B[:,2].std():>9.1f}".replace(',','.'))
            fil.append(['error_monte_carlo', f'n={n}', et, round(p95-p5,1),
                        round((B[:,2]-B[:,0]).std(),1), None])

    print("\n(b) familia distribucional del factor de recuperación, escenario A")
    print(f"    {'familia':>28}{'E[Rg]':>8}{'indep':>8}{'veces':>7}{'corr':>8}{'veces':>7}")
    T = prop.escenario_A(prop.N)
    for fam, et, esp in (('triangular','Triangular declarada',MEDIA_T),
                         ('lognormal','Lognormal, iguales momentos',MEDIA_T),
                         ('uniforme','Uniforme sobre el rango',(LO+HI)/2)):
        ai = np.diff(np.percentile(correr(fam, False, T, prop.N), [5,95]))[0]
        ac = np.diff(np.percentile(correr(fam, True,  T, prop.N), [5,95]))[0]
        print(f"    {et:>28}{esp:>8.3f}{ai:>8.0f}{ai/ANCHO_OFIC:>7.1f}{ac:>8.0f}{ac/ANCHO_OFIC:>7.1f}")
        fil.append(['familia_distribucional', et, round(esp,3), round(ai), round(ai/ANCHO_OFIC,1),
                    round(ac)])
    pd.DataFrame(fil, columns=['control','caso','c1','c2','c3','c4']).to_csv(DER/'robustez.csv', index=False)

if __name__ == '__main__':
    main()
