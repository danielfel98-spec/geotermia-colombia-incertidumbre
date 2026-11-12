"""06 · Propagacion de las incertidumbres declaradas y agregacion nacional.

Tres escenarios, que se diferencian solo en el origen de la temperatura:
  A  temperatura fija, la que implica la tabla de eficiencia del informe
  B  12 areas del ajuste con 12 observaciones; las 9 restantes, predictiva
  C  las 21 del ajuste con 21 observaciones, 9 de ellas con silice propia
Dos estructuras de dependencia: independiente y correlacionada.

Salida: datos/derivados/escenarios_propagacion.csv
"""
import pathlib, sys
import numpy as np, pandas as pd
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
cadena = __import__('05_cadena_oficial')

RAIZ = pathlib.Path(__file__).resolve().parents[1]
DER  = RAIZ/'datos'/'derivados'
N          = 100_000
ANCHO_OFIC = 1201.58 - 1138.81          # 62,77 MWe

def hiperparametros(etiqueta, n, semilla=11):
    P  = np.load(DER/f'posterior_mu_tau_{etiqueta}.npy')
    mu = np.load(DER/'rejilla_mu.npy'); tau = np.load(DER/'rejilla_tau.npy')
    g  = np.random.default_rng(semilla); p = P.ravel()/P.sum()
    i  = g.choice(p.size, n, p=p)
    return mu[i//P.shape[1]], np.clip(tau[i % P.shape[1]], 1e-6, None), g

def escenario_A(n):
    return {a: np.full(n, cadena.temperatura_implicita(a)[0]) for a in cadena.AREAS}

def escenario_B(n, obs):
    MU, TAU, g = hiperparametros('12obs', n)
    lit = obs[obs.nivel_evidencia != 'geoquimica propia'].set_index('id_area')
    T = {}
    for a in cadena.AREAS:
        if a in lit.index:
            y, s = lit.loc[a, 'T_observada_C'], lit.loc[a, 'sigma_C']
            T[a] = g.normal((y/s**2 + MU/TAU**2)/(1/s**2 + 1/TAU**2),
                            np.sqrt(1/(1/s**2 + 1/TAU**2)))
        else:
            T[a] = g.normal(MU, TAU)        # predictiva del nivel poblacional
    return T

def escenario_C(n, obs):
    MU, TAU, g = hiperparametros('21obs', n)
    o = obs.set_index('id_area')
    return {a: g.normal((o.loc[a,'T_observada_C']/o.loc[a,'sigma_C']**2 + MU/TAU**2)
                        /(1/o.loc[a,'sigma_C']**2 + 1/TAU**2),
                        np.sqrt(1/(1/o.loc[a,'sigma_C']**2 + 1/TAU**2)))
            for a in cadena.AREAS}

def main():
    obs = pd.read_csv(DER/'observaciones_21_areas.csv')
    fil = []
    for nom, T in (('A. Temperatura fija del informe', escenario_A(N)),
                   ('B. Predictiva, sin geoquímica',   escenario_B(N, obs)),
                   ('C. Con geoquímica propia',        escenario_C(N, obs))):
        for corr, et in ((False, 'independiente'), (True, 'correlacionada')):
            x = cadena.correr(T, N, correlacionado=corr)
            p5, p50, p95 = np.percentile(x, [5, 50, 95])
            fil.append([nom, et, round(p5), round(p50), round(p95), round(x.mean()),
                        round(p95-p5), round((p95-p5)/ANCHO_OFIC, 1),
                        round((p95-p5)/x.mean(), 3)])
    d = pd.DataFrame(fil, columns=['escenario','dependencia','P5','P50','P95','media',
                                   'ancho','veces_oficial','ancho_sobre_media'])
    d.to_csv(DER/'escenarios_propagacion.csv', index=False)
    print(d.to_string(index=False))

if __name__ == '__main__':
    main()
