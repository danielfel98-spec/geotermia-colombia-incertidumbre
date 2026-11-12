"""02 · Geotermometros quimicos y clasificacion de madurez de Giggenbach.

Entrada : datos/derivados/manantiales_geoquimica_sgc.csv (de 01_ingesta.py)
Salida  : la misma tabla con t_qtz, t_cal, t_nak, t_nakca, t_kmg, f_mg y clase

Geotermometros:
  cuarzo sin perdida de vapor y calcedonia .... Fournier (1977)
  Na-K ........................................ Giggenbach (1988)
  K-Mg ........................................ Giggenbach (1988)
  Na-K-Ca ..................................... Fournier y Truesdell (1973)

Madurez en el diagrama ternario Na/1000 - K/100 - raiz(Mg) (Giggenbach, 1988):
  f_Mg >= 0,5 inmadura | 0,2 a 0,5 equilibrio parcial | < 0,2 equilibrio
"""
import pathlib
import numpy as np, pandas as pd

RAIZ = pathlib.Path(__file__).resolve().parents[1]
TABLA = RAIZ/'datos'/'derivados'/'manantiales_geoquimica_sgc.csv'

def t_cuarzo(sio2):      return 1309.0/(5.19 - np.log10(sio2)) - 273.15
def t_calcedonia(sio2):  return 1032.0/(4.69 - np.log10(sio2)) - 273.15
def t_na_k(na, k):       return 1390.0/(1.75 + np.log10(na/k)) - 273.15
def t_k_mg(k, mg):       return 4410.0/(14.0 - np.log10(k**2/mg)) - 273.15

def t_na_k_ca(na, k, ca, beta=4/3):
    """Fournier y Truesdell (1973). beta = 4/3 por debajo de 100 C, 1/3 encima."""
    r = np.log10(na/k) + beta*(np.log10(np.sqrt(ca)/na) + 2.06)
    t = 1647.0/(r + 2.24) - 273.15
    r2 = np.log10(na/k) + (1/3)*(np.log10(np.sqrt(ca)/na) + 2.06)
    t2 = 1647.0/(r2 + 2.24) - 273.15
    return np.where(t > 100, t2, t)

def fraccion_mg(na, k, mg):
    """Coordenada de raiz(Mg) en el ternario de Giggenbach, normalizada."""
    a, b, c = na/1000.0, k/100.0, np.sqrt(np.clip(mg, 0, None))
    s = a + b + c
    return np.where(s > 0, c/s, np.nan)

def clasificar(f):
    return np.select([np.isnan(f), f >= 0.5, f >= 0.2],
                     ['sin dato', 'inmadura', 'parcial'], default='equilibrio')

def main():
    d = pd.read_csv(TABLA)
    # el SGC reporta silicio; la silice es estequiometrica
    sio2 = d['sio2'].where(d['sio2'].notna(), d['sio2_estequio'])
    with np.errstate(all='ignore'):
        d['t_qtz']   = np.where(sio2 > 0, t_cuarzo(sio2), np.nan)
        d['t_cal']   = np.where(sio2 > 0, t_calcedonia(sio2), np.nan)
        d['t_nak']   = np.where((d.na > 0) & (d.k > 0), t_na_k(d.na, d.k), np.nan)
        d['t_kmg']   = np.where((d.k > 0) & (d.mg > 0), t_k_mg(d.k, d.mg), np.nan)
        d['t_nakca'] = np.where((d.na > 0) & (d.k > 0) & (d.ca > 0),
                                t_na_k_ca(d.na, d.k, d.ca), np.nan)
        d['f_mg']    = fraccion_mg(d.na, d.k, d.mg)
    d['clase'] = clasificar(d['f_mg'].values)
    d.to_csv(TABLA, index=False)

    n = len(d); cl = d.clase.value_counts()
    con_quimica = d[['co3','hco3','cl','so4','ca','mg','na','k','si','sio2','b','li']].notna().any(axis=1).sum()
    print(f"{n} manantiales del inventario, {con_quimica} con alguna determinacion")
    for k_, v in cl.items(): print(f"  {k_:11} {v:3d}  ({100*v/n:.1f} %)")
    apl = int((d.t_qtz.notna()).sum())
    print(f"silice aplicable a {apl}: mediana {d.t_qtz.median():.1f} °C, P90 {d.t_qtz.quantile(.9):.1f} °C")
    viol = int((d.t_qtz < d.t_sup).sum())
    print(f"violaciones T_reservorio < T_surgencia: {viol} ({100*viol/apl:.1f} % de los aplicables)")

if __name__ == '__main__':
    main()
