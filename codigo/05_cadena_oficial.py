"""Cadena reconstruida con los parametros declarados en Alfaro et al. (2020).
Ancla: el calor almacenado por area publicado en la Tabla 9, que ya incorpora
area, espesor, rho c, T_reservorio y T_referencia tal como el informe los calculo.
Lo unico que este trabajo cambia es (a) la fuente de la temperatura y (b) la
estructura de dependencia al agregar."""
import numpy as np

FPS=0.9; VIDA=30*365.25*24*3600.0
RG_TRI=(0.05,0.10,0.25)                      # Tabla 6 del informe
T_EC=np.array([120,140,160,180,200,220,240,260,280],float)
EC_V=np.array([0.01,0.02,0.04,0.05,0.08,0.11,0.13,0.15,0.17],float)
_cof=np.polyfit(T_EC,EC_V,1)                 # extrapolacion lineal, nota de la Tabla 7

def EC(t):
    t=np.asarray(t,float)
    return np.clip(np.where((t>=120)&(t<=280), np.interp(t,T_EC,EC_V), np.polyval(_cof,t)), 0.0, None)

# Tabla 9: calor almacenado (EJ) con su intervalo de confianza del 90 %
Q={'PP':(4.31,3.41,5.22),'IZA':(2.72,2.27,3.14),'SD':(12.51,11.45,13.60),'VCB':(7.94,6.96,8.92),
   'VT':(4.83,4.03,5.62),'NBL':(12.19,10.55,13.83),'HG':(11.57,10.39,12.76),'VSR':(10.66,9.27,12.05),
   'LO':(0.63,0.30,0.95),'VNT':(8.66,7.50,9.82),'VCM':(10.05,8.29,11.81),'VH':(0.76,0.37,1.14),
   'CGL':(5.15,4.55,5.75),'CP':(14.27,12.86,15.67),'VSS':(2.82,2.37,3.27),'VDJA':(5.30,4.62,5.99),
   'VGM':(4.87,4.22,5.51),'VS':(3.09,2.66,3.52),'VAZ':(9.60,8.69,10.52),'VC':(2.56,1.59,3.51),
   'CHCN':(4.14,3.50,4.80)}
MW_OF={'PP':21.50,'IZA':12.09,'SD':141.85,'VCB':79.73,'VT':38.50,'NBL':100.72,'HG':67.24,
   'VSR':137.24,'LO':0.08,'VNT':82.70,'VCM':129.94,'VH':0.10,'CGL':24.78,'CP':117.96,
   'VSS':17.43,'VDJA':37.84,'VGM':29.49,'VS':9.80,'VAZ':81.90,'VC':15.66,'CHCN':23.77}
N_CLUS={'PP':4,'IZA':3,'SD':6,'VCB':4,'VT':3,'NBL':5,'HG':9,'VSR':3,'LO':1,'VNT':4,'VCM':2,
   'VH':1,'CGL':4,'CP':8,'VSS':2,'VDJA':3,'VGM':4,'VS':3,'VAZ':6,'VC':2,'CHCN':3}
AREAS=list(Q)

def potencia(q_EJ, rg, t_c):
    """P (MWe) = Q (EJ) * Rg * EC(T) / (FPS * vida)."""
    return np.asarray(q_EJ,float)*1e18*np.asarray(rg,float)*EC(t_c)/(FPS*VIDA)/1e6

def correr(T_muestras, n=10000, correlacionado=True, semilla=7, escala_q=True):
    """T_muestras: dict area -> array de temperaturas posteriores."""
    g=np.random.default_rng(semilla)
    if correlacionado:
        u_rg=g.random(n); u_q=g.random(n)      # un solo sorteo comun a las 21 areas
    tot=np.zeros(n)
    for a in AREAS:
        if not correlacionado:
            u_rg=g.random(n); u_q=g.random(n)
        rg=_tri_ppf(u_rg,*RG_TRI)
        qe,qlo,qhi=Q[a]
        q=_tri_ppf(u_q,qlo,qe,qhi) if escala_q else qe
        tot+=potencia(q,rg,T_muestras[a])
    return tot

def _tri_ppf(u,lo,mo,hi):
    u=np.asarray(u,float); c=(mo-lo)/(hi-lo)
    return np.where(u<c, lo+np.sqrt(u*(hi-lo)*(mo-lo)), hi-np.sqrt((1-u)*(hi-lo)*(hi-mo)))


# --- inversion de la tabla de eficiencia: recupera la temperatura que el
# --- informe uso por area, dato que no aparece publicado en el documento
def temperatura_implicita(area):
    """EC implicita = P*FPS*VIDA/(Q*Rg_moda); se invierte sobre la tabla 7."""
    ec = MW_OF[area]*1e6*FPS*VIDA/(Q[area][0]*1e18*0.10)
    if ec <= EC_V[0]:
        return 120.0, ec          # fila "menor a 120 C" de la tabla del informe
    return float(np.interp(ec, EC_V, T_EC)), ec


if __name__ == '__main__':
    tot = sum(potencia(Q[a][0], 0.10, temperatura_implicita(a)[0]) for a in AREAS)
    print(f"reconstruccion en la moda de Rg: {tot:.2f} MWe  (informe: 1170,20)")
    print(f"error: {100*(tot/1170.20-1):+.2f} %")
    print(f"\n{'area':6}{'T implicita':>13}{'EC implicita':>14}")
    for a in AREAS:
        t, ec = temperatura_implicita(a); print(f"{a:6}{t:>13.1f}{ec:>14.4f}")
