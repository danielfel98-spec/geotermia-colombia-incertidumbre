"""01 · Ingesta de la base geoquimica del Inventario Nacional de Manifestaciones
Hidrotermales entregada por el SGC (radicado 2026-310-011732-2).

Entrada : datos/crudos/Manantiales/Ultima_Medicion_Manantial.dbf.xlsx (quimica)
          datos/crudos/Manantiales/Manantial_Termal.dbf              (ubicacion)
Salida  : datos/derivados/manantiales_geoquimica_sgc.csv

Tratamiento de censura: los valores reportados como "<x" quedan en x/2 y se
marcan con una bandera cen_*. Los "N.A." quedan nulos.
"""
import re, struct, sys, pathlib
import numpy as np, pandas as pd

RAIZ = pathlib.Path(__file__).resolve().parents[1]
CRUDOS = RAIZ/'datos'/'crudos'/'Manantiales'
SALIDA = RAIZ/'datos'/'derivados'

# limites de cuantificacion declarados en el documento de metodos del SGC
LC = {'co3':6.0,'hco3':12.2,'cl':0.2,'so4':0.2,'ca':0.1,'mg':0.1,'na':0.1,
      'k':0.1,'si':0.5,'sio2':1.1,'b':0.02,'li':0.05}

COLS = {'Carbonato':'co3','Bicarbonato':'hco3','Cloruros':'cl','Sulfatos':'so4',
        'Calcio':'ca','Magnesio':'mg','Sodio':'na','Potasio':'k','Silicio':'si',
        'Sílice':'sio2','Boro':'b','Litio':'li'}


def leer_dbf(ruta, encoding='latin-1'):
    """Lector DBF minimo: no hay dependencia de dbfread ni de geopandas."""
    b = pathlib.Path(ruta).read_bytes()
    nreg, hlen, rlen = struct.unpack('<I H H', b[4:12])
    campos, pos = [], 32
    while b[pos] != 0x0D:
        d = b[pos:pos+32]
        campos.append((d[:11].split(b'\0')[0].decode(encoding).strip(),
                       chr(d[11]), d[16]))
        pos += 32
    filas = []
    for i in range(nreg):
        r = b[hlen + i*rlen: hlen + (i+1)*rlen]
        if r[:1] == b'*':            # registro borrado
            continue
        o, fila = 1, {}
        for nom, tipo, ancho in campos:
            v = r[o:o+ancho].decode(encoding, 'replace').strip(); o += ancho
            if tipo in 'NF':
                try: v = float(v) if v else None
                except ValueError: v = None
            elif tipo == 'L':
                v = v.upper() in ('T','Y')
            fila[nom] = v
        filas.append(fila)
    return pd.DataFrame(filas)


def valor_censurado(x, lc):
    """'<0,05' -> (lc/2, True). Numero -> (numero, False). 'N.A.' -> (nan, False)."""
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return np.nan, False
    s = str(x).strip().replace(',', '.')
    if s.upper() in ('N.A.', 'NA', 'ND', 'N.D.', ''):
        return np.nan, False
    if s.startswith('<'):
        # se usa el limite de cuantificacion declarado por el SGC cuando existe;
        # si no, el numero que acompaña al signo
        if lc: return lc/2.0, True
        try:    return float(re.sub(r'[^0-9.]', '', s))/2.0, True
        except ValueError: return np.nan, True
    try:    return float(s), False
    except ValueError: return np.nan, False


def main():
    q = pd.read_excel(CRUDOS/'Ultima_Medicion_Manantial.dbf.xlsx')
    q.columns = [str(c).split('\n')[0].strip().strip('()') for c in q.columns]
    # la ultima fila del libro es una nota al pie, no un manantial
    q = q[q['ID_MANANTIAL'].notna()].reset_index(drop=True)
    d = pd.DataFrame({'id_man': q['ID_MANANTIAL'], 'nombre': q['NOMBRE'],
                      'fecha': pd.to_datetime(q['FECHA'], errors='coerce')})
    for crudo, corto in [('Temperatura','t_sup'), ('Caudal','caudal'),
                         ('pH','ph'), ('Conductividad eléctrica','cond')]:
        d[corto] = [valor_censurado(v, None)[0] for v in q[crudo]]
    for crudo, corto in COLS.items():
        par = [valor_censurado(v, LC.get(corto)) for v in q[crudo]]
        d[corto]          = [p[0] for p in par]
        d['cen_'+corto]   = [p[1] for p in par]

    # la silice del SGC es estequiometrica a partir del silicio; se recalcula
    # para dejar constancia de la diferencia
    d['sio2_estequio'] = d['si']*(60.08/28.0855)

    # El universo de manantiales es el shapefile (328 registros), no el libro de
    # quimica (327): el manantial 1522 tiene ubicacion pero no analitica.
    loc = leer_dbf(CRUDOS/'Manantial_Termal.dbf')
    col_id = next(c for c in loc.columns if 'MANANTIAL' in c.upper() or c.upper().startswith('ID'))
    for cand in ('LATITUD','LAT','Y'):
        if cand in loc.columns: col_lat = cand; break
    for cand in ('LONGITUD','LON','X'):
        if cand in loc.columns: col_lon = cand; break
    loc = loc[[col_id, col_lat, col_lon]].rename(
        columns={col_id:'id_man', col_lat:'lat', col_lon:'lon'})
    d = loc.merge(d, on='id_man', how='left')          # left sobre el shapefile

    SALIDA.mkdir(parents=True, exist_ok=True)
    d.to_csv(SALIDA/'manantiales_geoquimica_sgc.csv', index=False)
    cen = {c.replace('cen_',''): int(d[c].fillna(False).sum()) for c in d.columns if c.startswith('cen_')}
    print(f"{len(d)} manantiales del inventario; {int(d.ph.notna().sum())} con registro de analitica")
    print("censurados por analito:", {k:v for k,v in cen.items() if v})


if __name__ == '__main__':
    main()
