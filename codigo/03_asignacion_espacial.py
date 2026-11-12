"""03 · Asignacion de cada manantial a su area geotermica por interseccion
punto en poligono. El nombre de area viene nulo en los 328 registros del SGC,
de modo que la asignacion hay que calcularla.

Punto en poligono por ray casting (regla de Jordan), solo con numpy, sin
geopandas ni shapely.

Entrada : datos/derivados/manantiales_geoquimica_sgc.csv
          datos/areas_21.geojson   (capa de areas geotermicas del SIGT, EPSG:4326)
Salida  : la misma tabla con la columna area_asignada
"""
import json, pathlib
import numpy as np, pandas as pd

RAIZ = pathlib.Path(__file__).resolve().parents[1]
TABLA = RAIZ/'datos'/'derivados'/'manantiales_geoquimica_sgc.csv'
CAPA  = RAIZ/'datos'/'areas_21.geojson'

def punto_en_poligono(x, y, poly):
    """Ray casting. poly: lista de (lon, lat). Devuelve True si (x,y) esta dentro."""
    dentro = False; n = len(poly); j = n - 1
    for i in range(n):
        xi, yi = poly[i]; xj, yj = poly[j]
        if ((yi > y) != (yj > y)) and (x < (xj - xi)*(y - yi)/(yj - yi + 1e-300) + xi):
            dentro = not dentro
        j = i
    return dentro

def anillos(geom):
    """Devuelve los anillos exteriores de un Polygon o MultiPolygon."""
    t, c = geom['type'], geom['coordinates']
    return [c[0]] if t == 'Polygon' else [p[0] for p in c]

def main():
    d = pd.read_csv(TABLA)
    capa = json.loads(CAPA.read_text(encoding='utf-8'))
    areas = [(f['properties'].get('NOMBRE') or f['properties'].get('nombre')
              or f['properties'].get('AREA_GEOTERMICA'), anillos(f['geometry']))
             for f in capa['features']]
    asign = []
    for _, r in d.iterrows():
        if pd.isna(r.lat) or pd.isna(r.lon): asign.append(None); continue
        hit = None
        for nombre, rings in areas:
            if any(punto_en_poligono(r.lon, r.lat, ring) for ring in rings):
                hit = nombre; break
        asign.append(hit)
    d['area_asignada'] = asign
    d.to_csv(TABLA, index=False)
    n, dentro = len(d), int(d.area_asignada.notna().sum())
    print(f"{dentro} de {n} manantiales ({100*dentro/n:.1f} %) caen dentro de alguna de las 21 areas")
    print(f"areas con al menos un manantial: {d.area_asignada.nunique()}")

if __name__ == '__main__':
    main()
