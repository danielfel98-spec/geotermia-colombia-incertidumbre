# Incertidumbre de la estimación oficial del potencial geotérmico de Colombia

Código y datos del trabajo de grado **«Estimación probabilística del potencial
geotérmico eléctrico de Colombia: reconstrucción de la cadena oficial y
propagación de las incertidumbres declaradas»**, Maestría en Estadística
Aplicada y Ciencia de Datos, Universidad El Bosque, 2026.

El trabajo reconstruye la cadena de cálculo con que el Servicio Geológico
Colombiano estima el potencial geotérmico del país en 1.170,20 MWe, propaga las
incertidumbres que el propio informe declara y contrasta el resultado con el
intervalo de confianza publicado.

## Qué reproduce este repositorio

| Resultado | Valor | Script |
|---|---|---|
| Reconstrucción de la cifra oficial | 1.171,77 MWe, error +0,13 % | `05_cadena_oficial.py` |
| Madurez de las aguas (Giggenbach) | 308 de 328 inmaduras (94 %) | `02_geotermometria.py` |
| Asignación espacial | 198 de 328 dentro de las 21 áreas | `03_asignacion_espacial.py` |
| Modelo jerárquico | μ = 192,0 °C, τ = 19,0 °C | `04_modelo_jerarquico.py` |
| Intervalo nacional propagado | 471 a 1.641 MWe de ancho, 7,5 a 26 veces el oficial | `06_propagacion.py` |
| Sensibilidad | *S*~T~ 0,778 para el factor de recuperación | `07_sobol.py` |

## Cómo ejecutarlo

```bash
pip install -r requirements.txt
python codigo/01_ingesta.py            # lee los archivos del SGC
python codigo/02_geotermometria.py     # geotermómetros y madurez
python codigo/03_asignacion_espacial.py  # punto en polígono
python codigo/04_modelo_jerarquico.py  # posterior de (mu, tau) por cuadratura
python codigo/05_cadena_oficial.py     # reconstrucción de la cadena del informe
python codigo/06_propagacion.py        # los tres escenarios de temperatura
python codigo/07_sobol.py              # Spearman e índices de Sobol
python codigo/08_figuras.py            # figuras 1 y 2
```

Los scripts van en orden y cada uno escribe en `datos/derivados/`. No hay
dependencias de geopandas, shapely, dbfread ni PyMC: el lector de DBF, el punto
en polígono y el ajuste bayesiano están implementados aquí con numpy y scipy.

Las semillas están fijas. El posterior del modelo jerárquico se obtiene por
cuadratura sobre una rejilla y no por MCMC, de modo que es exactamente
reproducible y no depende de la convergencia de un muestreador.

## Estructura

```
datos/
  crudos/             archivos tal como los entregó el SGC
    Manantiales/      Manantial_Termal.{dbf,shp,shx,prj,cpg}
                      Ultima_Medicion_Manantial.dbf.xlsx  (química, 327 registros)
                      Métodos analíticos y límites de cuantificación.docx
    Fumarolas/        Fumarolas.{dbf,shp,...}, Medicion_Fumarola.dbf
  informe_sgc/        transcripción de las tablas del informe de 2020
  derivados/          todo lo que producen los scripts
  areas_21.geojson    capa de áreas geotérmicas del SIGT (EPSG:4326)
codigo/               los ocho scripts, en orden de ejecución
figuras/              figuras 1 y 2 del artículo
```

## Fuentes de datos

**Base geoquímica.** Inventario Nacional de Manifestaciones Hidrotermales,
entregado por el Servicio Geológico Colombiano en respuesta al derecho de
petición con radicado **2026-310-011732-2**. Contiene 328 manantiales termales
georreferenciados, de los cuales 327 tienen registro de analítica con once
determinaciones de laboratorio, y 15 fumarolas. No estaba disponible en las
capas de acceso abierto del SIGT.

**Capas geoespaciales.** Sistema de Información de Geotermia (SIGT) 2025 del
Servicio Geológico Colombiano, vía ArcGIS REST:
`https://srvags.sgc.gov.co/arcgis/rest/services/Geotermia/Sistema_Informacion_Geotermia_2025/MapServer`

**Informe auditado.** Alfaro, C. M., Rueda-Gutiérrez, J. B., Casallas, Y. P.,
Rodríguez, G. Z. y Malo, J. E. (2020). *Estimación preliminar del potencial
geotérmico de Colombia*. Bogotá: Servicio Geológico Colombiano. Las tablas
transcritas en `datos/informe_sgc/` corresponden a sus páginas 46, 51 y 58.

## Nota sobre los datos crudos

Los archivos de `datos/crudos/` se publican tal como los entregó la entidad, sin
corrección ni limpieza, para que cualquiera pueda verificar el procesamiento
desde el principio. Dos cosas que conviene saber al usarlos:

- El campo de nombre de área viene nulo en los 328 manantiales, por eso la
  asignación a áreas geotérmicas se calcula por intersección espacial.
- La última fila del libro de química es una nota al pie, no un manantial, y el
  manantial con identificador 1522 tiene ubicación pero no registro de analítica.

## Licencia

Código bajo licencia MIT (ver `LICENSE`).

Los datos de `datos/crudos/` son del Servicio Geológico Colombiano. La respuesta
al derecho de petición no impone restricción de uso, reserva ni confidencialidad
sobre ellos, y la propia entidad señala que su información está disponible para
consulta pública. Se republican como información pública en los términos de la
Ley 1712 de 2014; cite al Servicio Geológico Colombiano como fuente y el radicado
2026-310-011732-2. Ver `datos/crudos/LEEME.md`.

## Cómo citar

> Ávila Pedraza, D. F. (2026). *Estimación probabilística del potencial
> geotérmico eléctrico de Colombia: reconstrucción de la cadena oficial y
> propagación de las incertidumbres declaradas* \[Trabajo de grado, Maestría en
> Estadística Aplicada y Ciencia de Datos\]. Universidad El Bosque.

Para consultas sobre este repositorio, use el sistema de incidencias de GitHub.
