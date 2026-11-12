# Datos crudos entregados por el Servicio Geológico Colombiano

Esta carpeta contiene los archivos **tal como los entregó la entidad**, sin
corrección ni limpieza, para que el procesamiento pueda verificarse desde el
principio.

## Procedencia

Respuesta del Servicio Geológico Colombiano al derecho de petición presentado en
ejercicio del artículo 23 de la Constitución Política y regulado por la Ley 1755
de 2015.

- **Radicado de respuesta:** 2026-310-011732-2, del 21 de septiembre de 2026
- **Radicado padre:** 2026-290-008496-1
- **Dependencia:** Dirección de Geociencias Básicas, Grupo de Investigación Geotérmica

El punto 1 de la respuesta entrega «la base de datos disponible correspondiente a
la información geoquímica asociada al aplicativo de Hidrotermales del Servicio
Geológico Colombiano», en el formato disponible y conforme a los registros del
aplicativo.

## Condiciones de uso

La respuesta **no impone restricción de uso, reserva ni confidencialidad** sobre
la información entregada. La entidad describe su propia función como de
«generación, integración y divulgación de conocimiento geocientífico» y señala
que su información «se encuentra disponible para consulta pública a través del
Motor de Integración de Información Geocientífica (MIIG)».

Se republica aquí como información pública en los términos de la Ley 1712 de
2014. Al reutilizarla, cite al Servicio Geológico Colombiano como fuente y el
radicado de la respuesta.

## Dos omisiones deliberadas

El **oficio de respuesta** no se incluye porque su encabezado contiene datos
personales del peticionario.

Los **metadatos de autoría** de los tres archivos de Office (`dc:creator` y
`cp:lastModifiedBy`) se vaciaron antes de publicar, porque nombraban a personas
del equipo que preparó los archivos. El contenido de los archivos no se tocó: las
tablas, los registros y los valores son los entregados por la entidad.

## Contenido

```
Manantiales/
  Manantial_Termal.{dbf,shp,shx,prj,cpg,sbn,sbx,shp.xml}   328 manantiales georreferenciados
  Ultima_Medicion_Manantial.dbf.xlsx                       327 registros de analítica
  Métodos analíticos de referencia y límites de cuantificación.docx
Fumarolas/
  Fumarolas.{dbf,shp,shx,prj,cpg,sbn,sbx,shp.xml}          15 fumarolas
  Medicion_Fumarola.dbf / .xlsx                            serie de mediciones
  Ult_Medicion_Fumarola.dbf / .xml                         última medición por fumarola
```

## Dos cosas que conviene saber al usarlos

- El campo de nombre de área viene nulo en los 328 manantiales, de modo que la
  asignación a áreas geotérmicas hay que calcularla por intersección espacial.
- La última fila del libro de química es una nota al pie, no un manantial, y el
  manantial con identificador 1522 tiene ubicación pero no registro de analítica.
