# bone-electrical-conductivity

# Conductividad eléctrica del tejido óseo

Este repositorio contiene modelos numéricos, datos y scripts desarrollados para el estudio de las propiedades eléctricas del tejido óseo mediante un enfoque multiescala.

El trabajo se organiza en tres niveles de modelado: trabecular, cortical y macroescala.

## Hueso trabecular
Modelos computacionales de hueso trabecular con microestructuras inspiradas en imágenes reales de microtomografía computarizada (μCT). 
La respuesta eléctrica se obtiene mediante el método de elementos finitos (FEM). Los resultados se complementan con un modelo de medio efectivo basado en la formulación de Bruggeman para analizar la relación entre la microestructura ósea y la conductividad eléctrica efectiva.

## Hueso cortical

Modelos computacionales de hueso cortical con microestructuras paramétricas que representan los principales componentes del sistema vascular óseo. La respuesta eléctrica se obtiene mediante el método de elementos finitos (FEM), permitiendo analizar la relación entre la microestructura, la porosidad y la conductividad eléctrica en las direcciones axial y transversal.
os resultados se complementan con un modelo de medio efectivo basado en la formulación de Bruggeman para analizar la relación entre la porosidad y la conductividad eléctrica en las direcciones axial y transversal.

## Macroescala

Modelo computacional a macroescala de hueso que integra las propiedades eléctricas de los tejidos cortical, trabecular y medular. La respuesta eléctrica se obtiene mediante el método de elementos finitos (FEM), considerando diferentes condiciones del tejido.

## Estructura del repositorio

- `trabecular/` – Modelos, datos y scripts correspondientes al hueso trabecular.
- `cortical/` – Modelos, datos y scripts correspondientes al hueso cortical.
- `macroscale/` – Modelos, datos y scripts correspondientes al modelo macroscópico.

## Frecuencia

Los resultados eléctricos presentados corresponden a una frecuencia de **100 kHz**.
