# JoanCLl
Engineering projects
# Análisis numérico de flujos magnetohidrodinámicos (MHD)

Trabajo de Fin de Grado — Ingeniería Mecánica, UPC (EEBE).

Análisis numérico de flujos magnetohidrodinámicos (fluidos conductores
bajo campos magnéticos) aplicados a los sistemas de refrigeración por
metales líquidos de los reactores de fusión nuclear. El trabajo estudia
dos configuraciones clásicas —el flujo en canales rectangulares
(problema de Hunt) y el flujo en expansiones bruscas— y valida los
resultados frente a datos experimentales de referencia.

## Qué hace este proyecto

A partir de simulaciones numéricas, se desarrolló un conjunto de
herramientas de postproceso en Python para extraer, analizar y
visualizar las variables físicas del flujo (velocidad, presión,
potencial eléctrico y densidad de corriente), y compararlas con datos
publicados en la literatura.

## Herramientas

- **Python** — desarrollo de todas las herramientas de postproceso.
- **ParaView (pvpython)** — lectura de mallas y extracción de datos de
  las simulaciones.
- **VTK / numpy_support** — conversión de datos científicos a arrays de NumPy.
- **NumPy** — cálculo numérico y tratamiento de datos.
- **Matplotlib** — generación automática de gráficas.
- **Julia** — ejecución de las simulaciones numéricas.

## Contenido

- `Script_Menu_Postprocess.py` — herramienta principal con menú
  interactivo de postproceso para el problema de Hunt: perfiles de
  velocidad y de densidad de corriente, isolíneas, comparativa de la
  fuerza de Lorentz frente a la solución teórica, cálculo de caudales
  (Q, K) con exportación a Excel, perfiles de velocidad en 3D y
  streamlines del campo de corriente.
- `ComparativaCaidaPresiones.py` — comparación automática de la caída de
  presión en la expansión entre la simulación y los datos experimentales
  de referencia (archivos .agr), para distintos números de Hartmann y Reynolds.
- `PerfilesV3D_Expansiones.py` — extracción y visualización 3D de los
  perfiles de velocidad en distintas secciones transversales de la expansión.

## Conceptos clave

Números de Hartmann y Reynolds, capas de Hartmann, fuerza de Lorentz,
refrigeración por metales líquidos, reactores de fusión (ITER / breeding
blankets).
