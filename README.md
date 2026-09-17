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
# Optimización dinámica de una plataforma Stewart

Proyecto de la asignatura de Sistemas Mecánicos — Ingeniería Mecánica, UPC.

Optimización del comportamiento dinámico de una plataforma Stewart, un
mecanismo robótico paralelo de seis grados de libertad (6 DoF) formado
por seis actuadores que mueven una plataforma respecto a una base fija.
El objetivo es llevar la plataforma de una configuración inicial a un
estado final deseado minimizando criterios como el esfuerzo de los
actuadores o el consumo de energía.

## Qué hace este proyecto

El problema se formula como un problema de control óptimo y se resuelve
mediante optimización dinámica, respetando las restricciones reales del
sistema: límites de los actuadores, evitación de colisiones y las
ecuaciones dinámicas del mecanismo. Incluye también el cálculo de las
fuerzas en cada actuador mediante equilibrio de fuerzas y momentos
(dinámica inversa).

## Herramientas

- **MATLAB** — desarrollo del modelo y de los scripts de cálculo.
- **CasADi** — optimización numérica y control óptimo (Direct Collocation).
- **OpenSim** — modelado del mecanismo (`.osim`) y simulación del movimiento.

## Contenido

- `TreballMS1.m` — script principal del proyecto.
- `calculateUnitVectors.m` — cálculo de los vectores unitarios de los
  actuadores y resolución del sistema de fuerzas y momentos (dinámica
  inversa) mediante CasADi.
- `calculateUnitVectors_val.m` — versión de validación del cálculo anterior.
- `stewartplatform.osim` — modelo OpenSim de la plataforma.
- `motion2.mot`, `motion3.mot` — archivos de movimiento de la simulación.
- `Project_MS_JComabella.pdf` — memoria del proyecto.

## Conceptos clave

Mecanismo paralelo de 6 grados de libertad, control óptimo, Direct
Collocation, dinámica inversa, ángulos de Euler, optimización de
trayectorias.
# Diseño y análisis estructural de un kart

Proyecto académico de diseño mecánico — Ingeniería Mecánica, UPC.

Modelización completa de un kart para analizar los aspectos
fundamentales de su diseño estructural y dinámico. El diseño está
pensado con un enfoque práctico: componentes fabricables con materiales
accesibles y piezas comerciales, para que pueda construirse en un
entorno no industrial.

## Qué hace este proyecto

Definición de la geometría del chasis, selección de materiales, diseño
de cada componente (dirección, transmisión, asiento, soporte de motor)
y análisis del comportamiento del kart bajo distintas condiciones
mediante simulación estructural, evaluando resistencia, tensiones (von
Mises) y deformaciones de los componentes críticos.

## Herramientas

- **SolidWorks** — modelado 3D del chasis y de todos los componentes.
- **SolidWorks Simulation** — análisis estructural por elementos finitos
  (FEA): tensiones de von Mises, deformaciones y factor de seguridad.

## Materiales

- Chasis en acero AISI 1020 (acero al carbono de bajo contenido, por su
  resistencia y ductilidad).
- ~90% del kart en aluminio 6061 (ligereza y disponibilidad).
- Ruedas en caucho SBR.

## Contenido

- `MemoriaDiseñoKart.pdf` — memoria completa del diseño (geometría,
  materiales, componentes y análisis).
- `CDIMProyecto.pdf` — documentación del proyecto.

## Conceptos clave

Diseño CAD, análisis estructural (FEA), tensión de von Mises, diseño
para fabricación, selección de materiales, sistemas de dirección y
transmisión.
