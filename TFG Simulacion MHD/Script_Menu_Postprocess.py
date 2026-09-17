# coding: utf-8

# In[1]:


from paraview.simple import *
from vtk.util import numpy_support
import numpy as np
import matplotlib.pyplot as plt
import os
import re

paraview.simple._DisableFirstRenderCameraReset()

# =========================
# SELECCIÓN DE CARPETA
# =========================
base_path = r"C:\Users\joanc\Desktop\TFG\Modelos"

print("Selecciona carpeta:")
print("2 → fe_comparison_hunt_2")
print("3 → fe_comparison_hunt_3")

choice = input("Opción (2 o 3): ").strip()

if choice == "2":
    folder = "fe_comparison_hunt_2"
elif choice == "3":
    folder = "fe_comparison_hunt_3"
else:
    print(" Opción no válida, usando carpeta 2 por defecto")
    folder = "fe_comparison_hunt_2"

datadir = os.path.join(base_path, folder)
outdir = os.path.join(datadir, "postprocess")

os.makedirs(outdir, exist_ok=True)

print(f"\n Usando: {datadir}\n")

# =========================
# MENÚ PRINCIPAL
# =========================
print("\n=== MENÚ DE POSTPROCESO ===")
print("1 → Perfiles de velocidad y de corriente")
print("2 → Isolíneas")
print("3 → Comparativa Perfil de velocidad en funcion del numero de Ha")
print("4 → Comparativa fuerza de Lorentz")
print("5 → Calcular Q, K y exportar a Excel")
print("6 → Perfil de velocidad 3D ")
print("7 → Streamlines (campo de corriente j)")
print("10 → Ejecutar TODO")

# =========================
# FUNCIÓN Ha DESDE NOMBRE
# =========================
def extract_Ha(file):
    match = re.search(r"ha(\d+)", file)
    return int(match.group(1)) if match else None


# =========================
# FUERZA TEÓRICA
# =========================
def lorentz_force(y, Ha, u_hat=1.0):
    return -1 + u_hat * Ha**2 * np.cosh(Ha * y) / np.cosh(Ha)

option = input("Selecciona opción: ").strip()

if option == "1" or option == "3" or option == "10":
    print("\nDefine la línea:")

    x1 = float(input("x1: "))
    y1 = float(input("y1: "))
    z1 = float(input("z1: "))

    x2 = float(input("x2: "))
    y2 = float(input("y2: "))
    z2 = float(input("z2: "))

if option == "1" or option == "10":

    # =========================
    # LISTA DE ARCHIVOS
    # =========================
    files = [f for f in os.listdir(datadir) if f.endswith(".pvtu")]
    # =========================
    # ORDEN DE ARCHIVOS
    # =========================
    def sort_key(file):
        match = re.search(r"ha(\d+)", file)
        return int(match.group(1)) if match else 9999

    files = sorted(files, key=sort_key)

    # =========================
    # LOOP PRINCIPAL
    # =========================
    for file in files:

        print(f">> {file}")

        filename = os.path.join(datadir, file)

        # =========================
        # LEER MESH
        # =========================
        mesh = XMLPartitionedUnstructuredGridReader(
            registrationName=file,
            FileName=[filename]
        )

        mesh.TimeArray = 'None'
        SetActiveSource(mesh)

        # =========================
        # PLOT OVER LINE
        # =========================
        line = PlotOverLine(Input=mesh)
        line.Point1 = [x1, y1, z1]
        line.Point2 = [x2, y2, z2]
        line.SamplingPattern = 'Sample Uniformly'
        line.Resolution = 500

        data = servermanager.Fetch(line)
        point_data = data.GetPointData()

        array_names = [point_data.GetArrayName(i) for i in range(point_data.GetNumberOfArrays())]

        if not all(name in array_names for name in ['u', 'uh', 'u_ref']):
            print(f" Variables no encontradas en {file}")
            print("   Disponibles:", array_names)
            Delete(line)
            Delete(mesh)
            continue

        # =========================
        # EXTRAER DATOS
        # =========================
        s = numpy_support.vtk_to_numpy(point_data.GetArray('arc_length'))

        u_vec = numpy_support.vtk_to_numpy(point_data.GetArray('u'))
        uh_vec = numpy_support.vtk_to_numpy(point_data.GetArray('uh'))
        u_ref_vec = numpy_support.vtk_to_numpy(point_data.GetArray('u_ref'))

        #  componente correcta (Z)
        u = u_vec[:, 2]
        uh = uh_vec[:, 2]
        u_ref = u_ref_vec[:, 2]

        j_vec = numpy_support.vtk_to_numpy(point_data.GetArray('j'))
        jh_vec = numpy_support.vtk_to_numpy(point_data.GetArray('jh'))
        j_ref_vec = numpy_support.vtk_to_numpy(point_data.GetArray('j_ref'))

        j = np.linalg.norm(j_vec, axis=1)
        jh = np.linalg.norm(jh_vec, axis=1)
        j_ref = np.linalg.norm(j_ref_vec, axis=1)

        # =========================
        # PLOT VELOCIDAD
        # =========================
        plt.figure(figsize=(8,5))
        plt.plot(s, u, label='u', linewidth=2)
        plt.plot(s, uh, label='uh', linewidth=2)
        plt.plot(s, u_ref, '--', label='u_ref', linewidth=2)

        plt.xlabel("arc_length")
        plt.ylabel("velocity")
        plt.title(file.replace(".pvtu",""))
        plt.legend()
        plt.grid()

        name = file.replace(".pvtu", "")
        name_base = name.split("_ou")[0]

        output_vel = os.path.join(outdir, name_base + "_PerfilVelocidad.png")
        plt.savefig(output_vel, dpi=150)
        plt.close()

        print(f" Guardado: {output_vel}")

        # =========================
        # PLOT CORRIENTE
        # =========================
        plt.figure(figsize=(8,5))
        plt.plot(s, j, label='j', linewidth=2)
        plt.plot(s, jh, label='jh', linewidth=2)
        plt.plot(s, j_ref, '--', label='j_ref', linewidth=2)

        plt.xlabel("arc_length")
        plt.ylabel("current density")
        plt.title(file.replace(".pvtu","") + " - corriente")
        plt.legend()
        plt.grid()

        output_j = os.path.join(outdir, name_base + "_PerfilCorriente.png")
        plt.savefig(output_j, dpi=150)
        plt.close()

        print(f" Guardado: {output_j}")

        # =========================
        # LIMPIEZA
        # =========================
        Delete(line)
        Delete(mesh)

    print("\n PROCESO COMPLETADO")

if option == "2" or option == "10":
    # ISOLÍNEAS
    # =========================
    # PARÁMETROS
    # =========================
    ny = 50   # número de líneas en Y
    y_vals = np.linspace(-1, 1, ny)

    # =========================
    # LOOP ARCHIVOS
    # =========================
    for file in os.listdir(datadir):

        if not file.endswith(".pvtu"):
            continue

        print(f">> {file}")

        filename = os.path.join(datadir, file)

        mesh = XMLPartitionedUnstructuredGridReader(
            registrationName=file,
            FileName=[filename]
        )

        mesh.PointArrayStatus = ['u']
        mesh.TimeArray = 'None'

        SetActiveSource(mesh)

        # =========================
        # MATRIZ PARA DATOS
        # =========================
        U = []
        X = None

        # =========================
        # BUCLE EN Y
        # =========================
        for y in y_vals:

            line = PlotOverLine(Input=mesh)
            line.Point1 = [-1.0, y, 0.0]
            line.Point2 = [ 1.0, y, 0.0]
            line.SamplingPattern = 'Sample At Segment Centers'

            data = servermanager.Fetch(line)

            s_array = data.GetPointData().GetArray('arc_length')
            if s_array is None:
                print(" No se encontró arc_length")
                continue

            s = numpy_support.vtk_to_numpy(s_array)

            # velocidad
            u_array = data.GetPointData().GetArray('u')
            if u_array is None:
                print(" No se encontró 'u'")
                continue

            u_vec = numpy_support.vtk_to_numpy(u_array)
            u = u_vec[:, 2]   # componente X

            # guardar X solo una vez
            if X is None:
                X = s - 1.0

            U.append(u)

            Delete(line)

        # =========================
        # CONVERTIR A MATRIZ
        # =========================
        U = np.array(U)

        # crear grid
        X_grid, Y_grid = np.meshgrid(X, y_vals)

        # =========================
        # PLOT ISOLÍNEAS
        # =========================
        plt.figure(figsize=(6,6))

        U = np.nan_to_num(U)  # limpia NaN

        u_min = np.min(U)
        u_max = np.max(U)

        # evitar caso degenerado
        if abs(u_max - u_min) < 1e-12:
            print(" Campo constante, se omite")
            print(u_vec[:10])
            continue

        levels = np.linspace(u_min, u_max, 20)
        levels = np.sort(levels)

        plt.contour(X_grid, Y_grid, U, levels=levels, colors='black')

        plt.xlabel("x")
        plt.ylabel("y")
        plt.title(file.replace(".pvtu",""))

        plt.grid()

        # guardar
        name = file.replace(".pvtu","")
        output = os.path.join(outdir, name + "_Isolines_fromLines.png")

        plt.savefig(output, dpi=150)
        plt.close()

        print(f" Guardado: {output}")

        Delete(mesh)

    print("Proceso Completado")


if option == "3" or option == "10":

    # COMPARATIVA
    # =========================
    # FIGURA
    # =========================
    plt.figure(figsize=(6,6))

    # =========================
    # ORDENAR POR Ha
    # =========================
    def get_Ha(file):
        match = re.search(r"ha(\d+)", file)
        return int(match.group(1)) if match else 9999

    files = sorted([f for f in os.listdir(datadir) if f.endswith(".pvtu")], key=get_Ha)

    # =========================
    # LOOP
    # =========================
    for file in files:

        filename = os.path.join(datadir, file)
        print(f">> {file}")

        mesh = XMLPartitionedUnstructuredGridReader(
            registrationName=file,
            FileName=[filename]
        )

        mesh.PointArrayStatus = ['u']
        mesh.TimeArray = 'None'

        SetActiveSource(mesh)

        # =========================
        # LÍNEA VERTICAL
        # =========================
        line = PlotOverLine(Input=mesh)
        line.Point1 = [x1, y1, z1]
        line.Point2 = [x2, y2, z2]
        line.SamplingPattern = 'Sample At Segment Centers'

        data = servermanager.Fetch(line)

        # coordenada
        s = numpy_support.vtk_to_numpy(
            data.GetPointData().GetArray('arc_length')
        )

        # velocidad
        u_vec = numpy_support.vtk_to_numpy(
            data.GetPointData().GetArray('u')
        )

        # componente correcta (Z en tu caso)
        u = u_vec[:, 2]

        # normalizar eje y
        y = s

        # extraer Ha
        match = re.search(r"ha(\d+)", file)
        Ha = match.group(1) if match else "?"

        # =========================
        # PLOT
        # =========================
        plt.plot(y, u, label=f"Ha={Ha}")

        Delete(line)
        Delete(mesh)

    # =========================
    # ESTILO
    # =========================
    plt.xlabel("y")
    plt.ylabel("u")
    plt.title("Perfil de velocidad")

    plt.legend()
    plt.grid()

    # guardar
    output = os.path.join(outdir, "Perfiles_velocidad.png")
    plt.savefig(output, dpi=150)
    plt.close()

    print(f"\n Guardado: {output}")

if option == "4" or option == "10":

    print("\nCalculando fuerza de Lorentz ADIMENSIONAL (una gráfica por archivo)\n")

    # =========================
    # FUNCIÓN Ha
    # =========================
    def extract_Ha(file):
        match = re.search(r"ha(\d+)", file)
        return int(match.group(1)) if match else None

    # =========================
    # FUNCIÓN TEÓRICA
    # =========================
    def lorentz_force(y, Ha, u_hat=1.0):
        y = np.array(y)

        # forma estable numéricamente
        exp1 = np.exp(Ha * (y - 1))
        exp2 = np.exp(-Ha * (y + 1))

        term = (exp1 + exp2) / (1 + np.exp(-2 * Ha))

        return -1 + u_hat * Ha**2 * term

    # =========================
    # LISTA DE ARCHIVOS
    # =========================
    files = [f for f in os.listdir(datadir) if f.endswith(".pvtu")]
    files = sorted(files, key=extract_Ha)

    # =========================
    # LOOP
    # =========================
    for file in files:

        Ha = extract_Ha(file)

        if Ha is None:
            print(f"No se pudo extraer Ha de {file}")
            continue

        print(f">> {file}  (Ha={Ha})")

        # dominio
        # malla refinada en paredes
        y_core = np.linspace(-0.95, 0.95, 2000)
        y_wall1 = np.linspace(-1, -0.95, 4000)
        y_wall2 = np.linspace(0.95, 1, 4000)

        y = np.concatenate([y_wall1, y_core, y_wall2])

        # fuerza
        f = lorentz_force(y, Ha)

        # CLIPPING para ver signos
        f = np.clip(f, -1, 1)

        # =========================
        # FIGURA NUEVA
        # =========================
        plt.figure(figsize=(6,6))

        plt.plot(y, f, color='red',linewidth=2, label=f"Ha={Ha}")

        plt.xlabel("y")
        plt.ylabel("f(y)")
        plt.title(f"Fuerza de Lorentz (Ha={Ha})")

        plt.xlim([-1,1])
        plt.ylim([-1.2,1.2])

        plt.legend()
        plt.grid()

        # =========================
        # GUARDAR
        # =========================
        name = file.replace(".pvtu", "")
        output = os.path.join(outdir, name + "_Lorentz.png")

        plt.savefig(output, dpi=150)
        plt.close()

        print(f" Guardado: {output}")

    print("\n PROCESO COMPLETADO")
if option == "5" or option == "10":

    print("\nCalculando Q y K para todos los casos...\n")

    import pandas as pd

    # =========================
    # FUNCIÓN Ha
    # =========================
    def extract_Ha(file):
        match = re.search(r"ha(\d+)", file)
        return int(match.group(1)) if match else None

    # =========================
    # LISTA DE ARCHIVOS
    # =========================
    files = [f for f in os.listdir(datadir) if f.endswith(".pvtu")]
    files = sorted(files, key=extract_Ha)

    results = []

    # =========================
    # LOOP
    # =========================
    for file in files:

        filename = os.path.join(datadir, file)
        Ha = extract_Ha(file)

        if Ha is None:
            print(f"No se pudo extraer Ha de {file}")
            continue

        print(f">> {file} (Ha={Ha})")

        # =========================
        # LEER MESH
        # =========================
        mesh = XMLPartitionedUnstructuredGridReader(
            registrationName=file,
            FileName=[filename]
        )

        mesh.PointArrayStatus = ['u']
        mesh.TimeArray = 'None'
        SetActiveSource(mesh)

        # =========================
        # PERFIL EN Y (CENTRO)
        # =========================
        line = PlotOverLine(Input=mesh)
        line.Point1 = [0.0, -1.0, 0.0]
        line.Point2 = [0.0,  1.0, 0.0]
        line.SamplingPattern = 'Sample Uniformly'
        line.Resolution = 500

        data = servermanager.Fetch(line)

        # =========================
        # EXTRAER DATOS
        # =========================
        points = numpy_support.vtk_to_numpy(data.GetPoints().GetData())
        y = points[:, 1]   # eje y

        u_vec = numpy_support.vtk_to_numpy(
            data.GetPointData().GetArray('u')
        )

        u = u_vec[:, 2]  # componente del flujo

        # ordenar
        idx = np.argsort(y)
        y = y[idx]
        u = u[idx]

        # =========================
        # INTEGRAL DE LA VELOCIDAD · dy → Q
        # =========================
        Q = np.trapz(u, y)

        # =========================
        # CALCULAR K
        # =========================
        if abs(Q) < 1e-12:
            print(" Q ~ 0, se omite")
            continue

        K = 2.0 / Q
        HaL= 1/Ha
        # guardar resultados
        results.append({
            "file": file,
            "Ha": Ha,
            "Q": Q,
            "K": K,
            "Hartmannlayer": HaL,
        })

        Delete(line)
        Delete(mesh)

    # =========================
    # CONVERTIR A DATAFRAME
    # =========================
    df = pd.DataFrame(results)

    # ordenar por Ha
    df = df.sort_values(by="Ha")

    print("\nResultados:")
    print(df)

    # =========================
    # GUARDAR A CSV(se puede abrir con excel)
    # =========================
    output_csv = os.path.join(outdir, "Resultados_Q_K.csv")
    df.to_csv(output_csv, index=False)

    print(f"\n CSV guardado en: {output_csv}")

    # =========================
    # PLOT K vs Ha
    # =========================
    plt.figure(figsize=(6,6))

    plt.plot(df["Ha"], df["K"], 'o-', linewidth=2)

    plt.xscale("log")
    plt.yscale("log")

    plt.xlabel("Ha")
    plt.ylabel("K")
    plt.title("K vs Ha")

    plt.grid()

    output_plot = os.path.join(outdir, "K_vs_Ha.png")
    plt.savefig(output_plot, dpi=150)
    plt.close()

    print(f" Gráfica guardada en: {output_plot}")

if option == "6" or option == "10":

    print("\nGenerando perfil 3D \n")

    # =========================
    # FUNCIÓN Ha
    # =========================
    def extract_Ha(file):
        match = re.search(r"ha(\d+)", file)
        return int(match.group(1)) if match else None

    # =========================
    # PARÁMETROS MALLA
    # =========================
    nx = 20   # bajar si peta
    ny = 20

    x_vals = np.linspace(-1, 1, nx)
    y_vals = np.linspace(-1, 1, ny)

    # =========================
    # ARCHIVOS
    # =========================
    files = [f for f in os.listdir(datadir) if f.endswith(".pvtu")]
    files = sorted(files, key=extract_Ha)

    # =========================
    # LOOP ARCHIVOS
    # =========================
    for file in files:

        print(f">> {file}")

        filename = os.path.join(datadir, file)

        mesh = XMLPartitionedUnstructuredGridReader(
            registrationName=file,
            FileName=[filename]
        )

        mesh.PointArrayStatus = ['u']
        mesh.TimeArray = 'None'
        SetActiveSource(mesh)

        # =========================
        # MATRIZ U(x,y)
        # =========================
        U = np.zeros((ny, nx))

        # =========================
        # MUESTREO PUNTO A PUNTO
        # =========================
        for i, y in enumerate(y_vals):
            for j, x in enumerate(x_vals):

                probe = ProbeLocation(Input=mesh)
                probe.ProbeType = 'Fixed Radius Point Source'
                probe.ProbeType.Center = [x, y, 0.0]

                data = servermanager.Fetch(probe)

                u_array = data.GetPointData().GetArray('u')

                if u_array is None:
                    U[i, j] = 0.0
                else:
                    u_vec = numpy_support.vtk_to_numpy(u_array)
                    U[i, j] = u_vec[0, 2]  # componente z

                Delete(probe)

        Delete(mesh)

        # =========================
        # GRID
        # =========================
        X, Y = np.meshgrid(x_vals, y_vals)

        # =========================
        # PLOT 3D
        # =========================
        fig = plt.figure(figsize=(8,6))
        ax = fig.add_subplot(111, projection='3d')

        ax.plot_surface(X, Y, U, cmap='viridis')

        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('u')

        Ha = extract_Ha(file)
        ax.set_title(f"Perfil 3D de velocidad (Ha={Ha})")

        output = os.path.join(outdir, file.replace(".pvtu","_3D_xy.png"))
        plt.savefig(output, dpi=150)
        plt.close()

        print(f" Guardado: {output}")
if option == "7" or option == "10":

    print("\nGenerando streamlines...\n")

    files = [f for f in os.listdir(datadir) if f.endswith(".pvtu")]
    files = sorted(files, key=extract_Ha)

    for file in files:

        print(f">> {file}")

        filename = os.path.join(datadir, file)

        mesh = XMLPartitionedUnstructuredGridReader(
            registrationName=file,
            FileName=[filename]
        )

        mesh.PointArrayStatus = ['j']   # campo vectorial
        mesh.TimeArray = 'None'
        SetActiveSource(mesh)

        # =========================
        # STREAM TRACER
        # =========================
        stream = StreamTracer(Input=mesh,
                             SeedType='Point Cloud')

        # campo vectorial
        stream.Vectors = ['POINTS', 'j']


        stream.SeedType.Center = [0.0, 0.0, 0.0]
        stream.SeedType.Radius = 1
        stream.SeedType.NumberOfPoints = 300

        # longitud de integración
        stream.MaximumStreamlineLength = 5.0

        # =========================
        # FETCH
        # =========================
        data = servermanager.Fetch(stream)

        points = numpy_support.vtk_to_numpy(
            data.GetPoints().GetData()
        )
        j_vec = numpy_support.vtk_to_numpy(
            data.GetPointData().GetArray('j')
        )

        j_mag = np.linalg.norm(j_vec, axis=1)

        # =========================
        # PLOT (2D proyección)
        # =========================
        plt.figure(figsize=(6,6))

        x = points[:,0]
        y = points[:,1]

        plt.scatter(x, y, c=j_mag, s=2)
        plt.colorbar(label='|j|')

        plt.xlabel("x")
        plt.ylabel("y")
        plt.title(f"Streamlines j (Ha={extract_Ha(file)})")

        plt.grid()

        output = os.path.join(outdir, file.replace(".pvtu","_streamlines.png"))
        plt.savefig(output, dpi=150)
        plt.close()

        Delete(stream)
        Delete(mesh)

        print(f" Guardado: {output}")


# In[ ]:




