#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from paraview.simple import *
from vtk.util import numpy_support
import matplotlib.pyplot as plt
import numpy as np
import os
import re

paraview.simple._DisableFirstRenderCameraReset()

N_values = {
    1: "33.6",
    2: "34",
    3: "65.5",
    4: "127",
    5: "249",
    6: "250",
    7: "504",
    8: "917",
    9: "989",
    10: "1076",
    11: "1888",
    12: "3740"
}

# =========================================================
# RUTAS
# =========================================================
datadir = r"C:\Users\joanc\Desktop\TFG\Modelos\data\data"
agrdir  = r"C:\Users\joanc\Desktop\TFG\Modelos\Expansion Experiment Figures 19-23"
outdir  = os.path.join(datadir, "postprocess_comparison")

os.makedirs(outdir, exist_ok=True)

# =========================================================
# EXTRAER Ha Y Re
# =========================================================
def get_Ha_Re(file):

    match = re.search(r"Ha_(\d+)_Re_(\d+)", file)

    if match:
        return float(match.group(1)), float(match.group(2))

    return None, None

# =========================================================
# MAPA AGR
# =========================================================
def build_agr_map():

    agr_map = {}

    for f in os.listdir(agrdir):

        if not f.endswith(".agr"):
            continue

        match = re.search(r"Ha[_ ]?(\d+)", f)

        if match:

            Ha_file = float(match.group(1))

            agr_map[Ha_file] = os.path.join(
                agrdir,
                f
            )

    return agr_map

agr_map = build_agr_map()

print("\nAGR detectados:")

for k in agr_map:
    print(f"Ha={k} -> {agr_map[k]}")

# =========================================================
# LEER MALLA
# =========================================================
def read_mesh(filename):

    mesh = XMLPartitionedUnstructuredGridReader(
        registrationName=filename,
        FileName=[filename]
    )

    mesh.PointArrayStatus = ['uh', 'ph']
    mesh.TimeArray = 'None'

    SetActiveSource(mesh)

    return mesh

# =========================================================
# LINEA CENTRAL
# =========================================================
def read_line(mesh):

    line = PlotOverLine(Input=mesh)

    line.Point1 = [-8.0, 0.0, 0.0]
    line.Point2 = [ 8.0, 0.0, 0.0]

    line.SamplingPattern = 'Sample Uniformly'
    line.Resolution = 1000

    return line

# =========================================================
# EXTRAER ARRAYS
# =========================================================
def get_array(data, name):

    return numpy_support.vtk_to_numpy(
        data.GetPointData().GetArray(name)
    )

# =========================================================
# LEER S1 ... S13
# =========================================================
def read_agr_all(filename):

    datasets = {}

    current_set = None
    reading = False

    with open(filename, 'r') as f:

        for line in f:

            line = line.strip()

            match = re.match(
                r"@target G0\.S(\d+)",
                line
            )

            if match:

                s_num = int(match.group(1))

                # ahora también S13
                if 1 <= s_num <= 13:

                    current_set = s_num

                    datasets[current_set] = {
                        "x": [],
                        "y": []
                    }

                    reading = True

                else:

                    reading = False

                continue

            if reading:

                if line == "&":

                    reading = False
                    current_set = None
                    continue

                if line.startswith("@") or line.startswith("#"):
                    continue

                parts = line.split()

                if len(parts) >= 2:

                    try:

                        datasets[current_set]["x"].append(
                            float(parts[0])
                        )

                        datasets[current_set]["y"].append(
                            float(parts[1])
                        )

                    except:
                        pass

    return datasets

# =========================================================
# LOOP PRINCIPAL
# =========================================================
files = [
    f for f in os.listdir(datadir)
    if f.endswith(".pvtu")
]

for file in files:

    print("\n====================================")
    print(f">> {file}")

    Ha, Re = get_Ha_Re(file)

    if Ha is None:
        continue

    print(f"Ha={Ha}, Re={Re}")

    # =====================================================
    # LEER SIMULACIÓN
    # =====================================================
    filename = os.path.join(
        datadir,
        file
    )

    mesh = read_mesh(filename)

    line = read_line(mesh)

    data = servermanager.Fetch(line)

    s = get_array(data, 'arc_length')
    p = get_array(data, 'ph')

    x = s - 8.0

    # =====================================================
    # REFERENCIA DE PRESIÓN
    # =====================================================
    idx_ref = np.argmin(
        np.abs(x + 6.38)
    )

    p_ref = p[idx_ref]

    p_scaled = p - p_ref

    # =====================================================
    # PLOT
    # =====================================================
    plt.figure(figsize=(8,5))

    plt.plot(
        x,
        p_scaled,
        linewidth=2,
        color='black',
        label='Simulation'
    )

    # =====================================================
    # AGR DEL MISMO Ha
    # =====================================================
    agr_file = agr_map.get(
        Ha,
        None
    )

    if agr_file is not None:

        print(
            f"AGR encontrado: {agr_file}"
        )

        datasets = read_agr_all(
            agr_file
        )
        first_scatter = True
        for s_num, data_ref in datasets.items():

            x_ref = np.array(
                data_ref["x"]
            )

            p_ref_curve = np.array(
                data_ref["y"]
            )

            if len(x_ref) == 0:
                continue

            idx_ref_agr = np.argmin(
                np.abs(
                    x_ref + 6.38
                )
            )

            p_ref_curve = (
                p_ref_curve
                - p_ref_curve[idx_ref_agr]
            )

            # ======================================
            # S13 = curva N infinito
            # ======================================
            if s_num == 13:

                plt.plot(
                    x_ref,
                    p_ref_curve,
                    '--',
                    color='red',
                    linewidth=2,
                    label='N = ∞'
                )

            # ======================================
            # S1-S12 = datos experimentales
            # ======================================
            else:

                colors = [
                    'blue',
                    'green',
                    'orange',
                    'purple',
                    'brown',
                    'magenta',
                    'cyan',
                    'olive',
                    'gold',
                    'gray',
                    'pink',
                    'teal'
                ]

                idx = s_num - 1

                plt.scatter(
                    x_ref,
                    p_ref_curve,
                    s=25,
                    color=colors[idx],
                    label=f"N={N_values[s_num]}"
                )

        # =====================================================
        # ESTÉTICA
        # =====================================================
        plt.xlabel("x")
        plt.ylabel("Pressure")

        plt.title(
            f"Ha={Ha}, Re={Re}"
        )

        plt.grid()

        plt.legend(
            ncol=3,
            fontsize=8,
            loc='best'
        )

    # =====================================================
    # GUARDAR
    # =====================================================
    name = file.replace(
        ".pvtu",
        ""
    )

    output = os.path.join(
        outdir,
        name + "_comparison.png"
    )

    plt.savefig(
        output,
        dpi=150,
        bbox_inches='tight'
    )

    plt.close()

    print(
        f"Guardado: {output}"
    )

    Delete(mesh)
    Delete(line)

print("\nDONE")

