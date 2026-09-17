#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from paraview.simple import *
from vtk.util import numpy_support

import numpy as np
import matplotlib.pyplot as plt
import os
from mpl_toolkits.mplot3d import Axes3D

# =====================================================
# CONFIGURACIÓN
# =====================================================

datadir = os.getcwd()

sections = [-0.5, 0.0, 0.5]

tol = 0.01

# =====================================================
# BUSCAR ARCHIVOS
# =====================================================

files = sorted(
    [f for f in os.listdir(datadir)
     if f.endswith(".pvtu")]
)

# =====================================================
# LOOP PRINCIPAL
# =====================================================

for filename in files:

    print("Procesando:", filename)

    mesh = XMLPartitionedUnstructuredGridReader(
        registrationName=filename,
        FileName=[os.path.join(datadir, filename)]
    )

    data = servermanager.Fetch(mesh)

    # ==========================================
    # COORDENADAS
    # ==========================================

    points = numpy_support.vtk_to_numpy(
        data.GetPoints().GetData()
    )

    x = points[:,0]
    y = points[:,1]
    z = points[:,2]

    # ==========================================
    # VELOCIDAD
    # ==========================================

    u_vec = numpy_support.vtk_to_numpy(
        data.GetPointData().GetArray("u")
    )

    ux = u_vec[:,0]

    # ==========================================
    # CARPETA SALIDA
    # ==========================================

    outdir = os.path.join(
        datadir,
        "postprocess_3D"
    )

    os.makedirs(
        outdir,
        exist_ok=True
    )

    # ==========================================
    # SECCIONES
    # ==========================================

    for xsec in sections:

        mask = np.abs(x - xsec) < tol

        if np.sum(mask) == 0:

            print(
                f"No hay puntos para x={xsec}"
            )
            continue

        y_sec = y[mask]
        z_sec = z[mask]
        ux_sec = ux[mask]

        # ======================================
        # FIGURA
        # ======================================

        fig = plt.figure(
            figsize=(8,6)
        )

        ax = fig.add_subplot(
            111,
            projection='3d'
        )

        sc = ax.scatter(
            y_sec,
            z_sec,
            ux_sec,
            c=ux_sec,
            cmap='jet',
            s=8
        )

        plt.colorbar(
            sc,
            ax=ax,
            label='Ux'
        )

        ax.set_xlabel("Y")
        ax.set_ylabel("Z")
        ax.set_zlabel("Ux")

        ax.set_title(
            f"{filename}\n x={xsec}"
        )

        plt.tight_layout()

        output = os.path.join(
            outdir,
            filename.replace(
                ".pvtu",
                f"_x{xsec}.png"
            )
        )

        plt.savefig(
            output,
            dpi=300
        )

        plt.close()

        print(
            "Guardado:",
            output
        )

print("Terminado")

