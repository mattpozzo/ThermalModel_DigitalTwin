# -*- coding: utf-8 -*-
############## Thermal Model Block ##########################
  # Coded by Universitwin - Digital Twin Team
  # April 2020
  # Started by: Nicolás Valentín Conde
  # Version: V02
  # Latest version by:
  # Latest version date:

# Libraries...
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# Functions...

import funciones_termico as f

# ToDo actual:  - Agregar nodos de contorno a .csv
# ToDo          - Ver cómo determinar la carga solar para cada distancia sat-sol; implementar en cálculo de flujo incidente

def mainT(t: float, dt: float, alb: float, nad_rad: float, nad_sun: float, onoff: bool, umbra: bool, T=None, matrix_e=None, matrix_G=None, matrix_q=None, mass_Cp=None, Nodo=None):
    # t:            Tiempo actual de simulación [sec]
    # dt:           Paso del integrador [sec]
    # iterations:   Cantidad de puntos en una órbita [#]
    # orb:          Cantidad de órbitas a simular [#]
    # alb:          Factor de albedo del paso actual [ ]
    # nad_rad:      Ángulo entre nadir y vector radial a la Tierra [rad]
    # nad_sun:      Ángulo entre nadir y vector solar [rad]
    # onoff:        Lista de bools de largo num, donde True indica un nodo operativo, y False uno no operativo
    # umbra:        True dentro de la umbra de la órbita, False fuera de la umbra
    # T:            Array de temperaturas nodales del paso actual [K]
    # matrix_e:     Lista de emisividades nodales [ ]
    # matrix_G:     Matriz de conductancias entre nodos [W/mK]
    # matrix_q:     Lista de absortividades nodales [ ]
    # mass_Cp:      Lista de masas térmicas nodales [J/K]
    # Nodo:         Lista de nodos (con sus características)

    # Input

    # Input 2: Constants
    mu_earth = 398600.4415  # [km^3/sec^2] Parametro gravitacional
    R_earth = 6378.1363  # [km] Radio de la Tierra
    T_earth = 288.05  # [K] Temperatura promedio de la Tierra
    params = (mu_earth,)
    stef_boltz = 5.670373e-8  # [W m^-2 K^-4] Stefan-Boltzman

    # Ver en funciones_termico cómo es la clase Nodos. Cada nodo además de tener un ID,
    # tiene asignada una parte, por ejemplo "Acople Al" o "PCB Nadir". Esto es para que
    # cuando expandamos la cantidad de nodos, podamos seleccionar todos los nodos de una
    # pieza para algo en particular (no se qué todavía).
    # También, algunos tienen asignados componentes, para luego en el concepto de operaciones
    # se puedan prender o apagar estos nodos según el componente que tengan.

    if Nodo is None:
        # Input 1: Thermal Properties of materials
        df = pd.read_csv('4.8 CSVs/physical_properties.csv')
        # df.drop(['Unnamed: 0'], axis=1, inplace=True)

        Nodo = []
        for i in range(0, df.index[-1] + 1):
            Nodo.append(f.Nodos(df[i:i + 1].values[0]))

    num = len(Nodo)  # Number of nodes

    if T is None:
        T = np.ones(num) * 273.15  # Initial temperature

    if None in [matrix_e,matrix_G,matrix_q,mass_Cp]:
        # Step 2: Define matrix from equation
        matrix_e, matrix_G, matrix_q, mass_Cp = f.matrices(Nodo)

    """
    # Input 3: Simulator Settings (These settings will be passed by the master algorithm)
    tf = per * 60  # [sec] Tiempo final de simulacion
    iterations = int((tf - t0) / dt) + 1  # Numero de puntos por orbita
    orbits = orb  # Cantidad de órbitas
    """

    # Algorithm

    # T_list = np.zeros(16)  # Lista de evolución de temperaturas para graficar.
    cargas_externas, conduccion, radiacion, carga_interna = 0, 0, 0, 0  # Inicializo las variables

    # Step 3: Convergence analysis
    conv = 0
    for i in range(num):
        if dt / (mass_Cp[i] * Nodo[i].A) > 0.6:     # Dato de Nahuel
            conv += 1
    if conv != 0:
        print('No converge')
    else:
        print('Converge')

    # La ecuación que se resuelve es:
    #     "calor acumulado en el nodo n = cargas_externas + conduccion + radiacion + carga_interna"
    # Dentro del for se va resolviendo cada uno de los términos y se itera hasta llegar a un valor de
    # temperatura estable que haga que el calor acumulado sea 0.

    for i in range(0, num):  # Se abre un loop que itera entre los N nodos
        if Nodo[i].type != 'C':     # La evolución de T de los nodos de contorno es un input

            # Step 1: Concepto de operaciones
            # La idea es que acá se vayan encendiendo y apagando los componentes según el modelo de
            # experimento (Concepto de operaciones) y vaya cambiando la carga interna.
            if onoff[i]:
                carga_interna = Nodo[i].Qop
            else:
                carga_interna = Nodo[i].Qnop

            # Step 2: Cargas externas
            # La idea es que acá esté el modelo de cargas externas, incluyendo el de albedo.
            # Hay que buscar la forma de automatizar el cálculo de sombras para una órbita.
            # Actualmente puse que en la mitad de la órbita haya sombra y en la otra mitad sol.
            # Y que la radiación solo le da a los nodos externos.
            if Nodo[i].ext:
                carga_inf = matrix_e[i] * Nodo[i].A / 1000 ** 2 * stef_boltz * T_earth ** 4 # ToDo: AGREGAR F VISTA
                carga_sun = 0
                carga_alb = 0
                if not umbra:
                    pass
                    # carga_sun = matrix_q[i] * Nodo[i].A / 1000 ** 2 * sun_q[j]        # ToDo: AGREGAR F VISTA e implementar sun_q
                    # carga_alb = matrix_q[i] * Nodo[i].A / 1000 ** 2 * sun_q[j] * alb  # ToDo: AGREGAR F VISTA e implementar sun_q
                cargas_externas = carga_sun + carga_alb + carga_inf

            # Step 3: Cálculo del calor por conducción
            # Se calcula el balance entre el calor transmitido a otros nodos y el recibido.
            conduccion = (matrix_G[i] * T).sum() - matrix_G[i].sum() * T[i]

            # Step 4: Cálculo del calor por radiación ToDo: Todo básicamente
            # Se calcula el balance entre la radiación emitida y la recibida del satélite y otros nodos.
            # Actualmente solo se libera radiación al ambiente como si estuviera cada nodo solo en el
            # espacio.
            # Hay que agregar un modelo de radiación entre nodos según geometría. Y agregar un nodo imaginario
            # que sea el satélite.
            radiacion = - stef_boltz * matrix_e[i] * Nodo[i].A * T[i] ** 4 / 1000 ** 2

            # Step 5: Calcular el acumulado.
            acumulado_n = cargas_externas + conduccion + radiacion + carga_interna

            """
            # Step 6: Condiciones de contorno
            # En nuestro caso, el acople de aluminio va a estar en contacto con el radiador, asique
            # va a  tener su temperatura. Como está modelado como un nodo solo, la temperatura va a
            # ser constante en toda la pieza.
            # Hay que ver la forma de automatizar esto, que puedas ingresar qué nodo es constante como
            # input, porque sino cada vez que haya que agregar una condición de contorno va a ser una paja.
            if Nodo[i].part == 'Acople SC':                         # ToDo. Hacer este if al principio del loop[i]
                T[i] = T_SC[j] + 273      # Temperatura del spacecraft    # ToDo. para no operar al pedo.
            elif Nodo[i].part == 'Acople Rad':
                T[i] = T_rad[j] + 273    # Temperatura del radiador
            else:
            """

            # Step 6: Cálculo del cambio de temperatura para el nodo.
            T[i] += acumulado_n * dt / mass_Cp[i]  # Temp del nodo i para una posición orbital j.

"""
        # Concateno la temperatura de los N nodos para tiempo j, en la lista.
        T_list = np.concatenate((T_list, T), axis=0)

    # Hago un gráfico de la lista de temperaturas en el tiempo.
    f.grafics_TvsNu((int(len(T_list) / 16), 16), T_list - 273)

    tupla = (int(len(T_list) / 16), 16)
    T_list = T_list.reshape(tupla)

    # import seaborn as sns
    # sns.set(style="ticks", context="talk")
    # plt.style.use("dark_background")
    # blue, = sns.color_palette("muted", 1)
    # fig, ax = plt.subplots(figsize=(20,10))
    # ax.plot(T_list[1:,8]-273, color=blue, lw=3)
    # plt.title('Temperatura PCB Payload', fontsize=20)
    # plt.xlabel("time [sec]", labelpad=15, fontsize=20, color="#333533")
    # plt.ylabel("T [°C]", labelpad=15, fontsize=20, color="#333533")

    np.linspace(0, 180, int(iterations / 2))

    nu = np.linspace(0, 360, int(iterations))
    plt.plot(nu, T_list[1:, :] - 273)
    plt.xlabel("v [deg]", labelpad=15, fontsize=12, color="#333533")
    plt.ylabel("T [°C]", labelpad=15, fontsize=12, color="#333533")
    plt.title('T of Nodes vs True Anomally')
    plt.grid(True, color="#93a1a1", alpha=0.3)

    plt.show()
"""

# Algorithm...

# Output...

# Visualization
  # Input 1 visualization

########## ##     ####     ######
    ##     ##    ##  ##    ##   ##
    ##     ##   ##    ##   ##   ##
    ##     ##   ########   ######
    ##     ##  ##      ##  ##   ## 
    ##     ##  ##      ##  ##   ##
    ##     ##  ##      ##  ######
