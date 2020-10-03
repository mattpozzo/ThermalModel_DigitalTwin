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

# Input

# Input 1: Thermal Properties of materials
df = pd.read_csv('../4.8 CSVs/physical_properties.csv')
df.drop(['Unnamed: 0'], axis=1, inplace=True)
df

# Input 2: Constants
mu_earth = 398600.4415              # [km^3/sec^2] Parametro gravitacional
R_earth = 6378.1363                 # [km] Radio de la Tierra
params = (mu_earth,)
stef_boltz = 5.670373e-8            # [W m^-2 K^-4] Stefan-Boltzman

# Input 3: Simulator Settings (These settings will be passed by the master algorithm)
t0 = 0                              # [sec] Tiempo inicial de simulacion
dt = 1                              # [sec] Paso del integrador
tf = 94.2 * 60                      # [sec] Tiempo final de simulacion  REEMPLAZAR 94.2 POR VAR PERIODO
iterations = int((tf - t0)/dt) + 1  # Numero de puntos por orbita
orbits = 1                          # Cantidad de órbitas

# Input 4: Initial conditions
T = np.ones(16)*270                 # Initial temperature

# Algorithm

# Step 1: Nodes structure definition
Nodo = []
for i in range(0,df.index[-1]+1):
    Nodo.append(f.Nodos( df[i:i+1].values[0] ) )

num = len(Nodo) # Number of nodes

# Ver en funciones_termico cómo es la clase Nodos. Cada nodo además de tener un ID,
# tiene asignada una parte, por ejemplo "Acople Al" o "PCB Nadir". Esto es para que
# cuando expandamos la cantidad de nodos, podamos seleccionar todos los nodos de una
# pieza para algo en particular (no se qué todavía).
# También, algunos tienen asignados componentes, para luego en el concepto de operaciones
# se puedan prender o apagar estos nodos según el componente que tengan.

# Step 2: Define matrix from equation
emisividad_node, matrix_G, matrix_q, mass_Cp = f.matrices(Nodo)


#Step 3: Convergence analysis
converge = 1
for i in range(num):
    if dt / (mass_Cp[i] * Nodo[i].A) > 0.6:
        converge = 0
    #print(D_t / (mass_Cp[i] * Nodo[i].A))
        
if converge == 0:
    print('No Converge')
else:
    print('Converge') 
del converge

T_list = np.zeros(16) # Lista de evolución de temperaturas para graficar.
cargas_externas, conduccion, radiacion, carga_interna = 0,0,0,0 # Inicializo las variables

# La ecuación que se resuelve es:
#     "calor acumulado en el nodo n = cargas_externas + conduccion + radiacion + carga_interna"
# Dentro del for se va resolviendo cada uno de los términos y se itera hasta llegar a un valor de
# temperatura estable que haga que el calor acumulado sea 0.

for j in range(0, iterations*orbits): # Loop para cada instante de tiempo.
    for i in range(0,int(num/2)): # Se abre otro loop para que itera entre los N nodos.
        
        # Step 1: Concepto de operaciones
            # La idea es que acá se vayan encendiendo y apagando los componentes según el modelo de
            # experimento (Concepto de operaciones) y vaya cambiando la carga interna.
        if j/iterations-int(j/iterations)<0.5:
            carga_interna = Nodo[i].Qop
        else:
            carga_interna = Nodo[i].Qnop
        
        # Step 2: Cargas externas
            # La idea es que acá esté el modelo de cargas externas, incluyendo el de albedo. 
            # Hay que buscar la forma de automatizar el cálculo de sombras para una órbita.
            # Actualmente puse que en la mitad de la órbita haya sombra y en la otra mitad sol.
            # Y que la radiación solo le da al nodo "Payload".
        if j/iterations-int(j/iterations)<0.5:
            if Nodo[i].part == 'Payload':
                cargas_externas = matrix_q[i] * 1007.1377763/2
            else:
                cargas_externas = 0
        else:
            cargas_externas = 0
        
        # Step 3: Cálculo del calor por conducción
            # Se calcula el balance entre el calor transmitido a otros nodos y el recibido.
        conduccion = (matrix_G[i] * T).sum() - matrix_G[i].sum() * T[i]
        
        # Step 4: Cálculo del calor por radiación
            # Se calcula el balance entre la radiación emitida y la recibida del satélite y otros nodos.
            # Actualmente solo se libera radiación al ambiente como si estuviera cada nodo solo en el 
            # espacio.
            # Hay que agregar un modelo de radiación entre nodos según geometría. Y agregar un nodo imaginario
            # que sea el satélite.
        radiacion = - stef_boltz * emisividad_node[i] * Nodo[i].A * T[i]**4/1000**2
        
        # Step 5: Calcular el acumulado.
        acumulado_n = cargas_externas + conduccion + radiacion + carga_interna
        
        # Step 6: Condiciones de contorno
            # En nuestro caso, el acople de aluminio va a estar en contacto con el radiador, asique
            # va a  tener su temperatura. Como está modelado como un nodo solo, la temperatura va a
            # ser constante en toda la pieza.
            # Hay que ver la forma de automatizar esto, que puedas ingresar qué nodo es constante como
            # input, porque sino cada vez que haya que agregar una condición de contorno va a ser una paja.
        if Nodo[i].part == 'Acople Al':
            T[i] = 3+273    # ESTá BIEN, ES LA T CONTORNO (DEL SAT)
        else:
            # Step 6: Cálculo del cambio de temperatura para el nodo.
            T[i] += acumulado_n * dt / mass_Cp[i] # Temp del nodo i para una posición orbital j.
    
    # Concateno la temperatura de los N nodos para tiempo j, en la lista.
    T_list = np.concatenate((T_list, T), axis=0)

# Hago un gráfico de la lista de temperaturas en el tiempo.
f.grafics_TvsNu((int(len(T_list)/16),16), T_list-273)

tupla = (int(len(T_list)/16),16)
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

np.linspace(0, 180, int(iterations/2))

nu = np.linspace(0, 360, int(iterations))
plt.plot(nu, T_list[1:,:]-273)
plt.xlabel("v [deg]", labelpad=15, fontsize=12, color="#333533")
plt.ylabel("T [°C]", labelpad=15, fontsize=12, color="#333533")
plt.title('T of Nodes vs True Anomally')
plt.grid(True, color="#93a1a1", alpha=0.3)

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