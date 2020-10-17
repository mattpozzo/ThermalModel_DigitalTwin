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

# Functions...
class Nodos:
  def __init__(self, node):    
      
    # Properties of node
    self.ID = node[0]                                 # Identification
    self.type = node[1]                               # Arithmetic or Diffusion
    self.x = float(node[2])                           # metros
    self.y = float(node[3])                           # metros
    self.z = float(node[4])                           # metros
    self.part = node[5]                               # Part name
    self.espesor = float(node[6])                     # thickness
    self.A = float(node[7])                           # Area
    self.component = node[8]                          # Component name
    self.Qnop = float(node[9])                        # Generated heat - Non Operative
    self.Qop = float(node[10])                        # Generated heat - Operative
    self.P = node[11]                                 # power domain
    self.contact = np.array(node[12].split(',')).astype(float)   # Nodes in contact with this node
    
    # Properties of material
    self.material = node[13]                          # Material's name
    self.temp = node[14]                              # °C
    self.density = float(node[15])                    # Density in kg/m3
    self.Q = float(node[16])                          # Heat Capacity J/kgK
    self.k = float(node[17])                          # Conductivity W/mK
    self.a = float(node[18])                          # Absorptivity
    self.e = float(node[19])                          # Emissivity
    
def matrices(Nodo):
    """
    Coded by Universitwin - Digital Twin Team
    September 2020
    Started by: Nicolás Valentín Conde
    Version: V02
    Latest version by:
    Latest version date:
        
    Function
    ---------
    Thermal coefficients matrix for each node.
    
    Parameters
    ---------
        Nodo: Nodes structure
    
    Returns
    ---------
        emisividad_node: Radiation emissivity for each node. Vector [1xn]
        G: Conductance matrix between nodes [nxn]
        matrix_q: Absorptivity matrix [1xn]
        mass_Cp: Heat Capacity por each node [1xn]
    """
    
    num = len(Nodo)
    
    emisividad_node = []
    matrix_G = []
    matrix_q = []
    mass_Cp = []
    
    for i in range(0,num):
        # Step 1: Emissivity vector. 
            # Select from node properties and append to vector.
        emisividad_node.append(Nodo[i].e)
        
        # Step 2: Absorptivity vector. 
            # Select from node properties absorption and Area in [m], multiply and append to vector.
        matrix_q.append(Nodo[i].a)
        
        # Step 3: Heat capacity vector. 
            # Select from node properties and solve: Mass [kg] * Specific heat [J kg^-1 K^-1]
        mass_Cp.append(Nodo[i].A*Nodo[i].density*Nodo[i].espesor/1000**3*Nodo[i].Q)
        
        # Step 4: Conductance Matrix. 
            # In node structure we have two types of node: Arithmetic or Diffusion
            # Arithmetic is used as intermediate step for calculation of conductance 
            # between two Diffusion nodes. It represents changes on material or change
            # in path.
            # The conductance matrix only show conductance between diffusion nodes.
        for j in range(0,num):
            if i == j:
                matrix_G.append(0)
            elif sum(Nodo[i].contact-1 == j) > 0:
                L = ((Nodo[i].x-Nodo[j].x)**2+(Nodo[i].y-Nodo[j].y)**2+(Nodo[i].z-Nodo[j].z)**2)**0.5/1000
                matrix_G.append((areas(Nodo, i, j, L)*max(Nodo[i].espesor,Nodo[j].espesor)*max(Nodo[i].k,Nodo[j].k))/(L*1000**3))
            else:
                matrix_G.append(0)
    
    matrix_G = np.array(matrix_G).reshape(num,num)
    G = np.zeros([16,16]) 
    for i in range(16,num):
        if Nodo[i].type == 'A':
            filter_arr = np.isfinite(Nodo[i].contact)
            virtual_node = i
            nodes_cont = np.sort(Nodo[i].contact[filter_arr]-1)
            G[int(nodes_cont[0]),int(nodes_cont[1])] = 1/(1/matrix_G[int(nodes_cont[1])][virtual_node] + 1/matrix_G[int(nodes_cont[0])][virtual_node])
            G[int(nodes_cont[1]), int(nodes_cont[0])] = G[int(nodes_cont[0]),int(nodes_cont[1])]
            if len(nodes_cont)>2:
                G[int(nodes_cont[1]),int(nodes_cont[2])] = 1/(1/matrix_G[int(nodes_cont[1])][virtual_node] + 1/matrix_G[int(nodes_cont[2])][virtual_node])
                G[int(nodes_cont[2]), int(nodes_cont[1])] = G[int(nodes_cont[1]),int(nodes_cont[2])]
    
    emisividad_node = np.array(emisividad_node)
    where_are_NaNs = np.isnan(emisividad_node)
    emisividad_node[where_are_NaNs] = 0
    
    matrix_q = np.array(matrix_q)
    where_are_NaNs = np.isnan(matrix_q)
    matrix_q[where_are_NaNs] = 0
    
    mass_Cp = np.array(mass_Cp)
    where_are_NaNs = np.isnan(mass_Cp)
    mass_Cp[where_are_NaNs] = 0
    
    return emisividad_node, G, matrix_q, mass_Cp

def areas(Nodo, i, j, L):
    A = max(Nodo[i].A, Nodo[j].A)
    e = max(Nodo[i].e, Nodo[j].e)
    if A > 100:
        return A*e/L
    else:
        return A

def grafics_TvsNu(tupla, T_list):
    
    num = tupla[1]
    
    T_list = T_list.reshape(tupla)

    #Gráfico de Evolución
    j = 0
    n = 0
    fig, ax = plt.subplots(int(num/5)+1, 5, sharex='col', sharey='row', figsize=(25,22))
    fig.subplots_adjust(hspace=0.2, wspace=0.2)
    for i in range(0,num):
        n = int(i/5)

        if j > 4:
            j = 0

        ax[n, j].plot(T_list[1:,i], color="#52BD42")

        ax[n, j].set_xlabel("v [deg]", fontsize=17, color="#333533")
        ax[n, j].set_ylabel("T [K]", fontsize=17, color="#333533")

        ax[n, j].set_title('Node_'+str(i+1), fontweight="bold", size=20)

        j += 1

        plt.grid(True, color="#93a1a1", alpha=0.3)
        
def F(X, t, params):
    # Ecuación de los dos cuerpos
    mu = params[0]  # km^3/s^2 Parámetro gravitacional
    r_vec = X[:3]   # [km] Posicion
    v_vec = X[3:]   # [km/s] Velocidad
    r = np.linalg.norm(r_vec)
    
    total_acc = - mu * r_vec / r**3
    
    dX = np.concatenate([v_vec, total_acc])  # concatenate() vuelve a crear un vector de 6 componentes a partir de dos de 3

    return dX

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