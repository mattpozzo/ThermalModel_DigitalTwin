import numpy as np
from sympy import Plane, Point3D
from shapely.geometry import Polygon

#constantes
radiador= [[0, 0, 0], [1, 0, 0], [1, 0, 0.5], [0, 0, 2]] # lo que produce sombra
placa = [[0,0,0],[1,0,0],[1,1,0],[0,1,0]] # elemento a analizar
ang=45 #deg
vector=np.array([0,-np.cos(ang*np.pi/180),np.sin(ang*np.pi/180)]) #rayo de sol

#inicializando vectores
base_plane=Plane(Point3D(0,0,0), normal_vector=(0,0,1))
N=len(radiador)
normal=list(np.zeros(N))
constr_plane=list(np.zeros(N))
inter_plane=list(np.zeros(N))
sombra=list(np.zeros(N))
sombra_lista=list(np.zeros(N))


#Creo las lineas de la sombra del radiador
for i in range(N):
    #Calculo el vector normal de un plano que contiene a cada arista del radiador y al vector del rayo de sol
    if i!=N-1:
        normal[i]=np.cross(vector,np.array(radiador[i])-np.array(radiador[i+1]))
    else:
        normal[i]=np.cross(vector,np.array(radiador[i])-np.array(radiador[0]))
    #Creo el plano para el vector normal que acabo de crear
    constr_plane[i]=Plane(Point3D(radiador[i]), normal_vector=list(normal[i]))
    #Creo las aristas de la sombra del radiador usando la interseccion entre los planos creados y el plano base
    inter_plane[i]=Plane.intersection(base_plane, constr_plane[i])

#Calculo las coordenadas que definen a la sombra del radiador en el plano de la placa
for i in range(N):
    if i!=N-1:
        sombra[i]=(inter_plane[i])[0].intersection((inter_plane[i+1])[0])
    else:
        sombra[i]=(inter_plane[i])[0].intersection((inter_plane[0])[0])
    sombra_lista[i]=(sombra[i])[0] 

#Calculo que porcentaje de la placa esta iluminada
a = Polygon(placa)
b = Polygon(sombra_lista)
c = a.intersection(b)
not_ilum=c.area/a.area
ilum=1-not_ilum

#chin pum
text="Iluminado: "+str("%.2f" %(ilum*100))+"%\n"+"En sombra: "+str("%.2f" %(not_ilum*100))+"%"    
print(text)
    
    
#PLOT    
from matplotlib import pyplot
from descartes import PolygonPatch
BLUE = '#6699cc'
DARKGRAY = '#333333'
BLACK = '#000000'
fig = pyplot.figure(1)
ax = fig.add_subplot(111)
patch1 = PolygonPatch(a, fc=BLUE, ec=DARKGRAY, alpha=0.3, zorder=1)
ax.add_patch(patch1)
patch2 = PolygonPatch(b, fc=DARKGRAY, ec=DARKGRAY, alpha=0.3, zorder=1)
ax.add_patch(patch2)
patchc = PolygonPatch(c, fc=BLACK, ec=BLACK, alpha=0.8, zorder=1)
ax.add_patch(patchc)
ax.set_title('sombra')
props = dict(boxstyle='round', facecolor=BLUE, alpha=0.3)
ax.set_xlim(-0.5, 3)
ax.set_ylim(-0.5, 2.5)
ax.set_aspect("equal")
# place a text box in upper left in axes coords
ax.text(0.95, 0.75, text, transform=ax.transAxes, fontsize=14,
        verticalalignment='top', horizontalalignment='right', bbox=props)
pyplot.show()






    
    