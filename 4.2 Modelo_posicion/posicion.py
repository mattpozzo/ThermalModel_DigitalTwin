############ Attitude and Position Model Block ################
  # Coded by Universitwin - Digital Twin Team
  # April 2020
  # Started by: Nicolás Valentín Conde
  # Version: V02
  # Latest version by:
  # Latest version date:

# Libraries...
# Math and data treatment
import numpy as np
import pandas as pd
# Graphic visualization
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
# Time and orbital propagation model
from datetime import datetime, timedelta
from sgp4.earth_gravity import wgs84  # Modelo gravitatorio
from sgp4.io import twoline2rv

import posicion_funciones as f

# Input

# Input 1: Earth parameters
mu_earth = 398600.4415 # [km^3/sec^2] Parametro gravitacional
R_earth = 6378.1363  # [km] Radio de la Tierra
params = (mu_earth,)

# Input 2: Starlink 53°
line1= ('1 45673U 20035S   20221.91667824  .00035327  00000-0  24310-2 0  9997')
line2 = ('2 45673  53.0047 218.9931 0006967 329.6241 290.2676 15.05573709 11055')
satellite = twoline2rv(line1, line2, wgs84)
## Objeto `datetime`

year = 2020
epoch_day = satellite.epochdays
t0 = datetime(year, 1, 1) + timedelta(epoch_day - 1)

# Input 3: Vector de tiempos
t = t0
tf = t0 + timedelta(minutes=94.6*230)
dt = timedelta(seconds=100)

n = int((tf - t0)/dt + 1)  ## Número de puntos
position = np.zeros([3, n])
velocity = np.zeros([3, n])

# Algorithm...

# Propagación de la órbita
i = 0
t = t0
while t < tf:
    #print(t)
    position[:, i], velocity[:, i] = satellite.propagate(t.year, t.month, t.day, t.hour, t.minute, t.second)
    i += 1
    t += dt

# Groundtrack
g = np.array([])
for i in range(position.shape[1]-1):
    g = np.concatenate([g, f.groundtrack(position[:,i], satellite, i)])

# Output...
output = pd.DataFrame (np.concatenate([position, velocity]).T, columns = ['X','Y','Z','U','V','W'])
output.index = output.index*dt

#output.to_csv(r'Input_CSVs\orbita.csv')

groundtrack_earth = pd.DataFrame()
groundtrack_earth = groundtrack_earth.assign(lat=g[0::3])
groundtrack_earth = groundtrack_earth.assign(lon=g[1::3])
groundtrack_earth = groundtrack_earth.assign(alt=g[2::3])

#groundtrack_earth.to_csv(r'Input_CSVs\groundtrack_earthrotation.csv')

# Visualization
  # Input 1 visualization
#%matplotlib notebook

u = g[1::3]
v = g[0::3]
        
u2 = u-360/(24*60)*95/2
for i in range(0,len(u2)):
    if abs(u2[i]) >=180:
        u2[i] +=360
        
fig = plt.figure(figsize=(8,6), edgecolor='w')

th =  np.arange(0, 360)*2*np.pi/360
r = 0.0015*180/np.pi

jey = 57*229
e=0
p = 0


ax = plt.axes(projection=ccrs.PlateCarree())
ax.stock_img()

for i in range(0,jey):
    if u[i+1]> 150 and u[i]<-160:
        ax.plot(u[e+1:i], v[e+1:i],'k', linewidth=r)
        e = i
    if u2[i+1]> 150 and u2[i]<-160:
        ax.plot(u2[p+1:i], v[p+1:i],'r', linewidth=r)
        p = i
plt.title('Starlink Orbit, i=53°')
plt.show()

########## ##     ####     ######
    ##     ##    ##  ##    ##   ##
    ##     ##   ##    ##   ##   ##
    ##     ##   ########   ######
    ##     ##  ##      ##  ##   ## 
    ##     ##  ##      ##  ##   ##
    ##     ##  ##      ##  ######
