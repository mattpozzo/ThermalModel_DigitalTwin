# -*- coding: utf-8 -*-
"""
Created on Thu Sep 24 16:17:20 2020

@author: nicol
"""
import numpy as np
import matplotlib.pyplot as plt

def j0(year,month,day):
    # =============================================================================
    #     % This function computes the Julian day number at 0 UT for any
    #     % year between 1900 and 2100 using Equation 5.48.
    #     %
    #     % j0 - Julian day at 0 hr UT (Universal Time)
    #     % year - range: 1901 - 2099
    #     % month - range: 1 - 12
    #     % day - range: 1 - 31
    # =============================================================================
    return 367*year - np.fix(7*(year + np.fix((month + 9)/12))/4) + np.fix(275*month/9) + day + 1721013.5;

def JD(j0,UT):
    return j0+UT/24

def sideral_century(j0):
    J2000 = 2451545
    return (j0-J2000)/36525

def GMT0(T0):
    return 100.4606184 + 36000.77004*T0 + 0.000387933*T0**2 - 2.583e-8*T0**3

def gmst02gmst(gmst0, ut):
    return gmst0+360.98564724*ut/24

def time2right_asention_declination(year, month, day, hour, minute, second, east_longitude):
    # Code from: http://www.agopax.it/Libri_astronomia/pdf/Astronomical%20Algorithms.pdf
    
    epsilon = -23.45*np.pi/180 #Angle between ecliptic and equator
    J0=j0(year,month,day) # Julian day at 0 UT
    ut = hour+minute/60+second/60/60 # UT orbit
    T=sideral_century(J0) # Sideral Centuries of j0
    
    # Mean longitude of the sun
    L0 = (280.46646+36000.76983*T+0.0003032*T**2-2.583e-8*T**3+360.98564724*ut/24)*np.pi/180
    L0-=np.fix(L0/np.pi/2)*np.pi*2 # 0 < g < 360
    # Mean anomaly of the sun
    M = (357.52911+35999.05029*T+0.0001537*T**2)*np.pi/180
    M-=np.fix(M/np.pi/2)*np.pi*2 # 0 < g < 360
    
    # Center
    C = (1.914602-0.004817*T-0.000014*T**2)*np.sin(M)+(0.019993-0.000101*T)*np.sin(2*M)+0.000289*np.sin(3*M)
    C *= np.pi/180
    C-=np.fix(C/np.pi/2)*np.pi*2 # 0 < g < 360
    
    # Sun Raan
    lon = L0 + C
    
    # Declination and right ascension
    dec = -np.arcsin(np.sin(epsilon)*np.sin(lon)) # declination of Sun
    ra = np.arctan2(np.cos(epsilon)*np.sin(lon), np.cos(lon)) 
    return ra, dec

def raan_dot(i, a, ex):
    J2 = 0.0010826
    mu = 398600
    Re = 6378
    return -(1.5*mu**.5*J2*Re**2/(1-ex**2)**2/a**3.5)*np.cos(i)*180*24*60*60/np.pi

def beta_sun(ra, dec, i, raan):
    termino1 = np.cos(ra)*np.sin(raan)*np.sin(i)
    termino2 = np.sin(ra)*np.cos(dec)*np.cos(raan)*np.sin(i)
    termino3 = np.sin(ra)*np.sin(dec)*np.cos(i)
    angle = termino1-termino2+termino3 
    return np.arcsin(angle)

def fraction_eclipse(beta, e, i, h, Re):
    beta_critico = np.arcsin(Re/(h+Re))
    if abs(beta) < beta_critico:
        angle = (h*Re*2+h**2)**0.5/(Re+h)/np.cos(beta)
        return np.arccos(angle)/np.pi
    else:
        return 0

# Órbita 
Re = 6378*1000 #m
h = 702*1000 #m
T = 99.8 # min
y,m,d,ho,mi,se = 2020,9,25,0,0,0 # year, month, day, UT
lon = 0 # local longitude
i = 98.2 * np.pi / 180 # rad
ex = 0.0
mu = 3.986e14
raan = (157+90+360.9856*(22+20/60)/24-360*2)*np.pi/180
#raan = 358.77*np.pi/180
ra_sun = []
dec_sun = []
beta = []
fe = []
r_dot=[]
t_eclipse = []
raan_=[]
epsilon = -23.45

q = 2
for j in range(0,366*q):
    ra, dec= time2right_asention_declination(y, m, j+d, ho, mi, se, lon)
    #dec, epsilon = right_asc_and_dec_Sun(y,m,d,ho+mi/60,j)
    ra_sun += [ra*180/np.pi]
    dec_sun += [dec*180/np.pi]
    raan += raan_dot(i, (Re+h)/1000, ex)*np.pi/180
    if raan < -np.pi:
        raan += 2*np.pi
    elif raan > np.pi:
        raan -= 2*np.pi
    raan_ += [raan*180/np.pi]
    b = beta_sun(ra, dec, i, raan)
    beta += [b*180/np.pi]
    fe += [fraction_eclipse(b, epsilon, i, h, Re)]
    t_eclipse += [T*fe[j]]
 
s = 1
plt.figure(1)
plt.scatter(np.arange(0,366*q),raan_,s=s)
plt.title("Raan angle - inclination = "+str(round(i*180/np.pi)))
plt.ylabel("Raan angle [deg]")
plt.xlabel("Days from start")

plt.figure(2)
plt.scatter(np.arange(0,366*q),dec_sun,s=s)
plt.title("Declination - inclination = "+str(round(i*180/np.pi)))
plt.ylabel("dec [deg]")
plt.xlabel("Days from start")

plt.figure(3)
plt.scatter(np.arange(0,366*q),ra_sun,s=s)
plt.title("right aasc - inclination = "+str(round(i*180/np.pi)))
plt.ylabel("right asc [deg]")
plt.xlabel("Days from start")

plt.figure(4)
plt.scatter(np.arange(0,366*q),beta,s=s)
plt.title("Beta angle - inclination = "+str(round(i*180/np.pi)))
plt.ylim([-90,90])
plt.ylabel("beta angle [deg]")
plt.xlabel("Days from start")

plt.figure(5)
plt.plot(np.arange(0,366*q),t_eclipse)
plt.ylim([0, 45])
plt.title("Eclipse Time - inclination = "+str(round(i*180/np.pi)))
plt.ylabel("time eclipse [min]")
plt.xlabel("Days from start")
