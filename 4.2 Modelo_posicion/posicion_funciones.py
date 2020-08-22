############ Attitude and Position Model Block ################
  # Coded by Universitwin - Digital Twin Team
  # April 2020
  # Started by: Nicolás Valentín Conde
  # Version: V02
  # Latest version by:
  # Latest version date:
from math import atan2, cos, pi, sin, sqrt
# Functions...
def groundtrack(vector, satellite, i):
    x = vector[0]
    y = vector[1]
    z = vector[2]

    # Constants for WGS-87 ellipsoid
    a = 6378.1363
    e = satellite.ecco

    # Groundtrack
    b = sqrt(pow(a, 2) * (1-pow(e, 2)))
    ep = sqrt((pow(a, 2)-pow(b, 2))/pow(b, 2))
    p = sqrt(pow(x, 2)+pow(y, 2))
    th = atan2(a*z, b*p)
    lon = atan2(y, x)
    lat = atan2((z+ep*ep*b*pow(sin(th), 3)), (p-e*e*a*pow(cos(th), 3)))
    n = a/sqrt(1-e*e*pow(sin(lat), 2))

    alt = p/cos(lat)-n
    lat = (lat*180)/pi
    gmra = (lon*180)/pi
    
    # Calculate GST https://casa.nrao.edu/casadocs/casa-5.1.0/reference-material/time-reference-frames
    d = satellite.jdsatepoch - 2451545.0
    T = d / 36525
    T0 = 24110.54841 + 8640184.812866 * T + 0.093104 * T**2 - 0.0000062 * T**3
    T0 = T0 % 24
    
    UT =  (satellite.epochdays * 24 + i*100/3600) * 1.002737909
    T0 += UT
    gmst = T0 % 24

    lon = (gmst * 15.0) - gmra
    lon = ((lon + 180.0) % 360.0) - 180.0

    return (lat, lon, alt)

def datetime2decHours(time):
    """
    Converts a datetime.time or datetime.datetime object into decimal time.
    Parameters
    ----------
    time : datetime.time or datetime.datetime
    Returns
    -------
    decTime : float
        A decimal number representing the input time
    """
    return time.hour + time.minute/60.0 + time.second/3600.0 + time.microsecond/3600000000.0