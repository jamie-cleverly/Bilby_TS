#    meteorologicalfunctions.py  Provides meteorological functions to analysis routines
#    Copyright (C) 2015  Dr James Cleverly, Dr Peter Isaac
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

import constants as c
import numpy

def absolutehumidity(Ta,e):
    # Calculate absolute humidity from vapour pressure and temperature
    # Ah [g/m3] = e [kPa] / ((Ta [C] + C2K) * RV [J/(mg K)]
    # RV = 0.0004615  Gas constant for water vapor
    Ah = e / ((Ta + 273.15) * 0.0004615)
    return Ah

def aerodynamicresistance(Uavg,Ce):
    # Calculate the aerodynamic resistance, Stull 1988
    #  U - wind speed
    #  Ce - bulk transfer coefficient
    # Returns
    #  rav, s/m
    rav = 1 / (Uavg * Ce)
    return rav

def bulktransfercoefficient(wA,Uavg,qa,qs):
    # Calculate the bulk transfer coefficient to be used in aerodynamic resistance, Stull 1988
    wq = wA * c.g2kg / c.rho_water 
    Ce = -wq / (Uavg * (qa - qs))
    return Ce

def delta(Ta):
    # Calculate the slope of the saturation vapour pressure curve
    #  Ta - air temperature
    # Returns
    #  delta, kPa/K
    delta = (2503 / ((Ta + 237.3) ** 2 )) * numpy.ma.exp((17.27 * Ta) / (Ta + 237.3))
    return delta

def densitydryair(Ta,ps,vp):
    # Calculate density of dry air from temperature, pressure and vapour pressure
    #  Ta - air temperature, C
    #  ps - pressure, kPa
    #  vp - vapour pressure, kPa
    # Returns
    #  rhod - dry air density, kg/m3
    rhod = 1000*(ps-vp)/((Ta+273.15)*c.Rd)
    return rhod

def densitymoistair(Ta,ps,Ah):
    # Calculate density of  air from absolute humidity and temperature
    #  Ta - air temperature, C
    #  ps - pressure, kPa
    #  Ah - absolute humidity, g/m3
    # Returns
    #  rhom - moist air density, kg/m3
    vp = vapourpressure(Ah,Ta)
    rhom = (ps-vp)*1000/((Ta+273.15)*c.Rd) + vp*1000/((Ta+273.15)*c.Rv)
    return rhom

def efromrh(RH, T):
    # Vapour pressure (kPa) from relative humidity and temperature
    #  RH is the relative humidity, %
    #  T is the air temperature, C
    eRH = 0.01 * RH * es(T)
    return eRH

def es(T):
    # Saturation vapour pressure.
    #  T is the air temperature, C
    #  es is the saturation vapour pressure in kPa
    es = 0.6106 * numpy.exp(17.27 * T / (T + 237.3))
    return es

def gamma(ps,Cpm,Lv):
    # Calculate the approximate psychrometric coefficient
    #  ps - atmopsheric pressure, kPa
    # Returns
    #  gamma, kPa/K
    gamma = (Cpm * ps) / (0.622 * Lv)
    return gamma

def Lv(Ta):
    # Calculate Lv as a function of temperature, from Stull 1988
    #  Ta - air temperature, C
    # Returns
    #  Lv - latent heat of vapourisation, J/kg
    Lv = 2500800 - (2366.8 * Ta)
    return Lv

def mixingratio(ps,vp):
    # Calculate mixing ratio from vapour pressure and pressure
    #  ps - presure, kPa
    #  vp - vapour pressure, kPa
    # Returns
    #  mr - mixing ratio, kg/kg
    mr = 0.622*vp/(ps- vp)
    return mr

def molen(T,Ah,p,ustar,heatflux,fluxtype='sensible'):
    # Calculate the Monin-Obukhov length
    # Stull (1988) Eq. 5.7c, uses virtual potential temperature
    ustar = numpy.ma.sqrt(numpy.square(ustar))    # force the sign of ustar to be positive
    vp = vapourpressure(Ah,T)       # calculate the vapour pressure
    mr = mixingratio(p,vp)          # calculate the mixing ratio
    Tv = theta(T,p)                 # calculate potential temperature
    Tvp = virtualtheta(Tv,mr)
    if fluxtype=='sensible':
        L = -Tvp*densitydryair(T, p, vp)*c.Cp*(ustar**3)/(c.g*c.k*heatflux)
    elif fluxtype=='kinematic':
        L = -Tvp*(ustar**3)/(c.g*c.k*heatflux)
    else:
        raise Exception(" meteorologicalfunctions.molen: unkown value for fluxtype (="+str(fluxtype)+") encountered")
    return L

def qfromrh(RH, T, p):
    # Specific humidity (kg/kg) from relative humidity, temperature and pressure
    #  RH is the relative humidity, %
    #  T is the air temperature, C
    #  p is the atmospheric pressure, kPa
    qRH = (c.Mv / c.Md) * (0.01 * RH * es(T) / p)
    return qRH

def qsat(esat,ps):
    qsat = 0.622 * (esat / ps)
    return qsat

def specificheatmoistair(q):
    # Calculate Cp of moist air, from Stull 1988
    #  Cp - specific heat of dry air at constant pressure, J/kg-K
    #  q - specific humidity
    # Returns
    #  Cpm - specific heat of moist air at constant pressure, J/kg-K
    Cpm = c.Cpd * (1 + 0.84 * q)
    return Cpm

def specifichumidity(mr):
    # Calculate specific humidity from mixing ratio
    #  mr - mixing ration, kg/kg
    # Returns
    #  q = specific humidity, kg/kg
    q = mr/(1+mr)
    return q

def tafromtv(Tv,q):
    # Calculate air temperature from virtual temperature using formula
    # from Campbell Scientific CSAT manual.
    # NOTE: this differs from the usual definition by using 0.51 not 0.61
    #  Tv - virtual temperature, C
    #  q - specific humidity, kg/kg
    # Returns
    #  Ta - air temperature, C
    Ta = ((Tv+273.15)/(1+0.51*q))-273.15
    return Ta

def theta(T,p):
    # Calculate potential temperature from air temperature and pressure
    #  T - air temperature, C
    #  p - pressure, kPa
    # Returns
    #  theta - potential temperature, K
    return (T+273.15)*(100/p)**0.286

def vapourpressure(Ah,Ta):
    # Calculate vapour pressure from absolute humidity and temperature
    #  Ah - absolute humidity, g/m3
    #  Ta - air temperature, C
    # Returns
    #  vp - vapour pressure, kPa
    vp = 0.000001*Ah*(Ta+273.15)*c.R/c.Mv
    return vp

def virtualtheta(theta,mr):
    # Calculate virtual potential temperature
    #  theta - potential temperature, K
    #  mr - mixing ratio, kg/kg
    # Returns
    #  Tvp - virtual potential temperature, K
    Tvp = theta * (1 + (0.61 * mr))
    return Tvp
