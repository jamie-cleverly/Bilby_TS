#    qcts.py  Performs most of the analysis and synthesis functions
#    Copyright (C) 2015  Dr James Cleverly, Dr Peter Isaac
#    
#    footprint functions do_footprint_2d and footprint_2d
#        Applies the Kljun two-dimensional approximate analytical footprint
#        model to OzFlux standard L3 data
#    Copyright (C) 2014  James Cleverly, Eva Van Gorsel, Natascha Kljun
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

"""
    QC Data Function Module
    Used to perform the tasks queued by qcls.py
    """

import sys
import ast
from calendar import isleap
import constants as c
import datetime
from matplotlib.dates import date2num
import meteorologicalfunctions as mf
import numpy
import os
import qcck
import qcio
import qcts
import qcutils
from scipy import interpolate, signal
import time
import xlrd
from matplotlib.mlab import griddata
import xlwt
import logging
import math
import csv
import pysolar
import pdb

log = logging.getLogger('qc.ts')

def albedo(cf,ds):
    """
        Filter albedo measurements to:
            high solar angle specified by periods between 10.00 and 14.00, inclusive
            and
            full sunlight in which Fsd > 290 W/m2
        
        Usage qcts.albedo(ds)
        ds: data structure
        """
    log.info(' Applying albedo constraints')
    if 'Fsd' in ds.series.keys() and 'Fsu' in ds.series.keys() and 'albedo' in ds.series.keys():
        nightindex = numpy.ma.where(((ds.series['Fsd']['Data'] < 5) & (ds.series['Fsd']['Data'] > -9000)) | ((ds.series['Fsu']['Data'] < 2) & (ds.series['Fsu']['Data'] > -9000)))[0]
        ds.series['Fsd']['Data'][nightindex] = 0.
        ds.series['Fsu']['Data'][nightindex] = 0.
        ds.series['albedo']['Data'][nightindex] = 0.
    
    if 'albedo' not in ds.series.keys():
        if 'Fsd' in ds.series.keys() and 'Fsu' in ds.series.keys():
            Fsd,f,a = qcutils.GetSeriesasMA(ds,'Fsd')
            Fsu,f,a = qcutils.GetSeriesasMA(ds,'Fsu')
            albedo = Fsu / Fsd
            nightindex = numpy.ma.where(((ds.series['Fsd']['Data'] < 5) & (ds.series['Fsd']['Data'] > -9000)) | ((ds.series['Fsu']['Data'] < 2) & (ds.series['Fsu']['Data'] > -9000)))[0]
            ds.series['Fsd']['Data'][nightindex] = 0.
            ds.series['Fsu']['Data'][nightindex] = 0.
            ds.series['albedo']['Data'][nightindex] = 0.
            attr = qcutils.MakeAttributeDictionary(long_name='solar albedo',units='none',standard_name='solar_albedo')
            qcutils.CreateSeries(ds,'albedo',albedo,FList=['Fsd','Fsu'],Attr=attr)
        else:
            log.warning('  Fsd or Fsu not in ds, albedo not calculated')
            return
    else:
        albedo,f,a = qcutils.GetSeriesasMA(ds,'albedo')
        if 'Fsd' in ds.series.keys():
            Fsd,f,a = qcutils.GetSeriesasMA(ds,'Fsd')
        else:
            Fsd,f,a = qcutils.GetSeriesasMA(ds,'Fn')
    
    if qcutils.cfkeycheck(cf,ThisOne='albedo',key='Threshold'):
        Fsdbase = numpy.float64(cf['Variables']['albedo']['Threshold']['Fsd'])
        ds.series['albedo']['Attr']['FsdCutoff'] = Fsdbase
    else:
        Fsdbase = 290.
    index = numpy.ma.where((Fsd < Fsdbase) | (ds.series['Hdh']['Data'] < 10) | (ds.series['Hdh']['Data'] > 14))[0]
    index1 = numpy.ma.where(Fsd < Fsdbase)[0]
    index2 = numpy.ma.where((ds.series['Hdh']['Data'] < 10) | (ds.series['Hdh']['Data'] > 14))[0]
    albedo[index] = numpy.float64(c.missing_value)
    ds.series['albedo']['Flag'][index1] = numpy.int32(8)     # bad Fsd flag only if bad time flag not set
    ds.series['albedo']['Flag'][index2] = numpy.int32(9)     # bad time flag
    ds.series['albedo']['Data']=numpy.ma.filled(albedo,numpy.float64(c.missing_value))

def ApplyLinear(cf,ds,ThisOne):
    """
        Applies a linear correction to variable passed from qcls. Time period
        to apply the correction, slope and offset are specified in the control
        file.
        
        Usage qcts.ApplyLinear(cf,ds,x)
        cf: control file
        ds: data structure
        x: input/output variable in ds.  Example: 'Cc_7500_Av'
        """
    if ThisOne not in ds.series.keys(): return
    log.info('  Applying linear correction to '+ThisOne)
    if qcutils.incf(cf,ThisOne) and qcutils.haskey(cf,ThisOne,'Linear'):
        data = numpy.ma.masked_where(ds.series[ThisOne]['Data']==numpy.float64(c.missing_value),ds.series[ThisOne]['Data'])
        flag = ds.series[ThisOne]['Flag'].copy()
        ldt = ds.series['DateTime']['Data']
        LinearList = cf['Variables'][ThisOne]['Linear'].keys()
        for i in range(len(LinearList)):
            LinearItemList = ast.literal_eval(cf['Variables'][ThisOne]['Linear'][str(i)])
            try:
                si = ldt.index(datetime.datetime.strptime(LinearItemList[0],'%Y-%m-%d %H:%M'))
            except ValueError:
                si = 0
            try:
                ei = ldt.index(datetime.datetime.strptime(LinearItemList[1],'%Y-%m-%d %H:%M')) + 1
            except ValueError:
                ei = -1
            Slope = numpy.float64(LinearItemList[2])
            Offset = numpy.float64(LinearItemList[3])
            data[si:ei] = Slope * data[si:ei] + Offset
            index = numpy.where(flag[si:ei]==0)[0]
            flag[si:ei][index] = numpy.int32(10)
            ds.series[ThisOne]['Data'] = numpy.ma.filled(data,numpy.float64(c.missing_value))
            ds.series[ThisOne]['Flag'] = flag

def ApplyLinearDrift(cf,ds,ThisOne):
    """
        Applies a linear correction to variable passed from qcls. The slope is
        interpolated for each 30-min period between the starting value at time 0
        and the ending value at time 1.  Slope0, Slope1 and Offset are defined
        in the control file.  This function applies to a dataset in which the
        start and end times in the control file are matched by the time period
        in the dataset.
        
        Usage qcts.ApplyLinearDrift(cf,ds,x)
        cf: control file
        ds: data structure
        x: input/output variable in ds.  Example: 'Cc_7500_Av'
        """
    log.info('  Applying linear drift correction to '+ThisOne)
    if qcutils.incf(cf,ThisOne) and qcutils.haskey(cf,ThisOne,'Drift'):
        data = numpy.ma.masked_where(ds.series[ThisOne]['Data']==numpy.float64(c.missing_value),ds.series[ThisOne]['Data'])
        flag = ds.series[ThisOne]['Flag']
        ldt = ds.series['DateTime']['Data']
        DriftList = cf['Variables'][ThisOne]['Drift'].keys()
        for i in range(len(DriftList)):
            DriftItemList = ast.literal_eval(cf['Variables'][ThisOne]['Drift'][str(i)])
            try:
                si = ldt.index(datetime.datetime.strptime(DriftItemList[0],'%Y-%m-%d %H:%M'))
            except ValueError:
                si = 0
            try:
                ei = ldt.index(datetime.datetime.strptime(DriftItemList[1],'%Y-%m-%d %H:%M')) + 1
            except ValueError:
                ei = -1
            Slope = numpy.zeros(len(data),dtype=numpy.float64)
            Slope0 = numpy.float64(DriftItemList[2])
            Slope1 = numpy.float64(DriftItemList[3])
            Offset = numpy.float64(DriftItemList[4])
            nRecs = len(Slope[si:ei])
            for i in range(nRecs):
                ssi = si + i
                Slope[ssi] = ((((Slope1 - Slope0) / nRecs) * i) + Slope0)
            data[si:ei] = Slope[si:ei] * data[si:ei] + Offset
            flag[si:ei] = numpy.int32(10)
            ds.series[ThisOne]['Data'] = numpy.ma.filled(data,numpy.float64(c.missing_value))
            ds.series[ThisOne]['Flag'] = flag

def ApplyLinearDriftLocal(cf,ds,ThisOne):
    """
        Applies a linear correction to variable passed from qcls. The slope is
        interpolated since the starting value at time 0 using a known 30-min
        increment.  Slope0, SlopeIncrement and Offset are defined in the control
        file.  This function applies to a dataset in which the start time in the
        control file is matched by dataset start time, but in which the end time
        in the control file extends beyond the dataset end.
        
        Usage qcts.ApplyLinearDriftLocal(cf,ds,x)
        cf: control file
        ds: data structure
        x: input/output variable in ds.  Example: 'Cc_7500_Av'
        """
    log.info('  Applying linear drift correction to '+ThisOne)
    if qcutils.incf(cf,ThisOne) and qcutils.haskey(cf,ThisOne,'LocalDrift'):
        data = numpy.ma.masked_where(ds.series[ThisOne]['Data']==numpy.float64(c.missing_value),ds.series[ThisOne]['Data'])
        flag = ds.series[ThisOne]['Flag']
        ldt = ds.series['DateTime']['Data']
        DriftList = cf['Variables'][ThisOne]['LocalDrift'].keys()
        for i in range(len(DriftList)):
            DriftItemList = ast.literal_eval(cf['Variables'][ThisOne]['LocalDrift'][str(i)])
            try:
                si = ldt.index(datetime.datetime.strptime(DriftItemList[0],'%Y-%m-%d %H:%M'))
            except ValueError:
                si = 0
            try:
                ei = ldt.index(datetime.datetime.strptime(DriftItemList[1],'%Y-%m-%d %H:%M')) + 1
            except ValueError:
                ei = -1
            Slope = numpy.zeros(len(data),dtype=numpy.float64)
            Slope0 = numpy.float64(DriftItemList[2])
            SlopeIncrement = numpy.float64(DriftItemList[3])
            Offset = numpy.float64(DriftItemList[4])
            nRecs = len(Slope[si:ei])
            for i in range(nRecs):
                ssi = si + i
                Slope[ssi] = (SlopeIncrement * i) + Slope0
            data[si:ei] = Slope[si:ei] * data[si:ei] + Offset
            flag[si:ei] = numpy.int32(10)
            ds.series[ThisOne]['Data'] = numpy.ma.filled(data,numpy.float64(c.missing_value))
            ds.series[ThisOne]['Flag'] = flag

def AverageSeriesByElements(cf,ds,Av_out):
    """
        Calculates the average of multiple time series.  Multiple time series
        are entered and a single time series representing the average at each
        observational period is returned.
        
        Usage qcts.AverageSeriesByElements(cf,ds,Av_out)
        cf: control file object (must contain an entry for Av_out)
        ds: data structure
        Av_out: output variable to ds.  Example: 'Fg'
        Series_in: input variable series in ds.  Example: ['Fg_8cma','Fg_8cmb']
        """
    if Av_out not in cf['Variables'].keys(): return
    if Av_out in ds.averageserieslist: return
    srclist, standardname = qcutils.GetAverageSeriesKeys(cf,Av_out)
    log.info(' Averaging series in '+str(srclist)+' into '+Av_out)
    nSeries = len(srclist)
    if nSeries==0:
        log.error('  AverageSeriesByElements: no input series specified for'+str(Av_out))
        return
    if nSeries==1:
        tmp_data = ds.series[srclist[0]]['Data'].copy()
        tmp_flag = ds.series[srclist[0]]['Flag'].copy()
        tmp_attr = ds.series[srclist[0]]['Attr'].copy()
        Av_data = numpy.ma.masked_where(tmp_data==numpy.float64(c.missing_value),tmp_data)
        Mn_flag = tmp_flag
        SeriesNameString = srclist[0]
    else:
        tmp_data = ds.series[srclist[0]]['Data'].copy()
        tmp_flag = ds.series[srclist[0]]['Flag'].copy()

        index = numpy.where(numpy.mod(tmp_flag,10)==0)    # find the elements with flag = 0, 10, 20 etc
        tmp_flag[index] = numpy.int32(0)                               # set them all to 0
        
        tmp_attr = ds.series[srclist[0]]['Attr'].copy()
        SeriesNameString = srclist[0]
        srclist.remove(srclist[0])
        for ThisOne in srclist:
            SeriesNameString = SeriesNameString+', '+ThisOne
            tmp_data = numpy.vstack((tmp_data,ds.series[ThisOne]['Data'].copy()))
            tmp_flag = numpy.vstack((tmp_flag,ds.series[ThisOne]['Flag'].copy()))
        tmp_data = numpy.ma.masked_where(tmp_data==numpy.float64(c.missing_value),tmp_data)
        Av_data = numpy.ma.average(tmp_data,axis=0)
        Mn_flag = numpy.min(tmp_flag,axis=0)
    ds.averageserieslist.append(Av_out)
    DStr = 'Element-wise average of series '+SeriesNameString
    UStr = ds.series[srclist[0]]['Attr']['units']
    SStr = ds.series[srclist[0]]['Attr']['standard_name']
    attr = qcutils.MakeAttributeDictionary(long_name=DStr,units=UStr,standard_name=SStr)
    qcutils.CreateSeries(ds,Av_out,Av_data,Flag=Mn_flag,Attr=attr)

def AverageSeriesByElementsI(cf,ds,Av_out):
    """
        Calculates the average of multiple time series.  Multiple time series
        are entered and a single time series representing the average at each
        observational period is returned.
        
        Usage qcts.AverageSeriesByElements(cf,ds,Av_out)
        cf: control file object (must contain an entry for Av_out)
        ds: data structure
        Av_out: output variable to ds.  Example: 'Fg'
        Series_in: input variable series in ds.  Example: ['Fg_8cma','Fg_8cmb']
        """
    if Av_out not in cf['Variables'].keys(): return
    if Av_out in ds.averageserieslist: return
    srclist, standardname = qcutils.GetAverageSeriesKeys(cf,Av_out)
    log.info(' Averaging series in '+str(srclist)+' into '+Av_out)
    nSeries = len(srclist)
    if nSeries==0:
        log.error('  AverageSeriesByElements: no input series specified for'+str(Av_out))
        return
    if nSeries==1:
        tmp_data = ds.series[srclist[0]]['Data'].copy()
        tmp_flag = ds.series[srclist[0]]['Flag'].copy()
        tmp_attr = ds.series[srclist[0]]['Attr'].copy()
        Av_data = numpy.ma.masked_where(tmp_data==numpy.float64(c.missing_value),tmp_data)
        Mx_flag = tmp_flag
        SeriesNameString = srclist[0]
    else:
        tmp_data = ds.series[srclist[0]]['Data'].copy()
        tmp_flag = ds.series[srclist[0]]['Flag'].copy()
        tmp_attr = ds.series[srclist[0]]['Attr'].copy()
        SeriesNameString = srclist[0]
        srclist.remove(srclist[0])
        for ThisOne in srclist:
            SeriesNameString = SeriesNameString+', '+ThisOne
            tmp_data = numpy.vstack((tmp_data,ds.series[ThisOne]['Data'].copy()))
            tmp_flag = numpy.vstack((tmp_flag,ds.series[ThisOne]['Flag'].copy()))
        tmp_data = numpy.ma.masked_where(tmp_data==numpy.float64(c.missing_value),tmp_data)
        Av_data = numpy.ma.average(tmp_data,axis=0)
        Mx_flag = numpy.max(tmp_flag,axis=0)
    ds.averageserieslist.append(Av_out)
    DStr = 'Element-wise average of series '+SeriesNameString
    UStr = ds.series[srclist[0]]['Attr']['units']
    SStr = ds.series[srclist[0]]['Attr']['standard_name']
    attr = qcutils.MakeAttributeDictionary(long_name=DStr,units=UStr,standard_name=SStr)
    qcutils.CreateSeries(ds,Av_out,Av_data,Flag=Mx_flag,Attr=attr)

def BypassTcorr(cf,ds):
    if qcutils.cfkeycheck(cf,Base='Soil',ThisOne='BypassTcorrList'):
        subkeys = ast.literal_eval(cf['Soil']['BypassTcorrList'])
        for i in range(len(subkeys)):
            if subkeys[i] in ds.series.keys():
                Sws,f,a = qcutils.GetSeriesasMA(ds,subkeys[i])
                Sws_bypass = -0.0663 + (-0.0063 * Sws) + (0.0007 * Sws ** 2)
                attr = qcutils.MakeAttributeDictionary(long_name=ds.series[subkeys[i]]['Attr']['long_name'],units='frac',standard_name='soil_moisture_content')
                qcutils.CreateSeries(ds,subkeys[i],Sws_bypass,FList=[subkeys[i]],Attr=attr)
    else:
        for ThisOne in ds.series.keys():
            if 'Sws' in ThisOne:
                Sws,f,a = qcutils.GetSeriesasMA(ds,ThisOne)
                Sws_bypass = -0.0663 + (-0.0063 * Sws) + (0.0007 * Sws ** 2)
                attr = qcutils.MakeAttributeDictionary(long_name=ds.series[ThisOne]['Attr']['long_name'],units='frac',standard_name='soil_moisture_content')
                qcutils.CreateSeries(ds,ThisOne,Sws_bypass,FList=[ThisOne],Attr=attr)

def CalculateAhHMP(cf,ds,e_name='e',Ta_name='Ta',Ah_name='Ah',RH_name='RH'):
    """
        Calculate the absolute humidity from vapour pressure and temperature.
        
        """
    log.info('  Calculating Ah from vapour pressure and air temperature')
    if qcutils.cfkeycheck(cf,Base='FunctionArgs',ThisOne='AhVars'):
        vars = ast.literal_eval(cf['FunctionArgs']['AhVars'])
        e_name = vars[0]
        Ta_name = vars[1]
        Ah_name = vars[2]
    Ta,f,a = qcutils.GetSeriesasMA(ds,Ta_name)
    RH,f,a = qcutils.GetSeriesasMA(ds,RH_name)
    if e_name in ds.series.keys():
        e,f,a = qcutils.GetSeriesasMA(ds,e_name)
    else:
        e = mf.efromrh(RH,Ta)
        attr = qcutils.MakeAttributeDictionary(long_name='Vapour pressure',units='kPa',standard_name='water_vapor_partial_pressure_in_air')
        qcutils.CreateSeries(ds,e_name,e,FList=[Ta_name,RH_name],Attr=attr)
    
    Ah = mf.absolutehumidity(Ta,e)
    attr = qcutils.MakeAttributeDictionary(long_name='Absolute humidity (HMP)',units='g/m3',standard_name='mass_concentration_of_water_vapor_in_air')
    qcutils.CreateSeries(ds,Ah_name,Ah,FList=[Ta_name,e_name],Attr=attr)

def CalculateAvailableEnergy(ds,Fa_out='Fa',Fn_in='Fn',Fg_in='Fg'):
    """
        Calculate the average energy as Fn - G.
        
        Usage qcts.CalculateAvailableEnergy(ds,Fa_out='Fa',Fn_in='Fn',Fg_in='Fg')
        ds: data structure
        Fa_out: output available energy variable to ds.  Example: 'Fa'
        Fn_in: input net radiation in ds.  Example: 'Fn'
        Fg_in: input ground heat flux in ds.  Example: 'Fg'
        """
    log.info(' Calculating available energy from Fn and Fg')
    Fn,f,a = qcutils.GetSeriesasMA(ds,Fn_in)
    Fg,f,a = qcutils.GetSeriesasMA(ds,Fg_in)
    Fa = Fn - Fg
#    pdb.set_trace()
    attr = qcutils.MakeAttributeDictionary(long_name='Available energy using '+Fn_in+','+Fg_in,units='W/m2')
    qcutils.CreateSeries(ds,Fa_out,Fa,FList=[Fn_in,Fg_in],Attr=attr)

def CalculateET(cf,ds,level,Fe_in='Fe'):
    if qcutils.cfkeycheck(cf,Base='CalculateET',ThisOne='Fe_in'):
        Fe_in = cf['Sums']['Fe_in']
    Fe,f,a = qcutils.GetSeriesasMA(ds,Fe_in)
    Lv,f,a = qcutils.GetSeriesasMA(ds,'Lv')
    ET = Fe * 60 * 30 * 1000 / (Lv * c.rho_water)  # mm/30min for summing
    if 'ET' not in ds.series.keys():
        attr = qcutils.MakeAttributeDictionary(long_name='Evapotranspiration Flux',units='mm/30min',standard_name='lwe_water_evaporation_rate')
        qcutils.CreateSeries(ds,'ET',ET,FList=[Fe_in],Attr=attr)
        ds.series['ET']['Attr']['Level'] = level
    elif ds.series['ET']['Attr']['Level'] == level:
        log.warn('   ET already in dataset at level '+level+': ET not re-computed')
    else:
        attr = qcutils.MakeAttributeDictionary(long_name='',units='mm/30min',standard_name='lwe_water_evaporation_rate')
        qcutils.CreateSeries(ds,'ET',ET,FList=[Fe_in],Attr=attr)
        ds.series['ET']['Attr']['Level'] = level

def CalculateFluxes(cf,ds,Ta_name='Ta',ps_name='ps',Ah_name='Ah',wT_in='wT',wA_in='wA',wC_in='wC',uw_in='uw',vw_in='vw',Fh_out='Fh',Fe_out='Fe',Fc_out='Fc',Fm_out='Fm',ustar_out='ustar'):
    """
        Calculate the fluxes from the rotated covariances.
        
        Usage qcts.CalculateFluxes(ds)
        ds: data structure
        
        Pre-requisite: CoordRotation2D
        
        Accepts meteorological constants or variables
        """
    if qcutils.cfkeycheck(cf,Base='FunctionArgs',ThisOne='CF'):
        args = ast.literal_eval(cf['FunctionArgs']['CF'])
        Ta_name = args[0]
        Ah_name = args[1]
        ps_name = args[2]
        wT_in = args[3]
        wA_in = args[4]
        wC_in = args[5]
        uw_in = args[6]
        vw_in = args[7]
        Fh_out = args[8]
        Fe_out = args[9]
        Fc_out = args[10]
        Fm_out = args[11]
        ustar_out = args[12]
    long_name = ''
    if 'Massman' in ds.globalattributes['L3Functions']:
        long_name = ', frequency response corrected'
    Ta,f,a = qcutils.GetSeriesasMA(ds,Ta_name)
    ps,f,a = qcutils.GetSeriesasMA(ds,ps_name)
    Ah,f,a = qcutils.GetSeriesasMA(ds,Ah_name)
    rhom,f,a = qcutils.GetSeriesasMA(ds,'rhom')
    log.info(' Calculating fluxes from covariances')
    if wT_in in ds.series.keys():
        wT,f,a = qcutils.GetSeriesasMA(ds,wT_in)
        Fh = rhom * c.Cpd * wT
        attr = qcutils.MakeAttributeDictionary(long_name='rotated to natural wind coordinates'+long_name,units='W/m2',standard_name='surface_upward_sensible_heat_flux')
        qcutils.CreateSeries(ds,Fh_out,Fh,FList=[wT_in],Attr=attr)
        flagindex = numpy.where(numpy.mod(ds.series[Fh_out]['Flag'],10)!=0)[0]
        ds.series[Fh_out]['Data'][flagindex] = numpy.float64(c.missing_value)
    else:
        log.error('  CalculateFluxes: '+wT_in+' not found in ds.series, Fh not calculated')
    if wA_in in ds.series.keys():
        wA,f,a = qcutils.GetSeriesasMA(ds,wA_in)
        if 'Lv' in ds.series.keys():
            Lv,f,a = qcutils.GetSeriesasMA(ds,'Lv')
            Fe = Lv * wA / numpy.float64(1000)
        else:
            Fe = c.Lv * wA / numpy.float64(1000)
        attr = qcutils.MakeAttributeDictionary(long_name='rotated to natural wind coordinates'+long_name,units='W/m2',standard_name='surface_upward_latent_heat_flux')
        qcutils.CreateSeries(ds,Fe_out,Fe,FList=[wA_in],Attr=attr)
        flagindex = numpy.where(numpy.mod(ds.series[Fe_out]['Flag'],10)!=0)[0]
        ds.series[Fe_out]['Data'][flagindex] = numpy.float64(c.missing_value)
    else:
        log.error('  CalculateFluxes: '+wA_in+' not found in ds.series, Fe not calculated')
    if wC_in in ds.series.keys():
        wC,f,a = qcutils.GetSeriesasMA(ds,wC_in)
        Fc = wC
        attr = qcutils.MakeAttributeDictionary(long_name='rotated to natural wind coordinates'+long_name,units='mg/m2/s')
        qcutils.CreateSeries(ds,Fc_out,Fc,FList=[wC_in],Attr=attr)
        flagindex = numpy.where(numpy.mod(ds.series[Fc_out]['Flag'],10)!=0)[0]
        ds.series[Fc_out]['Data'][flagindex] = numpy.float64(c.missing_value)
    else:
        log.error('  CalculateFluxes: '+wC_in+' not found in ds.series, Fc_raw not calculated')
    if uw_in in ds.series.keys():
        if vw_in in ds.series.keys():
            uw,f,a = qcutils.GetSeriesasMA(ds,uw_in)
            vw,f,a = qcutils.GetSeriesasMA(ds,vw_in)
            vs = uw*uw + vw*vw
            Fm = rhom * numpy.ma.sqrt(vs)
            us = numpy.ma.sqrt(numpy.ma.sqrt(vs))
            attr = qcutils.MakeAttributeDictionary(long_name='rotated to natural wind coordinates'+long_name,units='kg/m/s2')
            qcutils.CreateSeries(ds,Fm_out,Fm,FList=['rhom',uw_in,vw_in],Attr=attr)
            attr = qcutils.MakeAttributeDictionary(long_name='rotated to natural wind coordinates'+long_name,units='m/s')
            qcutils.CreateSeries(ds,ustar_out,us,FList=[uw_in,vw_in],Attr=attr)
            for ThisOne in [Fm_out,ustar_out]:
                flagindex = numpy.where(numpy.mod(ds.series[ThisOne]['Flag'],10)!=0)[0]
                ds.series[ThisOne]['Data'][flagindex] = numpy.float64(c.missing_value)
        else:
            log.error('  CalculateFluxes: wy not found in ds.series, Fm and ustar not calculated')
    else:
        log.error('  CalculateFluxes: wx not found in ds.series, Fm and ustar not calculated')

def CalculateLongwave(ds,Fl_out,Fl_in,Tbody_in):
    """
        Calculate the longwave radiation given the raw thermopile output and the
        sensor body temperature.
        
        Usage qcts.CalculateLongwave(ds,Fl_out,Fl_in,Tbody_in)
        ds: data structure
        Fl_out: output longwave variable to ds.  Example: 'Flu'
        Fl_in: input longwave in ds.  Example: 'Flu_raw'
        Tbody_in: input sensor body temperature in ds.  Example: 'Tbody'
        """
    log.info(' Calculating longwave radiation')
    Fl_raw,f,a = qcutils.GetSeriesasMA(ds,Fl_in)
    Tbody,f,a = qcutils.GetSeriesasMA(ds,Tbody_in)
    Fl = Fl_raw + c.sb*(Tbody + 273.15)**4
    attr = qcutils.MakeAttributeDictionary(long_name='Calculated longwave radiation using '+Fl_in+','+Tbody_in,units='W/m2')
    qcutils.CreateSeries(ds,Fl_out,Fl,FList=[Fl_in,Tbody_in],Attr=attr)

def CalculateMeteorologicalVariables(ds,Ta_name='Ta',Tv_name='Tv_CSAT',ps_name='ps',q_name="q",Ah_name='Ah',RH_name='RH',Cc_name='Cc'):
    """
        Add time series of meteorological variables based on fundamental
        relationships (Stull 1988)

        Usage qcts.CalculateMeteorologicalVariables(ds,Ta_name,ps_name,Ah_name)
        ds: data structure
        Ta_name: data series name for air temperature
        ps_name: data series name for pressure
        Ah_name: data series name for absolute humidity

        Variables added:
            rhom: density of moist air, mf.densitymoistair(Ta,ps,Ah)
            Lv: latent heat of vapourisation, mf.Lv(Ta)
            q: specific humidity, mf.specifichumidity(mr)
                where mr (mixing ratio) = mf.mixingratio(ps,vp)
            Cpm: specific heat of moist air, mf.specificheatmoistair(q)
            VPD: vapour pressure deficit, VPD = esat - e
        """
    Ta,f,a = qcutils.GetSeriesasMA(ds,Ta_name)
    ps,f,a = qcutils.GetSeriesasMA(ds,ps_name)
    Ah,f,a = qcutils.GetSeriesasMA(ds,Ah_name)
    
    if 'e' in ds.series.keys():
        e,f,a = qcutils.GetSeriesasMA(ds,'e')
    else:
        e = mf.vapourpressure(Ah,Ta)
        attr = qcutils.MakeAttributeDictionary(long_name='Vapour pressure',units='kPa',standard_name='water_vapor_partial_pressure_in_air')
        qcutils.CreateSeries(ds,'e',e,FList=[Ta_name,Ah_name],Attr=attr)
    
    if 'esat' in ds.series.keys():
        esat,f,a = qcutils.GetSeriesasMA(ds,'esat')
    else:
        esat = mf.es(Ta)
        attr = qcutils.MakeAttributeDictionary(long_name='Saturation vapour pressure (HMP)',units='kPa')
        qcutils.CreateSeries(ds,'esat',esat,FList=[Ta_name],Attr=attr)
    
    if 'rhod' not in ds.series.keys():
        rhod = mf.densitydryair(Ta,ps,e)
        attr = qcutils.MakeAttributeDictionary(long_name='Density of dry air',units='kg/m3')
        qcutils.CreateSeries(ds,'rhod',rhod,FList=[Ta_name,ps_name],Attr=attr)
    
    if 'rhom' not in ds.series.keys():
        rhom = mf.densitymoistair(Ta,ps,Ah)
        attr = qcutils.MakeAttributeDictionary(long_name='Density of moist air',units='kg/m3',standard_name='air_density')
        qcutils.CreateSeries(ds,'rhom',rhom,FList=[Ta_name,ps_name,Ah_name],Attr=attr)
    
    if 'Lv' not in ds.series.keys():
        Lv = mf.Lv(Ta)
        attr = qcutils.MakeAttributeDictionary(long_name='Latent heat of vapourisation',units='J/kg')
        qcutils.CreateSeries(ds,'Lv',Lv,FList=[Ta_name],Attr=attr)
    
    if 'mr' not in ds.series.keys():
        mr = mf.mixingratio(ps,e)
        attr = qcutils.MakeAttributeDictionary(long_name='Mixing ratio',units='kg/kg',standard_name='humidity_mixing_ratio')
        qcutils.CreateSeries(ds,'mr',mr,FList=[Ta_name,ps_name,Ah_name],Attr=attr)
    else:
        mr,f,a = qcutils.GetSeriesasMA(ds,'mr')
    
    mrsat = mf.mixingratio(ps,esat)
    if 'qsat' not in ds.series.keys():
        qsat = mf.specifichumidity(mrsat)
    else:
        qsat,f,a = qcutils.GetSeriesasMA(ds,'qsat')
    
    if 'q' not in ds.series.keys():
        q = mf.specifichumidity(mr)
        attr = qcutils.MakeAttributeDictionary(long_name='Specific humidity',units='kg/kg',standard_name='specific_humidity')
        qcutils.CreateSeries(ds,'q',q,FList=[Ta_name,ps_name,Ah_name],Attr=attr)
    else:
        q,f,a = qcutils.GetSeriesasMA(ds,'q')
    
    if 'VPD' not in ds.series.keys():
        VPD = esat - e
        attr = qcutils.MakeAttributeDictionary(long_name='Vapour pressure deficit',units='kPa',standard_name='water_vapor_saturation_deficit_in_air')
        qcutils.CreateSeries(ds,'VPD',VPD,FList=[Ta_name,Ah_name],Attr=attr)
    else:
        VPD,f,a = qcutils.GetSeriesasMA(ds,'VPD')
    
    if 'Cpm' not in ds.series.keys():
        Cpm = mf.specificheatmoistair(q)
        attr = qcutils.MakeAttributeDictionary(long_name='Specific heat of moist air',units='J/kg-K')
        qcutils.CreateSeries(ds,'Cpm',Cpm,FList=[Ta_name,ps_name,Ah_name],Attr=attr)
    
    if 'SHD' not in ds.series.keys():
        SHD = qsat - q
        attr = qcutils.MakeAttributeDictionary(long_name='Specific humidity deficit',units='kg/kg')
        qcutils.CreateSeries(ds,'SHD',SHD,FList=[Ta_name,Ah_name],Attr=attr)
    

def CalculateNDVI(cf,ds):
    red_in,f1,a = qcutils.GetSeriesasMA(ds,'red_in')
    red_ref,f2,a = qcutils.GetSeriesasMA(ds,'red_ref')
    nIR_in,f3,a = qcutils.GetSeriesasMA(ds,'nIR_in')
    nIR_ref,f4,a = qcutils.GetSeriesasMA(ds,'nIR_ref')
    timeindex = numpy.ma.where((ds.series['Hdh']['Data'] < 9) | (ds.series['Hdh']['Data'] > 15))[0]
    badindex = numpy.ma.where((numpy.mod(f1,10)!=0)|(numpy.mod(f2,10)!=0)|(numpy.mod(f3,10)!=0)|(numpy.mod(f4,10)!=0))[0]
    goodindex = numpy.ma.where((red_in!=0)&(red_ref!=0)&(nIR_in!=0)&(nIR_ref!=0))[0]
    NDVI = numpy.ma.zeros(len(red_in),dtype=numpy.float64) + numpy.float64(c.missing_value)
    NDVI[goodindex] = ((nIR_ref[goodindex]/nIR_in[goodindex])-(red_ref[goodindex]/red_in[goodindex]))/((nIR_ref[goodindex]/nIR_in[goodindex])+(red_ref[goodindex]/red_in[goodindex]))
    NDVI[badindex] = numpy.float64(c.missing_value)
    NDVI[timeindex] = numpy.float64(c.missing_value)
    attr = qcutils.MakeAttributeDictionary(long_name='normalised vegetation difference index (NDVI)',units='none',standard_name='normalized_difference_vegetation_index')
    qcutils.CreateSeries(ds,'NDVI',NDVI,FList=['red_in','red_ref','nIR_in','nIR_ref'],Attr=attr)
    ds.series['NDVI']['Flag'][timeindex] = numpy.int32(9)     # bad time flag

def CalculateNetRadiation(ds,Fn_out,Fsd_in,Fsu_in,Fld_in,Flu_in):
    """
        Calculate the net radiation from the 4 components of the surface
        radiation budget.
        
        Usage qcts.CalculateNetRadiation(ds,Fn_out,Fsd_in,Fsu_in,Fld_in,Flu_in)
        ds: data structure
        Fn_out: output net radiation variable to ds.  Example: 'Fn_KZ'
        Fsd_in: input downwelling solar radiation in ds.  Example: 'Fsd'
        Fsu_in: input upwelling solar radiation in ds.  Example: 'Fsu'
        Fld_in: input downwelling longwave radiation in ds.  Example: 'Fld'
        Flu_in: input upwelling longwave radiation in ds.  Example: 'Flu'
        """
    log.info(' Calculating net radiation from 4 components')
    if Fsd_in in ds.series.keys() and Fsu_in in ds.series.keys() and Fld_in in ds.series.keys() and Flu_in in ds.series.keys():
        Fsd,f,a = qcutils.GetSeriesasMA(ds,Fsd_in)
        Fsu,f,a = qcutils.GetSeriesasMA(ds,Fsu_in)
        Fld,f,a = qcutils.GetSeriesasMA(ds,Fld_in)
        Flu,f,a = qcutils.GetSeriesasMA(ds,Flu_in)
        Fn = (Fsd - Fsu) + (Fld - Flu)
        attr = qcutils.MakeAttributeDictionary(long_name='Calculated net radiation using '+Fsd_in+','+Fsu_in+','+Fld_in+','+Flu_in,units='W/m2',standard_name='surface_net_downwawrd_radiative_flux')
        qcutils.CreateSeries(ds,Fn_out,Fn,FList=[Fsd_in,Fsu_in,Fld_in,Flu_in],Attr=attr)
        if "Variables" in cf:
            if Fn_out in cf["Variables"]:
                if "Linear" in cf["Variables"][Fn_out]:
                    ApplyLinear(cf,ds,Fn_out)
    else:
        nRecs = numpy.int32(ds.globalattributes['nc_nrecs'])
        Fn = numpy.array([c.missing_value]*nRecs,dtype=numpy.float64)
        flag = numpy.ones(nRecs,dtype=numpy.int32)
        attr = qcutils.MakeAttributeDictionary(long_name='Calculated net radiation (one or more components missing)',
                             standard_name='surface_net_downwawrd_radiative_flux',units='W/m2')
        qcutils.CreateSeries(ds,Fn_out,Fn,Flag=flag,Attr=attr)

def CalculateSpecificHumidityProfile(cf,ds):
    if qcutils.cfkeycheck(cf,Base='qTprofile',ThisOne='ps_in'):
        ps_in = cf['qTprofile']['p_in']
    else:
        ps_in = 'ps'
    
    if qcutils.cfkeycheck(cf,Base='qTprofile',ThisOne='Ta_in'):
        Ta_in = ast.literal_eval(cf['qTprofile']['Ta_in'])
    else:
        log.error('  No input air temperature variables identified')
        return
    
    if qcutils.cfkeycheck(cf,Base='qTprofile',ThisOne='e_in'):
        e_vars = ast.literal_eval(cf['qTprofile']['e_in'])
    else:
        log.error('  No input vapour pressure variables identified')
        return
    
    if qcutils.cfkeycheck(cf,Base='qTprofile',ThisOne='esat_in'):
        esat_vars = ast.literal_eval(cf['qTprofile']['esat_in'])
    else:
        log.error('  No input saturation vapour pressure variables identified')
        return
    
    if qcutils.cfkeycheck(cf,Base='qTprofile',ThisOne='q_out'):
        q_vars = ast.literal_eval(cf['qTprofile']['q_out'])
    else:
        log.error('  No output specific humidity variables identified')
        return
    
    if qcutils.cfkeycheck(cf,Base='qTprofile',ThisOne='qsat_out'):
        qsat_vars = ast.literal_eval(cf['qTprofile']['qsat_out'])
    else:
        log.error('  No output saturation specific humidity variables identified')
        return
    
    if qcutils.cfkeycheck(cf,Base='qTprofile',ThisOne='VPD_out'):
        VPD_vars = ast.literal_eval(cf['qTprofile']['VPD_out'])
    else:
        log.error('  No output vapour pressure deficit variables identified')
        return
    
    if qcutils.cfkeycheck(cf,Base='qTprofile',ThisOne='mr_out'):
        mr_vars = ast.literal_eval(cf['qTprofile']['mr_out'])
    else:
        log.error('  No output mixing ratio variables identified')
        return
    
    if qcutils.cfkeycheck(cf,Base='qTprofile',ThisOne='Tv_out'):
        Tv_vars = ast.literal_eval(cf['qTprofile']['Tv_out'])
    else:
        log.error('  No output virtual temperature variables identified')
        return
    
    if qcutils.cfkeycheck(cf,Base='qTprofile',ThisOne='Tvp_out'):
        Tvp_vars = ast.literal_eval(cf['qTprofile']['Tvp_out'])
    else:
        log.error('  No output virtual potential temperature variables identified')
        return
    
    if qcutils.cfkeycheck(cf,Base='qTprofile',ThisOne='q_attr'):
        q_attrs = ast.literal_eval(cf['qTprofile']['q_attr'])
    else:
        log.error('  Specific humidity attributes not identified')
        return
    
    if qcutils.cfkeycheck(cf,Base='qTprofile',ThisOne='qsat_attr'):
        qsat_attrs = ast.literal_eval(cf['qTprofile']['qsat_attr'])
    else:
        log.error('  Saturated specific humidity attributes not identified')
        return
    
    if qcutils.cfkeycheck(cf,Base='qTprofile',ThisOne='VPD_attr'):
        VPD_attrs = ast.literal_eval(cf['qTprofile']['VPD_attr'])
    else:
        log.error('  Vapour pressure deficit attributes not identified')
        return
    
    if qcutils.cfkeycheck(cf,Base='qTprofile',ThisOne='mr_attr'):
        mr_attrs = ast.literal_eval(cf['qTprofile']['mr_attr'])
    else:
        log.error('  Mixing ratio attributes not identified')
        return
    
    if qcutils.cfkeycheck(cf,Base='qTprofile',ThisOne='Tv_attr'):
        Tv_attrs = ast.literal_eval(cf['qTprofile']['Tv_attr'])
    else:
        log.error('  Virtual temperature attributes not identified')
        return
    
    if qcutils.cfkeycheck(cf,Base='qTprofile',ThisOne='Tvp_attr'):
        Tvp_attrs = ast.literal_eval(cf['qTprofile']['Tvp_attr'])
    else:
        log.error('  Virtual potential temperature attributes not identified')
        return
    
    if (len(e_vars) != len(q_vars) != len(q_attrs) != len(VPD_vars) != len(mr_vars) != len(Tv_vars) != len(Tvp_vars) != len(esat_vars) != len(qsat_vars) != len(qsat_attrs) != len(VPD_attrs) != len(mr_attrs) != len(Tv_attrs) != len(Tvp_attrs)):
        log.error('  Input and output vectors not of equal length')
        return
    
    ps,f,a = qcutils.GetSeriesasMA(ds,ps_in)
    for i in range(len(e_vars)):
        Ta,f,a = qcutils.GetSeriesasMA(ds,Ta_in[i])
        e,f,a = qcutils.GetSeriesasMA(ds,e_vars[i])
        esat,f,a = qcutils.GetSeriesasMA(ds,esat_vars[i])
        mr = mf.mixingratio(ps,e)
        q = mf.specifichumidity(mr)
        mrsat = mf.mixingratio(ps,esat)
        qsat = mf.specifichumidity(mrsat)
        VPD = esat - e
        Tv = mf.theta(Ta,ps)
        Tvp = mf.virtualtheta(Tv,mr)
        attr = qcutils.MakeAttributeDictionary(long_name=q_attrs[i],units='kg/kg',standard_name='specific_humidity')
        qcutils.CreateSeries(ds,q_vars[i],q,FList=[ps_in,e_vars[i]],Attr=attr)
        attr = qcutils.MakeAttributeDictionary(long_name=qsat_attrs[i],units='kg/kg')
        qcutils.CreateSeries(ds,qsat_vars[i],qsat,FList=[ps_in,esat_vars[i]],Attr=attr)
        attr = qcutils.MakeAttributeDictionary(long_name=VPD_attrs[i],units='kPa',standard_name='water_vapor_saturation_deficit_in_air')
        qcutils.CreateSeries(ds,VPD_vars[i],VPD,FList=[ps_in,e_vars[i],esat_vars[i]],Attr=attr)
        attr = qcutils.MakeAttributeDictionary(long_name=mr_attrs[i],units='kg/kg',standard_name='humidity_mixing_ratio')
        qcutils.CreateSeries(ds,mr_vars[i],mr,FList=[ps_in,e_vars[i]],Attr=attr)
        attr = qcutils.MakeAttributeDictionary(long_name=Tv_attrs[i],units='C',standard_name='virtual_temperature')
        qcutils.CreateSeries(ds,Tv_vars[i],Tv,FList=[ps_in,Ta_in[i]],Attr=attr)
        attr = qcutils.MakeAttributeDictionary(long_name=Tvp_attrs[i],units='C',standard_name='air_potential_temperature')
        qcutils.CreateSeries(ds,Tvp_vars[i],Tvp,FList=[ps_in,e_vars[i],Ta_in[i]],Attr=attr)
    
    log.info(' q and T profile computed')
    return

def CalculateWUEfromSimilarity(cf,ds):
    Fsd,f,a = qcutils.GetSeriesasMA(ds,'Fsd')
    Ta,f,a = qcutils.GetSeriesasMA(ds,'Ta')
    Fe,f,a = qcutils.GetSeriesasMA(ds,'Fe')
    Fc,WUEflag,a = qcutils.GetSeriesasMA(ds,'Fc')
    q,f,a = qcutils.GetSeriesasMA(ds,'Ah')
    esat,f,a = qcutils.GetSeriesasMA(ds,'esat')
    Cc,f,a = qcutils.GetSeriesasMA(ds,'Cc')
    Lv,f,a = qcutils.GetSeriesasMA(ds,'Lv')
    ustar,f,a = qcutils.GetSeriesasMA(ds,'ustar')
    zeta,f,a = qcutils.GetSeriesasMA(ds,'zeta')
    cica,f,a = qcutils.GetSeriesasMA(ds,'ci_ca')
    index = numpy.ma.where(Fsd < 10)[0]
    if qcutils.cfkeycheck(cf,Base='Massman',ThisOne='zmd'):
        zmd = ast.literal_eval(cf['Massman']['zmd'])
    else:
        log.error('  zmd not defined in cf section Massman')
        return
    
    if qcutils.cfkeycheck(cf,Base='Massman',ThisOne='z0'):
        z0 = ast.literal_eval(cf['Massman']['z0'])
    else:
        log.error('  z0 not defined in cf section Massman')
        return
    
    z0v = z0*0.1
    zeta[index] = numpy.float64(c.missing_value)
    stable = numpy.ma.where(zeta > 1/16)[0]
    zeta[stable] = numpy.float64(c.missing_value)
    psyv = 2*numpy.log((1+numpy.sqrt(1-(16*zeta)))/2)
    qs = q + ((Fe/Lv)/(0.41*ustar))*(numpy.log(zmd/z0v)-psyv)
    cs = Cc + (Fc/(0.41*ustar))*(numpy.log(zmd/z0v)-psyv)
    ci = cica * cs
    qi = mf.absolutehumidity(Ta,esat)
    WUE = 0.7*(cs-ci)/(qs-qi)
    WUE[index] = numpy.float64(c.missing_value)
    WUE[stable] = numpy.float64(c.missing_value)
    WUEflag[index] = numpy.int32(8)
    WUEflag[stable] = numpy.int32(22)
    attr = qcutils.MakeAttributeDictionary(long_name='Water-use efficiency from flux-gradient similarity',units='mg/g')
    qcutils.CreateSeries(ds,'WUE',WUE,Flag=WUEflag,Attr=attr)

def ComputeClimatology(cf,ds,OList):
    if qcutils.cfkeycheck(cf,Base='Climatology',ThisOne='EF'):
        efflag = cf['Climatology']['EF']
    else:
        efflag = 'True'
    
    if qcutils.cfkeycheck(cf,Base='Climatology',ThisOne='BR'):
        brflag = cf['Climatology']['BR']
    else:
        brflag = 'True'
    
    if qcutils.cfkeycheck(cf,Base='Climatology',ThisOne='WUE'):
        wueflag = cf['Climatology']['WUE']
    else:
        wueflag = 'True'
    
    if ((efflag != 'True') & (brflag != 'True') & (wueflag != 'True') & (len(OList) == 0)):
        log.warn('  Climatology:  no dataset to generate')
        return
    
    log.info(' Computing climatology...')
    if qcutils.cfkeycheck(cf,Base='Input',ThisOne='Fn_in'):
        Fn_in = cf['Climatology']['Fn_in']
    else:
        Fn_in = 'Fn'
    
    if qcutils.cfkeycheck(cf,Base='Input',ThisOne='Fg_in'):
        Fg_in = cf['Climatology']['Fg_in']
    else:
        Fg_in = 'Fg'
    
    if qcutils.cfkeycheck(cf,Base='Input',ThisOne='Fe_in'):
        Fe_in = cf['Climatology']['Fe_in']
    else:
        Fe_in = 'Fe'
    
    if qcutils.cfkeycheck(cf,Base='Input',ThisOne='Fh_in'):
        Fh_in = cf['Climatology']['Fh_in']
    else:
        Fh_in = 'Fh'
    
    if qcutils.cfkeycheck(cf,Base='Input',ThisOne='Fc_in'):
        Fc_in = cf['Climatology']['Fc_in']
    else:
        Fc_in = 'Fc'
    
    if qcutils.cfkeycheck(cf,Base='Params',ThisOne='firstMonth'):
        M1st = numpy.int32(cf['Params']['firstMonth'])
    else:
        M1st = 1
    
    if qcutils.cfkeycheck(cf,Base='Params',ThisOne='secondMonth'):
        M2nd = numpy.int32(cf['Params']['secondMonth'])
    else:
        M2nd = 12
    
    dt = numpy.int32(ds.globalattributes['time_step'])
    xlFileName = cf['Files']['Climatology']['xlFilePath']+cf['Files']['Climatology']['xlFileName']
    xlFile = xlwt.Workbook()
    monthabr = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
    Hdh,f,a = qcutils.GetSeriesasMA(ds,'Hdh')
    for ThisOne in OList:
        log.info('  Doing climatology for '+str(ThisOne))
        xlSheet = xlFile.add_sheet(ThisOne)
        xlCol = 0
        data,f,a = qcutils.GetSeriesasMA(ds,ThisOne)
        for month in range(M1st,M2nd+1):
            mi = numpy.where(ds.series['Month']['Data']==month)[0]
            Num,Hr,Av,Sd,Mx,Mn = qcutils.get_diurnalstats(Hdh[mi],data[mi],dt)
            Num = numpy.ma.filled(Num,numpy.float64(c.missing_value))
            Hr = numpy.ma.filled(Hr,numpy.float64(c.missing_value))
            if month==1:
                xlSheet.write(1,xlCol,'Hour')
                for j in range(len(Hr)):
                    xlSheet.write(j+2,xlCol,Hr[j])
                
                xlCol = xlCol + 1
            
            xlSheet.write(0,xlCol,monthabr[month-1])
            xlSheet.write(1,xlCol,'Num')
            xlSheet.write(1,xlCol+1,'Av')
            xlSheet.write(1,xlCol+2,'Sd')
            xlSheet.write(1,xlCol+3,'Mx')
            xlSheet.write(1,xlCol+4,'Mn')
            Av = numpy.ma.filled(Av,numpy.float64(c.missing_value))
            Sd = numpy.ma.filled(Sd,numpy.float64(c.missing_value))
            Mx = numpy.ma.filled(Mx,numpy.float64(c.missing_value))
            Mn = numpy.ma.filled(Mn,numpy.float64(c.missing_value))
            for j in range(len(Hr)):
                xlSheet.write(j+2,xlCol,Num[j])
                xlSheet.write(j+2,xlCol+1,Av[j])
                xlSheet.write(j+2,xlCol+2,Sd[j])
                xlSheet.write(j+2,xlCol+3,Mx[j])
                xlSheet.write(j+2,xlCol+4,Mn[j])
            
            xlCol = xlCol + 5
    
    if efflag != 'False':
        # calculate the evaporative fraction
        xlSheet = xlFile.add_sheet('EF')
        xlCol = 0
        EF = numpy.ma.zeros([48,12]) + numpy.float64(c.missing_value)
        log.info('  Doing evaporative fraction')
        for month in range(M1st,M2nd+1):
            mi = numpy.where((ds.series['Month']['Data']==month))[0]
            Hdh = numpy.ma.masked_where(abs(ds.series['Hdh']['Data'][mi]-numpy.float64(c.missing_value))<c.eps,
                                        ds.series['Hdh']['Data'][mi])
            Fn = numpy.ma.masked_where(abs(ds.series[Fn_in]['Data'][mi]-numpy.float64(c.missing_value))<c.eps,
                                       ds.series[Fn_in]['Data'][mi])
            Fg = numpy.ma.masked_where(abs(ds.series[Fg_in]['Data'][mi]-numpy.float64(c.missing_value))<c.eps,
                                       ds.series[Fg_in]['Data'][mi])
            Fa = Fn - Fg
            Fe = numpy.ma.masked_where(abs(ds.series[Fe_in]['Data'][mi]-numpy.float64(c.missing_value))<c.eps,
                                       ds.series[Fe_in]['Data'][mi])
            Fa_Num,Hr,Fa_Av,Sd,Mx,Mn = qcutils.get_diurnalstats(Hdh,Fa,dt)
            Fe_Num,Hr,Fe_Av,Sd,Mx,Mn = qcutils.get_diurnalstats(Hdh,Fe,dt)
            index = numpy.ma.where((Fa_Num>4)&(Fe_Num>4))
            EF[:,month-1][index] = Fe_Av[index]/Fa_Av[index]
        
        # reject EF values greater than or less than 1.5
        EF = numpy.ma.masked_where(abs(EF)>1.5,EF)
        EF = numpy.ma.filled(EF,numpy.float64(c.missing_value))
        # write the EF to the Excel object
        xlSheet.write(0,xlCol,'Hour')
        for j in range(len(Hr)):
            xlSheet.write(j+1,xlCol,Hr[j])
        xlCol = xlCol + 1
        d_xf = xlwt.easyxf(num_format_str='0.00')
        for month in range(M1st,M2nd+1):
            xlSheet.write(0,xlCol,monthabr[month-1])
            for j in range(len(Hr)):
                xlSheet.write(j+1,xlCol,EF[j,month-1],d_xf)
            xlCol = xlCol + 1
        # do the 2D interpolation to fill missing EF values
        EF_3x3=numpy.tile(EF,(3,3))
        nmn=numpy.shape(EF_3x3)[1]
        mni=numpy.arange(0,nmn)
        nhr=numpy.shape(EF_3x3)[0]
        hri=numpy.arange(0,nhr)
        mn,hr=numpy.meshgrid(mni,hri)
        EF_3x3_1d=numpy.reshape(EF_3x3,numpy.shape(EF_3x3)[0]*numpy.shape(EF_3x3)[1])
        mn_1d=numpy.reshape(mn,numpy.shape(mn)[0]*numpy.shape(mn)[1])
        hr_1d=numpy.reshape(hr,numpy.shape(hr)[0]*numpy.shape(hr)[1])
        index=numpy.where(EF_3x3_1d!=c.missing_value)
        EF_3x3i=griddata(mn_1d[index],hr_1d[index],EF_3x3_1d[index],mni,hri)
        EFi=numpy.ma.filled(EF_3x3i[nhr/3:2*nhr/3,nmn/3:2*nmn/3],0)
        xlSheet = xlFile.add_sheet('EFi')
        xlCol = 0
        # write the interpolated EF values to the Excel object
        xlSheet.write(0,xlCol,'Hour')
        for j in range(len(Hr)):
            xlSheet.write(j+1,xlCol,Hr[j])
        xlCol = xlCol + 1
        d_xf = xlwt.easyxf(num_format_str='0.00')
        for month in range(M1st,M2nd+1):
            xlSheet.write(0,xlCol,monthabr[month-1])
            for j in range(len(Hr)):
                xlSheet.write(j+1,xlCol,EFi[j,month-1],d_xf)
            xlCol = xlCol + 1
    
    if brflag != 'False':
        # calculate the Bowen ratio
        xlSheet = xlFile.add_sheet('BR')
        xlCol = 0
        BR = numpy.ma.zeros([48,12]) + numpy.float64(c.missing_value)
        log.info('  Doing Bowen ratio')
        for month in range(M1st,M2nd+1):
            mi = numpy.where((ds.series['Month']['Data']==month))[0]
            Hdh = numpy.ma.masked_where(abs(ds.series['Hdh']['Data'][mi]-numpy.float64(c.missing_value))<c.eps,
                                        ds.series['Hdh']['Data'][mi])
            Fe = numpy.ma.masked_where(abs(ds.series[Fe_in]['Data'][mi]-numpy.float64(c.missing_value))<c.eps,
                                       ds.series[Fe_in]['Data'][mi])
            Fh = numpy.ma.masked_where(abs(ds.series[Fh_in]['Data'][mi]-numpy.float64(c.missing_value))<c.eps,
                                       ds.series[Fh_in]['Data'][mi])
            Fh_Num,Hr,Fh_Av,Sd,Mx,Mn = qcutils.get_diurnalstats(Hdh,Fh,dt)
            Fe_Num,Hr,Fe_Av,Sd,Mx,Mn = qcutils.get_diurnalstats(Hdh,Fe,dt)
            index = numpy.ma.where((Fh_Num>4)&(Fe_Num>4))
            BR[:,month-1][index] = Fh_Av[index]/Fe_Av[index]
        # reject BR values greater than or less than 5
        BR = numpy.ma.masked_where(abs(BR)>20.0,BR)
        BR = numpy.ma.filled(BR,numpy.float64(c.missing_value))
        # write the BR to the Excel object
        xlSheet.write(0,xlCol,'Hour')
        for j in range(len(Hr)):
            xlSheet.write(j+1,xlCol,Hr[j])
        xlCol = xlCol + 1
        d_xf = xlwt.easyxf(num_format_str='0.00')
        for month in range(M1st,M2nd+1):
            xlSheet.write(0,xlCol,monthabr[month-1])
            for j in range(len(Hr)):
                xlSheet.write(j+1,xlCol,BR[j,month-1],d_xf)
            xlCol = xlCol + 1
        # do the 2D interpolation to fill missing BR values
        # tile to 3,3 array so we have a patch in the centre, this helps
        # deal with edge effects
        BR_3x3=numpy.tile(BR,(3,3))
        nmn=numpy.shape(BR_3x3)[1]
        mni=numpy.arange(0,nmn)
        nhr=numpy.shape(BR_3x3)[0]
        hri=numpy.arange(0,nhr)
        mn,hr=numpy.meshgrid(mni,hri)
        BR_3x3_1d=numpy.reshape(BR_3x3,numpy.shape(BR_3x3)[0]*numpy.shape(BR_3x3)[1])
        mn_1d=numpy.reshape(mn,numpy.shape(mn)[0]*numpy.shape(mn)[1])
        hr_1d=numpy.reshape(hr,numpy.shape(hr)[0]*numpy.shape(hr)[1])
        index=numpy.where(BR_3x3_1d!=c.missing_value)
        BR_3x3i=griddata(mn_1d[index],hr_1d[index],BR_3x3_1d[index],mni,hri)
        BRi=numpy.ma.filled(BR_3x3i[nhr/3:2*nhr/3,nmn/3:2*nmn/3],0)
        xlSheet = xlFile.add_sheet('BRi')
        xlCol = 0
        # write the interpolated BR values to the Excel object
        xlSheet.write(0,xlCol,'Hour')
        for j in range(len(Hr)):
            xlSheet.write(j+1,xlCol,Hr[j])
        xlCol = xlCol + 1
        d_xf = xlwt.easyxf(num_format_str='0.00')
        for month in range(M1st,M2nd+1):
            xlSheet.write(0,xlCol,monthabr[month-1])
            for j in range(len(Hr)):
                xlSheet.write(j+1,xlCol,BRi[j,month-1],d_xf)
            xlCol = xlCol + 1
    
    if wueflag != 'False':
        # calculate the ecosystem water use efficiency
        xlSheet = xlFile.add_sheet('WUE')
        xlCol = 0
        WUE = numpy.ma.zeros([48,12]) + numpy.float64(c.missing_value)
        log.info('  Doing ecosystem WUE')
        for month in range(M1st,M2nd+1):
            mi = numpy.where((ds.series['Month']['Data']==month))[0]
            Hdh = numpy.ma.masked_where(abs(ds.series['Hdh']['Data'][mi]-numpy.float64(c.missing_value))<c.eps,
                                        ds.series['Hdh']['Data'][mi])
            Fe = numpy.ma.masked_where(abs(ds.series[Fe_in]['Data'][mi]-numpy.float64(c.missing_value))<c.eps,
                                       ds.series[Fe_in]['Data'][mi])
            Fc = numpy.ma.masked_where(abs(ds.series[Fc_in]['Data'][mi]-numpy.float64(c.missing_value))<c.eps,
                                       ds.series[Fc_in]['Data'][mi])
            Fc_Num,Hr,Fc_Av,Sd,Mx,Mn = qcutils.get_diurnalstats(Hdh,Fc,dt)
            Fe_Num,Hr,Fe_Av,Sd,Mx,Mn = qcutils.get_diurnalstats(Hdh,Fe,dt)
            index = numpy.ma.where((Fc_Num>4)&(Fe_Num>4))
            WUE[:,month-1][index] = Fc_Av[index]/Fe_Av[index]
        # reject WUE values greater than 0.04 or less than -0.004
        WUE = numpy.ma.masked_where((WUE>0.04)|(WUE<-0.004),WUE)
        WUE = numpy.ma.filled(WUE,numpy.float64(c.missing_value))
        # write the WUE to the Excel object
        xlSheet.write(0,xlCol,'Hour')
        for j in range(len(Hr)):
            xlSheet.write(j+1,xlCol,Hr[j])
        xlCol = xlCol + 1
        d_xf = xlwt.easyxf(num_format_str='0.00000')
        for month in range(M1st,M2nd+1):
            xlSheet.write(0,xlCol,monthabr[month-1])
            for j in range(len(Hr)):
                xlSheet.write(j+1,xlCol,WUE[j,month-1],d_xf)
            xlCol = xlCol + 1
        # do the 2D interpolation to fill missing WUE values
        WUE_3x3=numpy.tile(WUE,(3,3))
        nmn=numpy.shape(WUE_3x3)[1]
        mni=numpy.arange(0,nmn)
        nhr=numpy.shape(WUE_3x3)[0]
        hri=numpy.arange(0,nhr)
        mn,hr=numpy.meshgrid(mni,hri)
        WUE_3x3_1d=numpy.reshape(WUE_3x3,numpy.shape(WUE_3x3)[0]*numpy.shape(WUE_3x3)[1])
        mn_1d=numpy.reshape(mn,numpy.shape(mn)[0]*numpy.shape(mn)[1])
        hr_1d=numpy.reshape(hr,numpy.shape(hr)[0]*numpy.shape(hr)[1])
        index=numpy.where(WUE_3x3_1d!=c.missing_value)
        WUE_3x3i=griddata(mn_1d[index],hr_1d[index],WUE_3x3_1d[index],mni,hri)
        WUEi=numpy.ma.filled(WUE_3x3i[nhr/3:2*nhr/3,nmn/3:2*nmn/3],0)
        xlSheet = xlFile.add_sheet('WUEi')
        xlCol = 0
        # write the interpolated WUE values to the Excel object
        xlSheet.write(0,xlCol,'Hour')
        for j in range(len(Hr)):
            xlSheet.write(j+1,xlCol,Hr[j])
        xlCol = xlCol + 1
        d_xf = xlwt.easyxf(num_format_str='0.00000')
        for month in range(M1st,M2nd+1):
            xlSheet.write(0,xlCol,monthabr[month-1])
            for j in range(len(Hr)):
                xlSheet.write(j+1,xlCol,WUEi[j,month-1],d_xf)
            xlCol = xlCol + 1
    
    log.info('  Saving Excel file '+str(xlFileName))
    xlFile.save(xlFileName)
    
    log.info(' Climatology: All done')

def ComputeDailySums(cf,ds,SumList,SubSumList,MinMaxList,MeanList,SoilList):
    """
        Computes daily sums, mininima and maxima on a collection variables in
        the L4 dataset containing gap filled fluxes.  Sums are computed only
        when the number of daily 30-min observations is equal to 48 (i.e., no
        missing data) to avoid biasing.  Output to an excel file that specified
        in the control file.
        
        Usage qcts.ComputeDailySums(cf,ds)
        cf: control file
        ds: data structure
        
        Parameters loaded from control file:
            M1st: dataset start month
            M2nd: dataset end month
            SumList: list of variables to be summed
            SubSumList: list of variables to sum positive and negative observations separately
            MinMaxList: list of variables to compute daily min & max
            SoilList: list of soil moisture measurements groups
            SW0, SW10, etc: list of soil moisture sensors at a common level (e.g., surface, 10cm, etc)
        
        Default List of sums:
            Rain, ET, Fe_MJ, Fh_MJ, Fg_MJ, Fld_MJ, Flu_MJ, Fn_MJ, Fsd_MJ,
            Fsu_MJ, Fc_g, Fc_mmol
        Default List of sub-sums (sums split between positive and negative observations)
            Fe_MJ, Fh_MJ, Fg_MJ
        Default List of min/max:
            Ta_HMP, Vbat, Tpanel, Fc_mg, Fc_umol
        Default List of soil moisture measurements:
        """
    OutList = []
    SumOutList = []
    SubOutList = []
    MinMaxOutList = []
    MeanOutList = []
    SoilOutList = []
    
    for ThisOne in SubSumList:
        if ThisOne not in SumList:
            SumList.append(ThisOne)
    
    for ThisOne in SumList:
        if ThisOne == 'ET':
            if 'ET' not in ds.series.keys():
                if qcutils.cfkeycheck(cf,Base='CalculateET',ThisOne='Fe_in'):
                    Fe_in = cf['Sums']['Fe_in']
                else:
                    Fe_in = 'Fe'
                Fe,f,a = qcutils.GetSeriesasMA(ds,Fe_in)
                if 'Lv' in ds.series.keys():
                    Lv,f,a = qcutils.GetSeriesasMA(ds,'Lv')
                else:
                    Lv = c.Lv
                ET = Fe * 60 * 30 * 1000 / (Lv * c.rho_water)  # mm/30min for summing
                attr = qcutils.MakeAttributeDictionary(long_name='Evapotranspiration Flux',units='mm/30min',standard_name='lwe_water_evaporation_rate')
                qcutils.CreateSeries(ds,'ET',ET,FList=[Fe_in],Attr=attr)
            
            SumOutList.append('ET')
            OutList.append('ET')
            if ThisOne in SubSumList:
                SubOutList.append('ET')
        elif ThisOne == 'Energy':
            if qcutils.cfkeycheck(cf,Base='Sums',ThisOne='Energyin'):
                EnergyIn = ast.literal_eval(cf['Sums']['Energyin'])
            else:
                EnergyIn = ['Fe', 'Fh', 'Fg']
            Fe,f,a = qcutils.GetSeriesasMA(ds,EnergyIn[0])
            Fh,f,a = qcutils.GetSeriesasMA(ds,EnergyIn[1])
            Fg,f,a = qcutils.GetSeriesasMA(ds,EnergyIn[2])
            EnergyOut = ['Fe_MJ','Fh_MJ','Fg_MJ']
            for index in range(0,3):
                convert_energy(ds,EnergyIn[index],EnergyOut[index])
                OutList.append(EnergyOut[index])
                SumOutList.append(EnergyOut[index])
                if ThisOne in SubSumList:
                    SubOutList.append(EnergyOut[index])
        elif ThisOne == 'Radiation':
            if qcutils.cfkeycheck(cf,Base='Sums',ThisOne='Radin'):
                RadiationIn = ast.literal_eval(cf['Sums']['Radin'])
            else:
                RadiationIn = ['Fld','Flu','Fn','Fsd','Fsu']
            Fld,f,a = qcutils.GetSeriesasMA(ds,RadiationIn[0])
            Flu,f,a = qcutils.GetSeriesasMA(ds,RadiationIn[1])
            Fn,f,a = qcutils.GetSeriesasMA(ds,RadiationIn[2])
            Fsd,f,a = qcutils.GetSeriesasMA(ds,RadiationIn[3])
            Fsu,f,a = qcutils.GetSeriesasMA(ds,RadiationIn[4])
            RadiationOut = ['Fld_MJ','Flu_MJ','Fn_MJ','Fsd_MJ','Fsu_MJ']
            for index in range(0,5):
                convert_energy(ds,RadiationIn[index],RadiationOut[index])
                OutList.append(RadiationOut[index])
                SumOutList.append(RadiationOut[index])
                if ThisOne in SubSumList:
                    log.error('  Subsum: Negative radiation flux not defined')
        elif ThisOne == 'Carbon':
            Fc,f,a = qcutils.GetSeriesasMA(ds,'Fc_c')
            Fco2,f,a = qcutils.GetSeriesasMA(ds,'Fc_co2')
            Fc_umol,f,a = qcutils.GetSeriesasMA(ds,'NEE')
            Fcp_umol,f,a = qcutils.GetSeriesasMA(ds,'NEP')
            Fc_mmol = Fc_umol * 1800 / 1000                # mmol/m2-30min for summing
            Fcp_mmol = Fcp_umol * 1800 / 1000                # mmol/m2-30min for summing
            Fc_gC = Fc * 1800 / 1000                        # g/m2-30min for summing
            Fc_gCO2 = Fco2 * 1800 / 1000                        # g/m2-30min for summing
            attr = qcutils.MakeAttributeDictionary(long_name='Cumulative 30-min Flux',units='mmol/m2',standard_name='surface_upward_mole_flux_of_carbon_dioxide')
            qcutils.CreateSeries(ds,'NEE_mmol',Fc_mmol,FList=['NEE'],Attr=attr)
            attr = qcutils.MakeAttributeDictionary(long_name='Cumulative 30-min Flux',units='mmol/m2',standard_name='surface_downward_mole_flux_of_carbon_dioxide')
            qcutils.CreateSeries(ds,'NEP_mmol',Fcp_mmol,FList=['NEP'],Attr=attr)
            attr = qcutils.MakeAttributeDictionary(long_name='Cumulative 30-min Flux',units='gC/m2')
            qcutils.CreateSeries(ds,'Fc_g',Fc_gC,FList=['Fc_c'],Attr=attr)
            attr = qcutils.MakeAttributeDictionary(long_name='Cumulative 30-min Flux',units='gCO2/m2')
            qcutils.CreateSeries(ds,'Fco2_g',Fc_gCO2,FList=['Fc_co2'],Attr=attr)
            #COut = ['NEE_mmol','NEP_mmol','Fc_g','Fco2_g']
            COut = ['Fc_g']
            OutList.append(COut[0])
            SumOutList.append(COut[0])
            if ThisOne in SubSumList:
                SubOutList.append(COut[0])
            
            if 'AliceSpringsMulga' in ds.globalattributes['site_name']:
                if qcutils.cfkeycheck(cf,Base='Sums',ThisOne='GPPin'):
                    GPPIn = ast.literal_eval(cf['Sums']['GPPin'])
                    GPP,f,a = qcutils.GetSeriesasMA(ds,GPPIn[0])
                    CE,f,a = qcutils.GetSeriesasMA(ds,GPPIn[1])
                    ER_night,f,a = qcutils.GetSeriesasMA(ds,GPPIn[2])
                    ER_dark,f,a = qcutils.GetSeriesasMA(ds,GPPIn[3])
                    CE_day,f,a = qcutils.GetSeriesasMA(ds,GPPIn[4])
                    CE_NEEmax,f,a = qcutils.GetSeriesasMA(ds,GPPIn[5])
                    NEE_day,NEE_day_flag,a = qcutils.GetSeriesasMA(ds,GPPIn[6])
                    Fsd,f,a = qcutils.GetSeriesasMA(ds,'Fsd')
                    GPP_mmol = GPP * 1800 / 1000
                    GPP_gC = GPP_mmol  * c.Mc
                    GPP_gCO2 = GPP_mmol * c.Mco2
                    CE_mmol = CE * 1800 / 1000
                    CE_gC = CE_mmol * c.Mc
                    CE_gCO2 = CE_mmol * c.Mco2
                    ER_night_mmol = ER_night * 1800 / 1000
                    ER_night_gC = ER_night_mmol  * c.Mc
                    ER_night_gCO2 = ER_night_mmol * c.Mco2
                    ER_dark_mmol = ER_dark * 1800 / 1000
                    ER_dark_gC = ER_dark_mmol * c.Mc
                    ER_dark_gCO2 = ER_dark_mmol * c.Mco2
                    CE_day_mmol = CE_day * 1800 / 1000
                    CE_day_gC = CE_day_mmol * c.Mc
                    CE_day_gCO2 = CE_day_mmol * c.Mco2
                    CE_NEEmax_mmol = CE_NEEmax * 1800 / 1000
                    CE_NEEmax_gC = CE_NEEmax_mmol * c.Mc
                    CE_NEEmax_gCO2 = CE_NEEmax_mmol * c.Mco2
                    NEE_day_mmol = NEE_day * 1800 / 1000
                    NEE_day_gC = NEE_day_mmol * c.Mc
                    NEE_day_gCO2 = NEE_day_mmol * c.Mco2
                    attr = qcutils.MakeAttributeDictionary(long_name='Cumulative 30-min GPP',units='mmol/m2')
                    qcutils.CreateSeries(ds,'GPP_mmol',GPP_mmol,FList=[GPPIn[0]],Attr=attr)
                    attr = qcutils.MakeAttributeDictionary(long_name='Cumulative 30-min GPP',units='gC/m2',standard_name='gross_primary_productivity_of_carbon')
                    qcutils.CreateSeries(ds,'GPP_gC',GPP_gC,FList=[GPPIn[0]],Attr=attr)
                    attr = qcutils.MakeAttributeDictionary(long_name='Cumulative 30-min GPP',units='gCO2/m2')
                    qcutils.CreateSeries(ds,'GPP_gCO2',GPP_gCO2,FList=[GPPIn[0]],Attr=attr)
                    attr = qcutils.MakeAttributeDictionary(long_name='Cumulative 30-min ecosystem-scale CO2 efflux CE',units='mmol/m2')
                    qcutils.CreateSeries(ds,'CE_mmol',CE_mmol,FList=[GPPIn[1]],Attr=attr)
                    attr = qcutils.MakeAttributeDictionary(long_name='Cumulative 30-min ecosystem-scale CO2 efflux CE',units='gC/m2',standard_name='surface_upward_mass_flux_of_carbon_dioxide_expressed_as_carbon_due_to_emission_from_natural_sources')
                    qcutils.CreateSeries(ds,'CE_gC',CE_gC,FList=[GPPIn[1]],Attr=attr)
                    attr = qcutils.MakeAttributeDictionary(long_name='Cumulative 30-min ecosystem-scale CO2 efflux CE',units='gCO2/m2')
                    qcutils.CreateSeries(ds,'CE_gCO2',CE_gCO2,FList=[GPPIn[1]],Attr=attr)
                    attr = qcutils.MakeAttributeDictionary(long_name='Cumulative 30-min ER, nocturnal accumulation',units='mmol/m2')
                    qcutils.CreateSeries(ds,'ER_night_mmol',ER_night_mmol,FList=[GPPIn[2]],Attr=attr)
                    attr = qcutils.MakeAttributeDictionary(long_name='Cumulative 30-min ER, nocturnal accumulation',units='gC/m2',standard_name='surface_upward_mass_flux_of_carbon_dioxide_expressed_as_carbon_due_to_emission_from_natural_sources')
                    qcutils.CreateSeries(ds,'ER_night_gC',ER_night_gC,FList=[GPPIn[2]],Attr=attr)
                    attr = qcutils.MakeAttributeDictionary(long_name='Cumulative 30-min ER, nocturnal accumulation',units='gCO2/m2')
                    qcutils.CreateSeries(ds,'ER_night_gCO2',ER_night_gCO2,FList=[GPPIn[2]],Attr=attr)
                    attr = qcutils.MakeAttributeDictionary(long_name='Cumulative 30-min dark ER, intercept of light response function',units='mmol/m2')
                    qcutils.CreateSeries(ds,'ER_dark_mmol',ER_dark_mmol,FList=[GPPIn[3]],Attr=attr)
                    attr = qcutils.MakeAttributeDictionary(long_name='Cumulative 30-min dark ER, intercept of light response function',units='gC/m2',standard_name='surface_upward_mass_flux_of_carbon_dioxide_expressed_as_carbon_due_to_emission_from_natural_sources')
                    qcutils.CreateSeries(ds,'ER_dark_gC',ER_dark_gC,FList=[GPPIn[3]],Attr=attr)
                    attr = qcutils.MakeAttributeDictionary(long_name='Cumulative 30-min dark ER, intercept of light response function',units='gCO2/m2')
                    qcutils.CreateSeries(ds,'ER_dark_gCO2',ER_dark_gCO2,FList=[GPPIn[3]],Attr=attr)
                    attr = qcutils.MakeAttributeDictionary(long_name='Cumulative 30-min CE, daytime accumulation',units='mmol/m2')
                    qcutils.CreateSeries(ds,'CE_day_mmol',CE_day_mmol,FList=[GPPIn[5]],Attr=attr)
                    attr = qcutils.MakeAttributeDictionary(long_name='Cumulative 30-min CE, daytime accumulation',units='gC/m2',standard_name='surface_upward_mass_flux_of_carbon_dioxide_expressed_as_carbon_due_to_emission_from_natural_sources')
                    qcutils.CreateSeries(ds,'CE_day_gC',CE_day_gC,FList=[GPPIn[5]],Attr=attr)
                    attr = qcutils.MakeAttributeDictionary(long_name='Cumulative 30-min CE, daytime accumulation',units='gCO2/m2')
                    qcutils.CreateSeries(ds,'CE_day_gCO2',CE_day_gCO2,FList=[GPPIn[5]],Attr=attr)
                    attr = qcutils.MakeAttributeDictionary(long_name='Cumulative 30-min CE, maximal NEE rates in the absence of GPP as a function of temperature',units='mmol/m2')
                    qcutils.CreateSeries(ds,'CE_NEEmax_mmol',CE_NEEmax_mmol,FList=[GPPIn[4]],Attr=attr)
                    attr = qcutils.MakeAttributeDictionary(long_name='Cumulative 30-min CE, maximal NEE rates in the absence of GPP as a function of temperature',units='gC/m2',standard_name='surface_upward_mass_flux_of_carbon_dioxide_expressed_as_carbon_due_to_emission_from_natural_sources')
                    qcutils.CreateSeries(ds,'CE_NEEmax_gC',CE_NEEmax_gC,FList=[GPPIn[4]],Attr=attr)
                    attr = qcutils.MakeAttributeDictionary(long_name='Cumulative 30-min CE, maximal NEE rates in the absence of GPP as a function of temperature',units='gCO2/m2')
                    qcutils.CreateSeries(ds,'CE_NEEmax_gCO2',CE_NEEmax_gCO2,FList=[GPPIn[4]],Attr=attr)
                    attr = qcutils.MakeAttributeDictionary(long_name='Cumulative 30-min NEE, daytime accumulation',units='mmol/m2')
                    qcutils.CreateSeries(ds,'NEE_day_mmol',NEE_day_mmol,Flag=NEE_day_flag,Attr=attr)
                    attr = qcutils.MakeAttributeDictionary(long_name='Cumulative 30-min NEE, daytime accumulation',units='gC/m2',standard_name='surface_upward_mass_flux_of_carbon_dioxide_expressed_as_carbon_due_to_emission_from_natural_sources')
                    qcutils.CreateSeries(ds,'NEE_day_gC',NEE_day_gC,Flag=NEE_day_flag,Attr=attr)
                    attr = qcutils.MakeAttributeDictionary(long_name='Cumulative 30-min NEE, daytime accumulation',units='gCO2/m2')
                    qcutils.CreateSeries(ds,'NEE_day_gCO2',NEE_day_gCO2,Flag=NEE_day_flag,Attr=attr)
                    #GPPOut = ['GPP_mmol','GPP_gC','GPP_gCO2','CE_mmol','CE_gC','CE_gCO2','ER_night_mmol','ER_night_gC','ER_night_gCO2','ER_dark_mmol','ER_dark_gC','ER_dark_gCO2','CE_day_mmol','CE_day_gC','CE_day_gCO2','CE_NEEmax_mmol','CE_NEEmax_gC','CE_NEEmax_gCO2','NEE_day_mmol','NEE_day_gC','NEE_day_gCO2']
                    GPPOut = ['GPP_gC','CE_gC','ER_night_gC','ER_dark_gC','CE_day_gC','CE_NEEmax_gC','NEE_day_gC']
                    for listindex in range(0,6):
                        OutList.append(GPPOut[listindex])
                        SumOutList.append(GPPOut[listindex])
            elif 'TiTreeEast' in ds.globalattributes['site_name']:
                if qcutils.cfkeycheck(cf,Base='Sums',ThisOne='GPPin'):
                    GPPIn = ast.literal_eval(cf['Sums']['GPPin'])
                    GPP,f,a = qcutils.GetSeriesasMA(ds,GPPIn[0])
                    CE,f,a = qcutils.GetSeriesasMA(ds,GPPIn[1])
                    ER_night,f,a = qcutils.GetSeriesasMA(ds,GPPIn[2])
                    ER_dark,f,a = qcutils.GetSeriesasMA(ds,GPPIn[3])
                    ER_bio,f,a = qcutils.GetSeriesasMA(ds,GPPIn[4])
                    CE_day,f,a = qcutils.GetSeriesasMA(ds,GPPIn[5])
                    PD,f,a = qcutils.GetSeriesasMA(ds,GPPIn[6])
                    NEE_day,NEE_day_flag,a = qcutils.GetSeriesasMA(ds,GPPIn[7])
                    Fsd,f,a = qcutils.GetSeriesasMA(ds,'Fsd')
                    GPP_mmol = GPP * 1800 / 1000
                    GPP_gC = GPP_mmol  * c.Mc
                    GPP_gCO2 = GPP_mmol * c.Mco2
                    CE_mmol = CE * 1800 / 1000
                    CE_gC = CE_mmol * c.Mc
                    CE_gCO2 = CE_mmol * c.Mco2
                    ER_night_mmol = ER_night * 1800 / 1000
                    ER_night_gC = ER_night_mmol  * c.Mc
                    ER_night_gCO2 = ER_night_mmol * c.Mco2
                    ER_dark_mmol = ER_dark * 1800 / 1000
                    ER_dark_gC = ER_dark_mmol * c.Mc
                    ER_dark_gCO2 = ER_dark_mmol * c.Mco2
                    ER_bio_mmol = ER_bio * 1800 / 1000
                    ER_bio_gC = ER_bio_mmol * c.Mc
                    ER_bio_gCO2 = ER_bio_mmol * c.Mco2
                    CE_day_mmol = CE_day * 1800 / 1000
                    CE_day_gC = CE_day_mmol * c.Mc
                    CE_day_gCO2 = CE_day_mmol * c.Mco2
                    PD_mmol = PD * 1800 / 1000
                    PD_gC = PD_mmol * c.Mc
                    PD_gCO2 = PD_mmol * c.Mco2
                    NEE_day_mmol = NEE_day * 1800 / 1000
                    NEE_day_gC = NEE_day_mmol * c.Mc
                    NEE_day_gCO2 = NEE_day_mmol * c.Mco2
                    attr = qcutils.MakeAttributeDictionary(long_name='Cumulative 30-min GPP',units='mmol/m2')
                    qcutils.CreateSeries(ds,'GPP_mmol',GPP_mmol,FList=[GPPIn[0]],Attr=attr)
                    attr = qcutils.MakeAttributeDictionary(long_name='Cumulative 30-min GPP',units='gC/m2',standard_name='gross_primary_productivity_of_carbon')
                    qcutils.CreateSeries(ds,'GPP_gC',GPP_gC,FList=[GPPIn[0]],Attr=attr)
                    attr = qcutils.MakeAttributeDictionary(long_name='Cumulative 30-min GPP',units='gCO2/m2')
                    qcutils.CreateSeries(ds,'GPP_gCO2',GPP_gCO2,FList=[GPPIn[0]],Attr=attr)
                    attr = qcutils.MakeAttributeDictionary(long_name='Cumulative 30-min ecosystem CO2 efflux CE',units='mmol/m2')
                    qcutils.CreateSeries(ds,'CE_mmol',CE_mmol,FList=[GPPIn[1]],Attr=attr)
                    attr = qcutils.MakeAttributeDictionary(long_name='Cumulative 30-min ecosystem CO2 efflux CE',units='gC/m2',standard_name='surface_upward_mass_flux_of_carbon_dioxide_expressed_as_carbon_due_to_emission_from_natural_sources')
                    qcutils.CreateSeries(ds,'CE_gC',CE_gC,FList=[GPPIn[1]],Attr=attr)
                    attr = qcutils.MakeAttributeDictionary(long_name='Cumulative 30-min ecosystem CO2 efflux CE',units='gCO2/m2')
                    qcutils.CreateSeries(ds,'CE_gCO2',CE_gCO2,FList=[GPPIn[1]],Attr=attr)
                    attr = qcutils.MakeAttributeDictionary(long_name='Cumulative 30-min ER, nocturnal accumulation',units='mmol/m2')
                    qcutils.CreateSeries(ds,'ER_night_mmol',ER_night_mmol,FList=[GPPIn[2]],Attr=attr)
                    attr = qcutils.MakeAttributeDictionary(long_name='Cumulative 30-min ER, nocturnal accumulation',units='gC/m2',standard_name='surface_upward_mass_flux_of_carbon_dioxide_expressed_as_carbon_due_to_emission_from_natural_sources')
                    qcutils.CreateSeries(ds,'ER_night_gC',ER_night_gC,FList=[GPPIn[2]],Attr=attr)
                    attr = qcutils.MakeAttributeDictionary(long_name='Cumulative 30-min ER, nocturnal accumulation',units='gCO2/m2')
                    qcutils.CreateSeries(ds,'ER_night_gCO2',ER_night_gCO2,FList=[GPPIn[2]],Attr=attr)
                    attr = qcutils.MakeAttributeDictionary(long_name='Cumulative 30-min dark ER, intercept of light response function',units='mmol/m2')
                    qcutils.CreateSeries(ds,'ER_dark_mmol',ER_dark_mmol,FList=[GPPIn[3]],Attr=attr)
                    attr = qcutils.MakeAttributeDictionary(long_name='Cumulative 30-min dark ER, intercept of light response function',units='gC/m2',standard_name='surface_upward_mass_flux_of_carbon_dioxide_expressed_as_carbon_due_to_emission_from_natural_sources')
                    qcutils.CreateSeries(ds,'ER_dark_gC',ER_dark_gC,FList=[GPPIn[3]],Attr=attr)
                    attr = qcutils.MakeAttributeDictionary(long_name='Cumulative 30-min dark ER, intercept of light response function',units='gCO2/m2')
                    qcutils.CreateSeries(ds,'ER_dark_gCO2',ER_dark_gCO2,FList=[GPPIn[3]],Attr=attr)
                    attr = qcutils.MakeAttributeDictionary(long_name='Cumulative 30-min ER, biological sources (AR+HR)',units='mmol/m2')
                    qcutils.CreateSeries(ds,'ER_bio_mmol',ER_bio_mmol,FList=[GPPIn[4]],Attr=attr)
                    attr = qcutils.MakeAttributeDictionary(long_name='Cumulative 30-min ER, biological sources (AR+HR)',units='gC/m2',standard_name='surface_upward_mass_flux_of_carbon_dioxide_expressed_as_carbon_due_to_emission_from_natural_sources')
                    qcutils.CreateSeries(ds,'ER_bio_gC',ER_bio_gC,FList=[GPPIn[4]],Attr=attr)
                    attr = qcutils.MakeAttributeDictionary(long_name='Cumulative 30-min ER, biological sources (AR+HR)',units='gCO2/m2')
                    qcutils.CreateSeries(ds,'ER_bio_gCO2',ER_bio_gCO2,FList=[GPPIn[4]],Attr=attr)
                    attr = qcutils.MakeAttributeDictionary(long_name='Cumulative 30-min CE, daytime accumulation',units='mmol/m2')
                    qcutils.CreateSeries(ds,'CE_day_mmol',CE_day_mmol,FList=[GPPIn[5]],Attr=attr)
                    attr = qcutils.MakeAttributeDictionary(long_name='Cumulative 30-min CE, daytime accumulation',units='gC/m2',standard_name='surface_upward_mass_flux_of_carbon_dioxide_expressed_as_carbon_due_to_emission_from_natural_sources')
                    qcutils.CreateSeries(ds,'CE_day_gC',CE_day_gC,FList=[GPPIn[5]],Attr=attr)
                    attr = qcutils.MakeAttributeDictionary(long_name='Cumulative 30-min CE, daytime accumulation',units='gCO2/m2')
                    qcutils.CreateSeries(ds,'CE_day_gCO2',CE_day_gCO2,FList=[GPPIn[5]],Attr=attr)
                    attr = qcutils.MakeAttributeDictionary(long_name='Cumulative 30-min photo-degradation',units='mmol/m2')
                    qcutils.CreateSeries(ds,'PD_mmol',PD_mmol,FList=[GPPIn[6]],Attr=attr)
                    attr = qcutils.MakeAttributeDictionary(long_name='Cumulative 30-min photo-degradation',units='gC/m2',standard_name='surface_upward_mass_flux_of_carbon_dioxide_expressed_as_carbon_due_to_emission_from_natural_sources')
                    qcutils.CreateSeries(ds,'PD_gC',PD_gC,FList=[GPPIn[6]],Attr=attr)
                    attr = qcutils.MakeAttributeDictionary(long_name='Cumulative 30-min photo-degradation',units='gCO2/m2')
                    qcutils.CreateSeries(ds,'PD_gCO2',PD_gCO2,FList=[GPPIn[6]],Attr=attr)
                    attr = qcutils.MakeAttributeDictionary(long_name='Cumulative 30-min NEE, daytime accumulation',units='mmol/m2')
                    qcutils.CreateSeries(ds,'NEE_day_mmol',NEE_day_mmol,Flag=NEE_day_flag,Attr=attr)
                    attr = qcutils.MakeAttributeDictionary(long_name='Cumulative 30-min NEE, daytime accumulation',units='gC/m2',standard_name='surface_upward_mass_flux_of_carbon_dioxide_expressed_as_carbon_due_to_emission_from_natural_sources')
                    qcutils.CreateSeries(ds,'NEE_day_gC',NEE_day_gC,Flag=NEE_day_flag,Attr=attr)
                    attr = qcutils.MakeAttributeDictionary(long_name='Cumulative 30-min NEE, daytime accumulation',units='gCO2/m2')
                    qcutils.CreateSeries(ds,'NEE_day_gCO2',NEE_day_gCO2,Flag=NEE_day_flag,Attr=attr)
                    #GPPOut = ['GPP_mmol','GPP_gC','GPP_gCO2','CE_mmol','CE_gC','CE_gCO2','ER_night_mmol','ER_night_gC','ER_night_gCO2','ER_dark_mmol','ER_dark_gC','ER_dark_gCO2','ER_bio_mmol','ER_bio_gC','ER_bio_gCO2','CE_day_mmol','CE_day_gC','CE_day_gCO2','PD_mmol','PD_gC','PD_gCO2','NEE_day_mmol','NEE_day_gC','NEE_day_gCO2']
                    GPPOut = ['GPP_gC','CE_gC','ER_night_gC','ER_dark_gC','ER_bio_gC','CE_day_gC','PD_gC','NEE_day_gC']
                    for listindex in range(0,7):
                        OutList.append(GPPOut[listindex])
                        SumOutList.append(GPPOut[listindex])
        elif ThisOne == 'PM':
            if 'GSv_2layer' not in ds.series.keys() and 'GSv_1layer' not in ds.series.keys() and 'GSm' not in ds.series.keys():
                SumList.remove('PM')
                log.error('  Penman-Monteith Daily sum: input Gst or Gc not located')
            
            if qcutils.cfkeycheck(cf,Base='PenmanMonteith',ThisOne='Cemethod') and cf['PenmanMonteith']['Cemethod'] == 'True':
                if 'GSv_1layer' in ds.series.keys():
                    Gst_1layer_mmol,f,a = qcutils.GetSeriesasMA(ds,'GSv_1layer')   # mmol/m2-s
                    Gst_1layer_mol =  Gst_1layer_mmol * 1800 / 1000                 # mol/m2-30min for summing
                    attr = qcutils.MakeAttributeDictionary(long_name='Cumulative 30-min Bulk Stomatal Conductance, 1-layer Ce method Penman-Monteith',units='mol/m2')
                    qcutils.CreateSeries(ds,'GSv_1layer_mol',Gst_1layer_mol,FList=['GSv_1layer'],Attr=attr)
                    OutList.append('GSv_1layer_mol')
                    if 'GSv_1layer_mol' in SubSumList:
                        log.error('  Subsum: Negative bulk stomatal conductance not defined')
                    SumOutList.append('GSv_1layer_mol')
                else:
                    log.error('  Penman-Monteith Daily sum: input Gst_1layer not located')
            
            if qcutils.cfkeycheck(cf,Base='PenmanMonteith',ThisOne='Ce_2layer') and cf['PenmanMonteith']['Ce_2layer'] == 'True':
                if 'GSv_2layer' in ds.series.keys() and 'GSv_top' in ds.series.keys() and 'GSv_base' in ds.series.keys() and 'GSv_full' in ds.series.keys():
                    for ThisOne in ['GSv_2layer','GSv_top','GSv_base','GSv_full']:
                        Gst_2layer_mmol,f,a = qcutils.GetSeriesasMA(ds,ThisOne)   # mmol/m2-s
                        Gst_2layer_mol =  Gst_2layer_mmol * 1800 / 1000                 # mol/m2-30min for summing
                        newvar = ThisOne + '_mol'
                        attr = qcutils.MakeAttributeDictionary(long_name='Cumulative 30-min Bulk Stomatal Conductance, 2-layer Ce method Penman-Monteith',units='mol/m2')
                        qcutils.CreateSeries(ds,newvar,Gst_2layer_mol,FList=[ThisOne],Attr=attr)
                        OutList.append(newvar)
                        if newvar in SubSumList:
                            log.error('  Subsum: Negative bulk stomatal conductance not defined')
                        SumOutList.append(newvar)
                else:
                    log.error('  Penman-Monteith Daily sum: input Gst_2layer not located')
            
            if qcutils.cfkeycheck(cf,Base='PenmanMonteith',ThisOne='Cdmethod') and cf['PenmanMonteith']['Cdmethod'] == 'True':
                if 'GSm' in ds.series.keys():
                    Gc_mmol,f,a = qcutils.GetSeriesasMA(ds,'GSm')   # mmol/m2-s
                    Gc_mol =  Gc_mmol * 1800 / 1000                 # mol/m2-30min for summing
                    attr = qcutils.MakeAttributeDictionary(long_name='Cumulative 30-min Canopy Conductance',units='mol/m2')
                    qcutils.CreateSeries(ds,'GSm_mol',Gc_mol,FList=['GSm'],Attr=attr)
                    OutList.append('GSm_mol')
                    if 'GSm_mol' in SubSumList:
                        log.error('  Subsum: Negative bulk stomatal conductance not defined')
                    SumOutList.append('GSm_mol')
                else:
                    log.error('  Penman-Monteith Daily sum: input Gc not located')
        elif ThisOne == 'Rainhours':
            rain,f,a = qcutils.GetSeriesasMA(ds,'Precip')
            day,f,a = qcutils.GetSeriesasMA(ds,'Ddd')
            rainhr = numpy.zeros(len(day), dtype=numpy.float64)
            for i in range(numpy.int32(numpy.floor(day[0])),numpy.int32(numpy.floor(day[-1]))):
                index = numpy.where(((day - i) < 1) & (((day - i) > 0) | ((day - i) == 0)))[0]
                for j in range(len(index)):
                    if rain[index[j]] > 0:
                        rainhr[index[j]] = 0.5
            attr = qcutils.MakeAttributeDictionary(long_name='Hourly rainfall occurrence (1) or absence (0)',units='hr')
            qcutils.CreateSeries(ds,'rainhrs',rainhr,FList=['Precip'],Attr=attr)
            OutList.append('rainhrs')
            SumOutList.append('rainhrs')
        else:
            OutList.append(ThisOne)
            SumOutList.append(ThisOne)
    
    for ThisOne in MinMaxList:
        if ThisOne == 'Carbon':
            COut = ['NEE','NEP','Fc_c','Fc_co2']
            for listindex in range(0,4):
                OutList.append(COut[listindex])
                MinMaxOutList.append(COut[listindex])
        elif ThisOne == 'PM':
            if 'GSv_1layer' not in ds.series.keys() and 'GSv_2layer' not in ds.series.keys() and 'GSv_base' not in ds.series.keys() and 'GSv_top' not in ds.series.keys() and 'GSm' not in ds.series.keys() and 'rSv_1layer' not in ds.series.keys() and 'rSv_2layer' not in ds.series.keys() and 'rSv_base' not in ds.series.keys() and 'rSv_top' not in ds.series.keys() and 'rSm' not in ds.series.keys() and 'rav_1layer' not in ds.series.keys() and 'rav_2layer' not in ds.series.keys() and 'rav_base' not in ds.series.keys() and 'rav_top' not in ds.series.keys() and 'ram' not in ds.series.keys():
                MinMaxList.remove('PM')
                log.error('  Penman-Monteith Daily min/max: input Gst, rst, rc or Gc not located')
            
            if qcutils.cfkeycheck(cf,Base='PenmanMonteith',ThisOne='Cemethod') and cf['PenmanMonteith']['Cemethod'] == 'True':
                if 'GSv_1layer' in ds.series.keys() and 'rSv_1layer' in ds.series.keys():
                    PMout = ['rSv_1layer','GSv_1layer']
                    for listindex in range(0,2):
                        if PMout[listindex] not in OutList:
                            OutList.append(PMout[listindex])
                        MinMaxOutList.append(PMout[listindex])
                else:
                    log.error('  Penman-Monteith Daily min/max: input Gst_1layer or rst_1layer not located')
            
            if qcutils.cfkeycheck(cf,Base='PenmanMonteith',ThisOne='Ce_2layer') and cf['PenmanMonteith']['Ce_2layer'] == 'True':
                if 'rav_2layer' in ds.series.keys() and 'GSv_2layer' in ds.series.keys() and 'rSv_2layer' in ds.series.keys() and 'rav_base' in ds.series.keys() and 'GSv_base' in ds.series.keys() and 'rSv_base' in ds.series.keys() and 'rav_top' in ds.series.keys() and 'GSv_top' in ds.series.keys() and 'rSv_top' in ds.series.keys() and 'rav_full' in ds.series.keys() and 'GSv_full' in ds.series.keys() and 'rSv_full' in ds.series.keys():
                    PMout = ['rav_2layer','rSv_2layer','GSv_2layer','rav_base','rSv_base','GSv_base','rav_top','rSv_top','GSv_top','rav_full','rSv_full','GSv_full']
                    for listindex in range(0,12):
                        if PMout[listindex] not in OutList:
                            OutList.append(PMout[listindex])
                        MinMaxOutList.append(PMout[listindex])
                else:
                    log.error('  Penman-Monteith Daily mean: input rav_2layer, Gst_2layer, rst_2layer, rav_base, Gst_base, rst_base, rav_top, Gst_top, rst_full, rav_full, Gst_full or rst_full not located')
            
            if qcutils.cfkeycheck(cf,Base='PenmanMonteith',ThisOne='Cdmethod') and cf['PenmanMonteith']['Cdmethod'] == 'True':
                if 'GSm' in ds.series.keys() and 'rSm' in ds.series.keys():
                    PMout = ['ram','rSm','GSm']
                    for listindex in range(0,3):
                        if PMout[listindex] not in OutList:
                            OutList.append(PMout[listindex])
                        MinMaxOutList.append(PMout[listindex])
                else:
                    log.error('  Penman-Monteith Daily min/max: input Gc or rc not located')
        else:
            if ThisOne not in OutList:
                OutList.append(ThisOne)
            MinMaxOutList.append(ThisOne)
    
    for ThisOne in MeanList:
        if ThisOne == 'Energy' or ThisOne == 'Carbon' or ThisOne == 'Radiation':
            log.error(' Mean error: '+ThisOne+' to be placed in SumList')
        elif ThisOne == 'PM':
            if 'GSv_1layer' not in ds.series.keys() and 'GSv_2layer' not in ds.series.keys() and 'GSv_base' not in ds.series.keys() and 'GSv_top' not in ds.series.keys() and 'GSm' not in ds.series.keys() and 'rSv_1layer' not in ds.series.keys() and 'rSv_2layer' not in ds.series.keys() and 'rSv_base' not in ds.series.keys() and 'rSv_top' not in ds.series.keys() and 'rSm' not in ds.series.keys() and 'rav_1layer' not in ds.series.keys() and 'rav_2layer' not in ds.series.keys() and 'rav_base' not in ds.series.keys() and 'rav_top' not in ds.series.keys() and 'ram' not in ds.series.keys():
                MeanList.remove('PM')
                log.error('  Penman-Monteith Daily mean: input Gst, rst, rc or Gc not located')
            
            if qcutils.cfkeycheck(cf,Base='PenmanMonteith',ThisOne='Cemethod') and cf['PenmanMonteith']['Cemethod'] == 'True':
                if 'GSv_1layer' in ds.series.keys() and 'rSv_1layer' in ds.series.keys():
                    PMout = ['rSv_1layer','GSv_1layer']
                    for listindex in range(0,2):
                        if PMout[listindex] not in OutList:
                            OutList.append(PMout[listindex])
                        MeanOutList.append(PMout[listindex])
                else:
                    log.error('  Penman-Monteith Daily mean: input Gst_1layer or rst_1layer not located')
                
            if qcutils.cfkeycheck(cf,Base='PenmanMonteith',ThisOne='Ce_2layer') and cf['PenmanMonteith']['Ce_2layer'] == 'True':
                if 'rav_2layer' in ds.series.keys() and 'GSv_2layer' in ds.series.keys() and 'rSv_2layer' in ds.series.keys() and 'rav_base' in ds.series.keys() and 'GSv_base' in ds.series.keys() and 'rSv_base' in ds.series.keys() and 'rav_top' in ds.series.keys() and 'GSv_top' in ds.series.keys() and 'rSv_top' in ds.series.keys() and 'rav_full' in ds.series.keys() and 'GSv_full' in ds.series.keys() and 'rSv_full' in ds.series.keys():
                    PMout = ['rav_2layer','rSv_2layer','GSv_2layer','rav_base','rSv_base','GSv_base','rav_top','rSv_top','GSv_top','rav_full','rSv_full','GSv_full']
                    for listindex in range(0,12):
                        if PMout[listindex] not in OutList:
                            OutList.append(PMout[listindex])
                        MeanOutList.append(PMout[listindex])
                else:
                    log.error('  Penman-Monteith Daily mean: input rav_2layer, Gst_2layer, rst_2layer, rav_base, Gst_base, rst_base, rav_top, Gst_top or rst_top not located')
            
            if qcutils.cfkeycheck(cf,Base='PenmanMonteith',ThisOne='Cdmethod') and cf['PenmanMonteith']['Cdmethod'] == 'True':
                if 'GSm' in ds.series.keys() and 'rSm' in ds.series.keys():
                    PMout = ['ram','rSm','GSm']
                    for listindex in range(0,3):
                        if PMout[listindex] not in OutList:
                            OutList.append(PMout[listindex])
                        MeanOutList.append(PMout[listindex])
                else:
                    log.error('  Penman-Monteith Daily mean: input Gc or rc not located')
        else:
            MeanOutList.append(ThisOne)
            if ThisOne not in OutList:
                OutList.append(ThisOne)
    
    if len(SoilList) > 0:
        for ThisOne in SoilList:
            if qcutils.cfkeycheck(cf,Base='Sums',ThisOne=ThisOne):
                vars = ast.literal_eval(cf['Sums'][ThisOne])
                for index in range(0,len(vars)):
                    SoilOutList.append(vars[index])
                OutList.append(ThisOne)
    
    outlevel = cf['Functions']['Sums']
    xlFileName = cf['Files'][outlevel]['xlSumFilePath']+cf['Files'][outlevel]['xlSumFileName']
    xlFile = xlwt.Workbook()
    
    for ThisOne in OutList:
        xlSheet = xlFile.add_sheet(ThisOne)
        xlCol = 0
        if ThisOne in SumOutList:
            if ThisOne in SubOutList:
                write_sums(cf,ds,ThisOne,xlCol,xlSheet,DoSum='True',DoSubSum='True')
            else:
                write_sums(cf,ds,ThisOne,xlCol,xlSheet,DoSum='True')
        
        if ThisOne in MinMaxOutList:
            if ThisOne in MeanOutList:
                write_sums(cf,ds,ThisOne,xlCol,xlSheet,DoMinMax='True',DoMean='True')
            else:
                write_sums(cf,ds,ThisOne,xlCol,xlSheet,DoMinMax='True')
        
        if ThisOne in MeanOutList:
            if ThisOne not in MinMaxOutList:
                write_sums(cf,ds,ThisOne,xlCol,xlSheet,DoMean='True')
        
        if ThisOne in SoilList:
            soilvars = ast.literal_eval(cf['Sums'][ThisOne])
            for n in soilvars:
                if n == soilvars[0]:
                    xC,xS = write_sums(cf,ds,n,xlCol,xlSheet,DoSoil='True')
                else:
                    xC,xS = write_sums(cf,ds,n,xlCol,xS,DoSoil='True')
                xlCol = xC + 1
        
    log.info(' Saving Excel file '+xlFileName)
    xlFile.save(xlFileName)

    log.info(' Daily sums: All done')

def convert_energy(ds,InVar,OutVar):
    """
        Integrate energy flux over 30-min time period.
        Converts flux in W/m2 to MJ/(m2 30-min)
        
        Usage qcts.convert_energy(ds,InVar,OutVar)
        ds: data structure
        InVar: name of input variable.  Example: 'Fe_gapfilled'
        OutVar: name of output variable.  Example: 'Fe_MJ'
        """
    Wm2,f,a = qcutils.GetSeriesasMA(ds,InVar)
    MJ = Wm2 * 1800 / 1e6
    attr = qcutils.MakeAttributeDictionary(long_name=ds.series[InVar]['Attr']['long_name'],units='MJ/m2',standard_name=ds.series[InVar]['Attr']['standard_name'])
    qcutils.CreateSeries(ds,OutVar,MJ,FList=[InVar],Attr=attr)

def CoordRotation2D(cf,ds):
    """
        2D coordinate rotation to force v = w = 0.  Based on Lee et al, Chapter
        3 of Handbook of Micrometeorology.  This routine does not do the third
        rotation to force v'w' = 0.
        
        Usage qcts.CoordRotation2D(ds)
        ds: data structure
        """
    log.info(' Applying 2D coordinate rotation to wind components and covariances')
    # get the raw wind velocity components
    Ux,f,a = qcutils.GetSeriesasMA(ds,'Ux')          # longitudinal component in CSAT coordinate system
    Uy,f,a = qcutils.GetSeriesasMA(ds,'Uy')          # lateral component in CSAT coordinate system
    Uz,f,a = qcutils.GetSeriesasMA(ds,'Uz')          # vertical component in CSAT coordinate system
    # get the raw covariances
    UxUz,f,a = qcutils.GetSeriesasMA(ds,'UxUz')      # covariance(Ux,Uz)
    UyUz,f,a = qcutils.GetSeriesasMA(ds,'UyUz')      # covariance(Uy,Uz)
    UxUy,f,a = qcutils.GetSeriesasMA(ds,'UxUy')      # covariance(Ux,Uy)
    UyUy,f,a = qcutils.GetSeriesasMA(ds,'UyUy')      # variance(Uy)
    UxUx,f,a = qcutils.GetSeriesasMA(ds,'UxUx')      # variance(Ux)
    UzUz,f,a = qcutils.GetSeriesasMA(ds,'UzUz')      # variance(Uz)
    UzC,f,a = qcutils.GetSeriesasMA(ds,'UzC')        # covariance(Uz,C)
    UzA,f,a = qcutils.GetSeriesasMA(ds,'UzA')        # covariance(Uz,A)
    UzT,f,a = qcutils.GetSeriesasMA(ds,'UzT')        # covariance(Uz,T)
    UxC,f,a = qcutils.GetSeriesasMA(ds,'UxC')        # covariance(Ux,C)
    UyC,f,a = qcutils.GetSeriesasMA(ds,'UyC')        # covariance(Uy,C)
    UxA,f,a = qcutils.GetSeriesasMA(ds,'UxA')        # covariance(Ux,A)
    UyA,f,a = qcutils.GetSeriesasMA(ds,'UyA')        # covariance(Ux,A)
    UxT,f,a = qcutils.GetSeriesasMA(ds,'UxT')        # covariance(Ux,T)
    UyT,f,a = qcutils.GetSeriesasMA(ds,'UyT')        # covariance(Uy,T)
    nRecs = numpy.int32(ds.globalattributes['nc_nrecs'])   # number of records
    # get the 2D and 3D wind speeds
    ws2d = numpy.ma.sqrt(Ux**2 + Uy**2)
    ws3d = numpy.ma.sqrt(Ux**2 + Uy**2 + Uz**2)
    # get the sine and cosine of the angles through which to rotate
    #  - first we rotate about the Uz axis by eta to get v = 0
    #  - then we rotate about the v axis by theta to get w = 0
    ce = Ux/ws2d          # cos(eta)
    se = Uy/ws2d          # sin(eta)
    ct = ws2d/ws3d        # cos(theta)
    st = Uz/ws3d          # sin(theta)
    # get the rotation angles
    theta = numpy.rad2deg(numpy.arctan2(st,ct))
    eta = numpy.rad2deg(numpy.arctan2(se,ce))
    # do the wind velocity components first
    u = Ux*ct*ce + Uy*ct*se + Uz*st           # longitudinal component in natural wind coordinates
    v = Uy*ce - Ux*se                         # lateral component in natural wind coordinates
    w = Uz*ct - Ux*st*ce - Uy*st*se           # vertical component in natural wind coordinates
    # do the variances
    uu = UxUx*ct**2*ce**2 + UyUy*ct**2*se**2 + UzUz*st**2 + 2*UxUy*ct**2*ce*se + 2*UxUz*ct*st*ce + 2*UyUz*ct*st*se
    vv = UyUy*ce**2 + UxUx*se**2 - 2*UxUy*ce*se
    ww = UzUz*ct**2 + UxUx*st**2*ce**2 + UyUy*st**2*se**2 - 2*UxUz*ct*st*ce - 2*UyUz*ct*st*se + 2*UxUy*st**2*ce*se
    # now do the covariances
    wT = UzT*ct - UxT*st*ce - UyT*st*se       # covariance(w,T) in natural wind coordinate system
    wA = UzA*ct - UxA*st*ce - UyA*st*se       # covariance(w,A) in natural wind coordinate system
    wC = UzC*ct - UxC*st*ce - UyC*st*se       # covariance(w,C) in natural wind coordinate system
    uw = UxUz*ce*(ct**2-st**2) - 2*UxUy*ct*st*ce*se + UyUz*se*(ct**2-st**2) - UxUx*ct*st*ce**2 - UyUy*ct*st*se**2 + UzUz*ct*st    # covariance(w,x) in natural wind coordinate system
    uv = UxUy*ct*(ce*ce-se*se) + UyUz*st*ce - UxUz*st*se - UxUx*ct*ce*se + UyUy*ct*ce*se    # covariance(x,y) in natural wind coordinate system
    vw = UyUz*ct*ce - UxUz*ct*se - UxUy*st*(ce**2-se**2) + UxUx*st*ce*se - UyUy*st*ce*se    # covariance(w,y) in natural wind coordinate system
    # now update the standard deviations
    u_Sd = numpy.ma.sqrt(uu)
    v_Sd = numpy.ma.sqrt(vv)
    w_Sd = numpy.ma.sqrt(ww)
    # store the rotated quantities in the nc object
    attr = qcutils.MakeAttributeDictionary(long_name='Horizontal rotation angle',units='deg')
    qcutils.CreateSeries(ds,'eta',eta,FList=['Ux','Uy','Uz'],Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='Vertical rotation angle',units='deg')
    qcutils.CreateSeries(ds,'theta',theta,FList=['Ux','Uy','Uz'],Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='Longitudinal component of wind-speed in natural wind coordinates',units='m/s')
    qcutils.CreateSeries(ds,'u',u,FList=['Ux','Uy','Uz'],Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='Lateral component of wind-speed in natural wind coordinates',units='m/s')
    qcutils.CreateSeries(ds,'v',v,FList=['Ux','Uy','Uz'],Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='Vertical component of wind-speed in natural wind coordinates',units='m/s')
    qcutils.CreateSeries(ds,'w',w,FList=['Ux','Uy','Uz'],Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='Kinematic heat flux, rotated to natural wind coordinates',units='mC/s')
    qcutils.CreateSeries(ds,'wT',wT,FList=['Ux','Uy','Uz','UxT','UyT','UzT'],Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='Kinematic vapour flux, rotated to natural wind coordinates',units='g/m2/s')
    qcutils.CreateSeries(ds,'wA',wA,FList=['Ux','Uy','Uz','UxA','UyA','UzA'],Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='Kinematic CO2 flux, rotated to natural wind coordinates',units='mg/m2/s')
    qcutils.CreateSeries(ds,'wC',wC,FList=['Ux','Uy','Uz','UxC','UyC','UzC'],Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='Momentum flux X component, rotated to natural wind coordinates',units='m2/s2')
    qcutils.CreateSeries(ds,'uw',uw,FList=['Ux','Uy','Uz','UxUz','UxUx','UxUy'],Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='Horizontal streamwise-crosswind covariance, rotated to natural wind coordinates',units='m2/s2')
    qcutils.CreateSeries(ds,'uv',uv,FList=['Ux','Uy','Uz','UxUz','UxUx','UxUy'],Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='Momentum flux Y component, rotated to natural wind coordinates',units='m2/s2')
    qcutils.CreateSeries(ds,'vw',vw,FList=['Ux','Uy','Uz','UyUz','UxUy','UyUy'],Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='Variance of streamwise windspeed, rotated to natural wind coordinates',units='m2/s2')
    qcutils.CreateSeries(ds,'uu',uu,FList=['Ux','Uy','Uz','UxUx','UxUy','UxUz'],Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='Variance of crossstream windspeed, rotated to natural wind coordinates',units='m2/s2')
    qcutils.CreateSeries(ds,'vv',vv,FList=['Ux','Uy','Uz','UyUy','UxUy'],Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='Variance of vertical windspeed, rotated to natural wind coordinates',units='m2/s2')
    qcutils.CreateSeries(ds,'ww',ww,FList=['Ux','Uy','Uz','UzUz','UxUz','UyUz'],Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='Longitudinal velocity component from CSAT, standard deviation, rotated to natural wind coordinates',units='m/s')
    qcutils.CreateSeries(ds,'u_Sd',u_Sd,FList=['Ux','Uy','Uz','UxUx','UxUy','UxUz'],Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='Lateral velocity component from CSAT, standard deviation, rotated to natural wind coordinates',units='m/s')
    qcutils.CreateSeries(ds,'v_Sd',v_Sd,FList=['Ux','Uy','Uz','UyUy','UxUy'],Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='Vertical velocity component from CSAT, standard deviation, rotated to natural wind coordinates',units='m/s')
    qcutils.CreateSeries(ds,'w_Sd',w_Sd,FList=['Ux','Uy','Uz','UzUz','UxUz','UyUz'],Attr=attr)
    if qcutils.cfkeycheck(cf,Base='General',ThisOne='RotateFlag') and cf['General']['RotateFlag'] == 'True':
        keys = ['eta','theta','u','v','w','wT','wA','wC','uw','vw','uu','vv','ww','u_Sd','v_Sd','w_Sd']
        for ThisOne in keys:
            testseries,f,a = qcutils.GetSeriesasMA(ds,ThisOne)
            mask = numpy.ma.getmask(testseries)
            index = numpy.where(mask.astype(numpy.int32)==1)
            ds.series[ThisOne]['Flag'][index] = numpy.int32(11)
    else:
        keys = ['eta','theta','u','v','w','wT','wA','wC','uw','vw','uu','vv','ww','u_Sd','v_Sd','w_Sd']
        for ThisOne in keys:
            testseries,f,a = qcutils.GetSeriesasMA(ds,ThisOne)
            mask = numpy.ma.getmask(testseries)
            index = numpy.where((numpy.mod(f,10)==0) & (mask.astype(numpy.int32)==1))    # find the elements with flag = 0, 10, 20 etc and masked (check for masked data with good data flag)
            ds.series[ThisOne]['Flag'][index] = numpy.int32(11)

def ConvertFc(cf,ds,Fco2_in='Fc'):
    """
        Converts CO2 flux (originally Fc) [mgCO2 m-2 s-1] to metabolic units [umol m-2 s-1]
        Converts CO2 flux (Fco2) to C flux (Fc) [mgC m-2 s-1]
        Calculates NEE [umol m-2 s-1], NEP [umol m-2 s-1], Fc_co2 [mgCO2 m-2 s-1], Fc_c [mgC m-2 s-1]
        NEE = -NEP
        """
    log.info(' Converting Fc units')
    if qcutils.cfkeycheck(cf,Base='FunctionArgs',ThisOne='convertFc'):
        Fc_list = ast.literal_eval(cf['FunctionArgs']['convertFc'])
        Fco2_in = Fc_list[0]
    Fc_co2,f,attr = qcutils.GetSeriesasMA(ds,Fco2_in)
    # conversion factors:  mol2umol: 1e6; mg2kg: 1e-6
    # conversions cancel:  mgCO2, kgCO2/mol, umol
    if 'umol' in attr['units']:
        NEE = Fc_co2
        Fc_co2 = NEE * c.Mco2
        NEP = -NEE
        Fc_c = NEE * c.Mc
    else:
        NEE = Fc_co2 / c.Mco2
        NEP = -NEE
        Fc_c = NEE * c.Mc
    
    attr_hist = 'Flux of carbon, Raw uncorrected'
    if 'rotated' in ds.series[Fco2_in]['Attr']['long_name']:
        attr_hist = attr_hist + ', rotated to natural wind coordinates'
    
    if 'frequency' in ds.series[Fco2_in]['Attr']['long_name']:
        attr_hist = attr_hist + ', frequency response corrected'
    
    if 'WPL' in ds.series[Fco2_in]['Attr']['long_name']:
        attr_hist = attr_hist + ', WPL'
    
    if 'SOLO' in ds.series[Fco2_in]['Attr']['long_name']:
        attr_hist = attr_hist + ', SOLO ANN gapfilled'
    
    for ThisOne in ['Fc_co2','Fc_c','NEE','NEP']:
        if ThisOne in ds.series.keys():
            ds.series[ThisOne]['Attr']['standard_name'] = 'not defined'
    
    if 'Fc_co2' in ds.series.keys():
        attr = qcutils.MakeAttributeDictionary(long_name='',units='mgCO2/m2/s')
        qcutils.CreateSeries(ds,'Fc_co2',Fc_co2,FList=[Fco2_in],Attr=attr)
        ds.series['Fc_co2']['Attr']['long_name'] = 'Fc_co2 (Flux of carbon dioxide), '+attr_hist
    else:
        attr = qcutils.MakeAttributeDictionary(long_name='Fc_co2 (Flux of carbon dioxide), '+attr_hist,units='mgCO2/m2/s')
        qcutils.CreateSeries(ds,'Fc_co2',Fc_co2,FList=[Fco2_in],Attr=attr)
    
    if 'Fc_c' in ds.series.keys():
        attr = qcutils.MakeAttributeDictionary(long_name='',units='mgC/m2/s')
        qcutils.CreateSeries(ds,'Fc_c',Fc_c,FList=[Fco2_in],Attr=attr)
        ds.series['Fc_c']['Attr']['long_name'] = 'Fc_c (Flux of carbon), '+attr_hist
    else:
        attr = qcutils.MakeAttributeDictionary(long_name='Fc_c (Flux of carbon), '+attr_hist,units='mgC/m2/s')
        qcutils.CreateSeries(ds,'Fc_c',Fc_c,FList=[Fco2_in],Attr=attr)
    
    if 'NEE' in ds.series.keys():
        attr = qcutils.MakeAttributeDictionary(long_name='',units='umol/m2/s',standard_name='surface_upward_mole_flux_of_carbon_dioxide')
        qcutils.CreateSeries(ds,'NEE',NEE,FList=[Fco2_in],Attr=attr)
        ds.series['NEE']['Attr']['long_name'] = 'NEE (Net ecosystem exchange of carbon), '+attr_hist
    else:
        attr = qcutils.MakeAttributeDictionary(long_name='NEE (Net ecosystem exchange of carbon), '+attr_hist,units='umol/m2/s',standard_name='surface_upward_mole_flux_of_carbon_dioxide')
        qcutils.CreateSeries(ds,'NEE',NEE,FList=[Fco2_in],Attr=attr)
    
    if 'NEP' in ds.series.keys():
        attr = qcutils.MakeAttributeDictionary(long_name='',units='umol/m2/s',standard_name='surface_downward_mole_flux_of_carbon_dioxide')
        qcutils.CreateSeries(ds,'NEP',NEP,FList=[Fco2_in],Attr=attr)
        ds.series['NEP']['Attr']['long_name'] = 'NEP (Net ecosystem photosynthesis), '+attr_hist
    else:
        attr = qcutils.MakeAttributeDictionary(long_name='NEP (Net ecosystem photosynthesis), '+attr_hist,units='umol/m2/s',standard_name='surface_downward_mole_flux_of_carbon_dioxide')
        qcutils.CreateSeries(ds,'NEP',NEP,FList=[Fco2_in],Attr=attr)
    
    if 'Fc' in ds.series.keys():
        attr = qcutils.MakeAttributeDictionary(long_name='',units='mgCO2/m2/s')
        qcutils.CreateSeries(ds,'Fc',Fc_co2,FList=[Fco2_in],Attr=attr)
        ds.series['Fc']['Attr']['long_name'] = 'Fc (Flux of carbon dioxide), '+attr_hist
    else:
        attr = qcutils.MakeAttributeDictionary(long_name='Fc (Flux of carbon dioxide), '+attr_hist,units='mgCO2/m2/s')
        qcutils.CreateSeries(ds,'Fc',Fc_co2,FList=[Fco2_in],Attr=attr)
    

def ConvertFcJason(cf,ds,Fco2_in='Fc'):
    """
        Converts CO2 flux (originally Fc) [mgCO2 m-2 s-1] to metabolic units [umol m-2 s-1]
        Converts CO2 flux (Fco2) to C flux (Fc) [mgC m-2 s-1]
        Calculates NEE [umol m-2 s-1], NEP [umol m-2 s-1], Fco2 [mgCO2 m-2 s-1], Fc [mgC m-2 s-1]
        NEE = -NEP
        """
    log.info(' Converting Fc units')
    if qcutils.cfkeycheck(cf,Base='FunctionArgs',ThisOne='JasonFc'):
        Fc_list = ast.literal_eval(cf['FunctionArgs']['JasonFc'])
        Fco2_in = Fc_list[0]
    Fco2,f,a = qcutils.GetSeriesasMA(ds,Fco2_in)
    # conversion factors:  mol2umol: 1e6; mg2kg: 1e-6
    # conversions cancel:  mgCO2, kgCO2/mol, umol
    NEE = Fco2 / c.Mco2
    NEP = -NEE
    Fc = NEE * c.Mc
    
    attr_hist = 'Flux of carbon, Raw uncorrected'
    if 'rotated' in ds.series[Fco2_in]['Attr']['long_name']:
        attr_hist = attr_hist + ', rotated to natural wind coordinates'
    
    if 'frequency' in ds.series[Fco2_in]['Attr']['long_name']:
        attr_hist = attr_hist + ', frequency response corrected'
    
    if 'WPL' in ds.series[Fco2_in]['Attr']['long_name']:
        attr_hist = attr_hist + ', WPL'
    
    if 'SOLO' in ds.series[Fco2_in]['Attr']['long_name']:
        attr_hist = attr_hist + ', SOLO ANN gapfilled'
    
    if 'Fc' in ds.series.keys():
        attr = qcutils.MakeAttributeDictionary(long_name='',units='umol/m2/s')
        qcutils.CreateSeries(ds,'Fc',NEE,FList=[Fco2_in],Attr=attr)
        ds.series['Fc']['Attr']['long_name'] = 'Fc (Flux of carbon), '+attr_hist
    else:
        attr = qcutils.MakeAttributeDictionary(long_name='Fc (Flux of carbon), '+attr_hist,units='umol/m2/s')
        qcutils.CreateSeries(ds,'Fc',NEE,FList=[Fco2_in],Attr=attr)

def CorrectIndividualFgForStorage(cf,ds,level='standard'):
    log.info(' Correcting soil heat flux for storage')
    zzz = 0
    if qcutils.cfkeycheck(cf,Base='FunctionArgs',ThisOne='CFgArgs'):
        List = cf['FunctionArgs']['CFgArgs'].keys()
        for i in range(len(List)):
            CFgArgs = ast.literal_eval(cf['FunctionArgs']['CFgArgs'][str(i)])
            if len(CFgArgs) == 4:
                CorrectFgForStorage(cf,ds,Fg_out=CFgArgs[0],Fg_in=CFgArgs[1],Ts_in=CFgArgs[2],Sws_in=CFgArgs[3])
            elif len(CFgArgs) == 5:
                CorrectFgForStorage(cf,ds,Fg_out=CFgArgs[0],Fg_in=CFgArgs[1],Ts_in=CFgArgs[2],Sws_in=CFgArgs[3],dTs_in=CFgArgs[4])
                zzz = 1
            else:
                log.error('  CorrectFg: bad parameters, check controlfile'); return
        if zzz == 1: qcts.AverageSeriesByElementsI(cf,ds,'dTs')
        return
    
    log.warn('  Correct Individual Fg for Storage:  no arguments provided, running default pre-storage average')
    qcts.AverageSeriesByElements(cf,ds3,'Fg')
    CorrectFgForStorage(cf,ds)

def CorrectGroupFgForStorage(cf,ds,level='standard'):
    if qcutils.cfkeycheck(cf,Base='FunctionArgs',ThisOne='CFgArgs'):
        CFgArgs = ast.literal_eval(cf['FunctionArgs']['CFgArgs'])
        if len(CFgArgs) == 1:
            qcts.AverageSeriesByElements(cf,ds,'dTs')
            log.info(' Correcting soil heat flux for storage')
            qcts.CorrectFgForStorage(cf,ds,dTs_in=CFgArgs[0])
        elif len(CFgArgs) == 4:
            log.info(' Correcting soil heat flux for storage')
            CorrectFgForStorage(cf,ds,Fg_out=CFgArgs[0],Fg_in=CFgArgs[1],Ts_in=CFgArgs[2],Sws_in=CFgArgs[3])
        elif len(CFgArgs) == 5:
            qcts.AverageSeriesByElements(cf,ds,'dTs')
            log.info(' Correcting soil heat flux for storage')
            CorrectFgForStorage(cf,ds,Fg_out=CFgArgs[0],Fg_in=CFgArgs[1],Ts_in=CFgArgs[2],Sws_in=CFgArgs[3],dTs_in=CFgArgs[4])
        return
    
    log.info(' Correcting soil heat flux for storage')
    CorrectFgForStorage(cf,ds)

def CorrectFgForStorage(cf,ds,Fg_out='Fg',Fg_in='Fg',Ts_in='Ts',Sws_in='Sws',dTs_in=''):
    """
        Correct ground heat flux for storage in the layer above the heat flux plate
        
        Usage qcts.CorrectFgForStorage(cf,ds,Fg_out,Fg_in,Ts_in,Sws_in)
        ds: data structure
        Fg_out: output soil heat flux variable to ds.  Example: 'Fg'
        Fg_in: input soil heat flux in ds.  Example: 'Fg_Av'
        Ts_in: input soil temperature in ds.  Example: 'Ts'
        
        Parameters loaded from control file:
            FgDepth: Depth of soil heat flux plates, m
            BulkDensity: soil bulk density, kg/m3
            OrganicContent: soil organic content, fraction
            SwsDefault: default value of soil moisture content used when no sensors present
        """
    d = max(0.0,min(0.5,numpy.float64(cf['Soil']['FgDepth'])))
    bd = max(1200.0,min(2500.0,numpy.float64(cf['Soil']['BulkDensity'])))
    oc = max(0.0,min(1.0,numpy.float64(cf['Soil']['OrganicContent'])))
    mc = 1.0 - oc
    Fg,Fg_flag,Fg_attr = qcutils.GetSeriesasMA(ds,Fg_in)  # raw soil heat flux
    nRecs = len(Fg)                               # number of records in series
    Ts,f,a = qcutils.GetSeriesasMA(ds,Ts_in)        # soil temperature
    Sws_default = min(1.0,max(0.0,numpy.float64(cf['Soil']['SwsDefault'])))
    if len(Sws_in) == 0:
        slist = []
        if qcutils.cfkeycheck(cf,Base='Soil',ThisOne='SwsSeries'):
            slist = ast.literal_eval(cf['Soil']['SwsSeries'])
        if len(slist)==0:
            log.info('  CorrectFgForStorage: Sws_default used for whole series')
            Sws = numpy.ones(nRecs,dtype=numpy.float64)*Sws_default
            Sws_flag = numpy.zeros(nRecs,dtype=numpy.int32)
        elif len(slist)==1:
            Sws,Sws_flag,a = qcutils.GetSeriesasMA(ds,slist[0])
        else:
            MergeSeries(cf,ds,'Sws',slist,[0,10])
            Sws,Sws_flag,a = qcutils.GetSeriesasMA(ds,'Sws')
    else:
        slist = Sws_in
        Sws,Sws_flag,Sws_attr = qcutils.GetSeriesasMA(ds,Sws_in)
    log.info('  CorrectForStorage: Fg_out is '+str(Fg_out)+' = f('+str(Fg_in)+', '+str(Ts_in)+', '+str(Sws_in)+', '+str(dTs_in)+')')
    iom = numpy.where(numpy.mod(Sws_flag,10)!=0)[0]
    #if len(iom)!=0:
    #    log.info('  CorrectFgForStorage: Sws_default used for '+str(len(iom))+' values')
    #    Sws[iom] = Sws_default
    #    Sws_flag[iom] = numpy.int32(0)
    if not dTs_in:
        dTs = numpy.ma.zeros(nRecs)
        dTs[1:] = numpy.ma.diff(Ts)
        # write the temperature difference into the data structure so we can use its flag later
        flag = numpy.zeros(nRecs,dtype=numpy.int32)
        index = numpy.ma.where(dTs.mask==True)[0]
        flag[index] = numpy.int32(1)
        attr = qcutils.MakeAttributeDictionary(long_name='Change in soil temperature',units='C')
        qcutils.CreateSeries(ds,"dTs",dTs,Flag=flag,Attr=attr)
    else:
        dTs,dTs_flag,dTs_attr = qcutils.GetSeriesasMA(ds,dTs_in)
    # get the time difference
    dt = numpy.ma.zeros(nRecs)
    dt[1:] = numpy.diff(date2num(ds.series['DateTime']['Data']))*numpy.float64(86400)
    dt[0] = dt[1]
    Cs = mc*bd*c.Cd + oc*bd*c.Co + Sws*c.rho_water*c.Cw
    S = Cs*(dTs/dt)*d
    Fg_o = Fg + S
    attr = qcutils.MakeAttributeDictionary(long_name='Soil heat flux corrected for storage',units='W/m2',standard_name='downward_heat_flux_in_soil')
    qcutils.CreateSeries(ds,Fg_out,Fg_o,FList=[Fg_in],Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='Soil heat flux storage',units='W/m2')
    qcutils.CreateSeries(ds,'S',S,FList=[Fg_in],Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='Specific heat capacity',units='J/m3/K')
    qcutils.CreateSeries(ds,'Cs',Cs,FList=[Fg_in],Attr=attr)

def CorrectSWC(cf,ds):
    """
        Correct soil moisture data using calibration curve developed from
        collected soil samples.  To avoid unrealistic or unphysical extreme
        values upon extrapolation, exponential and logarithmic using ln
        functions are applied to small and large values, respectively.
        Threshold values where one model replaces the other is determined where
        the functions cross.  The logarithmic curve is constrained at with a
        point at which the soil measurement = field porosity and the sensor
        measurement is maximised under saturation at field capacity.
        
        Usage qcts.CorrectSWC(cf,ds)
        cf: control file
        ds: data structure
        
        Parameters loaded from control file:
            SWCempList: list of raw CS616 variables
            SWCoutList: list of corrected CS616 variables
            SWCattr:  list of meta-data attributes for corrected CS616 variables
            SWC_a0: parameter in logarithmic model, actual = a1 * ln(sensor) + a0
            SWC_a1: parameter in logarithmic model, actual = a1 * ln(sensor) + a0
            SWC_b0: parameter in exponential model, actual = b0 * exp(b1 * sensor)
            SWC_b1: parameter in exponential model, actual = b0 * exp(b1 * sensor)
            SWC_t: threshold parameter for switching from exponential to logarithmic model
            TDRempList: list of raw CS610 variables
            TDRoutList: list of corrected CS610 variables
            TDRattr:  list of meta-data attributes for corrected CS610 variables
            TDRlinList: list of deep TDR probes requiring post-hoc linear correction to match empirical samples
            TDR_a0: parameter in logarithmic model, actual = a1 * ln(sensor) + a0
            TDR_a1: parameter in logarithmic model, actual = a1 * ln(sensor) + a0
            TDR_b0: parameter in exponential model, actual = b0 * exp(b1 * sensor)
            TDR_b1: parameter in exponential model, actual = b0 * exp(b1 * sensor)
            TDR_t: threshold parameter for switching from exponential to logarithmic model
        """
    SWCempList = ast.literal_eval(cf['Soil']['empSWCin'])
    SWCoutList = ast.literal_eval(cf['Soil']['empSWCout'])
    SWCattr = ast.literal_eval(cf['Soil']['SWCattr'])
    if cf['Soil']['TDR']=='Yes':
        TDRempList = ast.literal_eval(cf['Soil']['empTDRin'])
        TDRoutList = ast.literal_eval(cf['Soil']['empTDRout'])
        TDRlinList = ast.literal_eval(cf['Soil']['linTDRin'])
        TDRattr = ast.literal_eval(cf['Soil']['TDRattr'])
        TDR_a0 = numpy.float64(cf['Soil']['TDR_a0'])
        TDR_a1 = numpy.float64(cf['Soil']['TDR_a1'])
        TDR_b0 = numpy.float64(cf['Soil']['TDR_b0'])
        TDR_b1 = numpy.float64(cf['Soil']['TDR_b1'])
        TDR_t = numpy.float64(cf['Soil']['TDR_t'])
    SWC_a0 = numpy.float64(cf['Soil']['SWC_a0'])
    SWC_a1 = numpy.float64(cf['Soil']['SWC_a1'])
    SWC_b0 = numpy.float64(cf['Soil']['SWC_b0'])
    SWC_b1 = numpy.float64(cf['Soil']['SWC_b1'])
    SWC_t = numpy.float64(cf['Soil']['SWC_t'])
    
    for i in range(len(SWCempList)):
        log.info('  Applying empirical correction to '+SWCempList[i])
        invar = SWCempList[i]
        outvar = SWCoutList[i]
        attr = SWCattr[i]
        Sws,f,a = qcutils.GetSeriesasMA(ds,invar)
        
        nRecs = len(Sws)
        
        Sws_out = numpy.ma.empty(nRecs,numpy.float64)
        Sws_out.fill(c.missing_value)
        Sws_out.mask = numpy.ma.empty(nRecs,bool)
        Sws_out.mask.fill(True)
        
        index_high = numpy.ma.where((Sws.mask == False) & (Sws > SWC_t))[0]
        index_low = numpy.ma.where((Sws.mask == False) & (Sws < SWC_t))[0]
        
        Sws_out[index_low] = SWC_b0 * numpy.exp(SWC_b1 * Sws[index_low])
        Sws_out[index_high] = (SWC_a1 * numpy.log(Sws[index_high])) + SWC_a0
        
        attr = qcutils.MakeAttributeDictionary(long_name=attr,units='m3/m3',standard_name='soil_moisture_content')
        qcutils.CreateSeries(ds,outvar,Sws_out,FList=[invar],Attr=attr)
    if cf['Soil']['TDR']=='Yes':
        for i in range(len(TDRempList)):
            log.info('  Applying empirical correction to '+TDRempList[i])
            invar = TDRempList[i]
            outvar = TDRoutList[i]
            attr = TDRattr[i]
            Sws,f,a = qcutils.GetSeriesasMA(ds,invar)
            
            nRecs = len(Sws)
            
            Sws_out = numpy.ma.empty(nRecs,numpy.float64)
            Sws_out.fill(c.missing_value)
            Sws_out.mask = numpy.ma.empty(nRecs,bool)
            Sws_out.mask.fill(True)
            
            index_high = numpy.ma.where((Sws.mask == False) & (Sws > TDR_t))[0]
            index_low = numpy.ma.where((Sws.mask == False) & (Sws < TDR_t))[0]
            
            Sws_out[index_low] = TDR_b0 * numpy.exp(TDR_b1 * Sws[index_low])
            Sws_out[index_high] = (TDR_a1 * numpy.log(Sws[index_high])) + TDR_a0
            
            attr = qcutils.MakeAttributeDictionary(long_name=attr,units='m3/m3',standard_name='soil_moisture_content')
            qcutils.CreateSeries(ds,outvar,Sws_out,FList=[invar],Attr=attr)

def CorrectWindDirection(cf,ds,Wd_in):
    """
        Correct wind direction for mis-aligned sensor direction.
        
        Usage qcts.CorrectWindDirection(cf,ds,Wd_in)
        cf: control file
        ds: data structure
        Wd_in: input/output wind direction variable in ds.  Example: 'Wd_CSAT'
        """
    log.info(' Correcting wind direction')
    Wd,f,a = qcutils.GetSeriesasMA(ds,Wd_in)
    ldt = ds.series['DateTime']['Data']
    KeyList = cf['Variables'][Wd_in]['Correction'].keys()
    for i in range(len(KeyList)):
        ItemList = ast.literal_eval(cf['Variables'][Wd_in]['Correction'][str(i)])
        try:
            si = ldt.index(datetime.datetime.strptime(ItemList[0],'%Y-%m-%d %H:%M'))
        except ValueError:
            si = 0
        try:
            ei = ldt.index(datetime.datetime.strptime(ItemList[1],'%Y-%m-%d %H:%M')) + 1
        except ValueError:
            ei = -1
        Correction = numpy.float64(ItemList[2])
        Wd[si:ei] = Wd[si:ei] + Correction
    Wd = numpy.mod(Wd,numpy.float64(360))
    ds.series[Wd_in]['Data'] = numpy.ma.filled(Wd,numpy.float64(c.missing_value))

def do_attributes(cf,ds):
    """
        Import attriubes in xl2nc control file to netCDF dataset.  Included
        global and variable attributes.  Also attach flag definitions to global
        meta-data for reference.
        
        Usage qcts.do_attributes(cf,ds)
        cf: control file
        ds: data structure
        """
    log.info(' Getting the attributes given in control file')
    if 'Global' in cf.keys():
        for gattr in cf['Global'].keys():
            ds.globalattributes[gattr] = cf['Global'][gattr]
        ds.globalattributes['Flag000'] = 'Good data'
        ds.globalattributes['Flag001'] = 'QA/QC: -9999 in level 1 dataset'
        ds.globalattributes['Flag002'] = 'QA/QC: L2 Range Check'
        ds.globalattributes['Flag003'] = 'QA/QC: L2 Diurnal SD Check'
        ds.globalattributes['Flag004'] = 'QA/QC: CSAT Diagnostic'
        ds.globalattributes['Flag005'] = 'QA/QC: LI7500 Diagnostic'
        ds.globalattributes['Flag006'] = 'QA/QC: Excluded Dates'
        ds.globalattributes['Flag007'] = 'QA/QC: Excluded Hours'
        ds.globalattributes['Flag008'] = 'albedo: bad Fsd < threshold (290 W/m2 default) only if bad time flag not set'
        ds.globalattributes['Flag009'] = 'albedo: bad time flag (not midday 10.00 to 14.00)'
        ds.globalattributes['Flag010'] = 'Corrections: Apply Linear'
        ds.globalattributes['Flag011'] = 'Corrections/Combinations: Coordinate Rotation (Ux, Uy, Uz, UxT, UyT, UzT, UxA, UyA, UzA, UxC, UyC, UzC, UxUz, UxUx, UxUy, UyUz, UxUy, UyUy)'
        ds.globalattributes['Flag012'] = 'Corrections/Combinations: Massman Frequency Attenuation Correction (Coord Rotation, Tv_CSAT, Ah_HMP, ps)'
        ds.globalattributes['Flag013'] = 'Corrections/Combinations: Virtual to Actual Fh (Coord Rotation, Massman, Ta_HMP)'
        ds.globalattributes['Flag014'] = 'Corrections/Combinations: WPL correction for flux effects on density measurements (Coord Rotation, Massman, Fhv to Fh, Cc_7500_Av)'
        ds.globalattributes['Flag015'] = 'Corrections/Combinations: Ta from Tv'
        ds.globalattributes['Flag016'] = 'Corrections/Combinations: L3 Range Check'
        ds.globalattributes['Flag017'] = 'Corrections/Combinations: L3 Diurnal SD Check'
        ds.globalattributes['Flag018'] = 'Corrections/Combinations: u* filter'
        ds.globalattributes['Flag019'] = 'Corrections/Combinations: Gap coordination'
        ds.globalattributes['Flag021'] = 'Corrections/Combinations: Data-flag mismatch'
        ds.globalattributes['Flag022'] = 'Corrections/Combinations: L=0 in zeta=z/L'
        ds.globalattributes['Flag030'] = 'GapFilling: Gap Filled by ANN (SOLO)'
        ds.globalattributes['Flag031'] = 'GapFilling: Gap not filled'
        ds.globalattributes['Flag040'] = 'GapFilling: Gap Filled from one-minute average'
        ds.globalattributes['Flag050'] = 'GapFilling: Gap Filled from alternate'
        ds.globalattributes['Flag060'] = 'GapFilling: Gap Filled by interpolation'
        ds.globalattributes['Flag070'] = 'GapFilling: Gap Filled by replacement from paired tower'
        ds.globalattributes['Flag080'] = 'GapFilling: u* from Fh'
        ds.globalattributes['Flag081'] = 'GapFilling: u* not from Fh'
        ds.globalattributes['Flag082'] = 'GapFilling: L4 Range Check'
        ds.globalattributes['Flag083'] = 'GapFilling: L4 Diurnal SD Check'
        ds.globalattributes['Flag084'] = 'GapFilling: L5 Range Check'
        ds.globalattributes['Flag085'] = 'GapFilling: L5 Diurnal SD Check'
        ds.globalattributes['Flag086'] = 'GapFilling: L6 Range Check'
        ds.globalattributes['Flag087'] = 'GapFilling: L6 Diurnal SD Check'
        ds.globalattributes['Flag090'] = 'Partitioning Night: ER computed from exponential temperature response curves'
        ds.globalattributes['Flag100'] = 'Partitioning Day: GPP/ER computed from light-response curves, GPP = ER - Fc'
        ds.globalattributes['Flag110'] = 'Partitioning Day: GPP night mask'
        ds.globalattributes['Flag120'] = 'Partitioning Day: Fc > CE, GPP = 0, CE = Fc'
        ds.globalattributes['Flag130'] = 'Partitioning Day: ER_dark from sink period (+NEP) light response curves'
        ds.globalattributes['Flag140'] = 'Partitioning Day: ER_dark from source  period (+NEE) light response curves'
        ds.globalattributes['Flag150'] = 'Partitioning Day: PD, CE_day & GPP from conditional correlation'
        ds.globalattributes['Flag151'] = 'Partitioning Day: no solution from conditional correlation'
        ds.globalattributes['Flag161'] = 'Footprint: Date filter'
        ds.globalattributes['Flag162'] = 'Footprint: no solution'
        ds.globalattributes['Flag171'] = 'Penman-Monteith: bad rav or rSm only if bad Uavg, bad Fe and bad Fsd flags not set'
        ds.globalattributes['Flag172'] = 'Penman-Monteith: bad Fe < threshold (0 W/m2 default) only if bad Fsd flag not set'
        ds.globalattributes['Flag173'] = 'Penman-Monteith: bad Fsd < threshold (10 W/m2 default)'
        ds.globalattributes['Flag174'] = 'Penman-Monteith: Uavg == 0 (undefined aerodynamic resistance under calm conditions) only if bad Fe and bad Fsd flags not set'
        ds.globalattributes['Flag180'] = 'Penman-Monteith 2-layer: rav_base short-circuit'
        ds.globalattributes['Flag190'] = 'Penman-Monteith 2-layer: rav_top short-circuit'
        ds.globalattributes['Flag191'] = 'Penman-Monteith 2-layer: rav_top not short-circuit (rav_base undefined)'
        ds.globalattributes['Flag200'] = 'Penman-Monteith 2-layer: parallel circuit'
        ds.globalattributes['Flag201'] = 'Penman-Monteith 2-layer: not parallel circuit (rav_full short-circuit)'
        ds.globalattributes['Flag211'] = 'Bulk Richardson number flags: delta_U=0 (RiB infinite)'
    for ThisOne in ds.series.keys():
        if ThisOne in cf['Variables']:
            if 'Attr' in cf['Variables'][ThisOne].keys():
                ds.series[ThisOne]['Attr'] = {}
                for attr in cf['Variables'][ThisOne]['Attr'].keys():
                    ds.series[ThisOne]['Attr'][attr] = cf['Variables'][ThisOne]['Attr'][attr]

def do_bulkRichardson(cf,ds):
    log.info('  Computing bulk Richardson values from qT profile')
    IncludeList = cf['RiB']['series'].keys()
    NumHeights = len(IncludeList)
    if 'sigmaw' not in ds.series.keys():
        varw,f,a = qcutils.GetSeriesasMA(ds,'ww')
        windex = numpy.where(varw > 0)[0]
        n = len(varw)
        sigmaw = numpy.zeros(n,dtype=numpy.float64)
        sigmaw[windex] = numpy.sqrt(varw[windex])
        attr = qcutils.MakeAttributeDictionary(long_name='standard deviation of vertical windspeed',units='m^.5/s^.5')
        qcutils.CreateSeries(ds,'sigmaw',sigmaw,FList=['ww'],Attr=attr)
    
    # calculate layer values
    for i in range(NumHeights-1):
        IncludeHeightList = ast.literal_eval(cf['RiB']['series'][str(i)])
        NextHeightList = ast.literal_eval(cf['RiB']['series'][str(i+1)])
        U_top,f,a = qcutils.GetSeriesasMA(ds,NextHeightList[1])
        U_bottom,f,a = qcutils.GetSeriesasMA(ds,IncludeHeightList[1])
        Tvp_top,f,a = qcutils.GetSeriesasMA(ds,NextHeightList[2])
        Tvp_bottom,f,a = qcutils.GetSeriesasMA(ds,IncludeHeightList[2])
        delta_z = numpy.float64(NextHeightList[0]) - numpy.float64(IncludeHeightList[0])
        delta_Tvp = Tvp_top - Tvp_bottom
        delta_U = U_top - U_bottom
        Uindex = numpy.where(delta_U == 0)[0]
        badUindex = numpy.where(delta_U != 0)[0]
        AverageSeriesByElements(cf,ds,IncludeHeightList[3])
        Tvp_bar,f,a = qcutils.GetSeriesasMA(ds,IncludeHeightList[3])
        RiB = numpy.zeros(len(delta_U),dtype=numpy.float64)
        RiB[badUindex] = (9.8 * delta_Tvp[badUindex] * delta_z) / (Tvp_bar[badUindex] * (delta_U[badUindex] ** 2))
        RiB[Uindex] = numpy.float64(c.missing_value)
        RiB_flag = numpy.zeros(len(RiB),dtype=numpy.float64)
        attr = qcutils.MakeAttributeDictionary(long_name='U_upper less U_lower',units='m/s')
        qcutils.CreateSeries(ds,IncludeHeightList[4],delta_U,FList=[IncludeHeightList[1],IncludeHeightList[2],IncludeHeightList[3],'Fc'],Attr=attr)
        attr = qcutils.MakeAttributeDictionary(long_name='Tvp_upper less Tvp_lower',units='m/s')
        qcutils.CreateSeries(ds,IncludeHeightList[5],delta_Tvp,FList=[IncludeHeightList[1],IncludeHeightList[2],IncludeHeightList[3],'Fc'],Attr=attr)
        attr = qcutils.MakeAttributeDictionary(long_name='Bulk Richardson number',units='none')
        qcutils.CreateSeries(ds,IncludeHeightList[6],RiB,FList=[IncludeHeightList[1],IncludeHeightList[2],IncludeHeightList[3],'Fc'],Attr=attr)
        ds.series[IncludeHeightList[6]]['Flag'][Uindex] = numpy.int32(211)
        flaglist = [IncludeHeightList[4],IncludeHeightList[5],IncludeHeightList[6]]
        for ThisOne in flaglist:
            flag = numpy.where(numpy.mod(ds.series[ThisOne]['Flag'],10)!=0)[0]
            ds.series[ThisOne]['Data'][flag] = numpy.float64(c.missing_value)

def do_climatology(cf,ds):
    if qcutils.cfkeycheck(cf,Base='Climatology',ThisOne='met'):
        MList = ast.literal_eval(cf['Climatology']['met'])
    else:
        MList = ['Ta','Ah','Cc_7500_Av','Ws_CSAT','Wd_CSAT','ps']
    
    if qcutils.cfkeycheck(cf,Base='Climatology',ThisOne='rad'):
        RList = ast.literal_eval(cf['Climatology']['rad'])
    else:
        RList = ['Fld','Flu','Fn','Fsd','Fsu']
    
    if qcutils.cfkeycheck(cf,Base='Climatology',ThisOne='soil'):
        SList = ast.literal_eval(cf['Climatology']['soil'])
    else:
        SList = ['Ts','Sws','Fg']
    
    if qcutils.cfkeycheck(cf,Base='Climatology',ThisOne='flux'):
        FList = ast.literal_eval(cf['Climatology']['flux'])
    else:
        FList = ['Fc','Fe','Fh','Fm','ustar']
    
    OList = MList+RList+SList+FList                                   # output series
    qcts.ComputeClimatology(cf,ds,OList)

def do_footprint_2d(cf,ds,datalevel='L3',footprintlevel='V3'):
    if 'Footprint' not in cf.keys(): log.error('   no Footprint section in cf'); return
    if not qcutils.cfkeycheck(cf,Base='Footprint',ThisOne='zm'): log.error('   no measurement height zm in cf'); return
    if not qcutils.cfkeycheck(cf,Base='Footprint',ThisOne='zc'): log.error('   no canopy height zc in cf'); return
    log.info(' Calculating 2D footprint')
    n = len(ds.series['xlDateTime']['Data'])
    if cf['Footprint']['zc'] == 'variable':
        zc,f,a = qcutils.GetSeriesasMA(ds,'zc')
    else:
        zc = numpy.zeros(n,dtype=numpy.float64) + numpy.float64(cf['Footprint']['zc'])
    
    zm = numpy.zeros(n,dtype=numpy.float64) + numpy.float64(cf['Footprint']['zm'])
    if cf['Footprint']['zm']<1.0: log.error('   zm needs to be larger than 1 m'); return
    # roughness length, z0 = zc / 3e (Brutseart 1982, eqn 5.7), m
    denominator = 1 / (3 * numpy.exp(1))
    znot = zc * denominator
    d = 2 / 3 * zc
    zmd = zm - d
    if qcutils.cfkeycheck(cf,Base='Footprint',ThisOne='r'):
        r = numpy.int32(cf['Footprint']['r'])
    else:
        r = numpy.int32(96)
    
    if r>96:
        log.error('   r needs to be smaller than 96')
        return
    
    if qcutils.cfkeycheck(cf,Base='Footprint',ThisOne='ExcludeHours') and cf['Footprint']['ExcludeHours'] == 'True':
        qcck.do_excludehours(cf,ds,'L')
    
    varw,f,a = qcutils.GetSeriesasMA(ds,'ww')
    varv,f,a = qcutils.GetSeriesasMA(ds,'vv')
    ustar,f,a = qcutils.GetSeriesasMA(ds,'ustar')
    Ta,f,a = qcutils.GetSeriesasMA(ds,'Ta')
    Ah,f,a = qcutils.GetSeriesasMA(ds,'Ah')
    ps,f,a = qcutils.GetSeriesasMA(ds,'ps')
    Fh,f,a = qcutils.GetSeriesasMA(ds,'Fh')
    if 'L' not in ds.series.keys():
        L_calc = mf.molen(Ta, Ah, ps, ustar, Fh, fluxtype='sensible')
        attr = qcutils.MakeAttributeDictionary(long_name='Obukhov Length: Raw uncorrected',units='m')
        qcutils.CreateSeries(ds,'L',L_calc,FList=['Ta','Ah','ps','ustar','Fh'],Attr=attr)
    L,Lf,a = qcutils.GetSeriesasMA(ds,'L')
    Fsd,f,a = qcutils.GetSeriesasMA(ds,'Fsd')
    n = len(L)
    windex = numpy.where(varw > 0)[0]
    badw = numpy.where((varw < c.eps) & (numpy.mod(Lf,10)==0))[0]
    vindex = numpy.where(varv > 0)[0]
    badv = numpy.where((varv < c.eps) & (numpy.mod(Lf,10)==0))[0]
    includedate = len(L)
    excludehours = 0
    sigmaw = numpy.zeros(n,dtype=numpy.float64)
    sigmav = numpy.zeros(n,dtype=numpy.float64)
    sigmaw[windex] = numpy.sqrt(varw[windex])
    sigmav[vindex] = numpy.sqrt(varv[vindex])
    if 'zeta' not in ds.series.keys():
        zetaindex = numpy.where(L != 0)[0]
        zetaflagindex = numpy.where((L == 0)&(numpy.mod(Lf,10)==0))[0]
        zeta = numpy.ones(len(L),dtype=numpy.float64)*c.missing_value
        zeta[zetaindex] = zmd[zetaindex] / L[zetaindex]
        zf = numpy.zeros(n,dtype=numpy.int32) + Lf
        zf[zetaflagindex] = numpy.int32(22) 
        attr = qcutils.MakeAttributeDictionary(long_name='stability coefficient, (z-d)/L',units='unitless')
        qcutils.CreateSeries(ds,'zeta',zeta,Flag=zf,Attr=attr)
    else:
        zeta,zf,zeta_attr = qcutils.GetSeriesasMA(ds,'zeta')
    #Lmask = numpy.where(zf == 22)[0]
    #Lf[Lmask] = numpy.int(22)
    Lf[badw] = 9999
    Lf[badv] = 9999
    attr = qcutils.MakeAttributeDictionary(long_name='standard deviation of vertical windspeed',units='m^.5/s^.5')
    qcutils.CreateSeries(ds,'sigmaw',sigmaw,FList=['ww'],Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='standard deviation of crosswind speed (after rotation to natural wind coordinates)',units='m^.5/s^.5')
    qcutils.CreateSeries(ds,'sigmav',sigmav,FList=['vv'],Attr=attr)
    if qcutils.cfkeycheck(cf,Base='Footprint',ThisOne='AnalysisDates'):
        ldt = ds.series['DateTime']['Data']
        IncludeList = cf['Footprint']['AnalysisDates'].keys()
        NumDates = len(IncludeList)
        analysisflag = numpy.zeros(n,dtype=numpy.int32) + 1
        for i in range(NumDates):
            IncludeDateList = ast.literal_eval(cf['Footprint']['AnalysisDates'][str(i)])
            try:
                si = ldt.index(datetime.datetime.strptime(IncludeDateList[0],'%Y-%m-%d %H:%M'))
            except ValueError:
                si = 0
            
            try:
                ei = ldt.index(datetime.datetime.strptime(IncludeDateList[1],'%Y-%m-%d %H:%M')) + 1
            except ValueError:
                ei = -1
            
            analysisflag[si:ei] = numpy.int32(0)
        
        index = numpy.where(analysisflag == 1)[0]
        Lf[index] = 9999
        includedate = len(numpy.where(analysisflag == 0)[0])
    
    if qcutils.cfkeycheck(cf,Base='Footprint',ThisOne='Fluxes'):
        Fluxesin = ast.literal_eval(cf['Footprint']['Fluxes'])
        if len(Fluxesin) == 3:
            Fc,fFc,a = qcutils.GetSeriesasMA(ds,Fluxesin[0])
            Fe,fFe,a = qcutils.GetSeriesasMA(ds,Fluxesin[1])
            Fh,fFh,a = qcutils.GetSeriesasMA(ds,Fluxesin[2])
            footprintlevel = 'V3'
        elif len(Fluxesin) == 2:
            Fc,fFc,a = qcutils.GetSeriesasMA(ds,Fluxesin[0])
            Fe,fFe,a = qcutils.GetSeriesasMA(ds,Fluxesin[1])
            Fh = fFh = numpy.zeros(len(fFe),dtype=numpy.int32)
            footprintlevel = 'V2'
        else:
            log.error('  Footprint:  only 2 or 3 flux variables accepted.  Correct controlfile')
            return
    else:
        Fc,fFc,a = qcutils.GetSeriesasMA(ds,'Fc')
        Fe,fFe,a = qcutils.GetSeriesasMA(ds,'Fe')
        Fh,fFh,a = qcutils.GetSeriesasMA(ds,'Fh')
        footprintlevel = 'V3'
    
    if qcutils.cfkeycheck(cf,Base='Footprint',ThisOne='ExcludeDay') and cf['Footprint']['ExcludeDay'] == 'True':
        Dayindex = numpy.where((Fsd > 10) & (numpy.mod(Lf,10)==0))[0]
        excludehours = len(Dayindex)
        Lf[Dayindex] = 7
    
    if qcutils.cfkeycheck(cf,Base='Footprint',ThisOne='ExcludeNight') and cf['Footprint']['ExcludeNight'] == 'True':
        Nightindex = numpy.where((Fsd < 10) & (numpy.mod(Lf,10)==0))[0]
        excludehours = len(Nightindex)
        Lf[Nightindex] = 7
    
    ustarindex = numpy.where((ustar < 0.2) & (numpy.mod(Lf,10)==0))[0]
    Lf[ustarindex] = 18
    index_noFlux = numpy.ma.where((numpy.mod(Lf,10)==0) & ((numpy.mod(fFc,10)!=0) | (numpy.mod(fFe,10)!=0) | (numpy.mod(fFh,10)!=0)))[0]
    Lf[index_noFlux] = 31
    if 'zi' in ds.series.keys():
        h,zif,a = qcutils.GetSeriesasMA(ds,'zi')
    else:
        h = numpy.ma.zeros(n,dtype=numpy.float64)
        h[index_very_unstable] = numpy.float64(2000)
        h[index_unstable] = numpy.float64(1500)
        h[index_slight_unstable] = numpy.float64(1200)
        h[index_neutral] = numpy.float64(1000)
        h[index_slight_stable] = numpy.float64(800)
        h[index_stable] = numpy.float64(250)
        h[index_very_stable] = numpy.float64(200)
    
    ziindex = numpy.where((numpy.mod(zif,10)!=0) & (numpy.mod(Lf,10)==0))[0]
    Lf[ziindex] = 99999999
    zetaxindex = numpy.where(numpy.mod(Lf,10)==0)[0]
    Lfindex = numpy.where(numpy.mod(Lf,10)!=0)[0]
    zeta[Lfindex] = 9999999
    index_zero_L = numpy.ma.where(zeta > 9999998)[0]
    index_neutral = numpy.ma.where((zeta > -0.1) & (zeta < 0.1))[0]
    index_slight_stable = numpy.ma.where((zeta > 0.1) & (zeta < 1))[0]
    index_stable = numpy.ma.where((zeta > 1) & (zeta < 2))[0]
    index_very_stable = numpy.ma.where((zeta > 2) & (zeta < 9999999))[0]
    index_slight_unstable = numpy.ma.where((zeta < -0.1) & (zeta > -1))[0]
    index_unstable = numpy.ma.where((zeta < -1) & (zeta > -2))[0]
    index_very_unstable = numpy.ma.where(zeta < -2)[0]
    h[index_zero_L] = numpy.float64(0)
    log.info('  Total measurements:  '+str(includedate))
    log.info('  Masked measurements (not included):  '+str(includedate-len(zetaxindex)))
    log.info('   '+str(excludehours)+':  excluded times')
    log.info('   '+str(includedate-len(zetaxindex)-len(ustarindex)-len(index_noFlux)-excludehours)+':  turbulence gaps (sigma-w, L, u*)')
    log.info('   '+str(len(index_noFlux))+':  flux gaps')
    log.info('   '+str(len(ustarindex))+':  u* < 0.2 m/s')
    log.info('   '+str(len(ziindex))+':  zi gaps')
    log.info('  Included measurements:  '+str(len(zetaxindex)))
    log.info('   '+str(len(index_neutral))+':  neutral, -0.1 < zeta < 0.1')
    log.info('   '+str(len(index_slight_stable))+':  slightly stable, 0.1 < zeta < 1')
    log.info('   '+str(len(index_stable))+':  stable, 1 < zeta < 2')
    log.info('   '+str(len(index_very_stable))+':  very stable, zeta > 2')
    log.info('   '+str(len(index_slight_unstable))+':  slightly unstable, -0.1 > zeta > -1')
    log.info('   '+str(len(index_unstable))+':  unstable, -1 > zeta > -2')
    log.info('   '+str(len(index_very_unstable))+':  very unstable, zeta < -2')
    if (qcutils.cfkeycheck(cf,Base='Output',ThisOne='FootprintSummaryFile') and cf['Output']['FootprintSummaryFile'] == 'True') or (qcutils.cfkeycheck(cf,Base='Output',ThisOne='FootprintDataFile') and cf['Output']['FootprintDataFile'] == 'True'):
        filename = cf['Files']['Footprint']['FootprintFilePath']+'stats_'+datalevel+'.txt'
        out = open(filename,'w')
        out.write('Total measurements:  '+str(includedate)+'\n')
        out.write('Masked measurements (not included):  '+str(includedate-len(zetaxindex))+'\n')
        out.write(' '+str(excludehours)+':  excluded times\n')
        out.write(' '+str(includedate-len(zetaxindex)-len(ustarindex)-len(index_noFlux)-excludehours)+':  turbulence gaps (sigma-w, L, u*)\n')
        out.write(' '+str(len(index_noFlux))+':  flux gaps\n')
        out.write(' '+str(len(ustarindex))+':  u* < 0.2 m/s\n')
        out.write(' '+str(len(ziindex))+':  zi gaps\n')
        out.write('Included measurements:  '+str(len(zetaxindex))+'\n')
        out.write(' '+str(len(index_neutral))+':  neutral, -0.1 < zeta < 0.1\n')
        out.write(' '+str(len(index_slight_stable))+':  slightly stable, 0.1 < zeta < 1\n')
        out.write(' '+str(len(index_stable))+':  stable, 1 < zeta < 2\n')
        out.write(' '+str(len(index_very_stable))+':  very stable, zeta > 2\n')
        out.write(' '+str(len(index_slight_unstable))+':  slightly unstable, -0.1 > zeta > -1\n')
        out.write(' '+str(len(index_unstable))+':  unstable, -1 > zeta > -2\n')
        out.write(' '+str(len(index_very_unstable))+':  very unstable, zeta < -2\n')
        out.close()
    
    if len(zetaxindex) == 0:
        log.warn('   Footprint:  no observations passed filtering')
        return
    
    if qcutils.cfkeycheck(cf,Base='Footprint',ThisOne='wd'):
        wdin = cf['Footprint']['wd']
    else:
        wdin = 'Wd_CSAT'
    
    wd,f,a = qcutils.GetSeriesasMA(ds,wdin)
    eta,f,a = qcutils.GetSeriesasMA(ds,'eta')
    xr = numpy.ma.zeros(n,dtype=numpy.float64)
    if qcutils.cfkeycheck(cf,Base='Output',ThisOne='FootprintDataFileType') and cf['Output']['FootprintDataFileType'] == 'Climatology':
        if qcutils.cfkeycheck(cf,Base='Footprint',ThisOne='ClimateXmin') and qcutils.cfkeycheck(cf,Base='Footprint',ThisOne='ClimateXmax') and qcutils.cfkeycheck(cf,Base='Footprint',ThisOne='ClimateYmin') and qcutils.cfkeycheck(cf,Base='Footprint',ThisOne='ClimateYmax') and qcutils.cfkeycheck(cf,Base='Footprint',ThisOne='ClimatePixel'):
            n = ((numpy.float64(cf['Footprint']['ClimateXmax'])) - (numpy.float64(cf['Footprint']['ClimateXmin']))) / (numpy.float64(cf['Footprint']['ClimatePixel']))
            m = ((numpy.float64(cf['Footprint']['ClimateYmax'])) - (numpy.float64(cf['Footprint']['ClimateYmin']))) / (numpy.float64(cf['Footprint']['ClimatePixel']))
            xmin = (numpy.float64(cf['Footprint']['ClimateXmin']))
            xmax = (numpy.float64(cf['Footprint']['ClimateXmax']))
            ymin = (numpy.float64(cf['Footprint']['ClimateYmin']))
            ymax = (numpy.float64(cf['Footprint']['ClimateYmax']))
            p = (numpy.float64(cf['Footprint']['ClimatePixel']))
            fc_c1_2d = numpy.ma.zeros((m,n),dtype=numpy.float64)
            fc_e_2d = numpy.ma.zeros((m,n),dtype=numpy.float64)
            fc_h_2d = numpy.ma.zeros((m,n),dtype=numpy.float64)
            fc_c2_2d = numpy.ma.zeros((m,n),dtype=numpy.float64)
            fc_c3_2d = numpy.ma.zeros((m,n),dtype=numpy.float64)
        else:
            log.error('  Footprint climatology: Spatial parameters missing from controlfile')
            return
    
    hcheck = h - zm
    do_index = numpy.where((sigmaw > 0) & (sigmav > 0) & (ustar > 0.2) & (h > 1) & (hcheck > 0) & (numpy.mod(Lf,10)==0))[0]
    for i in range(len(do_index)):
        if qcutils.cfkeycheck(cf,Base='Footprint',ThisOne='loglist') and cf['Footprint']['loglist'] == 'True':
            log.info('    Footprint: '+str(ds.series['DateTime']['Data'][do_index[i]]))
        
        if qcutils.cfkeycheck(cf,Base='Output',ThisOne='FootprintDataFileType') and cf['Output']['FootprintDataFileType'] == 'Climatology' and qcutils.cfkeycheck(cf,Base='Output',ThisOne='FootprintDataFile') and cf['Output']['FootprintDataFile'] == 'True':
            if footprintlevel == 'V3':
                xr[do_index[i]],x_2d,y_2d,fw_c1_2d,fw_e_2d,fw_h_2d,fw_c2_2d,fw_c3_2d,filenames,labels = footprint_2d(cf,sigmaw[do_index[i]],sigmav[do_index[i]],ustar[do_index[i]],zm[do_index[i]],h[do_index[i]],znot[do_index[i]],r,wd[do_index[i]],zeta[do_index[i]],L[do_index[i]],zc[do_index[i]],ds.series['DateTime']['Data'][do_index[i]],eta[do_index[i]],Fc[do_index[i]],Fe[do_index[i]],Fh[do_index[i]],footprintlevel)
                fc_c1_2d,fc_e_2d,fc_h_2d,fc_c2_2d,fc_c3_2d,fc_x2d,fc_y2d = footprint_climatology(fc_c1_2d,fc_e_2d,x_2d,y_2d,fw_c1_2d,fw_e_2d,n,m,xmin,xmax,ymin,ymax,p,fc_h=fc_h_2d,fc_c2=fc_c2_2d,fc_c3=fc_c3_2d,fw_h=fw_h_2d,fw_c2=fw_c2_2d,fw_c3=fw_c3_2d)
            
            if footprintlevel == 'V2':
                xr[do_index[i]],x_2d,y_2d,fw_c1_2d,fw_e_2d,filenames,labels = footprint_2d(cf,sigmaw[do_index[i]],sigmav[do_index[i]],ustar[do_index[i]],zm[do_index[i]],h[do_index[i]],znot[do_index[i]],r,wd[do_index[i]],zeta[do_index[i]],L[do_index[i]],zc[do_index[i]],ds.series['DateTime']['Data'][do_index[i]],eta[do_index[i]],Fc[do_index[i]],Fe[do_index[i]],Fh[do_index[i]],footprintlevel)
                fc_c1_2d,fc_e_2d,fc_x2d,fc_y2d = footprint_climatology(fc_c1_2d,fc_e_2d,x_2d,y_2d,fw_c1_2d,fw_e_2d,n,m,xmin,xmax,ymin,ymax,p)
        else:
            xr[do_index[i]] = footprint_2d(cf,sigmaw[do_index[i]],sigmav[do_index[i]],ustar[do_index[i]],zm[do_index[i]],h[do_index[i]],znot[do_index[i]],r,wd[do_index[i]],zeta[do_index[i]],L[do_index[i]],zc[do_index[i]],ds.series['DateTime']['Data'][do_index[i]],eta[do_index[i]],Fc[do_index[i]],Fe[do_index[i]],Fh[do_index[i]],footprintlevel)
        
    if qcutils.cfkeycheck(cf,Base='Output',ThisOne='FootprintDataFileType') and cf['Output']['FootprintDataFileType'] == 'Climatology' and qcutils.cfkeycheck(cf,Base='Output',ThisOne='FootprintDataFile') and cf['Output']['FootprintDataFile'] == 'True':
        if ((not qcutils.cfkeycheck(cf, Base='Output', ThisOne='FootprintDataType')) or ((qcutils.cfkeycheck(cf, Base='Output', ThisOne='FootprintDataType')) and ((cf['Output']['FootprintDataType'] != 'Vector') and (cf['Output']['FootprintDataType'] != 'Matrix') and (cf['Output']['FootprintDataType'] != 'Both')))):
            log.error('  Footprint climatology:  FootprintFileType (Vector or Matrix) not defined in controlfile')
        
        if ((qcutils.cfkeycheck(cf, Base='Output', ThisOne='FootprintDataType')) and ((cf['Output']['FootprintDataType'] == 'Vector') or (cf['Output']['FootprintDataType'] == 'Both'))):
            for i in range(len(labels)):
                if labels[i] == 'fc_GPP' or labels[i] == 'fc_Fc':
                    footprint_vector_out(filenames[0][i],fc_x2d,fc_y2d,fc_c1_2d,labels[i])
                if labels[i] == 'fc_ER' or labels[i] == 'fc_Fe':
                    footprint_vector_out(filenames[0][i],fc_x2d,fc_y2d,fc_e_2d,labels[i])
                if labels[i] == 'fc_Fh':
                    footprint_vector_out(filenames[0][i],fc_x2d,fc_y2d,fc_h_2d,labels[i])
                if labels[i] == 'fc_PosFc':
                    footprint_vector_out(filenames[0][i],fc_x2d,fc_y2d,fc_c2_2d,labels[i])
                if labels[i] == 'fc_NegFc':
                    footprint_vector_out(filenames[0][i],fc_x2d,fc_y2d,fc_c3_2d,labels[i])
        
        if ((qcutils.cfkeycheck(cf, Base='Output', ThisOne='FootprintDataType')) and ((cf['Output']['FootprintDataType'] == 'Matrix') or (cf['Output']['FootprintDataType'] == 'Both'))):
            for i in range(len(labels)):
                if labels[i] == 'fc_GPP' or labels[i] == 'fc_Fc':
                    footprint_matrix_out(filenames[1][i],fc_x2d,fc_y2d,fc_c1_2d)
                if labels[i] == 'fc_CE' or labels[i] == 'fc_Fe':
                    footprint_matrix_out(filenames[1][i],fc_x2d,fc_y2d,fc_e_2d)
                if labels[i] == 'fc_Fh':
                    footprint_matrix_out(filenames[1][i],fc_x2d,fc_y2d,fc_h_2d)
                if labels[i] == 'fc_PosFc':
                    footprint_matrix_out(filenames[1][i],fc_x2d,fc_y2d,fc_c2_2d)
                if labels[i] == 'fc_NegFc':
                    footprint_matrix_out(filenames[1][i],fc_x2d,fc_y2d,fc_c3_2d)
        
    attr = qcutils.MakeAttributeDictionary(long_name='integrated footprint in the direction of the wind',units='m')
    qcutils.CreateSeries(ds,'xr',xr,FList=['L','ww','vv','ustar'],Attr=attr)
    flag_index = numpy.ma.where((xr == 0) & (numpy.mod(ds.series['xr']['Flag'],10)==0))[0]
    ustar_index = numpy.ma.where(ustar < 0.2)[0]
    date_index = numpy.ma.where(Lf == 9999)[0]
    ds.series['xr']['Flag'][flag_index] = numpy.int32(162)
    ds.series['xr']['Flag'][ustar_index] = numpy.int32(18)
    ds.series['xr']['Flag'][date_index] = numpy.int32(201)
    index = numpy.where((numpy.mod(ds.series['xr']['Flag'],10)!=0))[0]    # find the elements with flag != 0, 10, 20 etc
    ds.series['xr']['Data'][index] = numpy.float64(c.missing_value)

def do_functions(cf,ds):
    log.info(' Getting variances from standard deviations & vice versa')
    if 'AhAh' in ds.series.keys() and 'Ah_7500_Sd' not in ds.series.keys():
        AhAh,flag,attr = qcutils.GetSeriesasMA(ds,'AhAh')
        Ah_7500_Sd = numpy.ma.sqrt(AhAh)
        attr = qcutils.MakeAttributeDictionary(long_name='Absolute humidity from Li-7500, standard deviation',units='g/m3')
        qcutils.CreateSeries(ds,'Ah_7500_Sd',Ah_7500_Sd,Flag=flag,Attr=attr)
    if 'Ah_7500_Sd' in ds.series.keys() and 'AhAh' not in ds.series.keys():
        Ah_7500_Sd,flag,attr = qcutils.GetSeriesasMA(ds,'Ah_7500_Sd')
        AhAh = Ah_7500_Sd*Ah_7500_Sd
        attr = qcutils.MakeAttributeDictionary(long_name='Absolute humidity from Li-7500, variance',units='(g/m3)2')
        qcutils.CreateSeries(ds,'AhAh',AhAh,Flag=flag,Attr=attr)
    if 'CcCc' in ds.series.keys() and 'Cc_7500_Sd' not in ds.series.keys():
        CcCc,flag,attr = qcutils.GetSeriesasMA(ds,'CcCc')
        Cc_7500_Sd = numpy.ma.sqrt(CcCc)
        attr = qcutils.MakeAttributeDictionary(long_name='CO2 concentration from Li-7500, standard deviation',units='mg/m3')
        qcutils.CreateSeries(ds,'Cc_7500_Sd',Cc_7500_Sd,Flag=flag,Attr=attr)
    if 'Cc_7500_Sd' in ds.series.keys() and 'CcCc' not in ds.series.keys():
        Cc_7500_Sd,flag,attr = qcutils.GetSeriesasMA(ds,'Cc_7500_Sd')
        CcCc = Cc_7500_Sd*Cc_7500_Sd
        attr = qcutils.MakeAttributeDictionary(long_name='CO2 concentration from Li-7500, variance',units='(mg/m3)2')
        qcutils.CreateSeries(ds,'CcCc',CcCc,Flag=flag,Attr=attr)
    if 'Ux_Sd' in ds.series.keys() and 'UxUx' not in ds.series.keys():
        Ux_Sd,flag,attr = qcutils.GetSeriesasMA(ds,'Ux_Sd')
        UxUx = Ux_Sd*Ux_Sd
        attr = qcutils.MakeAttributeDictionary(long_name='Longitudinal velocity component from CSAT, variance',units='(m/s)2')
        qcutils.CreateSeries(ds,'UxUx',UxUx,Flag=flag,Attr=attr)
    if 'UxUx' in ds.series.keys() and 'Ux_Sd' not in ds.series.keys():
        UxUx,flag,attr = qcutils.GetSeriesasMA(ds,'UxUx')
        Ux_Sd = numpy.ma.sqrt(UxUx)
        attr = qcutils.MakeAttributeDictionary(long_name='Longitudinal velocity component from CSAT, standard deviation',units='m/s')
        qcutils.CreateSeries(ds,'Ux_Sd',Ux_Sd,Flag=flag,Attr=attr)
    if 'Uy_Sd' in ds.series.keys() and 'UyUy' not in ds.series.keys():
        Uy_Sd,flag,attr = qcutils.GetSeriesasMA(ds,'Uy_Sd')
        UyUy = Uy_Sd*Uy_Sd
        attr = qcutils.MakeAttributeDictionary(long_name='Lateral velocity component from CSAT, variance',units='(m/s)2')
        qcutils.CreateSeries(ds,'UyUy',UyUy,Flag=flag,Attr=attr)
    if 'UyUy' in ds.series.keys() and 'Uy_Sd' not in ds.series.keys():
        UyUy,flag,attr = qcutils.GetSeriesasMA(ds,'UyUy')
        Uy_Sd = numpy.ma.sqrt(UyUy)
        attr = qcutils.MakeAttributeDictionary(long_name='Lateral velocity component from CSAT, standard deviation',units='m/s')
        qcutils.CreateSeries(ds,'Uy_Sd',Uy_Sd,Flag=flag,Attr=attr)
    if 'Uz_Sd' in ds.series.keys() and 'UzUz' not in ds.series.keys():
        Uz_Sd,flag,attr = qcutils.GetSeriesasMA(ds,'Uz_Sd')
        UzUz = Uz_Sd*Uz_Sd
        attr = qcutils.MakeAttributeDictionary(long_name='Vertical velocity component from CSAT, variance',units='(m/s)2')
        qcutils.CreateSeries(ds,'UzUz',UzUz,Flag=flag,Attr=attr)
    if 'UzUz' in ds.series.keys() and 'Uz_Sd' not in ds.series.keys():
        UzUz,flag,attr = qcutils.GetSeriesasMA(ds,'UzUz')
        Uz_Sd = numpy.ma.sqrt(UzUz)
        attr = qcutils.MakeAttributeDictionary(long_name='Vertical velocity component from CSAT, standard deviation',units='m/s')
        qcutils.CreateSeries(ds,'Uz_Sd',Uz_Sd,Flag=flag,Attr=attr)
    if 'Tv_Sd' in ds.series.keys() and 'TvTv' not in ds.series.keys():
        Tv_Sd,flag,attr = qcutils.GetSeriesasMA(ds,'Tv_Sd')
        TvTv = Tv_Sd*Tv_Sd
        attr = qcutils.MakeAttributeDictionary(long_name='Virtual air temperature from CSAT, variance',units='C2')
        qcutils.CreateSeries(ds,'TvTv',TvTv,Flag=flag,Attr=attr)
    if 'TvTv' in ds.series.keys() and 'Tv_Sd' not in ds.series.keys():
        TvTv,flag,attr = qcutils.GetSeriesasMA(ds,'TvTv')
        Tv_Sd = numpy.ma.sqrt(TvTv)
        attr = qcutils.MakeAttributeDictionary(long_name='Virtual air temperature from CSAT, standard deviation',units='C')
        qcutils.CreateSeries(ds,'Tv_Sd',Tv_Sd,Flag=flag,Attr=attr)

def do_PenmanMonteith(cf,ds):
    if qcutils.cfkeycheck(cf,Base='PenmanMonteith',ThisOne='Cdmethod'):
        Cdmethod = cf['PenmanMonteith']['Cdmethod']
    else:
        Cdmethod = 'False'
    
    if qcutils.cfkeycheck(cf,Base='PenmanMonteith',ThisOne='Cemethod'):
        Cemethod = cf['PenmanMonteith']['Cemethod']
    else:
        Cemethod = 'False'
    
    if qcutils.cfkeycheck(cf,Base='PenmanMonteith',ThisOne='Ce_2layer'):
        Ce_2layer = cf['PenmanMonteith']['Ce_2layer']
    else:
        Ce_2layer = 'False'
    
    if Cdmethod != 'True' and Cemethod != 'True' and Ce_2layer != 'True':
        log.error(' PenmanMontieth:  no method selected')
        return
        
    prep_aerodynamicresistance(cf,ds,Cdmethod,Cemethod,Ce_2layer)
    return

def do_sums(cf,ds):
    # compute daily statistics
    if qcutils.cfkeycheck(cf,Base='Sums',ThisOne='SumList'):
        SumList = ast.literal_eval(cf['Sums']['SumList'])
    else:
        SumList = ['Precip','ET','Energy','Radiation','Carbon']
    
    if qcutils.cfkeycheck(cf,Base='Sums',ThisOne='SubSumList'):
        SubSumList = ast.literal_eval(cf['Sums']['SubSumList'])
    else:
        SubSumList = []
    
    if qcutils.cfkeycheck(cf,Base='Sums',ThisOne='MinMaxList'):
        MinMaxList = ast.literal_eval(cf['Sums']['MinMaxList'])
    else:
        MinMaxList = ['Ta','Vbat','Tpanel','Carbon']
    
    if qcutils.cfkeycheck(cf,Base='Sums',ThisOne='MeanList'):
        MeanList = ast.literal_eval(cf['Sums']['MeanList'])
    else:
        MeanList = ['Ta','Tpanel']
    
    if qcutils.cfkeycheck(cf,Base='Sums',ThisOne='SoilList'):
        SoilList = ast.literal_eval(cf['Sums']['SoilList'])
    else:
        SoilList = []
    
    StatsList = SumList + MinMaxList + MeanList + SoilList
    if len(StatsList) > 0:
        qcts.ComputeDailySums(cf,ds,SumList,SubSumList,MinMaxList,MeanList,SoilList)

def do_WPL(cf,ds,cov=''):
    if cov == 'True':
        Fe_WPLcov(cf,ds)
        Fc_WPLcov(cf,ds)
        return
    
    Fe_WPL(cf,ds)
    Fc_WPL(cf,ds)

def extrapolate_humidity(zq_ref,zq_high,zq_low,q_high,q_low,fqh,fql):
    rise = zq_high - zq_low
    run = q_high - q_low
    qref = numpy.zeros(len(run),dtype=numpy.float64) + numpy.float64(c.missing_value)
    slope = numpy.zeros(len(run),dtype=numpy.float64) + numpy.float64(c.missing_value)
    index = numpy.where((numpy.mod(fqh,10)==0) & (numpy.mod(fql,10)==0))[0]
    slope[index] = rise / run[index]
    qref[index] = q_high[index] - ((zq_high - zq_ref) / slope[index])
    return qref

def Fc_WPLcov(cf,ds,Fc_wpl_out='Fc',wC_out='wC',wC_in='wC',Fh_in='Fh',wA_in='wA',Ta_in='Ta',Ah_in='Ah',Cc_in='Cc_7500_Av',ps_in='ps'):
    """
        Apply Webb, Pearman and Leuning correction to carbon flux using the
        original formulation (WPL80).  This correction is necessary to account
        for flux effects on density measurements.  This method uses the
        originally-published formulation using covariances rather than fluxes.
        The difference in the corrected fluxes using the two routines is minor
        and related to scaling the met variables.
        
        Usage qcts.Fc_WPLcov(ds,Fc_wpl_out,wC,Fh,wA,Ta,Ah,Cc,ps)
        ds: data structure
        Fc_wpl_out: output corrected carbon flux to ds.  Example: 'Fc_wpl'
        wC: input covariance(wC) in ds.  Example: 'wCM'
        Fh: input sensible heat flux in ds.  Example: 'Fh_rmv'
        wA: input covariance(wA) in ds.  Example: 'wAwpl'
        Ta: input air temperature in ds.  Example: 'Ta_HMP'
        Ah: input absolute humidity in ds.  Example: 'Ah_HMP'
        Cc: input co2 density in ds.  Example: 'Cc_7500_Av'
        ps: input atmospheric pressure in ds.  Example: 'ps'
        
        Pre-requisite: FhvtoFh
        Pre-requisite: Fe_WPLcov
        
        Accepts meteorological constants or variables
        """
    log.info(' Applying WPL correction to Fc')
    if qcutils.cfkeycheck(cf,Base='FunctionArgs',ThisOne='CWPL'):
        Cargs = ast.literal_eval(cf['FunctionArgs']['CWPL'])
        Fc_wpl_out = Cargs[0]
        wC_in = Cargs[1]
        Fh_in = Cargs[2]
        wA_in = Cargs[3]
        Ta_in = Cargs[4]
        Ah_in = Cargs[5]
        Cc_in = Cargs[6]
        ps_in = Cargs[7]
    wC,f,a = qcutils.GetSeriesasMA(ds,wC_in)
    Fh,f,a = qcutils.GetSeriesasMA(ds,Fh_in)
    wA,f,a = qcutils.GetSeriesasMA(ds,wA_in)
    Ta,f,a = qcutils.GetSeriesasMA(ds,Ta_in)
    Ah,f,a = qcutils.GetSeriesasMA(ds,Ah_in)
    Cc,f,a = qcutils.GetSeriesasMA(ds,Cc_in)
    ps,f,a = qcutils.GetSeriesasMA(ds,ps_in)
    Cpm,f,a = qcutils.GetSeriesasMA(ds,'Cpm')
    rhom,f,a = qcutils.GetSeriesasMA(ds,'rhom')
    rhod,f,a = qcutils.GetSeriesasMA(ds,'rhod')
    nRecs = numpy.size(wC)
    TaK = Ta + 273.15
    Ah = Ah/numpy.float64(1000)                       # Absolute humidity from g/m3 to kg/m3
    Ccg = Cc/numpy.float64(1000)                  # CO2 from mg/m3 to g/m3
    sigma_wpl = Ah/rhod
    wT = Fh / (rhom * Cpm)
    # unit check:  TERM1 wC: mg/m2/s; TERM2 Ccg: g/m3; rhod: kg/m3; wA: g/m2/s; TERM3 mu: unitless; sigma: unitless; Cc: mg/m3; Ta: K; wT: mC/s
    # unit check:  [mg/m2/s] + ([g/m3] / [kg/m3]) * [g/m2/s] + ([mg/m3] * [mC/s]) / K; [mg/m2/s] + [g2/kg/m2/s] + [mg/m2/s]
    # unit check:  [g2/kg/m2/s] * [kg/1000g]; [g/m2/s] * [1000mg/g]; [mg/m2/s]
    Fc_wpl_data = wC + (c.mu * (Ccg / rhod) * wA) + ((1 + (c.mu * sigma_wpl)) * (Cc / TaK) * wT)
    attr = qcutils.MakeAttributeDictionary(long_name='WPL',units='mg/m2/s')
    qcutils.CreateSeries(ds,Fc_wpl_out,Fc_wpl_data,FList=[wC_in,Fh_in,wA_in,Ta_in,Ah_in,Cc_in,ps_in],Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='WPL',units='mg/m2/s')
    qcutils.CreateSeries(ds,wC_out,Fc_wpl_data,FList=[wC_in,Fh_in,wA_in,Ta_in,Ah_in,Cc_in,ps_in],Attr=attr)
    if qcutils.cfkeycheck(cf,Base='General',ThisOne='WPLFlag') and cf['General']['WPLFlag'] == 'True':
        keys = [Fc_wpl_out,wC_out]
        for ThisOne in keys:
            testseries,f,a = qcutils.GetSeriesasMA(ds,ThisOne)
            mask = numpy.ma.getmask(testseries)
            index = numpy.where(mask.astype(numpy.int32)==1)
            ds.series[ThisOne]['Flag'][index] = numpy.int32(14)
    else:
        keys = [Fc_wpl_out,wC_out]
        for ThisOne in keys:
            testseries,f,a = qcutils.GetSeriesasMA(ds,ThisOne)
            mask = numpy.ma.getmask(testseries)
            index = numpy.where((numpy.mod(f,10)==0) & (mask.astype(numpy.int32)==1))    # find the elements with flag = 0, 10, 20 etc and masked (check for masked data with good data flag)
            ds.series[ThisOne]['Flag'][index] = numpy.int32(14)

def Fe_WPLcov(cf,ds,Fe_wpl_out='Fe',wA_in='wA',Fh_in='Fh',Ta_in='Ta',Ah_in='Ah',ps_in='ps',wA_out='wA'):
    """
        Apply Webb, Pearman and Leuning correction to vapour flux using the
        original formulation (WPL80).  This correction is necessary to account
        for flux effects on density measurements.  This method uses the
        originally-published formulation using covariances rather than fluxes.
        The difference in the corrected fluxes using the two routines is minor
        and related to scaling the met variables.
        
        Usage qcts.Fe_WPLcov(ds,Fe_wpl_out,wA,Fh,Ta,Ah,ps)
        ds: data structure
        Fe_wpl_out: output corrected water vapour flux to ds.  Example: 'Fe_wpl'
        wA: input covariance(wA) in ds.  Example: 'wAM'
        Fh: input sensible heat flux in ds.  Example: 'Fh_rmv'
        Ta: input air temperature in ds.  Example: 'Ta_HMP'
        Ah: input absolute humidity in ds.  Example: 'Ah_HMP'
        ps: input atmospheric pressure in ds.  Example: 'ps'
        
        Pre-requisite: FhvtoFh
        Pre-requisite: Fe_WPLcov
        
        Accepts meteorological constants or variables
        """
    log.info(' Applying WPL correction to Fe')
    if qcutils.cfkeycheck(cf,Base='FunctionArgs',ThisOne='EWPL'):
        Eargs = ast.literal_eval(cf['FunctionArgs']['EWPL'])
        Fe_wpl_out = Eargs[0]
        wA_in = Eargs[1]
        Fh_in = Eargs[2]
        Ta_in = Eargs[3]
        Ah_in = Eargs[4]
        ps_in = Eargs[5]
        wA_out = Eargs[6]
    wA,f,a = qcutils.GetSeriesasMA(ds,wA_in)
    Fh,f,a = qcutils.GetSeriesasMA(ds,Fh_in)
    Ta,f,a = qcutils.GetSeriesasMA(ds,Ta_in)
    Ah,f,a = qcutils.GetSeriesasMA(ds,Ah_in)
    ps,f,a = qcutils.GetSeriesasMA(ds,ps_in)
    rhom,f,a = qcutils.GetSeriesasMA(ds,'rhom')
    rhod,f,a = qcutils.GetSeriesasMA(ds,'rhod')
    Lv,f,a = qcutils.GetSeriesasMA(ds,'Lv')
    Cpm,f,a = qcutils.GetSeriesasMA(ds,'Cpm')
    nRecs = numpy.size(wA)
    TaK = Ta + 273.15
    Ah = Ah * c.g2kg                       # Absolute humidity from g/m3 to kg/m3
    sigma_wpl = Ah / rhod
    wT = Fh / (rhom * Cpm)
    # unit check:  TERM1 Lv J/kg; TERM2 c.mu: unitless; sigma: unitless; TERM3 wA: g/m2/s; Ah: kg/m3; Ta: K; wT: mC/s
    # unit check:  J/kg * ([g/m2/s * kg/1000g] + [kg/m3 * mC/K/s]); J/kg * ([kg/m2/s] + [kg/m2/s]); W/m2
    Fe_wpl_data = Lv * (1 + (c.mu * sigma_wpl)) * ((wA * c.g2kg) + ((Ah / TaK) * wT))
    # unit check:  Fe [J/m2/s] * Lv-1 [kg/J] * g2kg-1; g/m2/s
    wAwpl = Fe_wpl_data / (Lv * c.g2kg)
    attr = qcutils.MakeAttributeDictionary(long_name='WPL',units='W/m2',standard_name='surface_upward_latent_heat_flux')
    qcutils.CreateSeries(ds,Fe_wpl_out,Fe_wpl_data,FList=[wA_in,Fh_in,Ta_in,Ah_in,ps_in],Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='WPL',units='g/(m2 s)')
    qcutils.CreateSeries(ds,wA_out,wAwpl,FList=[wA_in,Fh_in,Ta_in,Ah_in,ps_in],Attr=attr)
    if qcutils.cfkeycheck(cf,Base='General',ThisOne='WPLFlag') and cf['General']['WPLFlag'] == 'True':
        keys = [Fe_wpl_out,wA_out]
        for ThisOne in keys:
            testseries,f,a = qcutils.GetSeriesasMA(ds,ThisOne)
            mask = numpy.ma.getmask(testseries)
            index = numpy.where(mask.astype(numpy.int32)==1)
            ds.series[ThisOne]['Flag'][index] = numpy.int32(14)
    else:
        keys = [Fe_wpl_out,wA_out]
        for ThisOne in keys:
            testseries,f,a = qcutils.GetSeriesasMA(ds,ThisOne)
            mask = numpy.ma.getmask(testseries)
            index = numpy.where((numpy.mod(f,10)==0) & (mask.astype(numpy.int32)==1))    # find the elements with flag = 0, 10, 20 etc and masked (check for masked data with good data flag)
            ds.series[ThisOne]['Flag'][index] = numpy.int32(14)

def FhvtoFh(cf,ds,Ta_in='Ta',Fh_in='Fh',Tv_in='Tv_CSAT',Fe_in='Fe',ps_in='ps',Ah_in='Ah',Fh_out='Fh',wT_out='wT'):
    """
        Corrects sensible heat flux calculated on the covariance between w'
        and theta', deviations in vertical windspeed and sonically-derived
        virtual temperature.  Uses the formulation developed by Ed Swiatek,
        Campbell Scientific and located in the open path eddy covariance manual.
        
        Usage qcts.FhvtoFh(ds,Ta_in,Fh_in,Tv_in,Fe_in,ps_in,Ah_in,Fh_out,attr)
        ds: data structure
        Ta_in: input air temperature in ds.  Example: 'Ta_HMP'
        Fh_in: input sensible heat flux in ds.  Example: 'Fh'
        Tv_in: input sonic virtual temperature in ds.  Example: 'Tv_CSAT'
        Fe_in: input water vapour flux in ds.  Example: 'Fe_raw'
        ps_in: input atmospheric pressure in ds.  Example: 'ps'
        Ah_in: input absolute pressure in ds.  Example: 'Ah_EC'
        Fh_out: output corrected sensible heat flux to ds.  Example: 'Fh_rv'
        attr: attribute field for variable meta-data in ds.  Example: 'Fh rotated and converted from virtual heat flux'
        
        Typically used following:
            CoordRotation, MassmanApprox, Massman, CalculateFluxesRM (recommended)
            or
            CoordRotation, CalculateFluxes
            or
            CalculateFluxes_Unrotated
        
        Accepts meteorological constants or variables
        """
    log.info(' Converting virtual Fh to Fh')
    if qcutils.cfkeycheck(cf,Base='FunctionArgs',ThisOne='FhvtoFhArgs'):
        args = ast.literal_eval(cf['FunctionArgs']['FhvtoFhArgs'])
        Ta_in = args[0]
        Fh_in = args[1]
        Tv_in = args[2]
        Fe_in = args[3]
        ps_in = args[4]
        Ah_in = args[5]
        Fh_out = args[6]
        wT_out = args[7]
    Ta,f,a = qcutils.GetSeriesasMA(ds,Ta_in)
    Fh,f,a = qcutils.GetSeriesasMA(ds,Fh_in)
    Tv,f,a = qcutils.GetSeriesasMA(ds,Tv_in)
    Fe,f,a = qcutils.GetSeriesasMA(ds,Fe_in)
    ps,f,a = qcutils.GetSeriesasMA(ds,ps_in)
    Ah,f,a = qcutils.GetSeriesasMA(ds,Ah_in)
    rhom,f,a = qcutils.GetSeriesasMA(ds,'rhom')
    Lv,f,a = qcutils.GetSeriesasMA(ds,'Lv')
    Cpm,f,a = qcutils.GetSeriesasMA(ds,'Cpm')
    nRecs = len(Fh)
    psPa = ps * 1000
    TaK = Ta + c.C2K
    TvK = Tv + c.C2K
    Fh_o = (TaK / TvK) * (Fh - (rhom * Cpm * ((0.51 * c.Rd * (TaK ** 2)) / psPa) * (Fe / Lv)))
    wT = Fh_o / (rhom * Cpm)
    
    attr = qcutils.MakeAttributeDictionary(long_name='converted from virtual heat flux',units='W/m2',standard_name='surface_upward_sensible_heat_flux')
    qcutils.CreateSeries(ds,Fh_out,Fh_o,FList=[Ta_in, Fh_in, Tv_in, Fe_in, ps_in, Ah_in],Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='converted from virtual heat flux',units='mC/s')
    qcutils.CreateSeries(ds,wT_out,wT,FList=[Ta_in, Fh_in, Tv_in, Fe_in, ps_in, Ah_in],Attr=attr)
    if qcutils.cfkeycheck(cf,Base='General',ThisOne='FhvtoFhFlag') and cf['General']['FhvtoFhFlag'] == 'True':
        for ThisOne in [Fh_out,wT_out]:
            testseries,f,a = qcutils.GetSeriesasMA(ds,ThisOne)
            mask = numpy.ma.getmask(testseries)
            index = numpy.where(mask.astype(numpy.int32)==1)
            ds.series[ThisOne]['Flag'][index] = numpy.int32(13)
            ds.series[ThisOne]['Flag'][index] = numpy.int32(13)
    else:
        keys = [Fh_out,wT_out]
        for ThisOne in keys:
            testseries,f,a = qcutils.GetSeriesasMA(ds,ThisOne)
            mask = numpy.ma.getmask(testseries)
            index = numpy.where((numpy.mod(f,10)==0) & (mask.astype(numpy.int32)==1))    # find the elements with flag = 0, 10, 20 etc and masked (check for masked data with good data flag)
            ds.series[ThisOne]['Flag'][index] = numpy.int32(13)

def FilterUstar(cf,ds,ustar_in='ustar',ustar_out='ustar_filtered'):
    """
        Filter ustar for low turbulence periods.  The filtering is done by checking the
        friction velocity for each time period.  If ustar is less than or equal to the
        threshold specified in the control file then ustar is set to missing.  If
        the ustar is greater than the threshold, no action is taken.  Filtering is not
        done "in place", a new series is created with the label given in the control file.
        The QC flag is set to 18 to indicate the missing low ustar values.
        
        Usage: qcts.FilterUstar(cf,ds)
        cf: control file object
        ds: data structure object
        """
    if ustar_out not in cf['Variables'].keys(): return
    if 'ustar_threshold' in cf['Variables'][ustar_out].keys():
        log.info(' Filtering ustar to remove values below threshold')
        ustar_threshold = numpy.float64(cf['Variables'][ustar_out]['ustar_threshold'])
        ustar,ustar_flag,ustar_attr = qcutils.GetSeriesasMA(ds,ustar_in)
        index = numpy.ma.where(ustar<=ustar_threshold)[0]
        ustar = numpy.ma.masked_where(ustar<=ustar_threshold,ustar)
        ustar_flag[index] = numpy.int32(18)
        descr = 'ustar filtered for low turbulence conditions (<'+str(ustar_threshold)+')'
        units = qcutils.GetUnitsFromds(ds, ustar_in)
        attr = qcutils.MakeAttributeDictionary(long_name=descr,units=units)
        qcutils.CreateSeries(ds,ustar_out,ustar,Flag=ustar_flag,Attr=attr)
    else:
        log.error(' ustar threshold (ustar_threshold) not found in '+ustar_out+' section of control file')

def footprint_2d(cf,sigmaw,sigmav,ustar,zm,h,znot,r,wd,zeta,L,zc,timestamp,eta,Fc,Fe,Fh,footprintlevel):
    """
        footprint_2d.py
        
        Derive a footprint estimate based on a simple parameterisation
        
        Details see Kljun, N., Calanca, P., Rotach, M.W., Schmid, H.P., 2004:
        Boundary-Layer Meteorology 112, 503-532.
         
        online version: http://footprint.kljun.net
        contact: n.kljun@swansea.ac.uk
        
        Usage: footprint_2d.py <measurement height[m]> <roughness length [m]> <planetary boundary
                              layer height [m]> <sigma_w [m s-1]> <sigma_v [m s-1]> <u* [m s-1]> <r [%]> 
         Output: crosswind integrated footprint: x,f
                 extent of footprint up to r%: xr
                 matrix of 2d footprint: x_2d, y_2d, f_2d
        
         Created: May 15 2012, natascha kljun
         Last change: May 15 2012, nk
        """
    
    STList = []
    for fmt in ['%Y','%m','%d','%H','%M']:
        STList.append('_')
        STList.append(timestamp.strftime(fmt))
    
    if qcutils.cfkeycheck(cf, Base='Output', ThisOne='FootprintSummaryFile') and cf['Output']['FootprintSummaryFile'] == 'True':
        summaryFileName = cf['Files']['Footprint']['FootprintFilePath']+'footprint_2d_summary'+''.join(STList)+'.xls'
    
    if qcutils.cfkeycheck(cf,Base='Output',ThisOne='FootprintDataFile') and (cf['Output']['FootprintDataFile'] == 'True'):
        vectorFileName = cf['Files']['Footprint']['FootprintFilePath']+'footprint_2d_vectors'+''.join(STList)+'.xls'
        matrixFileName = cf['Files']['Footprint']['FootprintFilePath']+'footprint_2d_matrix'+''.join(STList)+'.xls'
    
    # -----------------------------------------------
    # Initialize local variables
    n = 400
    nr = 96
    af = 0.175
    bb = 3.418
    ac = 4.277
    ad = 1.685
    b = 3.69895
    xstep = 0.5
    a3 = 0.08
    k  = 0.4
    sqrt2pi = numpy.sqrt(2*numpy.pi)
    a4 = 3.25
    
    lall  = (0.000000, 0.302000, 0.368000, 0.414000, 0.450000,
        0.482000, 0.510000, 0.536000, 0.560000, 0.579999,
        0.601999, 0.621999, 0.639999, 0.657998, 0.675998,
        0.691998, 0.709998, 0.725998, 0.741997, 0.755997,
        0.771997, 0.785997, 0.801997, 0.815996, 0.829996,
        0.843996, 0.857996, 0.871996, 0.885995, 0.899995,
        0.911995, 0.925995, 0.939995, 0.953995, 0.965994,
        0.979994, 0.993994, 1.005994, 1.019994, 1.033994,
        1.045993, 1.059993, 1.073993, 1.085993, 1.099993,
        1.113993, 1.127992, 1.141992, 1.155992, 1.169992,
        1.183992, 1.197991, 1.211991, 1.225991, 1.239991,
        1.253991, 1.269991, 1.283990, 1.299990, 1.315990,
        1.329990, 1.345990, 1.361989, 1.379989, 1.395989,
        1.411989, 1.429989, 1.447988, 1.465988, 1.483988,
        1.501988, 1.521987, 1.539987, 1.559987, 1.581987,
        1.601986, 1.623986, 1.647986, 1.669985, 1.693985,
        1.719985, 1.745984, 1.773984, 1.801984, 1.831983,
        1.863983, 1.895983, 1.931982, 1.969982, 2.009982,
        2.053984, 2.101986, 2.153988, 2.211991, 2.279994,
        2.355998, 2.450002, 2.566008, 2.724015, 2.978027, 
        3.864068)
    
    a = (af/(bb-numpy.log(znot)))
    c = (ac*(bb-numpy.log(znot)))
    d = (ad*(bb-numpy.log(znot)))
    axa = numpy.zeros(n,dtype=numpy.float64) + a
    bxb = numpy.zeros(n,dtype=numpy.float64) + b
    cxc = numpy.zeros(n,dtype=numpy.float64) + c
    dxd = numpy.zeros(n,dtype=numpy.float64) + d
    zmxzm = numpy.zeros(n,dtype=numpy.float64) + zm
    sigmawxsigmaw = numpy.zeros(n,dtype=numpy.float64) + sigmaw
    ustarxustar = numpy.zeros(n,dtype=numpy.float64) + ustar
    coeff8xcoeff = numpy.zeros(n,dtype=numpy.float64) + 0.8
    negcoeff8x = numpy.zeros(n,dtype=numpy.float64) - 0.8
    hxh = numpy.zeros(n,dtype=numpy.float64) + h
    
    xstar = numpy.ones(n, dtype=numpy.float64)
    fstar = numpy.ones(n, dtype=numpy.float64)
    x = numpy.ones(n, dtype=numpy.float64)
    f_ci = numpy.ones(n, dtype=numpy.float64)
    
    xstar[0] = -5
    while xstar[0] < -d:
        xstar[0] = xstar[0]+1
    
    for i in range(0, n):
        # Calculate X*
        if i>0:
            xstar[i] = xstar[i-1] + xstep
    
    # Calculate F*
    fstar = axa*((xstar+dxd)/cxc)**bxb * numpy.exp(bxb*(1-(xstar+dxd)/cxc))
    
    # Calculate x and f
    x = xstar * zmxzm * (sigmawxsigmaw/ustar)**(negcoeff8x)
    f_ci = fstar / zmxzm * (1-(zmxzm/hxh)) * (sigmawxsigmaw/ustar)**(coeff8xcoeff)
    
    # Calculate maximum location of influence (peak location)
    xstarmax = c-d
    xmax = xstarmax * zm *(sigmaw/ustar)**(-0.8)
    
    # Get l corresponding to r
    lr = lall[r]
    
    # Calculate distance including R percentage of the footprint
    xstarr = lr*c - d
    xr = xstarr * zm *(sigmaw/ustar)**(-0.8)
    
    if qcutils.cfkeycheck(cf, Base='Output', ThisOne='FootprintSummaryFile') and cf['Output']['FootprintSummaryFile'] == 'True':
        xlFile = xlwt.Workbook()
        xlSheet = xlFile.add_sheet('summary')
        xlCol = 1
        xlRow = 0
        xlSheet.write(xlRow,xlCol,'Measurement height (zm)')
        xlCol = xlCol - 1
        xlSheet.write(xlRow,xlCol,zm)
        xlRow = xlRow + 1
        xlCol = 1
        xlSheet.write(xlRow,xlCol,'Canopy height (zc)')
        xlCol = xlCol - 1
        xlSheet.write(xlRow,xlCol,zc)
        xlRow = xlRow + 1
        xlCol = 1
        xlSheet.write(xlRow,xlCol,'Roughness length (z0)')
        xlCol = xlCol - 1
        xlSheet.write(xlRow,xlCol,znot)
        xlRow = xlRow + 1
        xlCol = 1
        xlSheet.write(xlRow,xlCol,'PBL height (zi)')
        xlCol = xlCol - 1
        xlSheet.write(xlRow,xlCol,h)
        xlRow = xlRow + 1
        xlCol = 1
        xlSheet.write(xlRow,xlCol,'Monin-Obukhov length (L)')
        xlCol = xlCol - 1
        xlSheet.write(xlRow,xlCol,L)
        xlRow = xlRow + 1
        xlCol = 1
        xlSheet.write(xlRow,xlCol,'Stability coefficient (z-d/L, zeta)')
        xlCol = xlCol - 1
        xlSheet.write(xlRow,xlCol,zeta)
        xlRow = xlRow + 1
        xlCol = 1
        xlSheet.write(xlRow,xlCol,'Friction coefficient (ustar)')
        xlCol = xlCol - 1
        xlSheet.write(xlRow,xlCol,ustar)
        xlRow = xlRow + 1
        xlCol = 1
        xlSheet.write(xlRow,xlCol,'sd(w) (sigma_w)')
        xlCol = xlCol - 1
        xlSheet.write(xlRow,xlCol,sigmaw)
        xlRow = xlRow + 1
        xlCol = 1
        xlSheet.write(xlRow,xlCol,'sd(v) (sigma_v)')
        xlCol = xlCol - 1
        xlSheet.write(xlRow,xlCol,sigmav)
        xlRow = xlRow + 1
        xlCol = 1
        if footprintlevel == 'V3':
            xlSheet.write(xlRow,xlCol,'Fc')
            xlCol = xlCol - 1
            xlSheet.write(xlRow,xlCol,Fc)
            xlRow = xlRow + 1
            xlCol = 1
            xlSheet.write(xlRow,xlCol,'Fe')
            xlCol = xlCol - 1
            xlSheet.write(xlRow,xlCol,Fe)
            xlRow = xlRow + 1
            xlCol = 1
            xlSheet.write(xlRow,xlCol,'Fh')
            xlCol = xlCol - 1
            xlSheet.write(xlRow,xlCol,Fh)
        elif footprintlevel == 'V2':
            xlSheet.write(xlRow,xlCol,'GPP')
            xlCol = xlCol - 1
            xlSheet.write(xlRow,xlCol,Fc)
            xlRow = xlRow + 1
            xlCol = 1
            xlSheet.write(xlRow,xlCol,'CE')
            xlCol = xlCol - 1
            xlSheet.write(xlRow,xlCol,Fe)
        
        xlRow = xlRow + 1
        xlCol = 1
        
        xlCol = 3
        xlRow = 0
        xlSheet.write(xlRow,xlCol,'eta')
        xlCol = xlCol - 1
        xlSheet.write(xlRow,xlCol,eta)
        xlRow = xlRow + 1
        xlCol = 3
        xlSheet.write(xlRow,xlCol,'Wind direction (wd)')
        xlCol = xlCol - 1
        xlSheet.write(xlRow,xlCol,wd)
        xlRow = xlRow + 1
        xlCol = 3
        xlSheet.write(xlRow,xlCol,'% of flux footprint (r)')
        xlCol = xlCol - 1
        xlSheet.write(xlRow,xlCol,r)
        xlRow = xlRow + 1
        xlCol = 3
        xlSheet.write(xlRow,xlCol,'Extent of footprint up to r% (xr)')
        xlCol = xlCol - 1
        xlSheet.write(xlRow,xlCol,xr)
        xlRow = xlRow + 1
        
        xlSheet = xlFile.add_sheet('crosswind_integrated')
        for i in range(0,n):
            xlSheet.write(i,0,x[i])
            xlSheet.write(i,1,f_ci[i])
        
        xlFile.save(summaryFileName)
    
    if qcutils.cfkeycheck(cf, Base='Output', ThisOne='FootprintDataFile') and cf['Output']['FootprintDataFile'] == 'True':
        # Calculate lateral dispersion
        u = ustar/k *(numpy.log(zm/znot) - (zm-znot)/zm)
        tau = numpy.sqrt((x/u)**2 + (a4*(zm-znot)/sigmaw)**2)
        tly = a3*h**2 /((h-zm)*ustar)
        fy_disp = 1/(1 + numpy.sqrt(tau/(2*tly)))
        sigmay = sigmav * tau * fy_disp
        
        x_lim = numpy.max(x)
        y0 = (x[x>0])
        y1 = (y0[y0<=x_lim/2])
        y2 = -(y1[::-1])
        y  = numpy.concatenate((y2,[0],y1))
        
        m = len(y)
        
        y_rot = numpy.reshape(y,(m,1))
        xy = numpy.broadcast_arrays(x, f_ci, sigmay, y_rot)
        if not qcutils.cfkeycheck(cf, Base='Footprint', ThisOne='etaadd'):
            log.error('   Footprint:  CSAT azimuth not provided for coordinate rotation')
            return
        
        x_2d = (xy[0] * numpy.cos(numpy.deg2rad(-eta+numpy.float64(cf['Footprint']['etaadd'])))) + (xy[3] * numpy.sin(numpy.deg2rad(-eta+numpy.float64(cf['Footprint']['etaadd']))))   # longitudal rotation of the x,y plane
        y_2d = (xy[3] * numpy.cos(numpy.deg2rad(-eta+numpy.float64(cf['Footprint']['etaadd'])))) - (xy[0] * numpy.sin(numpy.deg2rad(-eta+numpy.float64(cf['Footprint']['etaadd']))))   # lateral rotation of the x,y plane
        f_2d = xy[1] * 1/(sqrt2pi*xy[2]) *  numpy.exp(-xy[3]**2 / (2*xy[2]**2))
        
        if qcutils.cfkeycheck(cf,Base='Output',ThisOne='FootprintDataFileType') and (cf['Output']['FootprintDataFileType'] == 'Weighted' or cf['Output']['FootprintDataFileType'] == 'Climatology'):
            check = numpy.zeros(5,dtype=numpy.int32)
            matrixfiles = []
            vectorfiles = []
            labels = []
            if footprintlevel == 'V3':
                if cf['Output']['FootprintDataFileType'] == 'Weighted':
                    labels.append('fw_Fc')
                    labels.append('fw_Fe')
                    labels.append('fw_Fh')
                    vectorfiles.append(cf['Files']['Footprint']['FootprintFilePath']+'Fc_w_footprint_2d_vectors'+''.join(STList)+'.xls')
                    matrixfiles.append(cf['Files']['Footprint']['FootprintFilePath']+'Fc_w_footprint_2d_matrix'+''.join(STList)+'.xls')
                    vectorfiles.append(cf['Files']['Footprint']['FootprintFilePath']+'Fe_w_footprint_2d_vectors'+''.join(STList)+'.xls')
                    matrixfiles.append(cf['Files']['Footprint']['FootprintFilePath']+'Fe_w_footprint_2d_matrix'+''.join(STList)+'.xls')
                    vectorfiles.append(cf['Files']['Footprint']['FootprintFilePath']+'Fh_w_footprint_2d_vectors'+''.join(STList)+'.xls')
                    matrixfiles.append(cf['Files']['Footprint']['FootprintFilePath']+'Fh_w_footprint_2d_matrix'+''.join(STList)+'.xls')
                    fw_c1_2d = f_2d * Fc
                    fw_e_2d = f_2d * Fe
                    fw_h_2d = f_2d * Fh
                
                if cf['Output']['FootprintDataFileType'] == 'Climatology':
                    fw_c1_2d = numpy.zeros((len(f_2d),len(f_2d[0])),dtype=numpy.float64)
                    fw_c2_2d = numpy.zeros((len(f_2d),len(f_2d[0])),dtype=numpy.float64)
                    fw_c3_2d = numpy.zeros((len(f_2d),len(f_2d[0])),dtype=numpy.float64)
                    fw_e_2d = numpy.zeros((len(f_2d),len(f_2d[0])),dtype=numpy.float64)
                    fw_h_2d = numpy.zeros((len(f_2d),len(f_2d[0])),dtype=numpy.float64)
                
                if cf['Output']['FootprintDataFileType'] == 'Climatology' and (qcutils.cfkeycheck(cf,Base='Footprint',ThisOne='FcBoth') and cf['Footprint']['FcBoth'] == 'True'):
                    labels.append('fc_Fc')
                    vectorfiles.append(cf['Files']['Footprint']['FootprintFilePath']+'Fc_c_footprint_2d_vectors'+''.join(STList)+'.xls')
                    matrixfiles.append(cf['Files']['Footprint']['FootprintFilePath']+'Fc_c_footprint_2d_matrix'+''.join(STList)+'.xls')
                    fw_c1_2d = f_2d * Fc
                
                if cf['Output']['FootprintDataFileType'] == 'Climatology' and (qcutils.cfkeycheck(cf,Base='Footprint',ThisOne='Fe') and cf['Footprint']['Fe'] == 'True'):
                    labels.append('fc_Fe')
                    vectorfiles.append(cf['Files']['Footprint']['FootprintFilePath']+'Fe_c_footprint_2d_vectors'+''.join(STList)+'.xls')
                    matrixfiles.append(cf['Files']['Footprint']['FootprintFilePath']+'Fe_c_footprint_2d_matrix'+''.join(STList)+'.xls')
                    fw_e_2d = f_2d * Fe
                
                if cf['Output']['FootprintDataFileType'] == 'Climatology' and (qcutils.cfkeycheck(cf,Base='Footprint',ThisOne='Fh') and cf['Footprint']['Fh'] == 'True'):
                    labels.append('fc_Fh')
                    vectorfiles.append(cf['Files']['Footprint']['FootprintFilePath']+'Fh_c_footprint_2d_vectors'+''.join(STList)+'.xls')
                    matrixfiles.append(cf['Files']['Footprint']['FootprintFilePath']+'Fh_c_footprint_2d_matrix'+''.join(STList)+'.xls')
                    fw_h_2d = f_2d * Fh
                
                if cf['Output']['FootprintDataFileType'] == 'Climatology' and (qcutils.cfkeycheck(cf,Base='Footprint',ThisOne='FcPos') and cf['Footprint']['FcPos'] == 'True'):
                    labels.append('fc_PosFc')
                    vectorfiles.append(cf['Files']['Footprint']['FootprintFilePath']+'PosFc_c_footprint_2d_vectors'+''.join(STList)+'.xls')
                    matrixfiles.append(cf['Files']['Footprint']['FootprintFilePath']+'PosFc_c_footprint_2d_matrix'+''.join(STList)+'.xls')
                    if Fc > 0:
                        fw_c2_2d = f_2d * Fc
                
                if cf['Output']['FootprintDataFileType'] == 'Climatology' and (qcutils.cfkeycheck(cf,Base='Footprint',ThisOne='FcNeg') and cf['Footprint']['FcNeg'] == 'True'):
                    labels.append('fc_NegFc')
                    vectorfiles.append(cf['Files']['Footprint']['FootprintFilePath']+'NegFc_c_footprint_2d_vectors'+''.join(STList)+'.xls')
                    matrixfiles.append(cf['Files']['Footprint']['FootprintFilePath']+'NegFc_c_footprint_2d_matrix'+''.join(STList)+'.xls')
                    if Fc < 0:
                        fw_c3_2d = f_2d * Fc
            
            if footprintlevel == 'V2':
                fw_c1_2d = f_2d * Fc
                fw_e_2d = f_2d * Fe
                if cf['Output']['FootprintDataFileType'] == 'Weighted':
                    labels.append('fw_GPP')
                    labels.append('fw_CE')
                    vectorfiles.append(cf['Files']['Footprint']['FootprintFilePath']+'GPP_w_footprint_2d_vectors'+''.join(STList)+'.xls')
                    matrixfiles.append(cf['Files']['Footprint']['FootprintFilePath']+'GPP_w_footprint_2d_matrix'+''.join(STList)+'.xls')
                    vectorfiles.append(cf['Files']['Footprint']['FootprintFilePath']+'CE_w_footprint_2d_vectors'+''.join(STList)+'.xls')
                    matrixfiles.append(cf['Files']['Footprint']['FootprintFilePath']+'CE_w_footprint_2d_matrix'+''.join(STList)+'.xls')
                
                if cf['Output']['FootprintDataFileType'] == 'Climatology':
                    labels.append('fc_GPP')
                    labels.append('fc_CE')
                    vectorfiles.append(cf['Files']['Footprint']['FootprintFilePath']+'GPP_c_footprint_2d_vectors'+''.join(STList)+'.xls')
                    matrixfiles.append(cf['Files']['Footprint']['FootprintFilePath']+'GPP_c_footprint_2d_matrix'+''.join(STList)+'.xls')
                    vectorfiles.append(cf['Files']['Footprint']['FootprintFilePath']+'CE_c_footprint_2d_vectors'+''.join(STList)+'.xls')
                    matrixfiles.append(cf['Files']['Footprint']['FootprintFilePath']+'CE_c_footprint_2d_matrix'+''.join(STList)+'.xls')
            
            filenames = [vectorfiles, matrixfiles]
        
        if qcutils.cfkeycheck(cf, Base='Output', ThisOne='FootprintDataFileType') and cf['Output']['FootprintDataFileType'] == 'Footprint':
            if ((qcutils.cfkeycheck(cf, Base='Output', ThisOne='FootprintDataType')) and ((cf['Output']['FootprintDataType'] == 'Vector') or (cf['Output']['FootprintDataType'] == 'Both'))):
                footprint_vector_out(vectorFileName,x_2d,y_2d,f_2d,'f')
            
            if ((qcutils.cfkeycheck(cf, Base='Output', ThisOne='FootprintDataType')) and ((cf['Output']['FootprintDataType'] == 'Matrix') or (cf['Output']['FootprintDataType'] == 'Both'))):
                footprint_matrix_out(matrixFileName,x_2d,y_2d,f_2d)
        
        if qcutils.cfkeycheck(cf, Base='Output', ThisOne='FootprintDataFileType') and cf['Output']['FootprintDataFileType'] == 'Weighted':
            if ((qcutils.cfkeycheck(cf, Base='Output', ThisOne='FootprintDataType')) and ((cf['Output']['FootprintDataType'] == 'Vector') or (cf['Output']['FootprintDataType'] == 'Both'))):
                footprint_vector_out(filenames[0][0],x_2d,y_2d,fw_c1_2d,labels[0])
                footprint_vector_out(filenames[0][1],x_2d,y_2d,fw_e_2d,labels[1])
                footprint_vector_out(filenames[0][2],x_2d,y_2d,fw_h_2d,labels[2])
            
            if ((qcutils.cfkeycheck(cf, Base='Output', ThisOne='FootprintDataType')) and ((cf['Output']['FootprintDataType'] == 'Matrix') or (cf['Output']['FootprintDataType'] == 'Both'))):
                footprint_matrix_out(filenames[1][0],x_2d,y_2d,fw_c1_2d)
                footprint_matrix_out(filenames[1][1],x_2d,y_2d,fw_e_2d)
                footprint_matrix_out(filenames[1][2],x_2d,y_2d,fw_h_2d)
            
        if qcutils.cfkeycheck(cf, Base='Output', ThisOne='FootprintDataFileType') and cf['Output']['FootprintDataFileType'] == 'Climatology':
            if footprintlevel == 'V3':
                return xr, x_2d, y_2d, fw_c1_2d, fw_e_2d, fw_h_2d, fw_c2_2d, fw_c3_2d, filenames, labels
            elif footprintlevel == 'V2':
                return xr, x_2d, y_2d, fw_c1_2d, fw_e_2d, filenames, labels

    return xr

def footprint_climatology(fc_c,fc_e,x_2d,y_2d,fw_c,fw_e,n,m,xmin,xmax,ymin,ymax,p,fc_h='',fc_c2='',fc_c3='',fw_h='',fw_c2='',fw_c3=''):
    
    x = numpy.arange(xmin+(0.5*p),xmax,p,dtype=numpy.float64)
    y = numpy.arange(ymin+(0.5*p),ymax,p,dtype=numpy.float64)
    y_rot = numpy.reshape(y,(len(y),1))
    xy = numpy.broadcast_arrays(x, y_rot)
    jj = 0
    for j in range(numpy.int32(xmin),numpy.int32(xmax),numpy.int32(p)):
        ii = 0
        for i in range(numpy.int32(ymin),numpy.int32(ymax),numpy.int32(p)):
            xindex = numpy.where((x_2d > j) & (x_2d < j + p))
            yindex = numpy.where((y_2d > i) & (y_2d < i + p))
            if (len(xindex[0]) != 0) & (len(yindex[0]) != 0):
                mask = numpy.zeros((len(x_2d),len(y_2d[0])),dtype=numpy.int32)
                mask[xindex] = 1
                mask[yindex] = mask[yindex] + 1
                maskindex = numpy.where(mask > 1)
                if len(maskindex[0]) != 0:
                    fc_c[ii,jj] = fc_c[ii,jj] + numpy.mean(fw_c[maskindex])
                    fc_e[ii,jj] = fc_e[ii,jj] + numpy.mean(fw_e[maskindex])
                    if fc_h != '':
                        fc_h[ii,jj] = fc_h[ii,jj] + numpy.mean(fw_h[maskindex])
                        fc_c2[ii,jj] = fc_c2[ii,jj] + numpy.mean(fw_c2[maskindex])
                        fc_c3[ii,jj] = fc_c3[ii,jj] + numpy.mean(fw_c3[maskindex])
            
            ii = ii + 1
        
        jj = jj + 1
    
    if fc_h == '':
        return fc_c, fc_e, xy[0], xy[1]
    else:
        return fc_c, fc_e, fc_h, fc_c2, fc_c3, xy[0], xy[1]

def footprint_matrix_out(filename,x_2d,y_2d,f_2d):
    n = len(x_2d)
    m = len(x_2d[0])
    out = numpy.zeros((n+1,m+1),dtype=numpy.float64)
    out[0,0] = c.missing_value
    out[0,1:] = x_2d[0,:]
    out[1:,0] = y_2d[:,0]
    out[1:,1:] = f_2d[:,:]
    numpy.savetxt(filename,out,delimiter='\t')
    
    return

def footprint_vector_out(filename,x_2d,y_2d,f_2d,text):
    log.info('     '+filename+' cols: x, y, '+text)
    x = numpy.reshape(x_2d,(1,(len(x_2d)*len(x_2d[0]))))
    y = numpy.reshape(y_2d,(1,(len(x_2d)*len(x_2d[0]))))
    f = numpy.reshape(f_2d,(1,(len(x_2d)*len(x_2d[0]))))
    out = numpy.vstack((x,y,f)).T
    numpy.savetxt(filename,out,delimiter='\t')
    
    return

def get_averages(Data):
    """
        Get daily averages on days when no 30-min observations are missing.
        Days with missing observations return a value of c.missing_value
        Values returned are sample size (Num) and average (Av)
        
        Usage qcts.get_averages(Data)
        Data: 1-day dataset
        """
    li = numpy.ma.where(abs(Data-numpy.float64(c.missing_value))>c.eps)
    Num = numpy.size(li)
    if Num == 0:
        Av = c.missing_value
    elif Num == 48:
        Av = numpy.ma.mean(Data[li])
    else:
        x = 0
        index = numpy.ma.where(Data.mask == True)[0]
        if len(index) == 1:
            x = 1
        elif len(index) > 1:
            for i in range(len(Data)):
                if Data.mask[i] == True:
                    x = x + 1
        
        if x == 0:
            Av = numpy.ma.mean(Data[li])
        else:
            Av = c.missing_value
    return Num, Av

def get_averages_WUE(Data):
    """
        Get daily averages for WUE with missing observations.
        Values returned are sample size (Num) and average (Av)
        
        Usage qcts.get_averages_WUE(Data)
        Data: 1-day dataset
        """
    li = numpy.ma.where(abs(Data-numpy.float64(c.missing_value))>c.eps)
    Num = numpy.size(li)
    if Num == 0:
        Av = c.missing_value
    else:
        Av = numpy.ma.mean(Data[li])
    
    return Num, Av

def get_canopyresistance(cf,ds,Uavg,uindex,PMin,Level,critFsd,critFe):
    if qcutils.cfkeycheck(cf,Base='PenmanMonteith',ThisOne='zm'):
        zm = numpy.float64(cf['PenmanMonteith']['zm'])
    
    if qcutils.cfkeycheck(cf,Base='PenmanMonteith',ThisOne='z0'):
        z0m = numpy.float64(cf['PenmanMonteith']['z0'])
    
    if qcutils.cfkeycheck(cf,Base='PenmanMonteith',ThisOne='zc'):
        zc = numpy.float64(cf['PenmanMonteith']['zc'])
    
    Fe,f,a = qcutils.GetSeriesasMA(ds,PMin[0])
    Ta,f,a = qcutils.GetSeriesasMA(ds,PMin[1])
    Ah,f,a = qcutils.GetSeriesasMA(ds,PMin[2])
    ps,f,a = qcutils.GetSeriesasMA(ds,PMin[3])
    Fn,f,a = qcutils.GetSeriesasMA(ds,PMin[5])
    Fsd,f,a = qcutils.GetSeriesasMA(ds,PMin[6])
    Fg,f,a = qcutils.GetSeriesasMA(ds,PMin[7])
    flagList = [PMin[0],PMin[1],PMin[2],PMin[3],PMin[4],PMin[5],PMin[6],PMin[7]]
    VPD,f,a = qcutils.GetSeriesasMA(ds,'VPD')
    Lv,f,a = qcutils.GetSeriesasMA(ds,'Lv')
    Cpm,f,a = qcutils.GetSeriesasMA(ds,'Cpm')
    rhom,f,a = qcutils.GetSeriesasMA(ds,'rhom')
    gamma = mf.gamma(ps,Cpm,Lv)
    delta = mf.delta(Ta)
    if 'gamma' not in ds.series.keys():
        attr = qcutils.MakeAttributeDictionary(long_name='Psychrometric coefficient',units='kPa/C')
        qcutils.CreateSeries(ds,'gamma',gamma,FList=[PMin[3],'Cpm','Lv'],Attr=attr)
    
    if 'delta' not in ds.series.keys():
        attr = qcutils.MakeAttributeDictionary(long_name='Slope of the saturation vapour pressure v temperature curve',units='kPa/C')
        qcutils.CreateSeries(ds,'delta',delta,FList=[PMin[1]],Attr=attr)
    
    Feindex = numpy.ma.where(Fe < critFe)[0]
    Fsdindex = numpy.ma.where(Fsd < critFsd)[0]
    #Fnindex = numpy.where(Fn < 0)[0]
    z0v = 0.1 * z0m
    d = (2 / 3) * zc
    ra = (numpy.log((zm - d) / z0m) * numpy.log((zm - d) / z0v)) / (((c.k) ** 2) * Uavg)
    rc = ((((((delta * (Fn - Fg) / (Lv)) + (rhom * Cpm * (VPD / ((Lv) * ra)))) / (Fe / (Lv))) - delta) / gamma) - 1) * ra
    rcindex = numpy.ma.where(rc < 0)[0]
    Gc = (1 / rc) * (rhom * 1000) * (1000 / 18)
    if 'ram' not in ds.series.keys():
        attr = qcutils.MakeAttributeDictionary(long_name='Aerodynamic resistance from drag coefficient, Allen/Jensen formulation, '+Level,units='s/m')
    else:
        attr = qcutils.MakeAttributeDictionary(long_name=Level,units='s/m')
    qcutils.CreateSeries(ds,'ram',ra,FList=flagList,Attr=attr)
    if 'rSm' not in ds.series.keys():
        attr = qcutils.MakeAttributeDictionary(long_name='Canopy resistance from Penman-Monteith inversion, Allen/Jensen formulation, '+Level,units='s/m')
    else:
        attr = qcutils.MakeAttributeDictionary(long_name=Level,units='s/m')
    qcutils.CreateSeries(ds,'rSm',rc,FList=flagList,Attr=attr)
    if 'GSm' not in ds.series.keys():
        attr = qcutils.MakeAttributeDictionary(long_name='Canopy conductance from Penman-Monteith inversion, Allen/Jensen formulation, '+Level,units='mmolH2O/(m2ground s)')
    else:
        attr = qcutils.MakeAttributeDictionary(long_name=Level,units='mmolH2O/(m2ground s)')
    qcutils.CreateSeries(ds,'GSm',Gc,FList=flagList,Attr=attr)
    
    Label = ['ram','rSm','GSm']
    for listindex in range(0,3):
        ds.series[Label[listindex]]['Attr']['InputSeries'] = PMin
        ds.series[Label[listindex]]['Attr']['FsdCutoff'] = critFsd
        ds.series[Label[listindex]]['Attr']['FeCutoff'] = critFe
        ds.series[Label[listindex]]['Flag'][rcindex] = numpy.int32(171)
        ds.series[Label[listindex]]['Flag'][uindex] = numpy.int32(174)
        ds.series[Label[listindex]]['Flag'][Feindex] = numpy.int32(172)
        ds.series[Label[listindex]]['Flag'][Fsdindex] = numpy.int32(173)
        ds.series[Label[listindex]]['Data'][rcindex] = numpy.float64(c.missing_value)
        ds.series[Label[listindex]]['Data'][uindex] = numpy.float64(c.missing_value)
        ds.series[Label[listindex]]['Data'][Feindex] = numpy.float64(c.missing_value)
        ds.series[Label[listindex]]['Data'][Fsdindex] = numpy.float64(c.missing_value)
    
    return

def get_leafresistance(cf,ds,rinverted):
    Fsd,Fsd_flag,a = qcutils.GetSeriesasMA(ds,'Fsd')
    Hdh,f,a = qcutils.GetSeriesasMA(ds,'Hdh')
    Day,f,a = qcutils.GetSeriesasMA(ds,'Day')
    Month,f,a = qcutils.GetSeriesasMA(ds,'Month')
    nRecs = len(Fsd)
    if 'LAI' not in ds.series.keys():
        log.info('  Penman-Monteith: integrating daily LAI file into data structure')
        InLevel = 'L4LAI'
        OutLevel = 'L4LAI'
        qcio.autoxl2nc(cf,InLevel,OutLevel)
        dsLAI = qcio.nc_read_series(cf,'L4LAI')
        LAI,fd,a = qcutils.GetSeriesasMA(dsLAI,'LAI')
        Day_LAI,fd,a = qcutils.GetSeriesasMA(dsLAI,'Day')
        Month_LAI,fd,a = qcutils.GetSeriesasMA(dsLAI,'Month')
        nDays = len(LAI)
        LAI_expanded = numpy.ma.zeros(nRecs,numpy.float64)
        Night = numpy.ma.zeros(nRecs)
        LAI_flag = numpy.ma.zeros(nRecs,dtype=numpy.int32)
        
        for i in range(nRecs):
            if Month[i] == 1 or Month[i] == 3 or Month[i] == 5 or Month[i] == 7 or Month[i] == 8 or Month[i] == 10 or Month[i] == 12:
                dRan = 31
            if Month[i] == 2:
                if ds.series['Year']['Data'][0] % 4 == 0:
                    dRan = 29
                else:
                    dRan = 28
            if Month[i] == 4 or Month[i] == 6 or Month[i] == 9 or Month[i] == 11:
                dRan = 30
                
            Night[i] = Day[i]
        
        log.info(' Penman-Monteith: filling LAI from daily LAI')
        for z in range(nDays):
            for i in range(nRecs):
                if Night[i] == Day_LAI[z]:
                    if Month[i] == Month_LAI[z]:
                        LAI_expanded[i] = LAI[z]
                        LAI_flag[i] = dsLAI.series['LAI']['Flag'][z]
        
        attr = qcutils.MakeAttributeDictionary(long_name='Leaf area index, spline-fit interpolation from MODIS product',units='m2/m2',standard_name='leaf_area_index')
        qcutils.CreateSeries(ds,'LAI',LAI_expanded,Flag=0,Attr=attr)
        ds.series['LAI']['Flag'] = LAI_flag
    else:
        LAI_expanded, LAI_flag,a = qcutils.GetSeriesasMA(ds,'LAI')
    
    log.info('  Penman-Monteith: computing leaf resistance from inversion surface resistance and LAI')
    rSmin, rc_flag,a = qcutils.GetSeriesasMA(ds,rinverted)
    rl = rSmin * 0.5 * LAI_expanded
    rl_flag = numpy.zeros(nRecs,numpy.int32)
    Fsd_index = numpy.where(Fsd < 600)[0]
    Fsd2_index = numpy.where(numpy.mod(Fsd_flag,10)!=0)[0]
    rSm_index = numpy.where(numpy.mod(rc_flag,10)!=0)[0]
    LAI_index = numpy.where(numpy.mod(LAI_flag,10)!=0)[0]
    rl[Fsd_index] = numpy.float64(c.missing_value)
    rl_flag[Fsd_index] = numpy.int32(173)
    rl[Fsd2_index] = numpy.float64(c.missing_value)
    rl_flag[Fsd2_index] = Fsd_flag[Fsd2_index]
    rl[rSm_index] = numpy.float64(c.missing_value)
    rl_flag[rSm_index] = rc_flag[rSm_index]
    rl[LAI_index] = numpy.float64(c.missing_value)
    rl_flag[LAI_index] = LAI_flag[LAI_index]
    return rl, rl_flag

def get_minmax(Data):
    """
        Get daily minima and maxima on days when no 30-min observations are missing.
        Days with missing observations return a value of c.missing_value
        Values returned are sample size (Num), minimum (Min) and maximum (Max)
        
        Usage qcts.get_minmax(Data)
        Data: 1-day dataset
        """
    li = numpy.ma.where(abs(Data-numpy.float64(c.missing_value))>c.eps)
    Num = numpy.size(li)
    if Num == 0:
        Min = c.missing_value
        Max = c.missing_value
    elif Num == 48:
        Min = numpy.ma.min(Data[li])
        Max = numpy.ma.max(Data[li])
    else:
        x = 0
        index = numpy.ma.where(Data.mask == True)[0]
        if len(index) == 1:
            x = 1
        elif len(index) > 1:
            for i in range(len(Data)):
                if Data.mask[i] == True:
                    x = x + 1
        
        if x == 0:
            Min = numpy.ma.min(Data[li])
            Max = numpy.ma.max(Data[li])
        else:
            Min = numpy.ma.min(Data[li])
            Max = numpy.ma.max(Data[li])
            #Min = c.missing_value
            #Max = c.missing_value
    return Num, Min, Max

def get_minmax_WUE(Data):
    """
        Get daily minima and maxima for WUE with missing obs.
        Values returned are sample size (Num), minimum (Min) and maximum (Max)
        
        Usage qcts.get_minmax(Data)
        Data: 1-day dataset
        """
    li = numpy.ma.where(abs(Data-numpy.float64(c.missing_value))>c.eps)
    Num = numpy.size(li)
    if Num == 0:
        Min = c.missing_value
        Max = c.missing_value
    else:
        Min = numpy.ma.min(Data[li])
        Max = numpy.ma.max(Data[li])
    
    return Num, Min, Max

def get_nightsums(Data):
    """
        Get nightly sums and averages on nights when no 30-min observations are missing.
        Nights with missing observations return a value of c.missing_value
        Values returned are sample size (Num), sums (Sum) and average (Av)
        
        Usage qcts.get_nightsums(Data)
        Data: 1-day dataset
        """
    li = numpy.ma.where(Data.mask == False)[0]
    Num = numpy.size(li)
    if Num == 0:
        Sum = c.missing_value
        Av = c.missing_value
    else:
        x = 0
        for i in range(len(Data)):
            if Data.mask[i] == True:
                x = x + 1
        
        if x == 0:
            Sum = numpy.ma.sum(Data[li])
            Av = numpy.ma.mean(Data[li])
        else:
            Sum = c.missing_value
            Av = c.missing_value
    
    return Num, Sum, Av

def get_rav(cf,ds,Uavg,PMin,qList,layer='',method='Cemethod'):
    if qcutils.cfkeycheck(cf,Base='PenmanMonteith',ThisOne=method):
        if method == 'Ce_2layer' and layer == '':
            bothlayers = 'True'
            layer = 'base_'
        else:
            bothlayers = 'False'
        
        if qcutils.cfkeycheck(cf,Base='PenmanMonteith',ThisOne=layer+'zq_low'):
            zq_low = numpy.float64(cf['PenmanMonteith'][layer+'zq_low'])
        else:
            log.error('  PenmanMonteith:  zq_low (height of lower q sensor) not given')
            return
        
        if qcutils.cfkeycheck(cf,Base='PenmanMonteith',ThisOne=layer+'zq_surface'):
            zq_surface = numpy.float64(cf['PenmanMonteith'][layer+'zq_surface'])
        else:
            log.error('  PenmanMonteith:  zq_surface (height of q_surface extrapolation target) not given')
            return
        
        if bothlayers == 'True':
            layer = 'top_'
        
        if qcutils.cfkeycheck(cf,Base='PenmanMonteith',ThisOne=layer+'zq_high'):
            zq_high = numpy.float64(cf['PenmanMonteith'][layer+'zq_high'])
        else:
            log.error('  PenmanMonteith:  zq_high (height of upper q sensor) not given')
            return
        
        if qcutils.cfkeycheck(cf,Base='PenmanMonteith',ThisOne=layer+'zq_air'):
            zq_air = numpy.float64(cf['PenmanMonteith'][layer+'zq_air'])
        else:
            log.error('  PenmanMonteith:  zq_air (height of q_air extrapolation target) not given')
            return
        
        if bothlayers == 'True':
            layer = 'full_'
        
        log.info('   Ce method:  '+layer+'zq_high: '+str(zq_high)+', '+layer+'zq_low: '+str(zq_low))
        log.info('   Ce method:  '+layer+'zq_air: '+str(zq_air)+', '+layer+'zq_surface: '+str(zq_surface))
        q_high,fqh,a = qcutils.GetSeriesasMA(ds,qList[0])
        q_low,fql,a = qcutils.GetSeriesasMA(ds,qList[1])
        if zq_low == 0:
            qsurface = numpy.zeros(len(q_low), dtype=numpy.float64) + q_low
        elif zq_low == zq_surface:
            qsurface = numpy.zeros(len(q_low), dtype=numpy.float64) + q_low
        else:
            qsurface = extrapolate_humidity(zq_surface,zq_high,zq_low,q_high,q_low,fqh,fql)
        
        if zq_high == zq_air:
            qair = numpy.zeros(len(q_high), dtype=numpy.float64) + q_high
        else:
            qair = extrapolate_humidity(zq_air,zq_high,zq_low,q_high,q_low,fqh,fql)
        
    Fe,f,a = qcutils.GetSeriesasMA(ds,PMin[0])
    Ta,f,a = qcutils.GetSeriesasMA(ds,PMin[1])
    Ah,f,a = qcutils.GetSeriesasMA(ds,PMin[2])
    ps,f,a = qcutils.GetSeriesasMA(ds,PMin[3])
    Fn,f,a = qcutils.GetSeriesasMA(ds,PMin[5])
    Fsd,f,a = qcutils.GetSeriesasMA(ds,PMin[6])
    Fg,f,a = qcutils.GetSeriesasMA(ds,PMin[7])
    wA,f,a = qcutils.GetSeriesasMA(ds,PMin[8])
    Lv,f,a = qcutils.GetSeriesasMA(ds,'Lv')
    #if 'WPLcov' not in ds.globalattributes['L3Functions']:
    #    wA = Fe / (Lv * c.g2kg)
    if 'wA' not in ds.series.keys():
        wA = Fe / (Lv * c.g2kg)
    Ce = mf.bulktransfercoefficient(wA,Uavg,qair,qsurface)
    rav = mf.aerodynamicresistance(Uavg,Ce)
    ravindex = numpy.ma.where(rav < 0)[0]
    rav[ravindex] = numpy.float64(c.missing_value)
    return Ce, rav, ravindex

def get_rstGst(cf,ds,PMin,rav):
    Fe,f,a = qcutils.GetSeriesasMA(ds,PMin[0])
    Ta,f,a = qcutils.GetSeriesasMA(ds,PMin[1])
    Ah,f,a = qcutils.GetSeriesasMA(ds,PMin[2])
    ps,f,a = qcutils.GetSeriesasMA(ds,PMin[3])
    Fn,f,a = qcutils.GetSeriesasMA(ds,PMin[5])
    Fsd,f,a = qcutils.GetSeriesasMA(ds,PMin[6])
    Fg,f,a = qcutils.GetSeriesasMA(ds,PMin[7])
    VPD,f,a = qcutils.GetSeriesasMA(ds,'VPD')
    Lv,f,a = qcutils.GetSeriesasMA(ds,'Lv')
    Cpm,f,a = qcutils.GetSeriesasMA(ds,'Cpm')
    rhom,f,a = qcutils.GetSeriesasMA(ds,'rhom')
    gamma = mf.gamma(ps,Cpm,Lv)
    delta = mf.delta(Ta)
    if 'gamma' not in ds.series.keys():
        attr = qcutils.MakeAttributeDictionary(long_name='Psychrometric coefficient',units='kPa/C')
        qcutils.CreateSeries(ds,'gamma',gamma,FList=[PMin[3],'Cpm','Lv'],Attr=attr)
    
    if 'delta' not in ds.series.keys():
        attr = qcutils.MakeAttributeDictionary(long_name='Slope of the saturation vapour pressure v temperature curve',units='kPa/C')
        qcutils.CreateSeries(ds,'delta',delta,FList=[PMin[1]],Attr=attr)
    
    rst = ((((((delta * (Fn - Fg) / (Lv)) + (rhom * Cpm * (VPD / ((Lv) * rav)))) / (Fe / (Lv))) - delta) / gamma) - 1) * rav
    rstindex = numpy.ma.where(rst < 0)[0]
    Gst = (1 / rst) * (rhom * 1000) * (1000 / 18)
    Gstindex = numpy.ma.where(Gst < 0)[0]
    return rst, rstindex, Gst, Gstindex

def get_soilaverages(Data):
    """
        Get daily averages of soil water content on days when 15 or fewer 30-min observations are missing.
        Days with 16 or more missing observations return a value of c.missing_value
        Values returned are sample size (Num) and average (Av)
        
        Usage qcts.get_soilaverages(Data)
        Data: 1-day dataset
        """
    li = numpy.ma.where(abs(Data-numpy.float64(c.missing_value))>c.eps)
    Num = numpy.size(li)
    if Num > 28:
        Av = numpy.ma.mean(Data[li])
    else:
        Av = c.missing_value
    return Num, Av

def get_subsums(Data):
    """
        Get separate daily sums of positive and negative fluxes when no 30-min observations are missing.
        Days with missing observations return a value of c.missing_value
        Values returned are positive and negative sample sizes (PosNum and NegNum) and sums (SumPos and SumNeg)
        
        Usage qcts.get_subsums(Data)
        Data: 1-day dataset
        """
    li = numpy.ma.where(abs(Data-numpy.float64(c.missing_value))>c.eps)
    Num = numpy.size(li)
    if Num == 48:
        pi = numpy.ma.where(Data[li]>0)
        ni = numpy.ma.where(Data[li]<0)
        PosNum = numpy.size(pi)
        NegNum = numpy.size(ni)
        if PosNum > 0:
            SumPos = numpy.ma.sum(Data[pi])
        else:
            SumPos = 0
        if NegNum > 0:
            SumNeg = numpy.ma.sum(Data[ni])
        else:
            SumNeg = 0
    else:
        pi = numpy.ma.where(Data[li]>0)
        ni = numpy.ma.where(Data[li]<0)
        PosNum = numpy.size(pi)
        NegNum = numpy.size(ni)
        SumPos = c.missing_value
        SumNeg = c.missing_value
    return PosNum, NegNum, SumPos, SumNeg

def get_sums(Data):
    """
        Get daily sums when no 30-min observations are missing.
        Days with missing observations return a value of c.missing_value
        Values returned are sample size (Num) and sum (Sum)
        
        Usage qcts.get_sums(Data)
        Data: 1-day dataset
        """
    li = numpy.ma.where(abs(Data-numpy.float64(c.missing_value))>c.eps)
    Num = numpy.size(li)
    if Num == 0:
        Sum = c.missing_value
    elif Num == 48:
        Sum = numpy.ma.sum(Data[li])
    else:
        x = 0
        index = numpy.ma.where(Data.mask == True)[0]
        if len(index) == 1:
            x = 1
        elif len(index) > 1:
            for i in range(len(Data)):
                if Data.mask[i] == True:
                    x = x + 1
        
        if x == 0:
            Sum = numpy.ma.sum(Data[li])
        else:
            Sum = c.missing_value
    return Num, Sum

def get_qcflag(ds):
    """
        Set up flags during ingest of L1 data.
        Identifies missing observations as c.missing_value and sets flag value 1
        
        Usage qcts.get_qcflag(ds)
        ds: data structure
        """
    log.info(' Setting up the QC flags')
    nRecs = len(ds.series['xlDateTime']['Data'])
    for ThisOne in ds.series.keys():
        ds.series[ThisOne]['Flag'] = numpy.zeros(nRecs,dtype=numpy.int32)
        index = numpy.where(ds.series[ThisOne]['Data']==c.missing_value)[0]
        ds.series[ThisOne]['Flag'][index] = numpy.int32(1)

def get_synthetic_fsd(ds):
    if "DateTime_UTC" not in ds.series.keys(): return
    log.info(' Calculating synthetic Fsd')
    lat = numpy.float64(ds.globalattributes["latitude"])
    lon = numpy.float64(ds.globalattributes["longitude"])
    ldt_UTC = ds.series["DateTime_UTC"]["Data"]
    alt_solar = [pysolar.GetAltitude(lat,lon,dt) for dt in ldt_UTC]
    Fsd_syn = [pysolar.GetRadiationDirect(dt,alt) for dt,alt in zip(ldt_UTC,alt_solar)]
    Fsd_syn = numpy.ma.array(Fsd_syn)
    nRecs = len(Fsd_syn)
    flag = numpy.zeros(nRecs,dtype=numpy.int32)
    attr = qcutils.MakeAttributeDictionary(long_name='Synthetic downwelling shortwave radiation',units='W/m2',standard_name='surface_downwelling_shortwave_flux_in_air')
    qcutils.CreateSeries(ds,"Fsd_syn",Fsd_syn,Flag=flag,Attr=attr)

def InvertSign(ds,ThisOne):
    log.info(' Inverting sign of '+ThisOne)
    index = numpy.where(abs(ds.series[ThisOne]['Data']-numpy.float64(c.missing_value))>c.eps)[0]
    ds.series[ThisOne]['Data'][index] = numpy.float64(-1)*ds.series[ThisOne]['Data'][index]

def InterpolateOverMissing(cf,ds,series='',maxlen=1000):
    # check that series is in the data structure
    if series not in ds.series.keys():
        log.error('  InterpolateOverMissing: series '+str(series)+' not found in data structure')
        return
    # convert the Python datetime to a number
    DateNum = date2num(ds.series['DateTime']['Data'])
    # get the data
    data_org,flag_org,attr_org = qcutils.GetSeries(ds,series)
    # number of records
    nRecs = len(data_org)
    # index of good values
    iog = numpy.where(abs(data_org-numpy.float64(c.missing_value))>c.eps)[0]
    # index of missing values
    iom = numpy.where(abs(data_org-numpy.float64(c.missing_value))<=c.eps)[0]
    # return if there is not enough data to use
    if len(iog)<2:
        log.warn('  Interpolate over missing:  '+str(series)+' has insufficient good observations; gap-filling skipped')
        return
    if len(iom)<1:
        log.info('  Interpolate over missing:  '+str(series)+' has no gaps; gap-filling skipped')
        return
    # linear interpolation function
    f = interpolate.interp1d(DateNum[iog],data_org[iog],bounds_error=False,fill_value=numpy.float64(c.missing_value))
    # interpolate over the whole time series
    data_int = f(DateNum).astype(numpy.float64)
    flag_int = numpy.ones(nRecs,dtype=numpy.int32)*60
    # restore the original good data
    data_int[iog] = data_org[iog]
    flag_int[iog] = flag_org[iog]
    # now replace data in contiguous blocks of length > min with missing data
    # first, a conditional index, 0 where data is good, 1 where it is missing
    cond_ind = numpy.zeros(nRecs,dtype=numpy.int32)
    cond_ind[iom] = 1
    cond_bool = (cond_ind==1)
    # start and stop indices of contiguous blocks
    for start, stop in qcutils.contiguous_regions(cond_bool):
        # code to handle minimum segment length goes here
        duration = stop - start
        if duration>maxlen:
            #data_int[start:stop+1] = numpy.float(c.missing_value)
            #flag_int[start:stop+1] = flag_org[start:stop+1]
            data_int[start:stop] = numpy.float64(c.missing_value)
            flag_int[start:stop] = flag_org[start:stop]
    # put data_int back into the data structure
    attr_int = dict(attr_org)
    attr_int['long_name'] = ''
    qcutils.CreateSeries(ds,series,data_int,FList=[series],Attr=attr_int)
    log.info('  Interpolate over missing:  '+str(series)+' filled')

def MassmanStandard(cf,ds,Ta_in='Ta',Ah_in='Ah',ps_in='ps',ustar_in='ustar',ustar_out='ustar',L_in='L',L_out ='L',uw_out='uw',vw_out='vw',wT_out='wT',wA_out='wA',wC_out='wC',u_in='u',uw_in='uw',vw_in='vw',wT_in='wT',wC_in='wC',wA_in='wA'):
    """
       Massman corrections.
       The steps involved are as follows:
        1) calculate ustar and L using rotated but otherwise uncorrected covariances
       """
    if 'Massman' not in cf:
        log.info(' Massman section not found in control file, no corrections applied')
        return
    if qcutils.cfkeycheck(cf,Base='FunctionArgs',ThisOne='MassmanVars'):
        MArgs = ast.literal_eval(cf['FunctionArgs']['MassmanVars'])
        Ta_in = MArgs[0]
        Ah_in = MArgs[1]
        ps_in = MArgs[2]
        ustar_in = MArgs[3]
        L_in = MArgs[4]
    if qcutils.cfkeycheck(cf,Base='FunctionArgs',ThisOne='MassmanOuts'):
        MOut = ast.literal_eval(cf['FunctionArgs']['MassmanOuts'])
        ustar_out = MOut[0]
        L_out = MOut[1]
        uw_out = MOut[2]
        vw_out = MOut[3]
        wT_out = MOut[4]
        wA_out = MOut[5]
        wC_out = MOut[6]
    log.info(' Correcting for flux loss from spectral attenuation')
    zmd = numpy.float64(cf['Massman']['zmd'])             # z-d for site
    angle = numpy.float64(cf['Massman']['angle'])         # CSAT3-IRGA separation angle
    CSATarm = numpy.float64(cf['Massman']['CSATarm'])     # CSAT3 mounting distance
    IRGAarm = numpy.float64(cf['Massman']['IRGAarm'])     # IRGA mounting distance
    lLat = numpy.ma.sin(numpy.deg2rad(angle)) * IRGAarm
    lLong = CSATarm - (numpy.ma.cos(numpy.deg2rad(angle)) * IRGAarm)
    # *** Massman_1stpass starts here ***
    #  The code for the first and second passes is very similar.  It would be useful to make them the
    #  same and put into a loop to reduce the nu,ber of lines in this function.
    # calculate ustar and Monin-Obukhov length from rotated but otherwise uncorrected covariances
    Ta,f,a = qcutils.GetSeriesasMA(ds,Ta_in)
    Ah,f,a = qcutils.GetSeriesasMA(ds,Ah_in)
    ps,f,a = qcutils.GetSeriesasMA(ds,ps_in)
    nRecs = numpy.size(Ta)
    u,f,a = qcutils.GetSeriesasMA(ds,u_in)
    uw,f,a = qcutils.GetSeriesasMA(ds,uw_in)
    vw,f,a = qcutils.GetSeriesasMA(ds,vw_in)
    wT,f,a = qcutils.GetSeriesasMA(ds,wT_in)
    wC,f,a = qcutils.GetSeriesasMA(ds,wC_in)
    wA,f,a = qcutils.GetSeriesasMA(ds,wA_in)
    if ustar_in not in ds.series.keys():
        ustarm = numpy.ma.sqrt(numpy.ma.sqrt(uw ** 2 + vw ** 2))
        attr = qcutils.MakeAttributeDictionary(long_name='Friction velocity: Raw uncorrected',units='m/s')
        qcutils.CreateSeries(ds,ustar_in,ustarm,FList=[uw_in,vw_in],Attr=attr)
        ustarm,f,a = qcutils.GetSeriesasMA(ds,ustar_in)
        ustar_attr = ''
    else:
        ustarm,f,a = qcutils.GetSeriesasMA(ds,ustar_in)
        ustar_attr = ''
    if L_in not in ds.series.keys():
        Lm = mf.molen(Ta, Ah, ps, ustarm, wT, fluxtype='kinematic')
        attr = qcutils.MakeAttributeDictionary(long_name='Obukhov Length: Raw uncorrected',units='m')
        qcutils.CreateSeries(ds,L_in,Lm,FList=[Ta_in,Ah_in,ps_in,ustar_in,wT_in],Attr=attr)
        Lm,Lflag,a = qcutils.GetSeriesasMA(ds,L_in)
        L_attr = ''
    else:
        Lm,Lflag,a = qcutils.GetSeriesasMA(ds,L_in)
        L_attr = ''
    # now calculate z on L
    zetaindex = numpy.where(Lm != 0)[0]
    zetaflagindex = numpy.where(Lm == 0)[0]
    zoLm = numpy.ones(len(Lm),dtype=numpy.float64)*c.missing_value
    zoLm[zetaindex] = zmd / Lm[zetaindex]
    zetaflag = numpy.zeros(len(Lm),dtype=numpy.int32) + Lflag
    newindex = numpy.where((zoLm==c.missing_value)&(numpy.mod(zetaflag,10)==0))[0]
    zetaflag[newindex] = numpy.int32(22)
    if 'zeta' not in ds.series.keys():
        attr = qcutils.MakeAttributeDictionary(long_name='stability coefficient, (z-d)/L',units='unitless')
        qcutils.CreateSeries(ds,'zeta',zoLm,Flag=zetaflag,Attr=attr)
    else:
        attr = qcutils.MakeAttributeDictionary(long_name='',units='unitless')
        qcutils.CreateSeries(ds,'zeta',zoLm,Flag=zetaflag,Attr=attr)
    # start calculating the correction coefficients for approximate corrections
    #  create nxMom, nxScalar and alpha series with their unstable values by default
    nxMom, nxScalar, alpha = qcutils.nxMom_nxScalar_alpha(zoLm)
    # now calculate the fxMom and fxScalar coefficients
    fxMom = nxMom * u / zmd
    fxScalar = nxScalar * u / zmd
    # compute spectral filters
    tao_eMom = numpy.ma.sqrt(((c.lwVert / (5.7 * u)) ** 2) + ((c.lwHor / (2.8 * u)) ** 2))
    tao_ewT = numpy.ma.sqrt(((c.lwVert / (8.4 * u)) ** 2) + ((c.lTv / (4.0 * u)) ** 2))
    tao_ewIRGA = numpy.ma.sqrt(((c.lwVert / (8.4 * u)) ** 2) + ((c.lIRGA / (4.0 * u)) ** 2) + ((lLat / (1.1 * u)) ** 2) + ((lLong / (1.05 * u)) ** 2))
    tao_b = c.Tb / 2.8
    # calculate coefficients
    bMom = qcutils.bp(fxMom,tao_b)
    bScalar = qcutils.bp(fxScalar,tao_b)
    pMom = qcutils.bp(fxMom,tao_eMom)
    pwT = qcutils.bp(fxScalar,tao_ewT)
    # calculate corrections for momentum and scalars
    rMom = qcutils.r(bMom, pMom, alpha)        # I suspect that rMom and rwT are the same functions
    rwT = qcutils.r(bScalar, pwT, alpha)
    # determine approximately-true Massman fluxes
    uwm = uw / rMom
    vwm = vw / rMom
    wTm = wT / rwT
    # *** Massman_1stpass ends here ***
    # *** Massman_2ndpass starts here ***
    # we have calculated the first pass corrected momentum and temperature covariances, now we use
    # these to calculate the final corrections
    #  first, get the 2nd pass corrected friction velocity and Monin-Obukhov length
    ustarm = numpy.ma.sqrt(numpy.ma.sqrt(uwm ** 2 + vwm ** 2))
    attr = qcutils.MakeAttributeDictionary(long_name='',units='m/s')
    qcutils.CreateSeries(ds,ustar_in,ustarm,FList=[uw_in,vw_in,ustar_in],Attr=attr)
    ustarm,f,a = qcutils.GetSeriesasMA(ds,ustar_in)
    Lm = mf.molen(Ta, Ah, ps, ustarm, wTm, fluxtype='kinematic')
    attr = qcutils.MakeAttributeDictionary(long_name='',units='m')
    qcutils.CreateSeries(ds,L_in,Lm,FList=[Ta_in,Ah_in,ps_in,ustar_in,wT_in],Attr=attr)
    Lm,Lflag,a = qcutils.GetSeriesasMA(ds,L_in)
    nxMom, nxScalar, alpha = qcutils.nxMom_nxScalar_alpha(zoLm)
    fxMom = nxMom * (u / zmd)
    fxScalar = nxScalar * (u / zmd)
    # calculate coefficients
    bMom = qcutils.bp(fxMom,tao_b)
    bScalar = qcutils.bp(fxScalar,tao_b)
    pMom = qcutils.bp(fxMom,tao_eMom)
    pwT = qcutils.bp(fxScalar,tao_ewT)
    pwIRGA = qcutils.bp(fxScalar,tao_ewIRGA)
    # calculate corrections for momentum and scalars
    rMom = qcutils.r(bMom, pMom, alpha)
    rwT = qcutils.r(bScalar, pwT, alpha)
    rwIRGA = qcutils.r(bScalar, pwIRGA, alpha)
    # determine true fluxes
    uwM = uw / rMom
    vwM = vw / rMom
    wTM = wT / rwT
    wCM = wC / rwIRGA
    wAM = wA / rwIRGA
    ustarM = numpy.ma.sqrt(numpy.ma.sqrt(uwM ** 2 + vwM ** 2))
    attr = qcutils.MakeAttributeDictionary(long_name='',units='m/s')
    qcutils.CreateSeries(ds,ustar_in,ustarM,FList=[uw_in,vw_in,ustar_in],Attr=attr)
    ustarM,f,a = qcutils.GetSeriesasMA(ds,ustar_in)
    LM = mf.molen(Ta, Ah, ps, ustarM, wTM, fluxtype='kinematic')
    attr = qcutils.MakeAttributeDictionary(long_name='',units='m')
    qcutils.CreateSeries(ds,L_in,LM,FList=[Ta_in,Ah_in,ps_in,ustar_in,wT_in],Attr=attr)
    LM,Lflag,a = qcutils.GetSeriesasMA(ds,L_in)
    zetaindex = numpy.where(LM != 0)[0]
    zetaflagindex = numpy.where(LM == 0)[0]
    zeta = numpy.ones(len(LM),dtype=numpy.float64)*c.missing_value
    zeta[zetaindex] = zmd / LM[zetaindex]
    zetaflag = numpy.zeros(len(LM),dtype=numpy.int32) + Lflag
    newindex = numpy.where((zeta==c.missing_value)&(numpy.mod(zetaflag,10)==0))[0]
    zetaflag[newindex] = numpy.int32(22)
    attr = qcutils.MakeAttributeDictionary(long_name='',units='unitless')
    qcutils.CreateSeries(ds,'zeta',zeta,Flag=zetaflag,Attr=attr)
    # write the 2nd pass Massman corrected covariances to the data structure
    attr = qcutils.MakeAttributeDictionary(long_name=ustar_attr,units='m/s')
    qcutils.CreateSeries(ds,ustar_out,ustarM,FList=[Ta_in,Ah_in,ps_in,u_in,uw_in,vw_in,wT_in,wC_in,wA_in],Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name=L_attr+', rotated to natural wind coordinates, frequency response corrected',units='m')
    qcutils.CreateSeries(ds,L_out,LM,FList=[Ta_in,Ah_in,ps_in,u_in,uw_in,vw_in,wT_in,wC_in,wA_in],Attr=attr)
    L_flagindex = numpy.where(numpy.mod(ds.series[L_out]['Flag'],10)!=0)[0]
    ds.series[L_out]['Data'][L_flagindex] = numpy.float64(c.missing_value)
    attr = qcutils.MakeAttributeDictionary(long_name='frequency response corrected',units='m2/s2')
    qcutils.CreateSeries(ds,uw_out,uwM,FList=[Ta_in,Ah_in,ps_in,u_in,uw_in,vw_in,wT_in,wC_in,wA_in],Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='frequency response corrected',units='m2/s2')
    qcutils.CreateSeries(ds,vw_out,vwM,FList=[Ta_in,Ah_in,ps_in,u_in,uw_in,vw_in,wT_in,wC_in,wA_in],Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='frequency response corrected',units='mC/s')
    qcutils.CreateSeries(ds,wT_out,wTM,FList=[Ta_in,Ah_in,ps_in,u_in,uw_in,vw_in,wT_in,wC_in,wA_in],Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='frequency response corrected',units='g/m2/s')
    qcutils.CreateSeries(ds,wA_out,wAM,FList=[Ta_in,Ah_in,ps_in,u_in,uw_in,vw_in,wT_in,wC_in,wA_in],Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='frequency response corrected',units='mg/m2/s')
    qcutils.CreateSeries(ds,wC_out,wCM,FList=[Ta_in,Ah_in,ps_in,u_in,uw_in,vw_in,wT_in,wC_in,wA_in],Attr=attr)
    # *** Massman_2ndpass ends here ***
    
    if qcutils.cfkeycheck(cf,Base='General',ThisOne='MassmanFlag') and cf['General']['MassmanFlag'] == 'True':
        keys = [ustar_out,L_out,uw_out,vw_out,wT_out,wA_out,wC_out]
        for ThisOne in keys:
            testseries,f,a = qcutils.GetSeriesasMA(ds,ThisOne)
            mask = numpy.ma.getmask(testseries)
            index = numpy.where(mask.astype(numpy.int32)==1)
            ds.series[ThisOne]['Flag'][index] = numpy.int32(12)
    else:
        keys = [ustar_out,L_out,uw_out,vw_out,wT_out,wA_out,wC_out]
        for ThisOne in keys:
            testseries,f,a = qcutils.GetSeriesasMA(ds,ThisOne)
            mask = numpy.ma.getmask(testseries)
            index = numpy.where((numpy.mod(f,10)==0) & (mask.astype(numpy.int32)==1))    # find the elements with flag = 0, 10, 20 etc and masked (check for masked data with good data flag)
            ds.series[ThisOne]['Flag'][index] = numpy.int32(12)
            
    index_neutral = numpy.ma.where((zeta > -0.1) & (zeta < 0.1))[0]
    index_slight_stable = numpy.ma.where((zeta > 0.1) & (zeta < 1))[0]
    index_stable = numpy.ma.where((zeta > 1) & (zeta < 2))[0]
    index_very_stable = numpy.ma.where((zeta > 2) & (zeta < 9999999))[0]
    index_slight_unstable = numpy.ma.where((zeta < -0.1) & (zeta > -1))[0]
    index_unstable = numpy.ma.where((zeta < -1) & (zeta > -2))[0]
    index_very_unstable = numpy.ma.where(zeta < -2)[0]
    index_total = len(index_neutral) + len(index_slight_stable) + len(index_stable) + len(index_very_stable) + len(index_slight_unstable) + len(index_unstable) + len(index_very_unstable)
    log.info('   Total measurements:  '+str(index_total))
    log.info('   '+str(len(index_very_unstable))+':  very unstable, zeta < -2')
    log.info('   '+str(len(index_unstable))+':  unstable, -1 > zeta > -2')
    log.info('   '+str(len(index_slight_unstable))+':  slightly unstable, -0.1 > zeta > -1')
    log.info('   '+str(len(index_neutral))+':  neutral, -0.1 < zeta < 0.1')
    log.info('   '+str(len(index_slight_stable))+':  slightly stable, 0.1 < zeta < 1')
    log.info('   '+str(len(index_stable))+':  stable, 1 < zeta < 2')
    log.info('   '+str(len(index_very_stable))+':  very stable, zeta > 2')

def MergeSeries(cf,ds,series,okflags):
    """
        Merge two series of data to produce one series containing the best data from both.
        Calling syntax is: MergeSeries(cf,ds,series,okflags)
         where ds is the data structure containing all series
               series (str) is the label of the destination series
               okflags (list) is a list of QC flag values for which the data is considered acceptable
        If the QC flag for Primary is in okflags, the value from Primary is placed in destination.
        If the QC flag for Primary is not in okflags but the QC flag for Secondary is, the value
        from Secondary is placed in Destination.
        """
    # check to see if the series is specified in the control file
    section = qcutils.get_cfsection(cf,series=series)
    if len(section)==0: return
    # check to see if the entry for series in the control file has the MergeSeries key
    if 'MergeSeries' not in cf[section][series].keys(): return
    # check to see if the series has already been merged
    if series in ds.mergeserieslist: return
    # now get the source list and the standard name
    srclist, standardname = qcutils.GetMergeSeriesKeys(cf,series,section=section)
    nSeries = len(srclist)
    if nSeries==0:
        log.info(' MergeSeries: no input series specified for '+str(series))
        return
    if nSeries==1:
        log.info(' Merging series '+str(srclist)+' into '+series)
        if srclist[0] not in ds.series.keys():
            log.error('  MergeSeries: primary input series'+srclist[0]+'not found for'+str(series))
            return
        data = ds.series[srclist[0]]['Data'].copy()
        flag = ds.series[srclist[0]]['Flag'].copy()
        attr = ds.series[srclist[0]]['Attr'].copy()
        SeriesNameString = srclist[0]
    else:
        log.info(' Merging series '+str(srclist)+' into '+series)
        if srclist[0] not in ds.series.keys():
            log.error('  MergeSeries: primary input series'+srclist[0]+'not found')
            return
        data = ds.series[srclist[0]]['Data'].copy()
        flag = ds.series[srclist[0]]['Flag'].copy()
        attr = ds.series[srclist[0]]['Attr'].copy()
        SeriesNameString = srclist[0]
        srclist.remove(srclist[0])
        for ThisOne in srclist:
            if ThisOne in ds.series.keys():
                SeriesNameString = SeriesNameString+', '+ThisOne
                indx1 = numpy.zeros(numpy.size(data),dtype=numpy.int32)
                indx2 = numpy.zeros(numpy.size(data),dtype=numpy.int32)
                for okflag in okflags:
                    index = numpy.where((flag==okflag))[0]                             # index of acceptable primary values
                    indx1[index] = 1                                                   # set primary index to 1 when primary good
                    index = numpy.where((ds.series[ThisOne]['Flag']==okflag))[0]       # same process for secondary
                    indx2[index] = 1
                index = numpy.where((indx1!=1)&(indx2==1))[0]           # index where primary bad but secondary good
                data[index] = ds.series[ThisOne]['Data'][index]         # replace bad primary with good secondary
                flag[index] = ds.series[ThisOne]['Flag'][index]
            else:
                log.error('  MergeSeries: secondary input series'+ThisOne+'not found')
    ds.mergeserieslist.append(series)
    #attr = qcutils.MakeAttributeDictionary(long_name='Merged from '+SeriesNameString,
                             #standard_name=standardname,units=SeriesUnitString)
    attr["long_name"] = attr["long_name"]+", merged from " + SeriesNameString
    qcutils.CreateSeries(ds,series,data,Flag=flag,Attr=attr)

def prep_aerodynamicresistance(cf,ds,Cdmethod,Cemethod,Ce_2layer):
    if qcutils.cfkeycheck(cf,Base='PenmanMonteith',ThisOne='PMin'):
        PMin = ast.literal_eval(cf['PenmanMonteith']['PMin'])
    else:
        PMin = ['Fe', 'Ta', 'Ah', 'ps', 'Ws', 'Fn', 'Fsd', 'Fg', 'VPD', 'wA']     # ***
    
    Level = ds.globalattributes['nc_level']
    log.info(' Computing Penman-Monteith bulk stomatal resistance at level '+Level)
    if qcutils.cfkeycheck(cf,Base='PenmanMonteith',ThisOne='PMcritFsd'):
        critFsd = numpy.float64(cf['PenmanMonteith']['PMcritFsd'])
    else:
        critFsd = numpy.floag64(10)
    
    if qcutils.cfkeycheck(cf,Base='PenmanMonteith',ThisOne='PMcritFe'):
        critFe = numpy.float64(cf['PenmanMonteith']['PMcritFe'])
    else:
        critFe = numpy.float64(0)
    
    Fe,f,a = qcutils.GetSeriesasMA(ds,PMin[0])
    Fsd,f,a = qcutils.GetSeriesasMA(ds,PMin[6])
    Uavg,f,a = qcutils.GetSeriesasMA(ds,PMin[4])
    uindex = numpy.ma.where(Uavg == 0)[0]
    Feindex = numpy.ma.where(Fe < critFe)[0]
    Fsdindex = numpy.ma.where(Fsd < critFsd)[0]
    Uavg[uindex] = 0.000000000000001
    # use bulk transfer coefficient method
    if Cemethod == 'True':
        # rav parameterised with q profile measurements
        if qcutils.cfkeycheck(cf,Base='PenmanMonteith',ThisOne='profileq'):
            qList = ast.literal_eval(cf['PenmanMonteith']['profileq'])
        else:
            log.error('  PenmanMonteith:  profileq not given')
            return
        
        outList = ['Ce_1layer','rav_1layer','rSv_1layer','GSv_1layer']
        attribute = 'Ce with q profile, '
        log.info('  Ce method (Brutseart 1982, Stull 1988) used to estimate aerodynamic resistance, rav')
        log.info('   Ce method: qList: '+str(qList))
        Ce, rav, ravindex = get_rav(cf,ds,Uavg,PMin,qList)
        rst, rstindex, Gst, Gstindex = get_rstGst(cf,ds,PMin,rav)
        flagList = [PMin[0],PMin[1],PMin[2],PMin[3],PMin[4],PMin[5],PMin[6],PMin[7],qList[0],qList[1]]
        if outList[0] in ds.series.keys():
            attr = qcutils.MakeAttributeDictionary(long_name=Level,units='s/m')
        else:
            attr = qcutils.MakeAttributeDictionary(long_name='Bulk transfer coefficient, '+attribute+Level,units='s/m')
        qcutils.CreateSeries(ds,outList[0],Ce,FList=flagList,Attr=attr)
        if outList[1] in ds.series.keys():
            attr = qcutils.MakeAttributeDictionary(long_name=Level,units='s/m')
        else:
            attr = qcutils.MakeAttributeDictionary(long_name='Aerodynamic resistance, '+attribute+Level,units='s/m')
        qcutils.CreateSeries(ds,outList[1],rav,FList=flagList,Attr=attr)
        if outList[2] in ds.series.keys():
            attr = qcutils.MakeAttributeDictionary(long_name=Level,units='s/m')
        else:
            attr = qcutils.MakeAttributeDictionary(long_name='Stomatal resistance from Penman-Monteith inversion, '+attribute+Level,units='s/m')
        qcutils.CreateSeries(ds,outList[2],rst,FList=flagList,Attr=attr)
        if outList[3] in ds.series.keys():
            attr = qcutils.MakeAttributeDictionary(long_name=Level,units='mmolH2O/(m2ground s)')
        else:
            attr = qcutils.MakeAttributeDictionary(long_name='Conductance from Penman-Monteith inversion, '+attribute+Level,units='mmolH2O/(m2ground s)')
        qcutils.CreateSeries(ds,outList[3],Gst,FList=flagList,Attr=attr)
        for listindex in range(0,4):
            ds.series[outList[listindex]]['Attr']['InputSeries'] = PMin
            ds.series[outList[listindex]]['Attr']['FsdCutoff'] = critFsd
            ds.series[outList[listindex]]['Attr']['FeCutoff'] = critFe
            ds.series[outList[listindex]]['Flag'][rstindex] = numpy.int32(171)
            ds.series[outList[listindex]]['Flag'][ravindex] = numpy.int32(171)
            ds.series[outList[listindex]]['Flag'][uindex] = numpy.int32(174)
            ds.series[outList[listindex]]['Flag'][Feindex] = numpy.int32(172)
            ds.series[outList[listindex]]['Flag'][Fsdindex] = numpy.int32(173)
            ds.series[outList[listindex]]['Data'][rstindex] = numpy.float64(c.missing_value)
            ds.series[outList[listindex]]['Data'][ravindex] = numpy.float64(c.missing_value)
            ds.series[outList[listindex]]['Data'][uindex] = numpy.float64(c.missing_value)
            ds.series[outList[listindex]]['Data'][Feindex] = numpy.float64(c.missing_value)
            ds.series[outList[listindex]]['Data'][Fsdindex] = numpy.float64(c.missing_value)
        
        if qcutils.cfkeycheck(cf,Base='PenmanMonteith',ThisOne='rlvCe') and cf['PenmanMonteith']['rlvCe'] == 'True':
            rlv, flag = get_leafresistance(cf,ds,outList[2])
            if 'rLv_1layer' in ds.series.keys():
                attr = qcutils.MakeAttributeDictionary(long_name='',units='s/m')
            else:
                attr = qcutils.MakeAttributeDictionary(long_name='leaf resistance from Penman-Monteith inversion, Ce-method, under well-illuminated (> 600 W m-2 Fsd) conditions',units='s/m')
            qcutils.CreateSeries(ds,'rLv_1layer',rlv,Flag=0,Attr=attr)
            ds.series['rLv_1layer']['Flag'] = flag
    
    # use 2-layer bulk transfer coefficient method
    if Ce_2layer == 'True':
        # rav parameterised with q profile measurements
        if qcutils.cfkeycheck(cf,Base='PenmanMonteith',ThisOne='base_profileq'):
            base_qList = ast.literal_eval(cf['PenmanMonteith']['base_profileq'])
        else:
            log.error('  PenmanMonteith:  base_profileq not given')
            return
        
        if qcutils.cfkeycheck(cf,Base='PenmanMonteith',ThisOne='top_profileq'):
            top_qList = ast.literal_eval(cf['PenmanMonteith']['top_profileq'])
        else:
            log.error('  PenmanMonteith:  top_profileq not given')
            return
        
        outList = ['Ce_base','rav_base','rSv_base','GSv_base','Ce_top','rav_top','rSv_top','GSv_top','Ce_full','rav_full','rSv_full','GSv_full','rSv_2layer','GSv_2layer','rav_2layer']
        attribute = 'Ce 2-layer with q profile, '
        log.info('  Ce method (Brutseart 1982, Stull 1988) used to estimate aerodynamic resistance, rav')
        log.info('   2-layer Ce method: lower layer qList: '+str(base_qList))
        log.info('   2-layer Ce method: upper layer qList: '+str(top_qList))
        log.info('   2-layer Ce method: full layer qList: '+str(top_qList[0])+', '+str(base_qList[1]))
        flagListTop = [PMin[0],PMin[1],PMin[2],PMin[3],PMin[4],PMin[5],PMin[6],PMin[7],top_qList[0],top_qList[1]]
        flagListBase = [PMin[0],PMin[1],PMin[2],PMin[3],PMin[4],PMin[5],PMin[6],PMin[7],base_qList[0],base_qList[1]]
        flagListFull = [PMin[0],PMin[1],PMin[2],PMin[3],PMin[4],PMin[5],PMin[6],PMin[7],top_qList[0],base_qList[1]]
        flagList = [PMin[0],PMin[1],PMin[2],PMin[3],PMin[4],PMin[5],PMin[6],PMin[7],base_qList[0],base_qList[1],top_qList[0],top_qList[1]]
        Ce_top, rav_top, ravindex_top = get_rav(cf,ds,Uavg,PMin,top_qList,layer='top_',method='Ce_2layer')
        Ce_base, rav_base, ravindex_base = get_rav(cf,ds,Uavg,PMin,base_qList,layer='base_',method='Ce_2layer')
        Ce_full, rav_full, ravindex_full = get_rav(cf,ds,Uavg,PMin,[top_qList[0],base_qList[1]],method='Ce_2layer')
        rav_diff1 = numpy.setdiff1d(ravindex_base,ravindex_top)
        rav_diff2 = numpy.setdiff1d(ravindex_top,ravindex_base)
        rav_diff3 = numpy.concatenate((rav_diff2,rav_diff1), axis=0)
        ravindex_neither = numpy.setdiff1d(ravindex_top,rav_diff3)
        
        #determine rav when rav = rav_base (rav_top <= 0)
        rav = (rav_full * rav_top) / (rav_full + rav_top)
        
        #determine rav when rav = rav_top (rav_base <= 0)
        rav[ravindex_base] = rav_top[ravindex_base]
        
        #determine rav when rav = rav_base (rav_top <= 0)
        rav[ravindex_top] = rav_base[ravindex_top]
        
        ravindex = numpy.where(rav < 0)[0]
        ravgoodindex_top = numpy.where(rav_top > 0)[0]
        ravgoodindex_base = numpy.where(rav_base > 0)[0]
        ravgoodindex_full = numpy.where(rav_full > 0)[0]
        ravgoodindex = numpy.where(rav > 0)[0]
        rst_top, rstindex_top, Gst_top, Gstindex_top = get_rstGst(cf,ds,PMin,rav_top)
        rst_base, rstindex_base, Gst_base, Gstindex_base = get_rstGst(cf,ds,PMin,rav_base)
        rst_full, rstindex_full, Gst_full, Gstindex_full = get_rstGst(cf,ds,PMin,rav_full)
        rst, rstindex, Gst, Gstindex = get_rstGst(cf,ds,PMin,rav)
        if outList[0] in ds.series.keys():
            attr = qcutils.MakeAttributeDictionary(long_name=Level,units='s/m')
        else:
            attr = qcutils.MakeAttributeDictionary(long_name='Bulk transfer coefficient, '+attribute+'base layer, '+Level,units='s/m')
        qcutils.CreateSeries(ds,outList[0],Ce_base,FList=flagListBase,Attr=attr)
        if outList[1] in ds.series.keys():
            attr = qcutils.MakeAttributeDictionary(long_name=Level,units='s/m')
        else:
            attr = qcutils.MakeAttributeDictionary(long_name='Aerodynamic resistance, '+attribute+'base layer, '+Level,units='s/m')
        qcutils.CreateSeries(ds,outList[1],rav_base,FList=flagListBase,Attr=attr)
        if outList[2] in ds.series.keys():
            attr = qcutils.MakeAttributeDictionary(long_name=Level,units='s/m')
        else:
            attr = qcutils.MakeAttributeDictionary(long_name='Stomatal resistance from Penman-Monteith inversion, '+attribute+'base layer, '+Level,units='s/m')
        qcutils.CreateSeries(ds,outList[2],rst_base,FList=flagListBase,Attr=attr)
        if outList[3] in ds.series.keys():
            attr = qcutils.MakeAttributeDictionary(long_name=Level,units='mmolH2O/(m2ground s)')
        else:
            attr = qcutils.MakeAttributeDictionary(long_name='Conductance from Penman-Monteith inversion, '+attribute+'base layer, '+Level,units='mmolH2O/(m2ground s)')
        qcutils.CreateSeries(ds,outList[3],Gst_base,FList=flagListBase,Attr=attr)
        if outList[4] in ds.series.keys():
            attr = qcutils.MakeAttributeDictionary(long_name=Level,units='s/m')
        else:
            attr = qcutils.MakeAttributeDictionary(long_name='Bulk transfer coefficient, '+attribute+'top layer, '+Level,units='s/m')
        qcutils.CreateSeries(ds,outList[4],Ce_top,FList=flagListTop,Attr=attr)
        if outList[5] in ds.series.keys():
            attr = qcutils.MakeAttributeDictionary(long_name=Level,units='s/m')
        else:
            attr = qcutils.MakeAttributeDictionary(long_name='Aerodynamic resistance, '+attribute+'top layer, '+Level,units='s/m')
        qcutils.CreateSeries(ds,outList[5],rav_top,FList=flagListTop,Attr=attr)
        if outList[6] in ds.series.keys():
            attr = qcutils.MakeAttributeDictionary(long_name=Level,units='s/m')
        else:
            attr = qcutils.MakeAttributeDictionary(long_name='Stomatal resistance from Penman-Monteith inversion, '+attribute+'top layer, '+Level,units='s/m')
        qcutils.CreateSeries(ds,outList[6],rst_top,FList=flagListTop,Attr=attr)
        if outList[7] in ds.series.keys():
            attr = qcutils.MakeAttributeDictionary(long_name=Level,units='mmolH2O/(m2ground s)')
        else:
            attr = qcutils.MakeAttributeDictionary(long_name='Conductance from Penman-Monteith inversion, '+attribute+'top layer, '+Level,units='mmolH2O/(m2ground s)')
        qcutils.CreateSeries(ds,outList[7],Gst_top,FList=flagListTop,Attr=attr)
        if outList[8] in ds.series.keys():
            attr = qcutils.MakeAttributeDictionary(long_name=Level,units='s/m')
        else:
            attr = qcutils.MakeAttributeDictionary(long_name='Bulk transfer coefficient, '+attribute+'full layer, '+Level,units='s/m')
        qcutils.CreateSeries(ds,outList[8],Ce_full,FList=flagListFull,Attr=attr)
        if outList[9] in ds.series.keys():
            attr = qcutils.MakeAttributeDictionary(long_name=Level,units='s/m')
        else:
            attr = qcutils.MakeAttributeDictionary(long_name='Aerodynamic resistance, '+attribute+'full layer, '+Level,units='s/m')
        qcutils.CreateSeries(ds,outList[9],rav_full,FList=flagListFull,Attr=attr)
        if outList[10] in ds.series.keys():
            attr = qcutils.MakeAttributeDictionary(long_name=Level,units='s/m')
        else:
            attr = qcutils.MakeAttributeDictionary(long_name='Stomatal resistance from Penman-Monteith inversion, '+attribute+'full layer, '+Level,units='s/m')
        qcutils.CreateSeries(ds,outList[10],rst_full,FList=flagListFull,Attr=attr)
        if outList[11] in ds.series.keys():
            attr = qcutils.MakeAttributeDictionary(long_name=Level,units='mmolH2O/(m2ground s)')
        else:
            attr = qcutils.MakeAttributeDictionary(long_name='Conductance from Penman-Monteith inversion, '+attribute+'full layer, '+Level,units='mmolH2O/(m2ground s)')
        qcutils.CreateSeries(ds,outList[11],Gst_full,FList=flagListFull,Attr=attr)
        if outList[12] in ds.series.keys():
            attr = qcutils.MakeAttributeDictionary(long_name=Level,units='s/m')
        else:
            attr = qcutils.MakeAttributeDictionary(long_name='Stomatal resistance from Penman-Monteith inversion, '+attribute+Level,units='s/m')
        qcutils.CreateSeries(ds,outList[12],rst,FList=flagList,Attr=attr)
        if outList[13] in ds.series.keys():
            attr = qcutils.MakeAttributeDictionary(long_name=Level,units='mmolH2O/(m2ground s)')
        else:
            attr = qcutils.MakeAttributeDictionary(long_name='Conductance from Penman-Monteith inversion, '+attribute+Level,units='mmolH2O/(m2ground s)')
        qcutils.CreateSeries(ds,outList[13],Gst,FList=flagList,Attr=attr)
        if outList[14] in ds.series.keys():
            attr = qcutils.MakeAttributeDictionary(long_name=Level,units='s/m')
        else:
            attr = qcutils.MakeAttributeDictionary(long_name='Aerodynamic resistance, '+attribute+Level,units='s/m')
        qcutils.CreateSeries(ds,outList[14],rav,FList=flagList,Attr=attr)
        for listindex in range(0,15):
            ds.series[outList[listindex]]['Attr']['InputSeries'] = PMin
            ds.series[outList[listindex]]['Attr']['FsdCutoff'] = critFsd
            ds.series[outList[listindex]]['Attr']['FeCutoff'] = critFe
            goodflagindex = numpy.where(numpy.mod(ds.series[outList[listindex]]['Flag'],10)==0)[0]
            badflagindex = numpy.where(numpy.mod(ds.series[outList[listindex]]['Flag'],10)!=0)[0]
            if '_base' in outList[listindex]:
                goodflagonlyindex = numpy.setdiff1d(ravindex_top,goodflagindex)
                goodravonlyindex = numpy.setdiff1d(goodflagindex,ravindex_top)
                goodnotbothindex = numpy.concatenate((goodflagonlyindex,goodravonlyindex), axis=0)
                goodbothindex = numpy.setdiff1d(ravindex_top,goodnotbothindex)
                ds.series[outList[listindex]]['Flag'][goodbothindex] = numpy.int32(190)
                ds.series[outList[listindex]]['Flag'][rstindex_base] = numpy.int32(171)
                ds.series[outList[listindex]]['Flag'][ravindex_base] = numpy.int32(171)
                goodflagonlyindex = numpy.setdiff1d(ravgoodindex_top,goodflagindex)
                badravonlyindex = numpy.setdiff1d(goodflagindex,ravgoodindex_top)
                badnotbothindex = numpy.concatenate((goodflagonlyindex,badravonlyindex), axis=0)
                badbothindex = numpy.setdiff1d(ravgoodindex_top,badnotbothindex)
                ds.series[outList[listindex]]['Flag'][badbothindex] = numpy.int32(191)
                ds.series[outList[listindex]]['Data'][rstindex_base] = numpy.float64(c.missing_value)
                ds.series[outList[listindex]]['Data'][ravindex_base] = numpy.float64(c.missing_value)
                ds.series[outList[listindex]]['Data'][ravgoodindex_top] = numpy.float64(c.missing_value)
            if '_top' in outList[listindex]:
                goodflagonlyindex140 = numpy.setdiff1d(ravindex_base,goodflagindex)
                goodflagonlyindex160 = numpy.setdiff1d(ravgoodindex_base,goodflagindex)
                goodravonlyindex140 = numpy.setdiff1d(goodflagindex,ravindex_base)
                goodravonlyindex160 = numpy.setdiff1d(goodflagindex,ravgoodindex_base)
                goodnotbothindex140 = numpy.concatenate((goodflagonlyindex140,goodravonlyindex140), axis=0)
                goodnotbothindex160 = numpy.concatenate((goodflagonlyindex160,goodravonlyindex160), axis=0)
                goodbothindex140 = numpy.setdiff1d(ravindex_base,goodnotbothindex140)
                goodbothindex160 = numpy.setdiff1d(ravgoodindex_base,goodnotbothindex160)
                ds.series[outList[listindex]]['Flag'][goodbothindex140] = numpy.int32(180)
                ds.series[outList[listindex]]['Flag'][goodbothindex160] = numpy.int32(200)
                ds.series[outList[listindex]]['Flag'][rstindex_top] = numpy.int32(171)
                ds.series[outList[listindex]]['Flag'][ravindex_top] = numpy.int32(171)
                ds.series[outList[listindex]]['Data'][rstindex_top] = numpy.float64(c.missing_value)
                ds.series[outList[listindex]]['Data'][ravindex_top] = numpy.float64(c.missing_value)
            if '_full' in outList[listindex]:
                goodflagonlyindex = numpy.setdiff1d(ravgoodindex_full,goodflagindex)
                goodravonlyindex = numpy.setdiff1d(goodflagindex,ravgoodindex_full)
                goodnotbothindex = numpy.concatenate((goodflagonlyindex,goodravonlyindex), axis=0)
                goodbothindex = numpy.setdiff1d(ravgoodindex_full,goodnotbothindex)
                ds.series[outList[listindex]]['Flag'][goodbothindex] = numpy.int32(200)
                ds.series[outList[listindex]]['Flag'][rstindex_full] = numpy.int32(171)
                ds.series[outList[listindex]]['Flag'][ravindex_full] = numpy.int32(171)
                goodflagonlyindex_top = numpy.setdiff1d(ravindex_top,goodflagindex)
                goodflagonlyindex_base = numpy.setdiff1d(ravindex_base,goodflagindex)
                badravonlyindex_top = numpy.setdiff1d(goodflagindex,ravindex_top)
                badravonlyindex_base = numpy.setdiff1d(goodflagindex,ravindex_base)
                badnotbothindex_top = numpy.concatenate((goodflagonlyindex_top,badravonlyindex_top), axis=0)
                badnotbothindex_base = numpy.concatenate((goodflagonlyindex_base,badravonlyindex_base), axis=0)
                badbothindex_top = numpy.setdiff1d(ravindex_top,badnotbothindex_top)
                badbothindex_base = numpy.setdiff1d(ravindex_base,badnotbothindex_base)
                ds.series[outList[listindex]]['Flag'][badbothindex_top] = numpy.int32(201)
                ds.series[outList[listindex]]['Flag'][badbothindex_base] = numpy.int32(201)
                ds.series[outList[listindex]]['Data'][rstindex_full] = numpy.float64(c.missing_value)
                ds.series[outList[listindex]]['Data'][ravindex_full] = numpy.float64(c.missing_value)
                ds.series[outList[listindex]]['Data'][ravindex_top] = numpy.float64(c.missing_value)
                ds.series[outList[listindex]]['Data'][ravindex_base] = numpy.float64(c.missing_value)
            if '_2layer' in outList[listindex]:
                goodflagonlyindex140 = numpy.setdiff1d(ravindex_base,goodflagindex)
                goodflagonlyindex150 = numpy.setdiff1d(ravindex_top,goodflagindex)
                goodflagonlyindex160 = numpy.setdiff1d(ravgoodindex_full,goodflagindex)
                goodravonlyindex140 = numpy.setdiff1d(goodflagindex,ravindex_base)
                goodravonlyindex150 = numpy.setdiff1d(goodflagindex,ravindex_top)
                goodravonlyindex160 = numpy.setdiff1d(goodflagindex,ravgoodindex_full)
                goodnotbothindex140 = numpy.concatenate((goodflagonlyindex140,goodravonlyindex140), axis=0)
                goodnotbothindex150 = numpy.concatenate((goodflagonlyindex150,goodravonlyindex150), axis=0)
                goodnotbothindex160 = numpy.concatenate((goodflagonlyindex160,goodravonlyindex160), axis=0)
                goodbothindex140 = numpy.setdiff1d(ravindex_base,goodnotbothindex140)
                goodbothindex150 = numpy.setdiff1d(ravindex_top,goodnotbothindex150)
                goodbothindex160 = numpy.setdiff1d(ravgoodindex_full,goodnotbothindex160)
                ds.series[outList[listindex]]['Flag'][goodbothindex160] = numpy.int32(200)
                ds.series[outList[listindex]]['Flag'][goodbothindex140] = numpy.int32(180)
                ds.series[outList[listindex]]['Flag'][goodbothindex150] = numpy.int32(190)
                ds.series[outList[listindex]]['Flag'][rstindex] = numpy.int32(171)
                ds.series[outList[listindex]]['Flag'][ravindex] = numpy.int32(171)
                ds.series[outList[listindex]]['Data'][rstindex] = numpy.float64(c.missing_value)
                ds.series[outList[listindex]]['Data'][ravindex] = numpy.float64(c.missing_value)
            
            ds.series[outList[listindex]]['Flag'][uindex] = numpy.int32(174)
            ds.series[outList[listindex]]['Flag'][Feindex] = numpy.int32(172)
            ds.series[outList[listindex]]['Flag'][Fsdindex] = numpy.int32(173)
            ds.series[outList[listindex]]['Data'][badflagindex] = numpy.float64(c.missing_value)
            ds.series[outList[listindex]]['Data'][uindex] = numpy.float64(c.missing_value)
            ds.series[outList[listindex]]['Data'][Feindex] = numpy.float64(c.missing_value)
            ds.series[outList[listindex]]['Data'][Fsdindex] = numpy.float64(c.missing_value)
        
        # calculate FAO56 using LAI input file
        if qcutils.cfkeycheck(cf,Base='PenmanMonteith',ThisOne='rlv2layer') and cf['PenmanMonteith']['rlv2layer'] == 'True':
            rlv, flag = get_leafresistance(cf,ds,outList[12])
            if 'rLv_2layer' in ds.series.keys():
                attr = qcutils.MakeAttributeDictionary(long_name='',units='s/m')
            else:
                attr = qcutils.MakeAttributeDictionary(long_name='leaf resistance from Penman-Monteith inversion, 2layer Ce-method, under well-illuminated (> 600 W m-2 Fsd) conditions',units='s/m')
            qcutils.CreateSeries(ds,'rLv_2layer',rlv,Flag=0,Attr=attr)
            ds.series['rLv_2layer']['Flag'] = flag
        
        # calculate 2-layer Penman potential ET
        if qcutils.cfkeycheck(cf,Base='PenmanMonteith',ThisOne='potential') and cf['PenmanMonteith']['potential'] == 'True':
            Fe,f,a = qcutils.GetSeriesasMA(ds,PMin[0])
            Ta,f,a = qcutils.GetSeriesasMA(ds,PMin[1])
            Ah,f,a = qcutils.GetSeriesasMA(ds,PMin[2])
            ps,f,a = qcutils.GetSeriesasMA(ds,PMin[3])
            Fn,f,a = qcutils.GetSeriesasMA(ds,PMin[5])
            Fsd,f,a = qcutils.GetSeriesasMA(ds,PMin[6])
            Fg,f,a = qcutils.GetSeriesasMA(ds,PMin[7])
            VPD,f,a = qcutils.GetSeriesasMA(ds,'VPD')
            Lv,f,a = qcutils.GetSeriesasMA(ds,'Lv')
            Cpm,f,a = qcutils.GetSeriesasMA(ds,'Cpm')
            rhom,f,a = qcutils.GetSeriesasMA(ds,'rhom')
            rav,f,a = qcutils.GetSeriesasMA(ds,'rav_2layer')
            flagList = [PMin[0],PMin[1],PMin[2],PMin[3],PMin[4],PMin[5],PMin[6],PMin[7],'rav_2layer']
            gamma = mf.gamma(ps,Cpm,Lv)
            delta = mf.delta(Ta)
            if 'gamma' not in ds.series.keys():
                attr = qcutils.MakeAttributeDictionary(long_name='Psychrometric coefficient',units='kPa/C')
                qcutils.CreateSeries(ds,'gamma',gamma,FList=[PMin[3],'Cpm','Lv'],Attr=attr)
            
            if 'delta' not in ds.series.keys():
                attr = qcutils.MakeAttributeDictionary(long_name='Slope of the saturation vapour pressure v temperature curve',units='kPa/C')
                qcutils.CreateSeries(ds,'delta',delta,FList=[PMin[1]],Attr=attr)
            
            Fa = Fn - Fg
            Fe_p = ((delta * Fa) + ((rhom * Cpm * VPD) / rav)) / (delta + gamma)
            ETp = Fe_p * 60 * 30 * 1000 / (Lv * c.rho_water)  # mm/30min for summing
            attr = qcutils.MakeAttributeDictionary(long_name='Penman potential ET',units='mm')
            qcutils.CreateSeries(ds,'ETp',ETp,FList=flagList,Attr=attr)
    
    # use drag coefficient method
    if Cdmethod == 'True':
        log.info('  Cd method (Jensen 1990, Leuning 2008) used to estimate aerodynamic resistance, ra')
        get_canopyresistance(cf,ds,Uavg,uindex,PMin,Level,critFsd,critFe)
        
        if qcutils.cfkeycheck(cf,Base='PenmanMonteith',ThisOne='rlm') and cf['PenmanMonteith']['rlm'] == 'True':
            rlm, flag = get_leafresistance(cf,ds,'rSm')
            if 'rLm' in ds.series.keys():
                attr = qcutils.MakeAttributeDictionary(long_name='',units='s/m')
            else:
                attr = qcutils.MakeAttributeDictionary(long_name='leaf resistance from Penman-Monteith inversion, Cd-method, under well-illuminated (> 600 W m-2 Fsd) conditions',units='s/m')
            qcutils.CreateSeries(ds,'rLm',rlm,Flag=0,Attr=attr)
            ds.series['rLm']['Flag'] = flag
    
    return

def PT100(ds,T_out,R_in,m):
    log.info(' Calculating temperature from PT100 resistance')
    R,f,a = qcutils.GetSeriesasMA(ds,R_in)
    R = m*R
    T = (-c.PT100_alpha+numpy.sqrt(c.PT100_alpha**2-4*c.PT100_beta*(-R/100+1)))/(2*c.PT100_beta)
    attr = qcutils.MakeAttributeDictionary(long_name='Calculated PT100 temperature using '+str(R_in),units='degC')
    qcutils.CreateSeries(ds,T_out,T,FList=[R_in],Attr=attr)

def ReplaceOnDiff(cf,ds,series=''):
    # Gap fill using data from alternate sites specified in the control file
    ts = ds.globalattributes['time_step']
    if len(series)!=0:
        ds_alt = {}                     # create a dictionary for the data from alternate sites
        open_ncfiles = []               # create an empty list of open netCDF files
        for ThisOne in series:          # loop over variables in the series list
            # has ReplaceOnDiff been specified for this series?
            if qcutils.incf(cf,ThisOne) and qcutils.haskey(cf,ThisOne,'ReplaceOnDiff'):
                # loop over all entries in the ReplaceOnDiff section
                for Alt in cf['Variables'][ThisOne]['ReplaceOnDiff'].keys():
                    if 'FileName' in cf['Variables'][ThisOne]['ReplaceOnDiff'][Alt].keys():
                        alt_filename = cf['Variables'][ThisOne]['ReplaceOnDiff'][Alt]['FileName']
                        if 'AltVarName' in cf['Variables'][ThisOne]['ReplaceOnDiff'][Alt].keys():
                            alt_varname = cf['Variables'][ThisOne]['ReplaceOnDiff'][Alt]['AltVarName']
                        else:
                            alt_varname = ThisOne
                        if alt_filename not in open_ncfiles:
                            n = len(open_ncfiles)
                            open_ncfiles.append(alt_filename)
                            ds_alt[n] = qcio.nc_read_series_file(alt_filename)
                        else:
                            n = open_ncfiles.index(alt_filename)
                        if 'Transform' in cf['Variables'][ThisOne]['ReplaceOnDiff'][Alt].keys():
                            AltDateTime = ds_alt[n].series['DateTime']['Data']
                            AltSeriesData = ds_alt[n].series[alt_varname]['Data']
                            TList = ast.literal_eval(cf['Variables'][ThisOne]['ReplaceOnDiff'][Alt]['Transform'])
                            for TListEntry in TList:
                                qcts.TransformAlternate(TListEntry,AltDateTime,AltSeriesData,ts=ts)
                        if 'Range' in cf['Variables'][ThisOne]['ReplaceOnDiff'][Alt].keys():
                            RList = ast.literal_eval(cf['Variables'][ThisOne]['ReplaceOnDiff'][Alt]['Range'])
                            for RListEntry in RList:
                                qcts.ReplaceWhenDiffExceedsRange(ds.series['DateTime']['Data'],ds.series[ThisOne],
                                                                 ds.series[ThisOne],ds_alt[n].series[alt_varname],
                                                                 RListEntry)
                    elif 'AltVarName' in cf['Variables'][ThisOne]['ReplaceOnDiff'][Alt].keys():
                        alt_varname = ThisOne
                        if 'Range' in cf['Variables'][ThisOne]['ReplaceOnDiff'][Alt].keys():
                            RList = ast.literal_eval(cf['Variables'][ThisOne]['ReplaceOnDiff'][Alt]['Range'])
                            for RListEntry in RList:
                                qcts.ReplaceWhenDiffExceedsRange(ds.series['DateTime']['Data'],ds.series[ThisOne],
                                                                 ds.series[ThisOne],ds.series[alt_varname],
                                                                 RListEntry)
                    else:
                        log.error('ReplaceOnDiff: Neither AltFileName nor AltVarName given in control file')
    else:
        log.error('ReplaceOnDiff: No input series specified')

def ReplaceWhereMissing(Destination,Primary,Secondary,FlagOffset=None,FlagValue=None):
    #print time.strftime('%X')+' Merging series '+Primary+' and '+Secondary+' into '+Destination
    p_data = Primary['Data'].copy()
    p_flag = Primary['Flag'].copy()
    s_data = Secondary['Data'].copy()
    s_flag = Secondary['Flag'].copy()
    if numpy.size(p_data)>numpy.size(s_data):
        p_data = p_data[0:numpy.size(s_data)]
    if numpy.size(s_data)>numpy.size(p_data):
        s_data = s_data[0:numpy.size(p_data)]
    index = numpy.where((abs(p_data-numpy.float64(c.missing_value))<c.eps)&
                        (abs(s_data-numpy.float64(c.missing_value))>c.eps))[0]
    p_data[index] = s_data[index]
    if FlagValue==None and FlagOffset!=None:
        p_flag[index] = s_flag[index] + numpy.int32(FlagOffset)
    elif FlagValue!=None and FlagOffset==None:
        p_flag[index] = numpy.int32(FlagValue)
    else:
        p_flag[index] = s_flag[index]
    Destination['Data'] = Primary['Data'].copy()
    Destination['Flag'] = Primary['Flag'].copy()
    Destination['Data'][0:len(p_data)] = p_data
    Destination['Flag'][0:len(p_flag)] = p_flag
    Destination['Attr']['long_name'] = 'Merged from original and alternate'
    Destination['Attr']['units'] = Primary['Attr']['units']

def ReplaceWhenDiffExceedsRange(DateTime,Destination,Primary,Secondary,RList):
    #print time.strftime('%X')+' Replacing '+Primary+' with '+Secondary+' when difference exceeds threshold'
    # get the primary data series
    p_data = numpy.ma.array(Primary['Data'])
    p_flag = Primary['Flag'].copy()
    # get the secondary data series
    s_data = numpy.ma.array(Secondary['Data'])
    s_flag = Secondary['Flag'].copy()
    # truncate the longest series if the sizes do not match
    if numpy.size(p_data)!=numpy.size(s_data):
        log.warning(' ReplaceWhenDiffExceedsRange: Series lengths differ, longest will be truncated')
        if numpy.size(p_data)>numpy.size(s_data):
            p_data = p_data[0:numpy.size(s_data)]
        if numpy.size(s_data)>numpy.size(p_data):
            s_data = s_data[0:numpy.size(p_data)]
    # get the difference between the two data series
    d_data = p_data-s_data
    # normalise the difference if requested
    if RList[3]=='s':
        d_data = (p_data-s_data)/s_data
    elif RList[3]=='p':
        d_data = (p_data-s_data)/p_data
    #si = qcutils.GetDateIndex(DateTime,RList[0],0)
    #ei = qcutils.GetDateIndex(DateTime,RList[1],0)
    Range = RList[2]
    Upper = numpy.float64(Range[0])
    Lower = numpy.float64(Range[1])
    index = numpy.ma.where((abs(d_data)<Lower)|(abs(d_data)>Upper))
    p_data[index] = s_data[index]
    p_flag[index] = numpy.int32(70)
    Destination['Data'] = numpy.ma.filled(p_data,numpy.float64(c.missing_value))
    Destination['Flag'] = p_flag.copy()
    Destination['Attr']['long_name'] = 'Replaced original with alternate when difference exceeded threshold'
    Destination['Attr']['units'] = Primary['Attr']['units']

def savitzky_golay(y, window_size, order, deriv=0):
    ''' Apply Savitsky-Golay low-pass filter to data.'''
    try:
        window_size = numpy.abs(numpy.int32(window_size))
        order = numpy.abs(numpy.int32(order))
    except ValueError, msg:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = numpy.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = numpy.linalg.pinv(b).A[deriv]
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - numpy.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + numpy.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = numpy.concatenate((firstvals, y, lastvals))
    return numpy.convolve( m, y, mode='valid')

def Square(Series):
    tmp = numpy.array([c.missing_value]*numpy.size(Series),Series.dtype)
    index = numpy.where(Series!=numpy.float64(c.missing_value))[0]
    tmp[index] = Series[index] ** 2
    return tmp

def SquareRoot(Series):
    tmp = numpy.array([c.missing_value]*numpy.size(Series),Series.dtype)
    index = numpy.where(Series!=numpy.float64(c.missing_value))[0]
    tmp[index] = Series[index] ** .5
    return tmp

def TaFromTv(cf,ds,Ta_out='Ta_CSAT',Tv_in='Tv_CSAT',Ah_in='Ah',ps_in='ps'):
    # Calculate the air temperature from the virtual temperature, the
    # absolute humidity and the pressure.
    # NOTE: the virtual temperature is used in place of the air temperature
    #       to calculate the vapour pressure from the absolute humidity, the
    #       approximation involved here is of the order of 1%.
    log.info(' Calculating Ta from Tv')
    Tv,f,a = qcutils.GetSeriesasMA(ds,Tv_in)
    Ah,f,a = qcutils.GetSeriesasMA(ds,Ah_in)
    ps,f,a = qcutils.GetSeriesasMA(ds,ps_in)
    nRecs = numpy.int32(ds.globalattributes['nc_nrecs'])
    Ta_flag = numpy.zeros(nRecs,numpy.int32)
    vp = mf.vapourpressure(Ah,Tv)
    mr = mf.mixingratio(ps,vp)
    q = mf.specifichumidity(mr)
    Ta_data = mf.tafromtv(Tv,q)
    mask = numpy.ma.getmask(Ta_data)
    index = numpy.where(mask.astype(numpy.int32)==1)
    Ta_flag[index] = numpy.int32(15)
    attr = qcutils.MakeAttributeDictionary(long_name='Ta calculated from Tv using '+Tv_in,units='C',standard_name='air_temperature')
    qcutils.CreateSeries(ds,Ta_out,Ta_data,Flag=Ta_flag,Attr=attr)

def TransformAlternate(TList,DateTime,Series,ts=30):
    # Apply polynomial transform to data series being used as replacement data for gap filling
    #print time.strftime('%X')+' Applying polynomial transform to '+ThisOne
    si = qcutils.GetDateIndex(DateTime,TList[0],ts=ts,default=0,match='exact')
    ei = qcutils.GetDateIndex(DateTime,TList[1],ts=ts,default=-1,match='exact')
    Series = numpy.ma.masked_where(abs(Series-numpy.float64(c.missing_value))<c.eps,Series)
    Series[si:ei] = qcutils.polyval(TList[2],Series[si:ei])
    Series = numpy.ma.filled(Series,numpy.float64(c.missing_value))

def UstarFromFh(cf,ds,us_out='ustar', L_out='L', Fm_out='Fm', T_in='Ta', Ah_in='Ah', p_in='ps', Fh_in='Fh', u_in='Ws', us_in='ustar', L_in='L'):
    # Calculate ustar from sensible heat flux, wind speed and
    # roughness length using Wegstein's iterative method.
    #  T is the air temperature, C
    #  p is the atmospheric pressure, kPa
    #  H is the sensible heat flux, W/m^2
    #  u is the wind speed, m/s
    #  z is the measurement height minus the displacement height, m
    #  z0 is the momentum roughness length, m
    log.info(' Calculating ustar from (Fh,Ta,Ah,p,u)')
    # get z-d (measurement height minus displacement height) and z0 from the control file
    if qcutils.cfkeycheck(cf,Base='Params',ThisOne='zmd') and qcutils.cfkeycheck(cf,Base='Params',ThisOne='z0'):
        zmd = numpy.float64(cf['Params']['zmd'])   # z-d for site
        z0 = numpy.float64(cf['Params']['z0'])     # z0 for site
    else:
        log.error('Parameters zmd or z0 not found in control file.  u* not determined from Fh')
        return
    if qcutils.cfkeycheck(cf,Base='FunctionArgs',ThisOne='ustarFh'):
        args = ast.literal_eval(cf['FunctionArgs']['ustarFh'])
        us_out = args[0]
        T_in = args[1]
        Ah_in = args[2]
        p_in = args[3]
        Fh_in = args[4]
        u_in = args[5]
        us_in = args[6]
    T,T_flag,T_attr = qcutils.GetSeriesasMA(ds,T_in)
    Ah,Ah_flag,Ah_attr = qcutils.GetSeriesasMA(ds,Ah_in)
    p,p_flag,p_attr = qcutils.GetSeriesasMA(ds,p_in)
    Fh,Fh_flag,Fh_attr = qcutils.GetSeriesasMA(ds,Fh_in)
    u,u_flag,u_attr = qcutils.GetSeriesasMA(ds,u_in)
    L,L_flag,L_attr = qcutils.GetSeriesasMA(ds,L_in)
    rhom,rhmo_flag,rhom_attr = qcutils.GetSeriesasMA(ds,'rhom')
    nRecs = numpy.size(Fh)
    us = numpy.zeros(nRecs,dtype=numpy.float64) + numpy.float64(c.missing_value)
    us_flag = numpy.ones(nRecs,dtype=numpy.int32)*ds.series[us_in]['Flag']
    i = numpy.where((abs(ds.series[T_in]['Data']-numpy.float64(c.missing_value))>c.eps)&(abs(ds.series[Ah_in]['Data']-numpy.float64(c.missing_value))>c.eps)&(abs(ds.series[p_in]['Data']-numpy.float64(c.missing_value))>c.eps)&(abs(ds.series[Fh_in]['Data']-numpy.float64(c.missing_value))>c.eps)&(abs(ds.series[u_in]['Data']-numpy.float64(c.missing_value))>c.eps))
    us[i] = qcutils.Wegstein(T[i], Ah[i], p[i], Fh[i], u[i], zmd, z0)    #for i in range(nRecs):
    badindex = numpy.where(abs(us-numpy.float64(c.missing_value))<c.eps)
    goodindex = numpy.where(abs(us-numpy.float64(c.missing_value))>c.eps)
    us_flag[goodindex] = numpy.int32(80)
        #else:
    us[badindex] = numpy.float64(c.missing_value)
    us_flag[badindex] = numpy.int32(81)
    attr = qcutils.MakeAttributeDictionary(long_name='ustar from (Fh,Ta,Ah,p,u)',units='m/s')
    qcutils.CreateSeries(ds,'uscalc',us,Flag=us_flag,Attr=attr)
    MergeSeries(cf,ds,us_out,[0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180])
    us_filled,us_filled_flag,us_filled_attr = qcutils.GetSeriesasMA(ds,us_out)
    Fm = numpy.square(us_filled) * rhom
    attr = qcutils.MakeAttributeDictionary(long_name='Fm from gapfilled ustar',units='kg/m/s2')
    qcutils.CreateSeries(ds,Fm_out,Fm,FList=['rhom',us_out],Attr=attr)
    i = numpy.where((abs(ds.series[T_in]['Data']-numpy.float64(c.missing_value))>c.eps)&(abs(ds.series[Ah_in]['Data']-numpy.float64(c.missing_value))>c.eps)&(abs(ds.series[p_in]['Data']-numpy.float64(c.missing_value))>c.eps)&(abs(ds.series[Fh_in]['Data']-numpy.float64(c.missing_value))>c.eps)&(abs(ds.series[u_in]['Data']-numpy.float64(c.missing_value))>c.eps))
    L[i] = mf.molen(T[i], Ah[i], p[i], ds.series[us_out]['Data'][i], Fh[i], fluxtype='sensible')
    attr = qcutils.MakeAttributeDictionary(long_name='',units='m')
    qcutils.CreateSeries(ds,L_in,L,FList=[T_in,Ah_in,p_in,us_out,Fh_in],Attr=attr)
    L,Lflag,a = qcutils.GetSeriesasMA(ds,L_in)
    zetaindex = numpy.where(L != 0)[0]
    zetaflagindex = numpy.where((L == 0)&(numpy.mod(Lflag,10)==0))[0]
    zeta = numpy.ones(len(L),dtype=numpy.float64)*c.missing_value
    zeta[zetaindex] = zmd / L[zetaindex]
    zf = numpy.zeros(nRecs,dtype=numpy.int32) + Lflag
    zf[zetaflagindex] = numpy.int32(22)
    attr = qcutils.MakeAttributeDictionary(long_name='L from gapfilled ustar',units='m/s')
    qcutils.CreateSeries(ds,L_out,L,FList=[T_in,Ah_in,p_in,Fh_in,us_out],Attr=attr)
    attr = qcutils.MakeAttributeDictionary(long_name='zeta from gapfilled L',units='m/s')
    qcutils.CreateSeries(ds,'zeta',zeta,Flag=zf,Attr=attr)
    #return us_in, us_out
    return

def write_sums(cf,ds,ThisOne,xlCol,xlSheet,DoSum='False',DoMinMax='False',DoMean='False',DoSubSum='False',DoSoil='False'):
    monthabr = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
    if qcutils.cfkeycheck(cf,Base='Params',ThisOne='firstMonth'):
        M1st = numpy.int32(cf['Params']['firstMonth'])
    else:
        M1st = 1
    if qcutils.cfkeycheck(cf,Base='Params',ThisOne='secondMonth'):
        M2nd = numpy.int32(cf['Params']['secondMonth'])
    else:
        M2nd = 12
    log.info(' Doing daily sums for '+ThisOne)
    Units = ds.series[ThisOne]['Attr']['units']
    if Units == 'mm/30min':
        Units = 'mm'
    
    xlRow = 1
    if xlCol == 0:
        xlSheet.write(xlRow,xlCol,'Year')
        xlCol = xlCol + 1
        xlSheet.write(xlRow,xlCol,'Month')
        xlCol = xlCol + 1
        xlSheet.write(xlRow,xlCol,'Day')
        xlCol = xlCol + 1
    xlSheet.write(xlRow,xlCol,'n')
    xlCol = xlCol + 1
    if DoMinMax == 'True':
        xlSheet.write(xlRow,xlCol,ThisOne+'_min')
        xlSheet.write(xlRow-1,xlCol,Units)
        xlCol = xlCol + 1
        xlSheet.write(xlRow,xlCol,ThisOne+'_max')
        if DoMean == 'True':
            xlSheet.write(xlRow-1,xlCol,Units)
            xlCol = xlCol + 1
            xlSheet.write(xlRow,xlCol,ThisOne+'_mean')
    elif DoMinMax == 'False' and DoMean == 'True':
        xlSheet.write(xlRow,xlCol,ThisOne+'_mean')
    elif DoMinMax == 'False' and DoMean == 'False':
        xlSheet.write(xlRow,xlCol,ThisOne)
        
    xlSheet.write(xlRow-1,xlCol,Units)

    if DoSubSum == 'True':
        xlCol = xlCol + 1
        xlSheet.write(xlRow,xlCol,'Pos n')
        xlCol = xlCol + 1
        xlSheet.write(xlRow,xlCol,ThisOne+'_pos')
        xlSheet.write(xlRow-1,xlCol,Units)
        xlCol = xlCol + 1
        xlSheet.write(xlRow,xlCol,'Neg n')
        xlCol = xlCol + 1
        xlSheet.write(xlRow,xlCol,ThisOne+'_neg')
        xlSheet.write(xlRow-1,xlCol,Units)
    
    data,f,a = qcutils.GetSeriesasMA(ds,ThisOne)
    year,ff,aa = qcutils.GetSeriesasMA(ds,'Year')
    #data0,f,a = qcutils.GetSeriesasMA(ds,ThisOne)
    #data = numpy.zeros(len(data0), dtype=numpy.float64) + data0
    #data = numpy.ma.masked_where(abs(ds.series[ThisOne]['Data']-numpy.float64(c.missing_value))<c.eps,ds.series[ThisOne]['Data'])
    #data.dtype = numpy.float64
    for AnotherOne in range(year[0],year[len(year)-1]+1):
        yearindex = numpy.where(year==AnotherOne)[0]
        for month in range(M1st,M2nd+1):
            if month == 1 or month == 3 or month == 5 or month == 7 or month == 8 or month == 10 or month == 12:
                dRan = 31
            if month == 2:
                if ds.series['Year']['Data'][0] % 4 == 0:
                    dRan = 29
                else:
                    dRan = 28
            if month == 4 or month == 6 or month == 9 or month == 11:
                dRan = 30
                
            for day in range(1,dRan+1):
                xlRow = xlRow + 1
                PMList = ['GSv_1layer_mol', 'GSv_2layer_mol', 'GSv_top_mol', 'GSv_base_mol', 'GSv_full_mol', 'GSm_mol', 'rav_1layer', 'rSv_1layer', 'rLv_1layer', 'GSv_1layer', 'rav_2layer', 'rSv_2layer', 'rLv_2layer', 'GSv_2layer', 'rav_base', 'rSv_base', 'GSv_base', 'rav_top', 'rSv_top', 'GSv_top', 'rav_full', 'rSv_full', 'GSv_full', 'ram', 'rSm', 'rLm', 'GSm']
                CList = ['CE_mmol','ER_dark_mmol','ER_night_mmol','CE_NEEmax_mmol','GPP','GPP_mmol','C_ppm']
                VarList = PMList + CList
                if ThisOne in VarList:
                    di = numpy.where((year==AnotherOne) & (ds.series['Month']['Data']==month) & (ds.series['Day']['Data']==day) & (numpy.mod(ds.series[ThisOne]['Flag'],10) == 0))[0]
                    ti = numpy.where((year==AnotherOne) & (ds.series['Month']['Data']==month) & (ds.series['Day']['Data']==day))[0]
                    nRecs = len(ti)
                    check = numpy.ma.empty(nRecs,str)
                    for i in range(nRecs):
                        index = ti[i]
                        check[i] = ds.series['Day']['Data'][index]
                    if len(check) < 48:
                        di = []
                else:
                    di = numpy.where((year==AnotherOne) & (ds.series['Month']['Data']==month) & (ds.series['Day']['Data']==day))[0]
                    nRecs = len(di)
                    check = numpy.ma.empty(nRecs,str)
                    for i in range(nRecs):
                        index = di[i]
                        check[i] = ds.series['Day']['Data'][index]
                    if len(check) < 47:
                        di = []
                
                if DoSoil == 'True':
                    Num,Av = get_soilaverages(data[di])
                    if xlCol == 4:
                        xlCol = 3
                        xlSheet.write(xlRow,xlCol-3,AnotherOne)
                        xlSheet.write(xlRow,xlCol-2,monthabr[month-1])
                        xlSheet.write(xlRow,xlCol-1,day)
                    else:
                        xlCol = xlCol - 1
                else:
                    if DoSum == 'True':
                        Num,Sum = get_sums(data[di])
                    if DoMinMax == 'True':
                        if ThisOne == 'WUE':
                            Num,Min,Max = get_minmax_WUE(data[di])
                        else:
                            Num,Min,Max = get_minmax(data[di])
                    if DoMean == 'True':
                        if DoMinMax == 'True':
                            if ThisOne == 'WUE':
                                Num2,Av = get_averages_WUE(data[di])
                            else:
                                Num2,Av = get_averages(data[di])
                        else:
                            if ThisOne == 'WUE':
                                Num,Av = get_averages_WUE(data[di])
                            else:
                                Num,Av = get_averages(data[di])
                    if DoSubSum == 'True':
                        PosNum,NegNum,SumPos,SumNeg = get_subsums(data[di])
                    xlCol = 3
                    xlSheet.write(xlRow,xlCol-3,AnotherOne)
                    xlSheet.write(xlRow,xlCol-2,monthabr[month-1])
                    xlSheet.write(xlRow,xlCol-1,day)
                
                xlSheet.write(xlRow,xlCol,Num)
                xlCol = xlCol + 1
                if DoSoil == 'True':
                    xlSheet.write(xlRow,xlCol,Av)
                elif DoMinMax == 'True':
                    xlSheet.write(xlRow,xlCol,Min)
                    xlCol = xlCol + 1
                    xlSheet.write(xlRow,xlCol,Max)
                    if DoMean == 'True':
                        xlCol = xlCol + 1
                        xlSheet.write(xlRow,xlCol,Av)
                elif DoMinMax == 'False' and DoMean == 'True':
                    xlSheet.write(xlRow,xlCol,Av)
                elif DoSum == 'True':
                    xlSheet.write(xlRow,xlCol,Sum)
                    if DoSubSum == 'True':
                        xlCol = xlCol + 1
                        xlSheet.write(xlRow,xlCol,PosNum)
                        xlCol = xlCol + 1
                        xlSheet.write(xlRow,xlCol,SumPos)
                        xlCol = xlCol + 1
                        xlSheet.write(xlRow,xlCol,NegNum)
                        xlCol = xlCol + 1
                        xlSheet.write(xlRow,xlCol,SumNeg)
    
    if DoSoil == 'True': 
        return xlCol,xlSheet
    else:
        return
