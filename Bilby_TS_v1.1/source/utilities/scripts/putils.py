#    putils.py
#    Copyright (C) 2015  Dr James Cleverly
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

import ast
import constants as c
import datetime
import dateutil
import math
import meteorologicalfunctions as mf
import numpy
import logging
import xlwt

log = logging.getLogger('partition.utils')

def cfkeycheck(cf,Base='Variables',ThisOne=[],key=[]):
    if len(ThisOne) == 0:
        return
    if len(key) == 0:
        if Base in cf.keys() and ThisOne in cf[Base].keys():
            return ThisOne in cf[Base].keys()
        else:
            return
    else:
        if Base in cf.keys() and ThisOne in cf[Base].keys():
            return key in cf[Base][ThisOne].keys()
        else:
            return

def CreateSeries(ds,Label,Data,FList=[''],Flag=None,Descr='',Units='',Standard='not defined'):
    """
    Create a series (1d array) of data in the data structure.
    
    If the series already exists in the data structure, data values and QC flags will be
    overwritten but attributes will be preserved.  However, the long_name and Units attributes
    are treated differently.  The existing long_name will have Descr appended to it.  The
    existing units will be overwritten with units.
    
    This utility is the prefered method for creating or updating a data series because
    it imnplements a consistent method for creating series in the data structure.  Direct
    writes to the contents of the data structure are discouraged (unless PRI wrote the code!).
    """
    ds.series['_tmp_'] = {}                       # create a temporary series to avoid premature overwrites
    # put the data into the temporary series
    ds.series['_tmp_']['Data'] = numpy.ma.filled(Data,float(-9999))
    # copy or make the QC flag
    if Flag == None:
        ds.series['_tmp_']['Flag'] = MakeQCFlag(ds,FList)
    else:
        ds.series['_tmp_']['Flag'] = Flag
    # do the attributes
    ds.series['_tmp_']['Attr'] = {}
    if Label in ds.series.keys():                 # check to see if the series already exists
        for attr in ds.series[Label]['Attr']:     # if it does, copy the existing attributes
            if attr=='long_name':               # append new description to existing one
                attr_value = ds.series[Label]['Attr']['long_name']
                attr_value = attr_value + ', ' + Descr
                ds.series['_tmp_']['Attr'][attr] = attr_value
            elif attr=='units':                   # overwrite existing units with new ones
                ds.series['_tmp_']['Attr']['units'] = Units
            else:                                 # otherwise do a straight copy
                ds.series['_tmp_']['Attr'][attr] = ds.series[Label]['Attr'][attr]
    else:
        ds.series['_tmp_']['Attr']['long_name'] = Descr
        ds.series['_tmp_']['Attr']['units'] = Units
        ds.series['_tmp_']['Attr']['standard_name'] = Standard
    ds.series[unicode(Label)] = ds.series['_tmp_']     # copt temporary series to new series
    del ds.series['_tmp_']                        # delete the temporary series

def Fm(z, z0, L):
    ''' Integral form of the adiabatic correction to the wind speed profile.'''
    Fm = math.log(z/z0)                 # Neutral case
    if L<0:                             # Unstable case
        R0 = (1-c.gamma*z0/L)**0.25
        R1 = (1-c.gamma*z/L)**0.25
        x = ((R0+1)/(R1+1))**2
        Y = (R0*R0+1)/(R1*R1+1)
        w = z/z0
        V = 2 * numpy.arctan((R1-R0)/(1+R0*R1))
        Fm = math.log(w*Y*x)+V
    elif ((L>-200)|(L>200)):            # Neutral case
        Fm = math.log(z/z0)
    elif (z/L<=1):                      # Stable case, z < L
        x = math.log(z/z0)
        Y = c.beta*z/L
        Fm = x+Y
    elif ((z/L>1)&(z0/L<1)):            # Stable case, z > L > z0
        x = math.log(L/z0)
        Y = (1+c.beta)*math.log(z/L)
        Fm = x+c.beta+Y
    elif (z0/L>1):                      # Stable, L < z0
        Fm = (1+c.beta)*math.log(z/z0)
    else:
        print 'Error in function Fm'
    return Fm

def Fustar(T, Ah, p, Fh, u, z, z0, ustar):
#' Function used in iteration method to solve for ustar.
#' The function used is:
#'  ustar = u*k/Fm(z/L,z0/L)
#' where Fm is the integral form of the PHIm, the adiabatic
#' correction to the logarithmic wind speed profile.
#' Evaluate the function for ustar with this value for L.
    MO = mf.molen(T, Ah, p, ustar, Fh)
    Fustar = u*c.k/(Fm(z, z0, MO))
    return Fustar
    
def GetAverageList(cf,ThisOne,default=""):
    if cfkeycheck(cf,ThisOne=ThisOne):
        if default == "":
            log.error(' GetAverageList: no "Source" series provided for '+ThisOne)
        
        alist = default
    else:
        log.error('  GetAverageList: '+ThisOne+' not in control file')
        alist = ""
    
    return alist

def GetDateIndex(datetimeseries,date,default):
    # return the index of a date/datetime string in an array of datetime objects
    #  datetimeseries - array of datetime objects
    #  date - a date or date/time string in a format dateutils can parse
    #  default - default value, integer
    try:
        dateindex = datetimeseries.index(dateutil.parser.parse(date))
    except ValueError:
        dateindex = default
    return dateindex

def GetMergeList(cf,ThisOne,default=""):
    if cfkeycheck(cf,ThisOne=ThisOne):
        if default == "":
            log.error(' GetMergeList: no "Source" series provided for '+ThisOne)
        
        mlist = default
    else:
        log.error('  GetMergeList: '+ThisOne+' not in control file')
        mlist = ""
    
    return mlist

def GetSeries(ds,ThisOne,si=0,ei=-1):
    Series = ds.series[ThisOne]['Data']
    if 'Flag' in ds.series[ThisOne].keys():
        Flag = ds.series[ThisOne]['Flag']
    else:
        nRecs = numpy.size(ds.series[ThisOne]['Data'])
        Flag = numpy.zeros(nRecs,dtype=int)
    if ei==-1:
        Series = Series[si:]
        Flag = Flag[si:]
    else:
        Series = Series[si:ei+1]
        Flag = Flag[si:ei+1]
    return Series,Flag

def GetSeriesasMA(ds,ThisOne,si=0,ei=-1):
    Series = numpy.ma.masked_where(ds.series[ThisOne]['Data']==float(-9999),ds.series[ThisOne]['Data'])
    if 'Flag' in ds.series[ThisOne].keys():
        Flag = ds.series[ThisOne]['Flag']
    else:
        nRecs = numpy.size(ds.series[ThisOne]['Data'])
        Flag = numpy.zeros(nRecs,dtype=int)
    if ei==-1:
        Series = Series[si:]
        Flag = Flag[si:]
    else:
        Series = Series[si:ei+1]
        Flag = Flag[si:ei+1]
    return Series,Flag

def GetSeriesStats(cf,ds):
    # open an Excel file for the flag statistics
    level = ds.globalattributes['Level']
    xlFileName = cf['Files'][level]['xlFlagFilePath']+cf['Files'][level]['xlFlagFileName']
    log.info(' Writing flag stats to Excel file '+xlFileName)
    xlFile = xlwt.Workbook()
    if cf['General']['Platform'] == 'Mac':
        xlFile.dates_1904 = True
    xlFlagSheet = xlFile.add_sheet('Flag')
    # get the flag statistics
    xlRow = 0
    xlCol = 0
    xlFlagSheet.write(xlRow,xlCol,'0:')
    xlFlagSheet.write(xlRow,xlCol+1,ds.globalattributes['Flag0'])
    xlFlagSheet.write(xlRow,xlCol+2,'1:')
    xlFlagSheet.write(xlRow,xlCol+3,ds.globalattributes['Flag1'])
    xlFlagSheet.write(xlRow,xlCol+4,'2:')
    xlFlagSheet.write(xlRow,xlCol+5,ds.globalattributes['Flag2'])
    xlFlagSheet.write(xlRow,xlCol+6,'3:')
    xlFlagSheet.write(xlRow,xlCol+7,ds.globalattributes['Flag3'])
    xlFlagSheet.write(xlRow,xlCol+8,'4:')
    xlFlagSheet.write(xlRow,xlCol+9,ds.globalattributes['Flag4'])
    xlFlagSheet.write(xlRow,xlCol+10,'5:')
    xlFlagSheet.write(xlRow,xlCol+11,ds.globalattributes['Flag5'])
    xlFlagSheet.write(xlRow,xlCol+12,'6:')
    xlFlagSheet.write(xlRow,xlCol+13,ds.globalattributes['Flag6'])
    xlFlagSheet.write(xlRow,xlCol+14,'7:')
    xlFlagSheet.write(xlRow,xlCol+15,ds.globalattributes['Flag7'])
    xlRow = xlRow + 1
    xlFlagSheet.write(xlRow,xlCol,'10:')
    xlFlagSheet.write(xlRow,xlCol+1,ds.globalattributes['Flag10'])
    xlFlagSheet.write(xlRow,xlCol+2,'11:')
    xlFlagSheet.write(xlRow,xlCol+3,ds.globalattributes['Flag11'])
    xlFlagSheet.write(xlRow,xlCol+4,'12:')
    xlFlagSheet.write(xlRow,xlCol+5,ds.globalattributes['Flag12'])
    xlFlagSheet.write(xlRow,xlCol+6,'13:')
    xlFlagSheet.write(xlRow,xlCol+7,ds.globalattributes['Flag13'])
    xlFlagSheet.write(xlRow,xlCol+8,'14:')
    xlFlagSheet.write(xlRow,xlCol+9,ds.globalattributes['Flag14'])
    xlFlagSheet.write(xlRow,xlCol+10,'16:')
    xlFlagSheet.write(xlRow,xlCol+11,ds.globalattributes['Flag16'])
    xlFlagSheet.write(xlRow,xlCol+12,'17:')
    xlFlagSheet.write(xlRow,xlCol+13,ds.globalattributes['Flag17'])
    xlRow = xlRow + 1
    xlFlagSheet.write(xlRow,xlCol,'20:')
    xlFlagSheet.write(xlRow,xlCol+1,ds.globalattributes['Flag20'])
    xlFlagSheet.write(xlRow,xlCol+2,'21:')
    xlFlagSheet.write(xlRow,xlCol+3,ds.globalattributes['Flag21'])
    xlFlagSheet.write(xlRow,xlCol+4,'22:')
    xlFlagSheet.write(xlRow,xlCol+5,ds.globalattributes['Flag22'])
    bins = numpy.arange(-0.5,23.5)
    xlRow = 5
    xlCol = 1
    for Value in bins[:len(bins)-1]:
        xlFlagSheet.write(xlRow,xlCol,int(Value+0.5))
        xlCol = xlCol + 1
    xlRow = xlRow + 1
    xlCol = 0
    dsVarNames = ds.series.keys()
    dsVarNames.sort(key=unicode.lower)
    for ThisOne in dsVarNames:
        data,flag = GetSeries(ds, ThisOne)
        hist, bin_edges = numpy.histogram(flag, bins=bins)
        xlFlagSheet.write(xlRow,xlCol,ThisOne)
        xlCol = xlCol + 1
        for Value in hist:
            xlFlagSheet.write(xlRow,xlCol,float(Value))
            xlCol = xlCol + 1
        xlCol = 0
        xlRow = xlRow + 1
    xlFile.save(xlFileName)

def haskey(cf,ThisOne,key):
    return key in cf['Variables'][ThisOne].keys()

def incf(cf,ThisOne):
    return ThisOne in cf['Variables'].keys()

def MakeQCFlag(ds,SeriesList):
    flag = []
    if len(SeriesList)<=0:
        #log.info('  MakeQCFlag: no series list specified')
        pass
    if len(SeriesList)==1:
        if SeriesList[0] in ds.series.keys():
            flag = ds.series[SeriesList[0]]['Flag'].copy()
        else:
            log.error('  MakeQCFlag: series '+str(SeriesList[0])+' not in ds.series')
    if len(SeriesList)>1:
        for ThisOne in SeriesList:
            if ThisOne in ds.series.keys():
                if len(flag)==0:
                    #flag = numpy.ones(numpy.size(ds.series[ThisOne]['Flag']))
                    flag = ds.series[ThisOne]['Flag'].copy()
                else:
                    tmp_flag = ds.series[ThisOne]['Flag'].copy()      # get a temporary copy of the flag
                    index = numpy.where(numpy.mod(tmp_flag,10)==0)    # find the elements with flag = 0, 10, 20 etc
                    tmp_flag[index] = 0                               # set them all to 0
                    flag = numpy.maximum(flag,tmp_flag)               # now take the maximum
            else:
                log.error('  MakeQCFlag: series '+ThisOne+' not in ds.series')
    return flag

def polyval(p,x):
    """
    Replacement for the polyval routine in numpy.  This version doesnt check the
    input variables to make sure they are array_like.  This means that when
    masked arrays are treated correctly when they are passed to this routine.
    Parameters
    ----------
     p : a 1D array of coefficients, highest order first
     x : a 1D array of points at which to evaluate the polynomial described by
         the coefficents in p
    Example
    -------
    >>> x = numpy.array([1,2,3])
    >>> p = numpy.array([2,0])
    >>> qcutils.polyval(p,x)
        array([2,4,6])
    >>> y = numpy.array([1,-9999,3])
    >>> y = numpy.ma.masked_where(y==-9999,y)
    >>> qcutils.polyval(p,y)
    masked_array(data = [2 -- 6],
                 mask = [False True False],
                 fill_value = 999999)
    """
    y = 0
    for i in range(len(p)):
        y = x*y + p[i]
    return y

def startlog(loggername,loggerfile):
    logger = logging.getLogger(loggername)
    logger.setLevel(logging.DEBUG)
    fh = logging.FileHandler(loggerfile)
    fh.setLevel(logging.DEBUG)
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s %(name)-15s %(levelname)-6s %(message)s', '%d-%m-%y %H:%M')
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)
    logger.addHandler(fh)
    logger.addHandler(ch)
    return logger

def Wegstein(T, Ah, p, Fh, u, z, z0):

    NumIters = 50
    SolveErr = numpy.float64(0.001)
 
    FirstEst =  u*c.k/math.log(z/z0)
    ustar = Fustar(T, Ah, p, Fh, u, z, z0, FirstEst)
    Inc = ustar-FirstEst
    IncDiv = -Inc
    Residual = ustar-Fustar(T, Ah, p, Fh, u, z, z0, ustar)
 
    i = 1
    while (i<NumIters)&(float(Residual)>float(SolveErr)):
        IncDiv = (IncDiv/Residual)-1
        if (IncDiv == 0):
            print 'Wegstein: IncDiv equals 0'
            ustar = u*c.k/math.log(z/z0)
            break
        Inc = Inc/IncDiv
        ustar = ustar+Inc
        IncDiv = Residual
        Residual = ustar-Fustar(T, Ah, p, Fh, u, z, z0, ustar)
        if (abs(ustar)<=1):
            RangeErr = SolveErr
        else:
            RangeErr = SolveErr*abs(ustar)
        if (abs(Inc)<=RangeErr):
            if (abs(Residual)<=10*RangeErr):
                break
        i = i + 1
    if (i==NumIters):
        print 'Wegstein: did not converge'
        ustar = u*c.k/math.log(z/z0)
    return ustar
def GetMergeSeriesKeys(cf,ThisOne,section=''):
    if len(section)==0: section = 'Variables'
    if 'Source' in cf[section][ThisOne]['MergeSeries'].keys():
        mlist = ast.literal_eval(cf[section][ThisOne]['MergeSeries']['Source'])
    else:
        log.error('  GetMergeSeriesKeys: key "Source" not in control file MergeSeries section for '+ThisOne)
        mlist = []
    if 'standard_name' in cf[section][ThisOne]['MergeSeries'].keys():
        standardname = str(cf[section][ThisOne]['MergeSeries']['standard_name'])
    else:
        standardname = 'not defined'
    return mlist, standardname

def get_cfsection(cf,series='',mode='quiet'):
    '''
    Find the section in the control file that contains an entry for the series "series".
    USEAGE:  section = qcutils.get_cfsection(cf,series=<series_name>)
    INPUT:   cf            - a control file object (from ConfigObj)
             <series_name> - the name of the series (string)
    RETURNS: section       - the name of the section containing an entry for <series_name> (string)
    Note that the returned section name is an empty string if there is no entry for <series_name> in
    the control file.
    '''
    section = ''
    sectionlist = ['Variables','Drivers','Fluxes','Respiration','Partition']
    if len(series)==0:
        msgtxt = ' get_cfsection: no input series specified'
        if mode!='quiet': log.info(msgtxt)
        return section
    for ThisSection in sectionlist:
        if ThisSection in cf.keys():
            if series in cf[ThisSection]: section = ThisSection
    if len(section)==0:
        msgtxt = ' get_cfsection: series '+str(series)+' not found in control file'
        if mode!='quiet': log.info(msgtxt)
    return section

