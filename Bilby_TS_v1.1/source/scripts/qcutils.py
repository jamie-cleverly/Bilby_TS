#    qcutils.py  Utilities for database management
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

import ast
import constants as c
import datetime
import dateutil
import logging
import math
import meteorologicalfunctions as mf
import numpy
import os
import pytz
import xlrd
import xlwt

log = logging.getLogger('qc.utils')

def bp(fx,tao):
    """
    Function to calculate the b and p coeficients of the Massman frequency correction.
    """
    bp = 2 * c.Pi * fx * tao
    return bp

def cfkeycheck(cf,Base='Variables',ThisOne=[],key=[]):
    if len(ThisOne) == 0:
        if Base in cf.keys():
            return Base in cf.keys()
        else:
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

def CheckQCFlags(ds):
    """
    Purpose:
     Make sure that all values of -9999 in a data series have a non-zero QC flag value.
    Usage:
     qcutils.CheckQCFlags(ds)
    Author: PRI
    Date: August 2014
    """
    for ThisOne in ds.series.keys():
        data = numpy.ma.masked_values(ds.series[ThisOne]["Data"],c.missing_value)
        flag = numpy.ma.masked_equal(ds.series[ThisOne]["Flag"],0)
        mask = data.mask&flag.mask
        index = numpy.ma.where(mask==True)[0]
        ds.series[ThisOne]["Flag"][index] = numpy.int32(21)

def CheckTimeStep(ds,fix='gaps'):
    """
    Purpose:
     Checks the datetime series in the data structure ds to see if there are
     any missing time stamps.  If missing time stamps are found, this routine
     will print a message to the screen and log file and then check to see
     if all gaps in the datetime series are multiples of the time step.  If
     not, the routine prints a message to the screen and the log file.
     This function returns a logical variable that is true if any gaps exist
     in the time stamp.
    Useage:
     has_gaps = CheckTimeSTep(ds)
     if has_gaps:
         <do something about missing time stamps>
    Author: PRI
    Date: April 2013
    """
    has_gaps = False
    ts = numpy.int32(ds.globalattributes['time_step'])
    xldt = ds.series['xlDateTime']['Data']
    dt = numpy.zeros(len(xldt),dtype=numpy.float64)
    dt[0] = ts
    dt[1:-1] = (1440*(xldt[1:-1]-xldt[0:-2])+0.5).astype(numpy.int32)
    dt[-1] = (1440*(xldt[-1]-xldt[-2])+0.5).astype(numpy.int32)
    index = numpy.where(dt!=ts)[0]
    if len(index)!=0:
        has_gaps = True
        log.error(' CheckTimeStep: '+str(len(index))+ ' gaps found in the time series')
        dtmin = numpy.min(dt)
        dtmax = numpy.max(dt)
        if dtmin < ts:
            print dtmin,ds.series['DateTime']['Data'][numpy.where(dt==dtmin)[0]-3],numpy.where(dt==dtmin)[0]-3
            print dtmin,ds.series['DateTime']['Data'][numpy.where(dt==dtmin)[0]-2],numpy.where(dt==dtmin)[0]-2
            print dtmin,ds.series['DateTime']['Data'][numpy.where(dt==dtmin)[0]-1],numpy.where(dt==dtmin)[0]-1
            print dtmin,ds.series['DateTime']['Data'][numpy.where(dt==dtmin)[0]],numpy.where(dt==dtmin)[0]
            print dtmin,ds.series['DateTime']['Data'][numpy.where(dt==dtmin)[0]+1],numpy.where(dt==dtmin)[0]+1
            log.critical(' CheckTimeStep: duplicate or overlapping times found, fix the L1 spreadsheet')
            sys.exit()
        if numpy.min(numpy.mod(dt,ts))!=0 or numpy.max(numpy.mod(dt,ts))!=0:
            log.critical(' CheckTimeStep: time gaps are not multiples of the time step ('+str(ts)+'), fix the L1 spreadsheet')
            sys.exit()
        if dtmax > ts:
            log.info(' CheckTimeStep: one or more time gaps found')
            if fix.lower()=='gaps':
                log.info(' CheckTimeStep: calling FixTimeGaps to fix time gaps')
                FixTimeGaps(ds)
    else:
        log.info(' CheckTimeStep: no time gaps found')
    return has_gaps

def contiguous_regions(condition):
    """
    Purpose:
     Finds contiguous True regions of the boolean array "condition". Returns
     a 2D array where the first column is the start index of the region and the
     second column is the end index.
    Author: Joe Kington (via StackOverflow)
    Date: September 2014
    """
    # Find the indicies of changes in "condition"
    d = numpy.diff(condition)
    idx, = d.nonzero() 
    # We need to start things after the change in "condition". Therefore, 
    # we'll shift the index by 1 to the right.
    idx += 1
    if condition[0]:
        # If the start of condition is True prepend a 0
        idx = numpy.r_[0, idx]
    if condition[-1]:
        # If the end of condition is True, append the length of the array
        idx = numpy.r_[idx, condition.size] # Edit
    # Reshape the result into two columns
    idx.shape = (-1,2)
    return idx

def CreateSeries(ds,Label,Data,FList=None,Flag=None,Attr=None):
    """
        Create a series (1d array) of data in the data structure.
        
        If the series already exists in the data structure, data values and QC flags will be
        overwritten but attributes will be preserved.  However, the long_name and units attributes
        are treated differently.  The existing long_name will have long_name appended to it.  The
        existing units will be overwritten with units.
        
        This utility is the prefered method for creating or updating a data series because
        it implements a consistent method for creating series in the data structure.  Direct
        writes to the contents of the data structure are discouraged (unless PRI wrote the code:=P).
        """
    ds.series['_tmp_'] = {}                       # create a temporary series to avoid premature overwrites
    # put the data into the temporary series
    if numpy.ma.isMA(Data):
        ds.series['_tmp_']['Data'] = numpy.ma.filled(Data,numpy.float64(c.missing_value))
    else:
        ds.series['_tmp_']['Data'] = numpy.array(Data)
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
                if Attr[attr] != '':
                    attr_value = attr_value + ', ' + Attr[attr]
                ds.series['_tmp_']['Attr'][attr] = attr_value
            elif attr=='ancillary_variables':
                attr_value = str(Label)+' QC flag'
                ds.series['_tmp_']['Attr'][attr] = attr_value
            elif attr in Attr:
                if ds.series[Label]['Attr'][attr]!=Attr[attr]:
                    ds.series['_tmp_']['Attr'][attr] = Attr[attr]
                else:
                    ds.series['_tmp_']['Attr'][attr] = ds.series[Label]['Attr'][attr]
            else:
                ds.series['_tmp_']['Attr'][attr] = ds.series[Label]['Attr'][attr]
    else:
        for item in Attr:
            if item=='ancillary_variables':
                attr_value = str(Label)+' QC flag'
                ds.series['_tmp_']['Attr'][item] = attr_value
            else:
                ds.series['_tmp_']['Attr'][item] = Attr[item]
    ds.series[unicode(Label)] = ds.series['_tmp_']     # copy temporary series to new series
    del ds.series['_tmp_']                        # delete the temporary series

def file_exists(filename,mode="verbose"):
    if not os.path.exists(filename):
        if mode=="verbose":
            log.error(' File '+filename+' not found')
        return False
    else:
        return True

def Fm(z, z0, L):
    ''' Integral form of the adiabatic correction to the wind speed profile.'''
    try:
        nRecs = len(L)
        xtract = numpy.where((L.mask==True)|(L==0))[0]
        Lx = numpy.zeros(nRecs,dtype=numpy.float64) + L
        Ly = numpy.zeros(nRecs,dtype=numpy.float64) + L
        Lx[xtract] = numpy.float64(c.missing_value)
        Ly[xtract] = numpy.float64(-1) / numpy.float64(c.missing_value)
        xstable = numpy.where(z0/Lx>1)[0]                           # Stable, L < z0
        vstable = numpy.where((z/Lx>1)&(z0/Ly<1))[0]                 # Stable case, z > L > z0
        stable = numpy.where(z/Ly<=1)[0]                            # Stable case, z < L
        calm = numpy.where(L==0)[0]                                # Very stable case, u = 0
        unstable = numpy.where(L<0)[0]                             # Unstable case
        neutral = numpy.where((L<-200)|(L>200))[0]                 # Neutral case
        R0 = numpy.zeros(nRecs,dtype=numpy.float64) + numpy.float64(c.missing_value)
        R1 = numpy.zeros(nRecs,dtype=numpy.float64) + numpy.float64(c.missing_value)
        x = numpy.zeros(nRecs,dtype=numpy.float64) + numpy.float64(c.missing_value)
        Y = numpy.zeros(nRecs,dtype=numpy.float64) + numpy.float64(c.missing_value)
        V = numpy.zeros(nRecs,dtype=numpy.float64) + numpy.float64(c.missing_value)
        R0[unstable] = (1-c.gamma*z0/L[unstable])**0.25
        R1[unstable] = (1-c.gamma*z/L[unstable])**0.25
        x[vstable] = numpy.log(L[vstable]/z0)
        x[stable] = numpy.log(z/z0)
        x[unstable] = ((R0[unstable]+1)/(R1[unstable]+1))**2
        Y[vstable] = (1+c.beta)*numpy.log(z/L[vstable])
        Y[stable] = c.beta*z/L[stable]
        Y[unstable] = (R0[unstable]*R0[unstable]+1)/(R1[unstable]*R1[unstable]+1)
        w = z/z0
        V[unstable] = 2 * numpy.arctan((R1[unstable]-R0[unstable])/(1+R0[unstable]*R1[unstable]))
        Fm = numpy.zeros(nRecs,dtype=numpy.float64) + numpy.log(z/z0)     # Neutral case               neutral
        Fm[xstable] = (1+c.beta)*numpy.log(z/z0)                           # Stable, L < z0             xstable
        Fm[vstable] = x[vstable]+c.beta+Y[vstable]                         # Stable case, z > L > z0    vstable
        Fm[stable] = x[stable]+Y[stable]                                   # Stable case, z < L         stable
        Fm[calm] = (1+c.beta)*numpy.log(z/z0)                              # Very stable case, u = 0    calm
        Fm[unstable] = numpy.log(w*Y[unstable]*x[unstable])+V[unstable]    # Unstable case              unstable
        Fm[neutral] = numpy.log(z/z0)                                      # Neutral case               neutral
    except:
        if ((L<-200)|(L>200)):            # Neutral case    neutral
            Fm = math.log(z/z0)
        elif L<0:                             # Unstable case    unstable
            R0 = (1-c.gamma*z0/L)**0.25
            R1 = (1-c.gamma*z/L)**0.25
            x = ((R0+1)/(R1+1))**2
            Y = (R0*R0+1)/(R1*R1+1)
            w = z/z0
            V = 2 * numpy.arctan((R1-R0)/(1+R0*R1))
            Fm = math.log(w*Y*x)+V
        elif L == 0:                        # Very stable case, u = 0    calm
            Fm = (1+c.beta)*math.log(z/z0)
        elif (z/L<=1):                      # Stable case, z < L    stable
            x = math.log(z/z0)
            Y = c.beta*z/L
            Fm = x+Y
        elif ((z/L>1)&(z0/L<1)):            # Stable case, z > L > z0    vstable
            x = math.log(L/z0)
            Y = (1+c.beta)*math.log(z/L)
            Fm = x+c.beta+Y
        elif (z0/L>1):                      # Stable, L < z0    xstable
            Fm = (1+c.beta)*math.log(z/z0)
        else:
            log.error('   Error in function Fm')
    return Fm

def Fustar(T, Ah, p, Fh, u, z, z0, ustar):
#' Function used in iteration method to solve for ustar.
#' The function used is:
#'  ustar = u*k/Fm(z/L,z0/L)
#' where Fm is the integral form of the PHIm, the adiabatic
#' correction to the logarithmic wind speed profile.
#' Evaluate the function for ustar with this value for L.
    MO = mf.molen(T, Ah, p, ustar, Fh, fluxtype='sensible')
    Fustar = u*c.k/(Fm(z, z0, MO))
    return Fustar

def GetAverageSeriesKeys(cf,ThisOne):
    if incf(cf,ThisOne) and haskey(cf,ThisOne,'AverageSeries'):
        if 'Source' in cf['Variables'][ThisOne]['AverageSeries'].keys():
            alist = ast.literal_eval(cf['Variables'][ThisOne]['AverageSeries']['Source'])
        else:
            log.error('  GetAverageSeriesKeys: key "Source" not in control file AverageSeries section for '+ThisOne)
            alist = []
        if 'standard_name' in cf['Variables'][ThisOne]['AverageSeries'].keys():
            standardname = str(cf['Variables'][ThisOne]['AverageSeries']['standard_name'])
        else:
            standardname = "not defined"
    else:
        standardname = "not defined"
        log.info('  GetAverageSeriesKeys: '+ThisOne+ ' not in control file or it does not have the "AverageSeries" key')
        alist = []
    return alist, standardname

def GetAltName(cf,ds,ThisOne):
    '''
    Check to see if the specified variable name is in the data structure (ds).
    If it is, return the variable name unchanged.
    If it isn't, check the control file to see if an alternate name has been specified
     and return the alternate name if one exists.
    '''
    if ThisOne not in ds.series.keys():
        if ThisOne in cf['Variables'].keys():
            ThisOne = cf['Variables'][ThisOne]['AltVarName']
            if ThisOne not in ds.series.keys():
                log.error('GetAltName: alternate variable name not in ds')
        else:
            log.error('GetAltName: cant find ',ThisOne,' in ds or control file')
    return ThisOne

def GetAltNameFromCF(cf,ThisOne):
    '''
    Get an alternate variable name from the control file.
    '''
    if ThisOne in cf['Variables'].keys():
        if 'AltVarName' in cf['Variables'][ThisOne].keys():
            ThisOne = str(cf['Variables'][ThisOne]['AltVarName'])
        else:
            print 'GetAltNameFromCF: AltVarName key not in control file for '+str(ThisOne)
    else:
        print 'GetAltNameFromCF: '+str(ThisOne)+' not in control file'
    return ThisOne

def GetAttributeDictionary(ds,ThisOne):
    attr = {}
    # if series ThisOne is in the data structure
    if ThisOne in ds.series.keys():
        attr = ds.series[ThisOne]['Attr']
    else:
        attr = MakeAttributeDictionary()
    return attr

def GetcbTicksFromCF(cf,ThisOne):
    '''
    Get colour bar tick labels from the control file.
    '''
    if ThisOne in cf['Variables'].keys():
        if 'Ticks' in cf['Variables'][ThisOne].keys():
            Ticks = eval(cf['Variables'][ThisOne]['Ticks'])
        else:
            print 'GetcbTicksFromCF: Ticks key not in control file for '+str(ThisOne)
    else:
        print 'GetcbTicksFromCF: '+str(ThisOne)+' not in control file'
    return Ticks

def GetRangesFromCF(cf,ThisOne):
    '''
    Get lower and upper range limits from the control file.
    '''
    if ThisOne in cf['Variables'].keys():
        if 'Lower' in cf['Variables'][ThisOne].keys():
            lower = numpy.float64(cf['Variables'][ThisOne]['Lower'])
        else:
            print 'GetRangesFromCF: Lower key not in control file for '+str(ThisOne)
            lower = None
        if 'Upper' in cf['Variables'][ThisOne].keys():
            upper = numpy.float64(cf['Variables'][ThisOne]['Upper'])
        else:
            print 'GetRangesFromCF: Upper key not in control file for '+str(ThisOne)
            upper = None
    else:
        print 'GetRangesFromCF: '+str(ThisOne)+' not in control file'
        lower, upper = None
    return lower, upper

def GetDateIndex(datetimeseries,date,ts=30,default=0,match='exact'):
    # return the index of a date/datetime string in an array of datetime objects
    #  datetimeseries - array of datetime objects
    #  date - a date or date/time string in a format dateutils can parse
    #  default - default value, integer
    try:
        if len(date)!=0:
            i = datetimeseries.index(dateutil.parser.parse(date))
        else:
            if default==-1:
                i = len(datetimeseries)-1
            else:
                i = default
    except ValueError:
        if default==-1:
            i = len(datetimeseries)-1
        else:
            i = default
    if match=='startnextday':
        while abs(datetimeseries[i].hour+numpy.float64(datetimeseries[i].minute)/60-numpy.float64(ts)/60)>c.eps:
            i = i + 1
    if match=='endpreviousday':
        while abs(datetimeseries[i].hour+numpy.float64(datetimeseries[i].minute)/60)>c.eps:
            i = i - 1
    return i

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

def GetPlotTitleFromCF(cf, nFig):
    if 'Plots' in cf:
        if str(nFig) in cf['Plots']:
            if 'Title' in cf['Plots'][str(nFig)]:
                Title = str(cf['Plots'][str(nFig)]['Title'])
            else:
                print 'GetPlotTitleFromCF: Variables key not in control file for plot '+str(nFig)
        else:
            print 'GetPlotTitleFromCF: '+str(nFig)+' key not in Plots section of control file'
    else:
        print 'GetPlotTitleFromCF: Plots key not in control file'
    return Title

def GetPlotVariableNamesFromCF(cf, n):
    if 'Plots' in cf:
        if str(n) in cf['Plots']:
            if 'Variables' in cf['Plots'][str(n)]:
                SeriesList = eval(cf['Plots'][str(n)]['Variables'])
            else:
                print 'GetPlotVariableNamesFromCF: Variables key not in control file for plot '+str(n)
        else:
            print 'GetPlotVariableNamesFromCF: '+str(n)+' key not in Plots section of control file'
    else:
        print 'GetPlotVariableNamesFromCF: Plots key not in control file'
    return SeriesList

def GetSeries(ds,ThisOne,si=0,ei=-1,mode="truncate"):
    """ Returns the data, QC flag and attributes of a series from the data structure."""
    # number of records
    nRecs = numpy.int32(ds.globalattributes['nc_nrecs'])
    # check the series requested is in the data structure
    if ThisOne in ds.series.keys():
        # series is in the data structure
        if isinstance(ds.series[ThisOne]['Data'],list):
            # return a list if the series is a list
            Series = list(ds.series[ThisOne]['Data'])
        elif isinstance(ds.series[ThisOne]['Data'],numpy.ndarray):
            # return a numpy array if series is an array
            Series = ds.series[ThisOne]['Data'].copy()
        # now get the QC flag
        if 'Flag' in ds.series[ThisOne].keys():
            # return the QC flag if it exists
            Flag = ds.series[ThisOne]['Flag'].copy()
            Flag = Flag.astype(numpy.int32)
        else:
            # create a QC flag if one does not exist
            Flag = numpy.zeros(nRecs,dtype=numpy.int32)
        # now get the attribute dictionary
        if "Attr" in ds.series[ThisOne].keys():
            Attr = GetAttributeDictionary(ds,ThisOne)
        else:
            Attr = MakeAttributeDictionary()
    else:
        # make an empty series if the requested series does not exist in the data structure
        Series,Flag,Attr = MakeEmptySeries(ds,ThisOne)
    # tidy up
    if ei==-1: ei = nRecs - 1
    if mode=="truncate":
        # truncate to the requested start and end indices
        si = max(0,si)                  # clip start index at 0
        ei = min(nRecs,ei)              # clip end index to nRecs
        Series = Series[si:ei+1]        # truncate the data
        Flag = Flag[si:ei+1]            # truncate the QC flag
    elif mode=="pad":
        # pad with missing data at the start and/or the end of the series
        if si<0 and ei>nRecs-1:
            # pad at the start
            Series = numpy.append(numpy.float64(c.missing_value)*numpy.ones(abs(si),dtype=numpy.float64),Series)
            Flag = numpy.append(numpy.ones(*abs(si),dtype=numpy.int32),Flag)
            # pad at the end
            Series = numpy.append(Series,numpy.float64(c.missing_value)*numpy.ones((ei-(nRecs-1)),dtype=numpy.float64))
            Flag = numpy.append(Flag,numpy.ones((ei-(nRecs-1)),dtype=numpy.int32))
        elif si<0 and ei<=nRecs-1:
            # pad at the start, truncate the end
            Series = numpy.append(numpy.float64(c.missing_value)*numpy.ones(abs(si),dtype=numpy.float64),Series[:ei+1])
            Flag = numpy.append(numpy.ones(abs(si),dtype=numpy.int32),Flag[:ei+1])
        elif si>=0 and ei>nRecs-1:
            # truncate at the start, pad at the end
            Series = numpy.append(Series[si:],numpy.float64(c.missing_value)*numpy.ones((ei-(nRecs-1)),numpy.float64))
            Flag = numpy.append(Flag[si:],numpy.ones((ei-(nRecs-1)),dtype=numpy.int32))
        elif si>=0 and ei<=nRecs-1:
            # truncate at the start and end
            Series = Series[si:ei+1]
            Flag = Flag[si:ei+1]
        else:
            msg = 'GetSeries: unrecognised combination of si ('+str(si)+') and ei ('+str(ei)+')'
            raise ValueError(msg)
    else:
        raise ValueError('GetSeries: unrecognised mode option'+str(mode))
    return Series,Flag,Attr

def MakeEmptySeries(ds,ThisOne):
    nRecs = numpy.int32(ds.globalattributes['nc_nrecs'])
    Series = numpy.float64(c.missing_value)*numpy.ones(nRecs,dtype=numpy.float64)
    Flag = numpy.ones(nRecs,dtype=numpy.int32)
    Attr = MakeAttributeDictionary()
    return Series,Flag,Attr

def GetSeriesasMA(ds,ThisOne,si=0,ei=-1,mode="truncate"):
    '''
    PURPOSE:
     Returns a data series and the QC flag series from the data structure.
    USAGE:
     data,flag,attr = qcutils.GetSeriesasMA(ds,label,si=0,ei=-1)
    where the arguments are;
      ds    - the data structure (dict)
      label - label of the data series in ds (string)
      si    - start index (integer), default 0
      ei    - end index (integer), default -1
    and the returned values are;
      data - values for the requested series in ds
             (numpy masked array, float64)
      flag - QC flag for the requested series in ds
             (numpy masked array, int32)
      attr - attribute dictionary for series
    EXAMPLE:
     The code snippet below will return the incoming shortwave data values
     (Fsd) and the associated QC flag (f) as numpy masked arrays;
      ds = qcio.nc_read_series("HowardSprings_2011_L3.nc")
      Fsd,f,a = qcutils.GetSeriesasMA(ds,"Fsd")
    AUTHOR: PRI
    '''
    Series,Flag,Attr = GetSeries(ds,ThisOne,si=si,ei=ei,mode=mode)
    Series,WasND = SeriestoMA(Series)
    return Series,Flag,Attr

def GetUnitsFromds(ds, ThisOne):
    units = ds.series[ThisOne]['Attr']['units']
    return units

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

def get_coverage_groups(ds,rad=None,met=None,flux=None,soil=None):
    level = str(ds.globalattributes['nc_level'])
    if 'L1' in level: level=str('L1')
    if 'L2' in level: level=str('L2')
    if 'L3' in level: level=str('L3')
    if 'L4' in level: level=str('L4')
    if 'L5' in level: level=str('L5')
    if 'L6' in level: level=str('L6')
    rad = ['Fsd','Fsu','Fld','Flu','Fn']
    met = ['q','Precip','ps','Ta','Ws','Wd']
    flux = ['Fh','Fe','Fc','ustar']
    soil = ['Fg','Ts','Sws']
    for ThisGroup, ThisLabel in zip([rad,met,flux,soil],['radiation','meteorology','flux','soil']):
        sum_coverage = numpy.float64(0); count = numpy.float64(0)
        for ThisOne in ThisGroup:
            if ThisOne in ds.series.keys():
                sum_coverage = sum_coverage + numpy.float64(ds.series[ThisOne]['Attr']['coverage_'+level])
                count = count + 1
        if count!=0:
            coverage_group = sum_coverage/count
        else:
            coverage_group = 0
        ds.globalattributes['coverage_'+ThisLabel+'_'+level] = str('%d'%coverage_group)

def get_coverage_individual(ds):
    level = str(ds.globalattributes['nc_level'])
    if 'L1' in level: level=str('L1')
    if 'L2' in level: level=str('L2')
    if 'L3' in level: level=str('L3')
    if 'L4' in level: level=str('L4')
    if 'L5' in level: level=str('L5')
    if 'L6' in level: level=str('L6')
    SeriesList = ds.series.keys()
    for ThisOne in ["DateTime","DateTime_UTC"]:
        if ThisOne in SeriesList: SeriesList.remove(ThisOne)
    for ThisOne in SeriesList:
        num_good = len(numpy.where(abs(ds.series[ThisOne]['Data']-numpy.float64(c.missing_value))>c.eps)[0])
        coverage = 100*numpy.float64(num_good)/numpy.float64(ds.globalattributes['nc_nrecs'])
        ds.series[ThisOne]['Attr']['coverage_'+level] = str('%d'%coverage)

def get_datetimefromxldate(ds):
    ''' Creates a series of Python datetime objects from the Excel date read from the Excel file.
        Thanks to John Machin for the quick and dirty code
         see http://stackoverflow.com/questions/1108428/how-do-i-read-a-date-in-excel-format-in-python'''

    log.info(' Getting the Python datetime series from the Excel datetime')
    xldate = ds.series['xlDateTime']['Data']
    nRecs = len(ds.series['xlDateTime']['Data'])
    datemode = numpy.int32(ds.globalattributes['xl_datemode'])
    ds.series[unicode('DateTime')] = {}
    ds.series['DateTime']['Data'] = [None]*nRecs
    basedate = datetime.datetime(1899, 12, 30)
    for i in range(nRecs):
        ds.series['DateTime']['Data'][i] = basedate + datetime.timedelta(days=xldate[i] + 1462 * datemode)
    ds.series['DateTime']['Flag'] = numpy.zeros(nRecs,dtype=numpy.int32)
    ds.series['DateTime']['Attr'] = {}
    ds.series['DateTime']['Attr']['long_name'] = 'Datetime in local timezone'
    ds.series['DateTime']['Attr']['units'] = 'None'

def get_datetimefromymdhms(ds):
    ''' Creates a series of Python datetime objects from the year, month,
    day, hour, minute and second series stored in the netCDF file.'''
    SeriesList = ds.series.keys()
    if 'Year' not in SeriesList or 'Month' not in SeriesList or 'Day' not in SeriesList or 'Hour' not in SeriesList or 'Minute' not in SeriesList or 'Second' not in SeriesList:
        log.info(' get_datetimefromymdhms: unable to find all datetime fields required')
        return
    log.info(' Getting the date and time series')
    nRecs = get_nrecs(ds)
    ds.series[unicode('DateTime')] = {}
    ds.series['DateTime']['Data'] = [None]*nRecs
    for i in range(nRecs):
        ds.series['DateTime']['Data'][i] = datetime.datetime(numpy.int32(ds.series['Year']['Data'][i]),
                                                       numpy.int32(ds.series['Month']['Data'][i]),
                                                       numpy.int32(ds.series['Day']['Data'][i]),
                                                       numpy.int32(ds.series['Hour']['Data'][i]),
                                                       numpy.int32(ds.series['Minute']['Data'][i]),
                                                       numpy.int32(ds.series['Second']['Data'][i]))
    ds.series['DateTime']['Flag'] = numpy.zeros(nRecs,dtype=numpy.int32)
    ds.series['DateTime']['Attr'] = {}
    ds.series['DateTime']['Attr']['long_name'] = 'Date-time object'
    ds.series['DateTime']['Attr']['units'] = 'None'

def get_nrecs(ds):
    if 'nc_nrecs' in ds.globalattributes.keys():
        nRecs = numpy.int32(ds.globalattributes['nc_nrecs'])
    elif 'NumRecs' in ds.globalattributes.keys():
        nRecs = numpy.int32(ds.globalattributes['NumRecs'])
    else:
        nRecs = len(ds.series[SeriesList[0]]['Data'])
    return nRecs

def get_UTCfromlocaltime(ds):
    '''
    PURPOSE:
     Creates a UTC datetime series in the data structure from the
     local datetime series.
    USAGE:
     qcutils.get_UTCfromlocaltime(ds)
    ASSUMPTIONS:
     No daylight savings used in the local datetime
    AUTHOR: PRI
    '''
    # check the time_zone global attribute is set, we cant continue without it
    if "time_zone" not in ds.globalattributes.keys():
        log.error("get_UTCfromlocaltime: time_zone not in global attributes")
        return
    log.info(' Getting the UTC datetime from the local datetime')
    # get the number of records
    nRecs = len(ds.series['xlDateTime']['Data'])
    # get the time zone
    tz = ds.globalattributes["time_zone"]
    # create a timezone object
    loc_tz = pytz.timezone(tz)
    # local pointer to the datetime series in ds
    ldt = ds.series["DateTime"]["Data"]
    # localise the datetime
    ldt_loc = [loc_tz.localize(dt) for dt in ldt]
    # convert to UTC
    ldt_utc = [dt.astimezone(pytz.utc) for dt in ldt_loc]
    # put the UTC datetime into the data structure
    ds.series[unicode("DateTime_UTC")] = {}
    ds.series["DateTime_UTC"]["Data"] = ldt_utc
    ds.series['DateTime_UTC']["Flag"] = numpy.zeros(nRecs,dtype=numpy.int32)
    ds.series['DateTime_UTC']["Attr"] = {}
    ds.series['DateTime_UTC']["Attr"]["long_name"] = "Datetime as UTC"
    ds.series['DateTime_UTC']["Attr"]["units"] = "None"

def get_diurnalstats(DecHour,Data,dt):
    nInts = 24*numpy.int32((60/dt)+0.5)
    Hr = numpy.zeros(nInts,dtype=numpy.float64) + numpy.float64(c.missing_value)
    Av = numpy.zeros(nInts,dtype=numpy.float64) + numpy.float64(c.missing_value)
    Sd = numpy.zeros(nInts,dtype=numpy.float64) + numpy.float64(c.missing_value)
    Mx = numpy.zeros(nInts,dtype=numpy.float64) + numpy.float64(c.missing_value)
    Mn = numpy.zeros(nInts,dtype=numpy.float64) + numpy.float64(c.missing_value)
    Num = numpy.zeros(nInts,dtype=numpy.float64) + numpy.float64(c.missing_value)
    for i in range(nInts):
        Hr[i] = numpy.float64(i)*dt/60.
        li = numpy.where((abs(DecHour-Hr[i])<c.eps)&(abs(Data-numpy.float64(c.missing_value))>c.eps))
        Num[i] = numpy.size(li)
        if numpy.size(li)!=0:
            Av[i] = numpy.mean(Data[li])
            Sd[i] = numpy.std(Data[li])
            Mx[i] = numpy.max(Data[li])
            Mn[i] = numpy.min(Data[li])
    return Num, Hr, Av, Sd, Mx, Mn

def get_xldate_from_datetime(dt,datemode=0):
    '''
    PURPOSE:
     Returns a list of xldatetime (floating point number represent decimal days
     since 00:00 1/1/1900) from a list of Python datetimes
    USAGE:
    datemode = ds.globalattributes["xl_datemode"]
     xl_date = qcutils.get_xldate_from_datetime(py_date,datemode=datemode)
     ds.series["xlDateTime"] = numpy.array(xl_date)
    ASSUMPTIONS:
     datemode is set to 0 (1900 date system), use on Windows assumed
    AUTHOR: PRI
    '''
    xldate = [xlrd.xldate.xldate_from_datetime_tuple((dt[i].year,
                                                      dt[i].month,
                                                      dt[i].day,
                                                      dt[i].hour,
                                                      dt[i].minute,
                                                      dt[i].second),
                                                      datemode) for i in range(0,len(dt))]
    return xldate

def get_ymdhms_from_datetime(ds):
    '''
    PURPOSE:
     Gets the year, month, day, hour, minute and second from a list of
     Python datetimes.  The Python datetime series is read from
     the input data structure and the results are written back to the
     data structure.
    USAGE:
     qcutils.get_ymdhms_from_datetime(ds)
    ASSUMPTIONS:
     None
    AUTHOR: PRI
    '''
    nRecs = numpy.int32(ds.globalattributes["nc_nrecs"])
    dt = ds.series["DateTime"]["Data"]
    flag = numpy.zeros(nRecs,dtype=numpy.int32)
    Year = [dt[i].year for i in range(0,nRecs)]
    Month = [dt[i].month for i in range(0,nRecs)]
    Day = [dt[i].day for i in range(0,nRecs)]
    Hour = [dt[i].hour for i in range(0,nRecs)]
    Minute = [dt[i].minute for i in range(0,nRecs)]
    Second = [dt[i].second for i in range(0,nRecs)]
    Hdh = [numpy.float64(Hour[i])+numpy.float64(Minute[i])/60. for i in range(0,nRecs)]
    Ddd = [(dt[i] - datetime.datetime(Year[i],1,1)).days+1+Hdh[i]/24. for i in range(0,nRecs)]
    CreateSeries(ds,'Year',Year,Flag=flag,Attr=MakeAttributeDictionary(long_name='Year',units='none'))
    CreateSeries(ds,'Month',Month,Flag=flag,Attr=MakeAttributeDictionary(long_name='Month',units='none'))
    CreateSeries(ds,'Day',Day,Flag=flag,Attr=MakeAttributeDictionary(long_name='Day',units='none'))
    CreateSeries(ds,'Hour',Hour,Flag=flag,Attr=MakeAttributeDictionary(long_name='Hour',units='none'))
    CreateSeries(ds,'Minute',Minute,Flag=flag,Attr=MakeAttributeDictionary(long_name='Minute',units='none'))
    CreateSeries(ds,'Second',Second,Flag=flag,Attr=MakeAttributeDictionary(long_name='Second',units='none'))
    CreateSeries(ds,'Hdh',Hdh,Flag=flag,Attr=MakeAttributeDictionary(long_name='Decimal hour of the day',units='none'))
    CreateSeries(ds,'Ddd',Ddd,Flag=flag,Attr=MakeAttributeDictionary(long_name='Decimal day of the year',units='none'))

def get_ymdhmsfromxldate(ds):
    """
        Gets year, month, day, hour, and if available seconds, from
        excel-formatted Timestamp
        
        Usage qcts.get_ymdhmsfromxldate(ds)
        cf: control file
        ds: data structure
        """
    log.info(' Getting date and time variables')
    # get the date mode of the original Excel datetime
    datemode = numpy.int32(ds.globalattributes['xl_datemode'])
    nRecs = len(ds.series['xlDateTime']['Data'])
    Year = numpy.array([c.missing_value]*nRecs,numpy.int32)
    Month = numpy.array([c.missing_value]*nRecs,numpy.int32)
    Day = numpy.array([c.missing_value]*nRecs,numpy.int32)
    Hour = numpy.array([c.missing_value]*nRecs,numpy.int32)
    Minute = numpy.array([c.missing_value]*nRecs,numpy.int32)
    Second = numpy.array([c.missing_value]*nRecs,numpy.int32)
    Hdh = numpy.array([c.missing_value]*nRecs,numpy.float64)
    Ddd = numpy.array([c.missing_value]*nRecs,numpy.float64)
    flag = numpy.zeros(nRecs,dtype=numpy.int32)
    for i in range(nRecs):
        DateTuple = xlrd.xldate_as_tuple(ds.series['xlDateTime']['Data'][i],datemode)
        Year[i] = numpy.int32(DateTuple[0])
        Month[i] = numpy.int32(DateTuple[1])
        Day[i] = numpy.int32(DateTuple[2])
        Hour[i] = numpy.int32(DateTuple[3])
        Minute[i] = numpy.int32(DateTuple[4])
        Second[i] = numpy.int32(DateTuple[5])
        Hdh[i] = numpy.float64(DateTuple[3])+numpy.float64(DateTuple[4])/60.
        Ddd[i] = ds.series['xlDateTime']['Data'][i] - xlrd.xldate.xldate_from_date_tuple((Year[i],1,1),datemode) + 1
    CreateSeries(ds,'Year',Year,Flag=flag,Attr=MakeAttributeDictionary(long_name='Year',units='none'))
    CreateSeries(ds,'Month',Month,Flag=flag,Attr=MakeAttributeDictionary(long_name='Month',units='none'))
    CreateSeries(ds,'Day',Day,Flag=flag,Attr=MakeAttributeDictionary(long_name='Day',units='none'))
    CreateSeries(ds,'Hour',Hour,Flag=flag,Attr=MakeAttributeDictionary(long_name='Hour',units='none'))
    CreateSeries(ds,'Minute',Minute,Flag=flag,Attr=MakeAttributeDictionary(long_name='Minute',units='none'))
    CreateSeries(ds,'Second',Second,Flag=flag,Attr=MakeAttributeDictionary(long_name='Second',units='none'))
    CreateSeries(ds,'Hdh',Hdh,Flag=flag,Attr=MakeAttributeDictionary(long_name='Decimal hour of the day',units='none'))
    CreateSeries(ds,'Ddd',Ddd,Flag=flag,Attr=MakeAttributeDictionary(long_name='Decimal day of the year',units='none'))

def haskey(cf,ThisOne,key):
    return key in cf['Variables'][ThisOne].keys()

def incf(cf,ThisOne):
    return ThisOne in cf['Variables'].keys()

def MakeAttributeDictionary(**kwargs):
    default_list = ['ancillary_variables','height','instrument','serial_number','standard_name','long_name','units']
    attr = {}
    for item in kwargs:
        attr[item] = kwargs.get(item,'not defined')
        if item in default_list: default_list.remove(item)
    if len(default_list)!=0:
        for item in default_list: attr[item] = 'not defined'
    attr["missing_value"] = c.missing_value
    return attr

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
                    goodindex = numpy.where((numpy.mod(tmp_flag,10)==0) & (numpy.mod(flag,10)==0))    # find the elements with flag = 0, 10, 20 etc
                    badindex1 = numpy.where((numpy.mod(tmp_flag,10)==0) & (numpy.mod(flag,10)!=0))
                    badindex2 = numpy.where((numpy.mod(tmp_flag,10)!=0) & (numpy.mod(flag,10)==0))
                    badindex3 = numpy.where((numpy.mod(tmp_flag,10)!=0) & (numpy.mod(flag,10)!=0))
                    tmp_flag1 = numpy.zeros(len(tmp_flag),dtype=numpy.int32) + tmp_flag
                    flag2 = numpy.zeros(len(flag),dtype=numpy.int32) + flag
                    tmp_flag[badindex1] = 0                               # set them all to 0
                    flag2[badindex2] = 0                               # set them all to 0
                    flag[badindex1] = numpy.maximum(tmp_flag[badindex1],flag[badindex1])               # now take the maximum
                    flag[badindex2] = numpy.maximum(tmp_flag1[badindex2],flag2[badindex2])
                    flag[badindex3] = numpy.maximum(tmp_flag1[badindex3],flag[badindex3])
                    flag[goodindex] = numpy.maximum(tmp_flag1[goodindex],flag[goodindex])
            else:
                log.error('  MakeQCFlag: series '+ThisOne+' not in ds.series')
    return flag

def MAtoSeries(Series):
    """
    Convert a masked array to a numpy ndarray with masked elements set to c.missing_value.
    Useage:
     Series, WasMA = MAtoSeries(Series)
     where:
      Series (input)    is the data series to be converted.
      WasMA  (returned) is a logical, True if the input series was a masked array.
      Series (output)   is the input series convered to an ndarray with c.missing_value values
                        for missing data.
    """
    WasMA = False
    if numpy.ma.isMA(Series):
        WasMA = True
        Series = numpy.ma.filled(Series,numpy.float64(c.missing_value))
    return Series, WasMA

def nxMom_nxScalar_alpha(zoL):
    nRecs = numpy.size(zoL)
    nxMom = numpy.ma.ones(nRecs) * 0.079
    nxScalar = numpy.ma.ones(nRecs) * 0.085
    alpha = numpy.ma.ones(nRecs) * 0.925
    #  get the index of stable conditions
    stable = numpy.ma.where(zoL>0)[0]
    #  now set the series to their stable values
    nxMom[stable] = 0.079 * (1 + 7.9 * zoL[stable]) ** 0.75
    nxScalar[stable] = 2.0 - 1.915 / (1 + 0.5 * zoL[stable])
    alpha[stable] = 1
    return nxMom, nxScalar, alpha

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
    >>> y = numpy.array([1,c.missing_value,3])
    >>> y = numpy.ma.masked_where(y==c.missing_value,y)
    >>> qcutils.polyval(p,y)
    masked_array(data = [2 -- 6],
                 mask = [False True False],
                 fill_value = 999999)
    """
    y = 0
    for i in range(len(p)):
        y = x*y + p[i]
    return y

def r(b, p, alpha):
    """
    Function to calculate the r coeficient of the Massman frequency correction.
    """
    r = ((b ** alpha) / (b ** alpha + 1)) * \
           ((b ** alpha) / (b ** alpha + p ** alpha)) * \
           (1 / (p ** alpha + 1))
    return r

def SeriestoMA(Series):
    """
    Convert a numpy ndarray to a masked array.
    Useage:
     Series, WasND = SeriestoMA(Series)
     where:
      Series (input)    is the data series to be converted.
      WasND  (returned) is a logical, True if the input series was an ndarray
      Series (output)   is the input series convered to a masked array.
    """
    WasND = False
    if not numpy.ma.isMA(Series):
        WasND = True
        Series = numpy.ma.masked_where(abs(Series-numpy.float64(c.missing_value))<c.eps,Series)
    return Series, WasND

def SetUnitsInds(ds, ThisOne, units):
    ds.series[ThisOne]['Attr']['units'] = units

def prepOzFluxVars(cf,ds):
    invars = ['Ux', 'Uy', 'Uz', 'UxUx', 'UxUy', 'UxA', 'UxC', 'UxT', 'UyA', 'UyC', 'UyT', 'UyUy', 'AhAh', 'CcCc', 'UxUz', 'UyUz', 'UzA', 'UzC', 'UzT', 'UzUz', 'CE_umol', 'CcAh', 'Cc_7500_Sd', 'Ah_7500_Sd', 'TvAh']
    outvars = ['u', 'v', 'w', 'uu', 'uv', 'uA', 'uC', 'uT', 'vA', 'vC', 'vT', 'vv', 'AhAh', 'CcCc', 'uw', 'vw', 'wA', 'wC', 'wT', 'ww', 'CE', 'CA', 'C_Sd', 'Ah_Sd', 'TA']
    for i in range(len(invars)):
        if invars[i] in ds.series.keys() and outvars[i] not in ds.series.keys():
            attr = MakeAttributeDictionary(long_name=ds.series[invars[i]]['Attr']['long_name'],units=ds.series[invars[i]]['Attr']['units'],standard_name=ds.series[invars[i]]['Attr']['standard_name'])
            CreateSeries(ds,outvars[i],ds.series[invars[i]]['Data'],Flag=ds.series[invars[i]]['Flag'],Attr=attr)

def startlog(loggername,loggerfile):
    logger = logging.getLogger(loggername)
    logger.setLevel(logging.DEBUG)
    fh = logging.FileHandler(loggerfile)
    fh.setLevel(logging.DEBUG)
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s %(name)-8s %(levelname)-6s %(message)s', '%d-%m-%y %H:%M')
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)
    logger.addHandler(fh)
    logger.addHandler(ch)
    return logger

def Wegstein(T, Ah, p, Fh, u, z, z0):

    NumIters = 50
    SolveErr = numpy.zeros(len(u),dtype=numpy.float64) + numpy.float64(0.001)
    RangeErr = numpy.zeros(len(u),dtype=numpy.float64)
 
    FirstEst =  u*c.k/math.log(z/z0)
    ustar = Fustar(T, Ah, p, Fh, u, z, z0, FirstEst)
    Inc = ustar-FirstEst
    IncDiv = -Inc
    Residual = ustar-Fustar(T, Ah, p, Fh, u, z, z0, ustar)
 
    for j in range(len(u)):
        i = 1
        while (i<NumIters)&(numpy.float64(Residual[j])>numpy.float64(SolveErr[j])):
            IncDiv[j] = (IncDiv[j]/Residual[j])-1
            if (IncDiv[j] == 0):
                log.info('  Wegstein: IncDiv equals 0')
                log.info('   T: '+str(T[j])+', Ah: '+str(Ah[j])+', p: '+str(p[j])+', Fh: '+str(Fh[j])+', u: '+str(u[j]))
                ustar[j] = u[j]*c.k/numpy.log(z/z0)
                break
            Inc[j] = Inc[j]/IncDiv[j]
            ustar[j] = ustar[j]+Inc[j]
            IncDiv[j] = Residual[j]
            Residual[j] = ustar[j]-Fustar(T[j], Ah[j], p[j], Fh[j], u[j], z, z0, ustar[j])
            if (abs(ustar[j])<=1):
                RangeErr[j] = SolveErr[j]
            else:
                RangeErr[j] = SolveErr[j]*abs(ustar[j])
            if (abs(Inc[j])<=RangeErr[j]):
                if (abs(Residual[j])<=10*RangeErr[j]):
                    break
            i = i + 1
        if (i==NumIters):
            log.info('  Wegstein: did not converge')
            log.info('   T: '+str(T[j])+', Ah: '+str(Ah[j])+', p: '+str(p[j])+', Fh: '+str(Fh[j])+', u: '+str(u[j]))
            ustar[j] = u[j]*c.k/numpy.log(z/z0)
    return ustar
