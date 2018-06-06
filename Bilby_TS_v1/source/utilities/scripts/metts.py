#    metts.py
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

"""
    Partition Data Function Module
    Used to perform the tasks queued by pls.py
    """

import sys
import ast
import constants as c
import numpy
import metutils
import time
import xlrd
import xlwt
import logging

log = logging.getLogger('envelope.ts')

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
        ds.globalattributes['Flag0'] = 'Good data'
        ds.globalattributes['Flag1'] = 'QA/QC: -9999 in level 1 dataset'
        ds.globalattributes['Flag2'] = 'QA/QC: L2 Range Check'
        ds.globalattributes['Flag3'] = 'QA/QC: CSAT Diagnostic'
        ds.globalattributes['Flag4'] = 'QA/QC: LI7500 Diagnostic'
        ds.globalattributes['Flag5'] = 'QA/QC: L2 Diurnal SD Check'
        ds.globalattributes['Flag6'] = 'QA/QC: Excluded Dates'
        ds.globalattributes['Flag7'] = 'QA/QC: Excluded Hours'
        ds.globalattributes['Flag10'] = 'Corrections: Apply Linear'
        ds.globalattributes['Flag11'] = 'Corrections/Combinations: Coordinate Rotation (Ux, Uy, Uz, UxT, UyT, UzT, UxA, UyA, UzA, UxC, UyC, UzC, UxUz, UxUx, UxUy, UyUz, UxUy, UyUy)'
        ds.globalattributes['Flag12'] = 'Corrections/Combinations: Massman Frequency Attenuation Correction (Coord Rotation, Tv_CSAT, Ah_HMP, ps)'
        ds.globalattributes['Flag13'] = 'Corrections/Combinations: Virtual to Actual Fh (Coord Rotation, Massman, Ta_HMP)'
        ds.globalattributes['Flag14'] = 'Corrections/Combinations: WPL correction for flux effects on density measurements (Coord Rotation, Massman, Fhv to Fh, Cc_7500_Av)'
        ds.globalattributes['Flag15'] = 'Corrections/Combinations: Ta from Tv'
        ds.globalattributes['Flag16'] = 'Corrections/Combinations: L3 Range Check'
        ds.globalattributes['Flag17'] = 'Corrections/Combinations: L3 Diurnal SD Check'
        ds.globalattributes['Flag18'] = 'Corrections/Combinations: u* filter'
        ds.globalattributes['Flag19'] = 'Corrections/Combinations: Gap coordination'
        ds.globalattributes['Flag30'] = 'GapFilling: Flux Gap Filled by ANN (SOLO)'
        ds.globalattributes['Flag31'] = 'GapFilling: Flux Gap not Filled by ANN'
        ds.globalattributes['Flag32'] = 'GapFilling: Met Gap Filled from Climatology'
        ds.globalattributes['Flag33'] = 'GapFilling: Gap Filled from Ratios'
        ds.globalattributes['Flag34'] = 'GapFilling: Gap Filled by Interpolation'
        ds.globalattributes['Flag35'] = 'GapFilling: Gap Filled by Replacement'
        ds.globalattributes['Flag36'] = 'GapFilling: u* from Fh'
        ds.globalattributes['Flag37'] = 'GapFilling: u* not from Fh'
        ds.globalattributes['Flag38'] = 'GapFilling: L4 Range Check'
        ds.globalattributes['Flag39'] = 'GapFilling: L4 Diurnal SD Check'
        ds.globalattributes['Flag51'] = 'albedo: bad Fsd < threshold (290 W/m2 default) only if bad time flag (31) not set'
        ds.globalattributes['Flag52'] = 'albedo: bad time flag (not midday 10.00 to 14.00)'
        ds.globalattributes['Flag61'] = 'Penman-Monteith: bad rst (rst < 0) only if bad Uavg (35), bad Fe (33) and bad Fsd (34) flags not set'
        ds.globalattributes['Flag62'] = 'Penman-Monteith: bad Fe < threshold (0 W/m2 default) only if bad Fsd (34) flag not set'
        ds.globalattributes['Flag63'] = 'Penman-Monteith: bad Fsd < threshold (10 W/m2 default)'
        ds.globalattributes['Flag64'] = 'Penman-Monteith: Uavg == 0 (undefined aerodynamic resistance under calm conditions) only if bad Fe (33) and bad Fsd (34) flags not set'
        ds.globalattributes['Flag70'] = 'Partitioning Night: Re computed from exponential temperature response curves'
        ds.globalattributes['Flag80'] = 'Partitioning Day: GPP/Re computed from light-response curves, GPP = Re - Fc'
        ds.globalattributes['Flag81'] = 'Partitioning Day: GPP night mask'
        ds.globalattributes['Flag82'] = 'Partitioning Day: Fc > Re, GPP = 0, Re = Fc'
    for ThisOne in ds.series.keys():
        if ThisOne in cf['Variables']:
            if 'Attr' in cf['Variables'][ThisOne].keys():
                ds.series[ThisOne]['Attr'] = {}
                for attr in cf['Variables'][ThisOne]['Attr'].keys():
                    ds.series[ThisOne]['Attr'][attr] = cf['Variables'][ThisOne]['Attr'][attr]

def do_functions(cf,ds):
    log.info(' Resolving functions given in control file')
    for ThisOne in cf['Variables'].keys():
        if 'Function' in cf['Variables'][ThisOne].keys():
            ds.series[ThisOne] = {}
            FunctionList = cf['Variables'][ThisOne]['Function'].keys()
            if len(FunctionList) == 1:
                i = 0
                if 'Square' in cf['Variables'][ThisOne]['Function'][str(i)].keys() and 'Parent' in cf['Variables'][ThisOne]['Function'][str(i)]['Square'].keys():
                    Parent = cf['Variables'][ThisOne]['Function'][str(i)]['Square']['Parent']
                    ds.series[ThisOne]['Data'] = qcts.Square(ds.series[Parent]['Data'])
                    nRecs = numpy.size(ds.series[ThisOne]['Data'])
                    if 'Flag' not in ds.series[ThisOne].keys():
                        ds.series[ThisOne]['Flag'] = numpy.zeros(nRecs,int)
                        if 'Flag' in ds.series[Parent]:
                            ds.series[ThisOne]['Flag'] = ds.series[Parent]['Flag']
                        else:
                            ds.series[ThisOne]['Flag'] = numpy.zeros(nRecs,int)
                elif 'SquareRoot' in cf['Variables'][ThisOne]['Function'][str(i)].keys() and 'Parent' in cf['Variables'][ThisOne]['Function'][str(i)]['Square'].keys():
                    Parent = cf['Variables'][ThisOne]['Function'][str(i)]['Square']['Parent']
                    ds.series[ThisOne]['Data'] = qcts.SquareRoot(ds.series[Parent]['Data'])
                    nRecs = numpy.size(ds.series[ThisOne]['Data'])
                    if 'Flag' not in ds.series[ThisOne].keys():
                        ds.series[ThisOne]['Flag'] = numpy.zeros(nRecs,int)
                        if 'Flag' in ds.series[Parent]:
                            ds.series[ThisOne]['Flag'] = ds.series[Parent]['Flag']
                        else:
                            ds.series[ThisOne]['Flag'] = numpy.zeros(nRecs,int)
                else:
                    log.error ('Function missing or unknown for variable'+ThisOne)
                    return
            else:
                for i in range(len(FunctionList)):
                    if 'Square' in cf['Variables'][ThisOne]['Function'][str(i)].keys() and 'Parent' in cf['Variables'][ThisOne]['Function'][str(i)]['Square'].keys():
                        Parent = cf['Variables'][ThisOne]['Function'][str(i)]['Square']['Parent']
                        ds.series[ThisOne]['Data'] = qcts.Square(ds.series[Parent]['Data'])
                        nRecs = numpy.size(ds.series[ThisOne]['Data'])
                        if 'Flag' not in ds.series[ThisOne].keys():
                            ds.series[ThisOne]['Flag'] = numpy.zeros(nRecs,int)
                            if 'Flag' in ds.series[Parent]:
                                ds.series[ThisOne]['Flag'] = ds.series[Parent]['Flag']
                            else:
                                ds.series[ThisOne]['Flag'] = numpy.zeros(nRecs,int)
                    elif 'SquareRoot' in cf['Variables'][ThisOne]['Function'][str(i)].keys() and 'Parent' in cf['Variables'][ThisOne]['Function'][str(i)]['Square'].keys():
                        Parent = cf['Variables'][ThisOne]['Function'][str(i)]['Square']['Parent']
                        ds.series[ThisOne]['Data'] = qcts.SquareRoot(ds.series[Parent]['Data'])
                        nRecs = numpy.size(ds.series[ThisOne]['Data'])
                        if 'Flag' not in ds.series[ThisOne].keys():
                            ds.series[ThisOne]['Flag'] = numpy.zeros(nRecs,int)
                            if 'Flag' in ds.series[Parent]:
                                ds.series[ThisOne]['Flag'] = ds.series[Parent]['Flag']
                            else:
                                ds.series[ThisOne]['Flag'] = numpy.zeros(nRecs,int)
                    else:
                        log.error ('Function missing or unknown for variable'+ThisOne)
                        return

def get_nightsums(Data):
    """
        Get nightly sums and averages on nights when no 30-min observations are missing.
        Nights with missing observations return a value of -9999
        Values returned are sample size (Num), sums (Sum) and average (Av)
        
        Usage qcts.get_nightsums(Data)
        Data: 1-day dataset
        """
    li = numpy.ma.where(Data.mask == False)[0]
    Num = numpy.size(li)
    if Num == 0:
        Sum = -9999
        Av = -9999
    else:
        x = 0
        for i in range(len(Data)):
            if Data.mask[i] == True:
                x = x + 1
        
        if x == 0:
            Sum = numpy.ma.sum(Data[li])
            Av = numpy.ma.mean(Data[li])
        else:
            Sum = -9999
            Av = -9999
    
    return Num, Sum, Av

def get_yearmonthdayhourminutesecond(cf,ds):
    """
        Gets year, month, day, hour, and if available seconds, from
        excel-formatted Timestamp
        
        Usage qcts.get_yearmonthdayhourminutesecond(cf,ds)
        cf: control file
        ds: data structure
        """
    log.info(' Getting date and time variables')
    nRecs = len(ds.series['xlDateTime']['Data'])
    for ThisOne in ['Year','Month','Day','Hour','Minute','Second']:
        ds.series[ThisOne] = {}
        ds.series[ThisOne]['Data'] = numpy.array([-9999]*nRecs,numpy.int32)
    ds.series['Hdh'] = {}
    ds.series['Hdh']['Data'] = numpy.array([-9999]*nRecs,numpy.float32)
    for i in range(nRecs):
        if cf['General']['Platform'] == 'Mac':
            DateTuple = xlrd.xldate_as_tuple(ds.series['xlDateTime']['Data'][i],1)
        else:
            DateTuple = xlrd.xldate_as_tuple(ds.series['xlDateTime']['Data'][i],0)
        ds.series['Year']['Data'][i] = int(DateTuple[0])
        ds.series['Month']['Data'][i] = int(DateTuple[1])
        ds.series['Day']['Data'][i] = int(DateTuple[2])
        ds.series['Hour']['Data'][i] = int(DateTuple[3])
        ds.series['Minute']['Data'][i] = int(DateTuple[4])
        ds.series['Second']['Data'][i] = int(DateTuple[5])
        ds.series['Hdh']['Data'][i] = float(DateTuple[3])+float(DateTuple[4])/60.

def DBinFilters(cf,ds):
    '''Prepare data into bins for meteorological envelope calculation
    '''
    log.info('Preparing data for meteorological envelope calculation')
    if metutils.cfkeycheck(cf,Base='Params',ThisOne='dependent'):
        InVar = cf['Params']['dependent']
    else:
        InVar = 'Fc'
    
    if (InVar == 'Fc_wpl' or InVar == 'Fc'):
        Fcmg,f = metutils.GetSeriesasMA(ds,InVar)
        Fc1 = ((0-Fcmg) * (10 ** 6)) / (1000 * 44)
    elif (InVar == 'GPP'):
        Fc1,f = metutils.GetSeriesasMA(ds,InVar)
    
    Ts,f = metutils.GetSeriesasMA(ds,'Ts')
    Sws,f = metutils.GetSeriesasMA(ds,'Sws')
    Fsd,f = metutils.GetSeriesasMA(ds,'Fsd')
    D,f = metutils.GetSeriesasMA(ds,'VPD')
    nRecs = len(Fc1)
    
    D0 = numpy.ma.masked_where(D > 0.65,Fc1)
    Sws2D0 = numpy.ma.masked_where(Sws > 0.05,D0)
    Sws5D0 = numpy.ma.masked_where((Sws < 0.05) | (Sws > 0.08),D0)
    Sws8D0 = numpy.ma.masked_where((Sws < 0.08) | (Sws > 0.11),D0)
    Sws11D0 = numpy.ma.masked_where((Sws < 0.11) | (Sws > 0.14),D0)
    Sws14D0 = numpy.ma.masked_where((Sws < 0.14) | (Sws > 0.17),D0)
    Sws17D0 = numpy.ma.masked_where((Sws < 0.17) | (Sws > 0.20),D0)
    Sws20D0 = numpy.ma.masked_where((Sws < 0.20) | (Sws > 0.23),D0)
    Sws23D0 = numpy.ma.masked_where((Sws < 0.23) | (Sws > 0.26),D0)
    Sws26D0 = numpy.ma.masked_where((Sws < 0.26) | (Sws > 0.29),D0)
    Sws29D0 = numpy.ma.masked_where(Sws < 0.29,D0)
    
    metutils.CreateSeries(ds,'Sws2D0',Sws2D0,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D0',Sws5D0,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws8D0',Sws8D0,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws11D0',Sws11D0,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws14D0',Sws14D0,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws17D0',Sws17D0,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws20D0',Sws20D0,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws23D0',Sws23D0,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws26D0',Sws26D0,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws29D0',Sws29D0,FList=[InVar],Descr='',Units='umol/(m2 s)')
        
    D_65 = numpy.ma.masked_where((D < 0.65) | (D > 1.3),Fc1)
    Sws2D_65 = numpy.ma.masked_where(Sws > 0.05,D_65)
    Sws5D_65 = numpy.ma.masked_where((Sws < 0.05) | (Sws > 0.08),D_65)
    Sws8D_65 = numpy.ma.masked_where((Sws < 0.08) | (Sws > 0.11),D_65)
    Sws11D_65 = numpy.ma.masked_where((Sws < 0.11) | (Sws > 0.14),D_65)
    Sws14D_65 = numpy.ma.masked_where((Sws < 0.14) | (Sws > 0.17),D_65)
    Sws17D_65 = numpy.ma.masked_where((Sws < 0.17) | (Sws > 0.20),D_65)
    Sws20D_65 = numpy.ma.masked_where((Sws < 0.20) | (Sws > 0.23),D_65)
    Sws23D_65 = numpy.ma.masked_where((Sws < 0.23) | (Sws > 0.26),D_65)
    Sws26D_65 = numpy.ma.masked_where((Sws < 0.26) | (Sws > 0.29),D_65)
    Sws29D_65 = numpy.ma.masked_where(Sws < 0.29,D_65)
    
    metutils.CreateSeries(ds,'Sws2D_65',Sws2D_65,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D_65',Sws5D_65,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws8D_65',Sws8D_65,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws11D_65',Sws11D_65,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws14D_65',Sws14D_65,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws17D_65',Sws17D_65,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws20D_65',Sws20D_65,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws23D_65',Sws23D_65,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws26D_65',Sws26D_65,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws29D_65',Sws29D_65,FList=[InVar],Descr='',Units='umol/(m2 s)')
        
    D1_3 = numpy.ma.masked_where((D < 1.3) | (D > 1.95),Fc1)
    Sws2D1_3 = numpy.ma.masked_where(Sws > 0.05,D1_3)
    Sws5D1_3 = numpy.ma.masked_where((Sws < 0.05) | (Sws > 0.08),D1_3)
    Sws8D1_3 = numpy.ma.masked_where((Sws < 0.08) | (Sws > 0.11),D1_3)
    Sws11D1_3 = numpy.ma.masked_where((Sws < 0.11) | (Sws > 0.14),D1_3)
    Sws14D1_3 = numpy.ma.masked_where((Sws < 0.14) | (Sws > 0.17),D1_3)
    Sws17D1_3 = numpy.ma.masked_where((Sws < 0.17) | (Sws > 0.20),D1_3)
    Sws20D1_3 = numpy.ma.masked_where((Sws < 0.20) | (Sws > 0.23),D1_3)
    Sws23D1_3 = numpy.ma.masked_where((Sws < 0.23) | (Sws > 0.26),D1_3)
    Sws26D1_3 = numpy.ma.masked_where((Sws < 0.26) | (Sws > 0.29),D1_3)
    Sws29D1_3 = numpy.ma.masked_where(Sws < 0.29,D1_3)
    
    metutils.CreateSeries(ds,'Sws2D1_3',Sws2D1_3,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D1_3',Sws5D1_3,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws8D1_3',Sws8D1_3,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws11D1_3',Sws11D1_3,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws14D1_3',Sws14D1_3,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws17D1_3',Sws17D1_3,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws20D1_3',Sws20D1_3,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws23D1_3',Sws23D1_3,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws26D1_3',Sws26D1_3,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws29D1_3',Sws29D1_3,FList=[InVar],Descr='',Units='umol/(m2 s)')
        
    D1_95 = numpy.ma.masked_where((D < 1.95) | (D > 2.6),Fc1)
    Sws2D1_95 = numpy.ma.masked_where(Sws > 0.05,D1_95)
    Sws5D1_95 = numpy.ma.masked_where((Sws < 0.05) | (Sws > 0.08),D1_95)
    Sws8D1_95 = numpy.ma.masked_where((Sws < 0.08) | (Sws > 0.11),D1_95)
    Sws11D1_95 = numpy.ma.masked_where((Sws < 0.11) | (Sws > 0.14),D1_95)
    Sws14D1_95 = numpy.ma.masked_where((Sws < 0.14) | (Sws > 0.17),D1_95)
    Sws17D1_95 = numpy.ma.masked_where((Sws < 0.17) | (Sws > 0.20),D1_95)
    Sws20D1_95 = numpy.ma.masked_where((Sws < 0.20) | (Sws > 0.23),D1_95)
    Sws23D1_95 = numpy.ma.masked_where((Sws < 0.23) | (Sws > 0.26),D1_95)
    Sws26D1_95 = numpy.ma.masked_where((Sws < 0.26) | (Sws > 0.29),D1_95)
    Sws29D1_95 = numpy.ma.masked_where(Sws < 0.29,D1_95)
    
    metutils.CreateSeries(ds,'Sws2D1_95',Sws2D1_95,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D1_95',Sws5D1_95,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws8D1_95',Sws8D1_95,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws11D1_95',Sws11D1_95,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws14D1_95',Sws14D1_95,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws17D1_95',Sws17D1_95,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws20D1_95',Sws20D1_95,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws23D1_95',Sws23D1_95,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws26D1_95',Sws26D1_95,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws29D1_95',Sws29D1_95,FList=[InVar],Descr='',Units='umol/(m2 s)')
        
    D2_6 = numpy.ma.masked_where((D < 2.6) | (D > 3.25),Fc1)
    Sws2D2_6 = numpy.ma.masked_where(Sws > 0.05,D2_6)
    Sws5D2_6 = numpy.ma.masked_where((Sws < 0.05) | (Sws > 0.08),D2_6)
    Sws8D2_6 = numpy.ma.masked_where((Sws < 0.08) | (Sws > 0.11),D2_6)
    Sws11D2_6 = numpy.ma.masked_where((Sws < 0.11) | (Sws > 0.14),D2_6)
    Sws14D2_6 = numpy.ma.masked_where((Sws < 0.14) | (Sws > 0.17),D2_6)
    Sws17D2_6 = numpy.ma.masked_where((Sws < 0.17) | (Sws > 0.20),D2_6)
    Sws20D2_6 = numpy.ma.masked_where((Sws < 0.20) | (Sws > 0.23),D2_6)
    Sws23D2_6 = numpy.ma.masked_where((Sws < 0.23) | (Sws > 0.26),D2_6)
    Sws26D2_6 = numpy.ma.masked_where((Sws < 0.26) | (Sws > 0.29),D2_6)
    Sws29D2_6 = numpy.ma.masked_where(Sws < 0.29,D2_6)
    
    metutils.CreateSeries(ds,'Sws2D2_6',Sws2D2_6,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D2_6',Sws5D2_6,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws8D2_6',Sws8D2_6,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws11D2_6',Sws11D2_6,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws14D2_6',Sws14D2_6,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws17D2_6',Sws17D2_6,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws20D2_6',Sws20D2_6,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws23D2_6',Sws23D2_6,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws26D2_6',Sws26D2_6,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws29D2_6',Sws29D2_6,FList=[InVar],Descr='',Units='umol/(m2 s)')
        
    D3_25 = numpy.ma.masked_where((D < 3.25) | (D > 3.9),Fc1)
    Sws2D3_25 = numpy.ma.masked_where(Sws > 0.05,D3_25)
    Sws5D3_25 = numpy.ma.masked_where((Sws < 0.05) | (Sws > 0.08),D3_25)
    Sws8D3_25 = numpy.ma.masked_where((Sws < 0.08) | (Sws > 0.11),D3_25)
    Sws11D3_25 = numpy.ma.masked_where((Sws < 0.11) | (Sws > 0.14),D3_25)
    Sws14D3_25 = numpy.ma.masked_where((Sws < 0.14) | (Sws > 0.17),D3_25)
    Sws17D3_25 = numpy.ma.masked_where((Sws < 0.17) | (Sws > 0.20),D3_25)
    Sws20D3_25 = numpy.ma.masked_where((Sws < 0.20) | (Sws > 0.23),D3_25)
    Sws23D3_25 = numpy.ma.masked_where((Sws < 0.23) | (Sws > 0.26),D3_25)
    Sws26D3_25 = numpy.ma.masked_where((Sws < 0.26) | (Sws > 0.29),D3_25)
    Sws29D3_25 = numpy.ma.masked_where(Sws < 0.29,D3_25)
    
    metutils.CreateSeries(ds,'Sws2D3_25',Sws2D3_25,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D3_25',Sws5D3_25,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws8D3_25',Sws8D3_25,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws11D3_25',Sws11D3_25,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws14D3_25',Sws14D3_25,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws17D3_25',Sws17D3_25,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws20D3_25',Sws20D3_25,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws23D3_25',Sws23D3_25,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws26D3_25',Sws26D3_25,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws29D3_25',Sws29D3_25,FList=[InVar],Descr='',Units='umol/(m2 s)')
        
    D3_9 = numpy.ma.masked_where((D < 3.9) | (D > 4.55),Fc1)
    Sws2D3_9 = numpy.ma.masked_where(Sws > 0.05,D3_9)
    Sws5D3_9 = numpy.ma.masked_where((Sws < 0.05) | (Sws > 0.08),D3_9)
    Sws8D3_9 = numpy.ma.masked_where((Sws < 0.08) | (Sws > 0.11),D3_9)
    Sws11D3_9 = numpy.ma.masked_where((Sws < 0.11) | (Sws > 0.14),D3_9)
    Sws14D3_9 = numpy.ma.masked_where((Sws < 0.14) | (Sws > 0.17),D3_9)
    Sws17D3_9 = numpy.ma.masked_where((Sws < 0.17) | (Sws > 0.20),D3_9)
    Sws20D3_9 = numpy.ma.masked_where((Sws < 0.20) | (Sws > 0.23),D3_9)
    Sws23D3_9 = numpy.ma.masked_where((Sws < 0.23) | (Sws > 0.26),D3_9)
    Sws26D3_9 = numpy.ma.masked_where((Sws < 0.26) | (Sws > 0.29),D3_9)
    Sws29D3_9 = numpy.ma.masked_where(Sws < 0.29,D3_9)
    
    metutils.CreateSeries(ds,'Sws2D3_9',Sws2D3_9,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D3_9',Sws5D3_9,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws8D3_9',Sws8D3_9,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws11D3_9',Sws11D3_9,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws14D3_9',Sws14D3_9,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws17D3_9',Sws17D3_9,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws20D3_9',Sws20D3_9,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws23D3_9',Sws23D3_9,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws26D3_9',Sws26D3_9,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws29D3_9',Sws29D3_9,FList=[InVar],Descr='',Units='umol/(m2 s)')
        
    D4_55 = numpy.ma.masked_where((D < 4.55) | (D > 5.2),Fc1)
    Sws2D4_55 = numpy.ma.masked_where(Sws > 0.05,D4_55)
    Sws5D4_55 = numpy.ma.masked_where((Sws < 0.05) | (Sws > 0.08),D4_55)
    Sws8D4_55 = numpy.ma.masked_where((Sws < 0.08) | (Sws > 0.11),D4_55)
    Sws11D4_55 = numpy.ma.masked_where((Sws < 0.11) | (Sws > 0.14),D4_55)
    Sws14D4_55 = numpy.ma.masked_where((Sws < 0.14) | (Sws > 0.17),D4_55)
    Sws17D4_55 = numpy.ma.masked_where((Sws < 0.17) | (Sws > 0.20),D4_55)
    Sws20D4_55 = numpy.ma.masked_where((Sws < 0.20) | (Sws > 0.23),D4_55)
    Sws23D4_55 = numpy.ma.masked_where((Sws < 0.23) | (Sws > 0.26),D4_55)
    Sws26D4_55 = numpy.ma.masked_where((Sws < 0.26) | (Sws > 0.29),D4_55)
    Sws29D4_55 = numpy.ma.masked_where(Sws < 0.29,D4_55)
    
    metutils.CreateSeries(ds,'Sws2D4_55',Sws2D4_55,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D4_55',Sws5D4_55,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws8D4_55',Sws8D4_55,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws11D4_55',Sws11D4_55,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws14D4_55',Sws14D4_55,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws17D4_55',Sws17D4_55,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws20D4_55',Sws20D4_55,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws23D4_55',Sws23D4_55,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws26D4_55',Sws26D4_55,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws29D4_55',Sws29D4_55,FList=[InVar],Descr='',Units='umol/(m2 s)')
        
    D5_2 = numpy.ma.masked_where((D < 5.2) | (D > 5.85),Fc1)
    Sws2D5_2 = numpy.ma.masked_where(Sws > 0.05,D5_2)
    Sws5D5_2 = numpy.ma.masked_where((Sws < 0.05) | (Sws > 0.08),D5_2)
    Sws8D5_2 = numpy.ma.masked_where((Sws < 0.08) | (Sws > 0.11),D5_2)
    Sws11D5_2 = numpy.ma.masked_where((Sws < 0.11) | (Sws > 0.14),D5_2)
    Sws14D5_2 = numpy.ma.masked_where((Sws < 0.14) | (Sws > 0.17),D5_2)
    Sws17D5_2 = numpy.ma.masked_where((Sws < 0.17) | (Sws > 0.20),D5_2)
    Sws20D5_2 = numpy.ma.masked_where((Sws < 0.20) | (Sws > 0.23),D5_2)
    Sws23D5_2 = numpy.ma.masked_where((Sws < 0.23) | (Sws > 0.26),D5_2)
    Sws26D5_2 = numpy.ma.masked_where((Sws < 0.26) | (Sws > 0.29),D5_2)
    Sws29D5_2 = numpy.ma.masked_where(Sws < 0.29,D5_2)
    
    metutils.CreateSeries(ds,'Sws2D5_2',Sws2D5_2,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D5_2',Sws5D5_2,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws8D5_2',Sws8D5_2,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws11D5_2',Sws11D5_2,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws14D5_2',Sws14D5_2,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws17D5_2',Sws17D5_2,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws20D5_2',Sws20D5_2,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws23D5_2',Sws23D5_2,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws26D5_2',Sws26D5_2,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws29D5_2',Sws29D5_2,FList=[InVar],Descr='',Units='umol/(m2 s)')
        
    D5_85 = numpy.ma.masked_where(D < 5.85,Fc1)
    Sws2D5_85 = numpy.ma.masked_where(Sws > 0.05,D5_85)
    Sws5D5_85 = numpy.ma.masked_where((Sws < 0.05) | (Sws > 0.08),D5_85)
    Sws8D5_85 = numpy.ma.masked_where((Sws < 0.08) | (Sws > 0.11),D5_85)
    Sws11D5_85 = numpy.ma.masked_where((Sws < 0.11) | (Sws > 0.14),D5_85)
    Sws14D5_85 = numpy.ma.masked_where((Sws < 0.14) | (Sws > 0.17),D5_85)
    Sws17D5_85 = numpy.ma.masked_where((Sws < 0.17) | (Sws > 0.20),D5_85)
    Sws20D5_85 = numpy.ma.masked_where((Sws < 0.20) | (Sws > 0.23),D5_85)
    Sws23D5_85 = numpy.ma.masked_where((Sws < 0.23) | (Sws > 0.26),D5_85)
    Sws26D5_85 = numpy.ma.masked_where((Sws < 0.26) | (Sws > 0.29),D5_85)
    Sws29D5_85 = numpy.ma.masked_where(Sws < 0.29,D5_85)
    
    metutils.CreateSeries(ds,'Sws2D5_85',Sws2D5_85,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D5_85',Sws5D5_85,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws8D5_85',Sws8D5_85,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws11D5_85',Sws11D5_85,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws14D5_85',Sws14D5_85,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws17D5_85',Sws17D5_85,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws20D5_85',Sws20D5_85,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws23D5_85',Sws23D5_85,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws26D5_85',Sws26D5_85,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws29D5_85',Sws29D5_85,FList=[InVar],Descr='',Units='umol/(m2 s)')
    
    log.info('NEP-Met Envelopes: All ready')

def EsBinFilters(cf,ds):
    '''Prepare data into bins for meteorological envelope calculation
    '''
    log.info('Preparing data for meteorological envelope calculation')
    if metutils.cfkeycheck(cf,Base='Params',ThisOne='dependent'):
        InVar = cf['Params']['dependent']
    else:
        InVar = 'Fc'
    
    if (InVar == 'Fc_wpl' or InVar == 'Fc'):
        Fcmg,f = metutils.GetSeriesasMA(ds,InVar)
        Fc1 = ((0-Fcmg) * (10 ** 6)) / (1000 * 44)
    elif (InVar == 'GPP'):
        Fc1,f = metutils.GetSeriesasMA(ds,InVar)
    
    Ts,f = metutils.GetSeriesasMA(ds,'Ts')
    Sws,f = metutils.GetSeriesasMA(ds,'Sws')
    Fsd,f = metutils.GetSeriesasMA(ds,'Fsd')
    D,f = metutils.GetSeriesasMA(ds,'VPD')
    nRecs = len(Fc1)
    
    Es0 = numpy.ma.masked_where(Fsd > 100,Fc1)
    Sws2Es0 = numpy.ma.masked_where(Sws > 0.05,Es0)
    Sws5Es0 = numpy.ma.masked_where((Sws < 0.05) | (Sws > 0.08),Es0)
    Sws8Es0 = numpy.ma.masked_where((Sws < 0.08) | (Sws > 0.11),Es0)
    Sws11Es0 = numpy.ma.masked_where((Sws < 0.11) | (Sws > 0.14),Es0)
    Sws14Es0 = numpy.ma.masked_where((Sws < 0.14) | (Sws > 0.17),Es0)
    Sws17Es0 = numpy.ma.masked_where((Sws < 0.17) | (Sws > 0.20),Es0)
    Sws20Es0 = numpy.ma.masked_where((Sws < 0.20) | (Sws > 0.23),Es0)
    Sws23Es0 = numpy.ma.masked_where((Sws < 0.23) | (Sws > 0.26),Es0)
    Sws26Es0 = numpy.ma.masked_where((Sws < 0.26) | (Sws > 0.29),Es0)
    Sws29Es0 = numpy.ma.masked_where(Sws < 0.29,Es0)
    
    metutils.CreateSeries(ds,'Sws2Es0',Sws2Es0,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5Es0',Sws5Es0,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws8Es0',Sws8Es0,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws11Es0',Sws11Es0,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws14Es0',Sws14Es0,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws17Es0',Sws17Es0,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws20Es0',Sws20Es0,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws23Es0',Sws23Es0,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws26Es0',Sws26Es0,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws29Es0',Sws29Es0,FList=[InVar],Descr='',Units='umol/(m2 s)')
        
    Es100 = numpy.ma.masked_where((Fsd < 100) | (Fsd > 200),Fc1)
    Sws2Es100 = numpy.ma.masked_where(Sws > 0.05,Es100)
    Sws5Es100 = numpy.ma.masked_where((Sws < 0.05) | (Sws > 0.08),Es100)
    Sws8Es100 = numpy.ma.masked_where((Sws < 0.08) | (Sws > 0.11),Es100)
    Sws11Es100 = numpy.ma.masked_where((Sws < 0.11) | (Sws > 0.14),Es100)
    Sws14Es100 = numpy.ma.masked_where((Sws < 0.14) | (Sws > 0.17),Es100)
    Sws17Es100 = numpy.ma.masked_where((Sws < 0.17) | (Sws > 0.20),Es100)
    Sws20Es100 = numpy.ma.masked_where((Sws < 0.20) | (Sws > 0.23),Es100)
    Sws23Es100 = numpy.ma.masked_where((Sws < 0.23) | (Sws > 0.26),Es100)
    Sws26Es100 = numpy.ma.masked_where((Sws < 0.26) | (Sws > 0.29),Es100)
    Sws29Es100 = numpy.ma.masked_where(Sws < 0.29,Es100)
    
    metutils.CreateSeries(ds,'Sws2Es100',Sws2Es100,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5Es100',Sws5Es100,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws8Es100',Sws8Es100,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws11Es100',Sws11Es100,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws14Es100',Sws14Es100,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws17Es100',Sws17Es100,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws20Es100',Sws20Es100,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws23Es100',Sws23Es100,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws26Es100',Sws26Es100,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws29Es100',Sws29Es100,FList=[InVar],Descr='',Units='umol/(m2 s)')
        
    Es200 = numpy.ma.masked_where((Fsd < 200) | (Fsd > 300),Fc1)
    Sws2Es200 = numpy.ma.masked_where(Sws > 0.05,Es200)
    Sws5Es200 = numpy.ma.masked_where((Sws < 0.05) | (Sws > 0.08),Es200)
    Sws8Es200 = numpy.ma.masked_where((Sws < 0.08) | (Sws > 0.11),Es200)
    Sws11Es200 = numpy.ma.masked_where((Sws < 0.11) | (Sws > 0.14),Es200)
    Sws14Es200 = numpy.ma.masked_where((Sws < 0.14) | (Sws > 0.17),Es200)
    Sws17Es200 = numpy.ma.masked_where((Sws < 0.17) | (Sws > 0.20),Es200)
    Sws20Es200 = numpy.ma.masked_where((Sws < 0.20) | (Sws > 0.23),Es200)
    Sws23Es200 = numpy.ma.masked_where((Sws < 0.23) | (Sws > 0.26),Es200)
    Sws26Es200 = numpy.ma.masked_where((Sws < 0.26) | (Sws > 0.29),Es200)
    Sws29Es200 = numpy.ma.masked_where(Sws < 0.29,Es200)
    
    metutils.CreateSeries(ds,'Sws2Es200',Sws2Es200,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5Es200',Sws5Es200,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws8Es200',Sws8Es200,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws11Es200',Sws11Es200,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws14Es200',Sws14Es200,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws17Es200',Sws17Es200,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws20Es200',Sws20Es200,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws23Es200',Sws23Es200,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws26Es200',Sws26Es200,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws29Es200',Sws29Es200,FList=[InVar],Descr='',Units='umol/(m2 s)')
        
    Es300 = numpy.ma.masked_where((Fsd < 300) | (Fsd > 400),Fc1)
    Sws2Es300 = numpy.ma.masked_where(Sws > 0.05,Es300)
    Sws5Es300 = numpy.ma.masked_where((Sws < 0.05) | (Sws > 0.08),Es300)
    Sws8Es300 = numpy.ma.masked_where((Sws < 0.08) | (Sws > 0.11),Es300)
    Sws11Es300 = numpy.ma.masked_where((Sws < 0.11) | (Sws > 0.14),Es300)
    Sws14Es300 = numpy.ma.masked_where((Sws < 0.14) | (Sws > 0.17),Es300)
    Sws17Es300 = numpy.ma.masked_where((Sws < 0.17) | (Sws > 0.20),Es300)
    Sws20Es300 = numpy.ma.masked_where((Sws < 0.20) | (Sws > 0.23),Es300)
    Sws23Es300 = numpy.ma.masked_where((Sws < 0.23) | (Sws > 0.26),Es300)
    Sws26Es300 = numpy.ma.masked_where((Sws < 0.26) | (Sws > 0.29),Es300)
    Sws29Es300 = numpy.ma.masked_where(Sws < 0.29,Es300)
    
    metutils.CreateSeries(ds,'Sws2Es300',Sws2Es300,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5Es300',Sws5Es300,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws8Es300',Sws8Es300,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws11Es300',Sws11Es300,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws14Es300',Sws14Es300,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws17Es300',Sws17Es300,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws20Es300',Sws20Es300,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws23Es300',Sws23Es300,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws26Es300',Sws26Es300,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws29Es300',Sws29Es300,FList=[InVar],Descr='',Units='umol/(m2 s)')
        
    Es400 = numpy.ma.masked_where((Fsd < 400) | (Fsd > 500),Fc1)
    Sws2Es400 = numpy.ma.masked_where(Sws > 0.05,Es400)
    Sws5Es400 = numpy.ma.masked_where((Sws < 0.05) | (Sws > 0.08),Es400)
    Sws8Es400 = numpy.ma.masked_where((Sws < 0.08) | (Sws > 0.11),Es400)
    Sws11Es400 = numpy.ma.masked_where((Sws < 0.11) | (Sws > 0.14),Es400)
    Sws14Es400 = numpy.ma.masked_where((Sws < 0.14) | (Sws > 0.17),Es400)
    Sws17Es400 = numpy.ma.masked_where((Sws < 0.17) | (Sws > 0.20),Es400)
    Sws20Es400 = numpy.ma.masked_where((Sws < 0.20) | (Sws > 0.23),Es400)
    Sws23Es400 = numpy.ma.masked_where((Sws < 0.23) | (Sws > 0.26),Es400)
    Sws26Es400 = numpy.ma.masked_where((Sws < 0.26) | (Sws > 0.29),Es400)
    Sws29Es400 = numpy.ma.masked_where(Sws < 0.29,Es400)
    
    metutils.CreateSeries(ds,'Sws2Es400',Sws2Es400,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5Es400',Sws5Es400,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws8Es400',Sws8Es400,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws11Es400',Sws11Es400,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws14Es400',Sws14Es400,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws17Es400',Sws17Es400,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws20Es400',Sws20Es400,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws23Es400',Sws23Es400,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws26Es400',Sws26Es400,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws29Es400',Sws29Es400,FList=[InVar],Descr='',Units='umol/(m2 s)')
        
    Es500 = numpy.ma.masked_where((Fsd < 500) | (Fsd > 600),Fc1)
    Sws2Es500 = numpy.ma.masked_where(Sws > 0.05,Es500)
    Sws5Es500 = numpy.ma.masked_where((Sws < 0.05) | (Sws > 0.08),Es500)
    Sws8Es500 = numpy.ma.masked_where((Sws < 0.08) | (Sws > 0.11),Es500)
    Sws11Es500 = numpy.ma.masked_where((Sws < 0.11) | (Sws > 0.14),Es500)
    Sws14Es500 = numpy.ma.masked_where((Sws < 0.14) | (Sws > 0.17),Es500)
    Sws17Es500 = numpy.ma.masked_where((Sws < 0.17) | (Sws > 0.20),Es500)
    Sws20Es500 = numpy.ma.masked_where((Sws < 0.20) | (Sws > 0.23),Es500)
    Sws23Es500 = numpy.ma.masked_where((Sws < 0.23) | (Sws > 0.26),Es500)
    Sws26Es500 = numpy.ma.masked_where((Sws < 0.26) | (Sws > 0.29),Es500)
    Sws29Es500 = numpy.ma.masked_where(Sws < 0.29,Es500)
    
    metutils.CreateSeries(ds,'Sws2Es500',Sws2Es500,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5Es500',Sws5Es500,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws8Es500',Sws8Es500,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws11Es500',Sws11Es500,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws14Es500',Sws14Es500,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws17Es500',Sws17Es500,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws20Es500',Sws20Es500,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws23Es500',Sws23Es500,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws26Es500',Sws26Es500,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws29Es500',Sws29Es500,FList=[InVar],Descr='',Units='umol/(m2 s)')
        
    Es600 = numpy.ma.masked_where((Fsd < 600) | (Fsd > 700),Fc1)
    Sws2Es600 = numpy.ma.masked_where(Sws > 0.05,Es600)
    Sws5Es600 = numpy.ma.masked_where((Sws < 0.05) | (Sws > 0.08),Es600)
    Sws8Es600 = numpy.ma.masked_where((Sws < 0.08) | (Sws > 0.11),Es600)
    Sws11Es600 = numpy.ma.masked_where((Sws < 0.11) | (Sws > 0.14),Es600)
    Sws14Es600 = numpy.ma.masked_where((Sws < 0.14) | (Sws > 0.17),Es600)
    Sws17Es600 = numpy.ma.masked_where((Sws < 0.17) | (Sws > 0.20),Es600)
    Sws20Es600 = numpy.ma.masked_where((Sws < 0.20) | (Sws > 0.23),Es600)
    Sws23Es600 = numpy.ma.masked_where((Sws < 0.23) | (Sws > 0.26),Es600)
    Sws26Es600 = numpy.ma.masked_where((Sws < 0.26) | (Sws > 0.29),Es600)
    Sws29Es600 = numpy.ma.masked_where(Sws < 0.29,Es600)
    
    metutils.CreateSeries(ds,'Sws2Es600',Sws2Es600,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5Es600',Sws5Es600,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws8Es600',Sws8Es600,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws11Es600',Sws11Es600,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws14Es600',Sws14Es600,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws17Es600',Sws17Es600,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws20Es600',Sws20Es600,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws23Es600',Sws23Es600,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws26Es600',Sws26Es600,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws29Es600',Sws29Es600,FList=[InVar],Descr='',Units='umol/(m2 s)')
        
    Es700 = numpy.ma.masked_where((Fsd < 700) | (Fsd > 800),Fc1)
    Sws2Es700 = numpy.ma.masked_where(Sws > 0.05,Es700)
    Sws5Es700 = numpy.ma.masked_where((Sws < 0.05) | (Sws > 0.08),Es700)
    Sws8Es700 = numpy.ma.masked_where((Sws < 0.08) | (Sws > 0.11),Es700)
    Sws11Es700 = numpy.ma.masked_where((Sws < 0.11) | (Sws > 0.14),Es700)
    Sws14Es700 = numpy.ma.masked_where((Sws < 0.14) | (Sws > 0.17),Es700)
    Sws17Es700 = numpy.ma.masked_where((Sws < 0.17) | (Sws > 0.20),Es700)
    Sws20Es700 = numpy.ma.masked_where((Sws < 0.20) | (Sws > 0.23),Es700)
    Sws23Es700 = numpy.ma.masked_where((Sws < 0.23) | (Sws > 0.26),Es700)
    Sws26Es700 = numpy.ma.masked_where((Sws < 0.26) | (Sws > 0.29),Es700)
    Sws29Es700 = numpy.ma.masked_where(Sws < 0.29,Es700)
    
    metutils.CreateSeries(ds,'Sws2Es700',Sws2Es700,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5Es700',Sws5Es700,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws8Es700',Sws8Es700,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws11Es700',Sws11Es700,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws14Es700',Sws14Es700,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws17Es700',Sws17Es700,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws20Es700',Sws20Es700,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws23Es700',Sws23Es700,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws26Es700',Sws26Es700,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws29Es700',Sws29Es700,FList=[InVar],Descr='',Units='umol/(m2 s)')
        
    Es800 = numpy.ma.masked_where((Fsd < 800) | (Fsd > 900),Fc1)
    Sws2Es800 = numpy.ma.masked_where(Sws > 0.05,Es800)
    Sws5Es800 = numpy.ma.masked_where((Sws < 0.05) | (Sws > 0.08),Es800)
    Sws8Es800 = numpy.ma.masked_where((Sws < 0.08) | (Sws > 0.11),Es800)
    Sws11Es800 = numpy.ma.masked_where((Sws < 0.11) | (Sws > 0.14),Es800)
    Sws14Es800 = numpy.ma.masked_where((Sws < 0.14) | (Sws > 0.17),Es800)
    Sws17Es800 = numpy.ma.masked_where((Sws < 0.17) | (Sws > 0.20),Es800)
    Sws20Es800 = numpy.ma.masked_where((Sws < 0.20) | (Sws > 0.23),Es800)
    Sws23Es800 = numpy.ma.masked_where((Sws < 0.23) | (Sws > 0.26),Es800)
    Sws26Es800 = numpy.ma.masked_where((Sws < 0.26) | (Sws > 0.29),Es800)
    Sws29Es800 = numpy.ma.masked_where(Sws < 0.29,Es800)
    
    metutils.CreateSeries(ds,'Sws2Es800',Sws2Es800,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5Es800',Sws5Es800,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws8Es800',Sws8Es800,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws11Es800',Sws11Es800,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws14Es800',Sws14Es800,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws17Es800',Sws17Es800,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws20Es800',Sws20Es800,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws23Es800',Sws23Es800,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws26Es800',Sws26Es800,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws29Es800',Sws29Es800,FList=[InVar],Descr='',Units='umol/(m2 s)')
        
    Es900 = numpy.ma.masked_where(Fsd < 900,Fc1)
    Sws2Es900 = numpy.ma.masked_where(Sws > 0.05,Es900)
    Sws5Es900 = numpy.ma.masked_where((Sws < 0.05) | (Sws > 0.08),Es900)
    Sws8Es900 = numpy.ma.masked_where((Sws < 0.08) | (Sws > 0.11),Es900)
    Sws11Es900 = numpy.ma.masked_where((Sws < 0.11) | (Sws > 0.14),Es900)
    Sws14Es900 = numpy.ma.masked_where((Sws < 0.14) | (Sws > 0.17),Es900)
    Sws17Es900 = numpy.ma.masked_where((Sws < 0.17) | (Sws > 0.20),Es900)
    Sws20Es900 = numpy.ma.masked_where((Sws < 0.20) | (Sws > 0.23),Es900)
    Sws23Es900 = numpy.ma.masked_where((Sws < 0.23) | (Sws > 0.26),Es900)
    Sws26Es900 = numpy.ma.masked_where((Sws < 0.26) | (Sws > 0.29),Es900)
    Sws29Es900 = numpy.ma.masked_where(Sws < 0.29,Es900)
    
    metutils.CreateSeries(ds,'Sws2Es900',Sws2Es900,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5Es900',Sws5Es900,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws8Es900',Sws8Es900,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws11Es900',Sws11Es900,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws14Es900',Sws14Es900,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws17Es900',Sws17Es900,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws20Es900',Sws20Es900,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws23Es900',Sws23Es900,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws26Es900',Sws26Es900,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws29Es900',Sws29Es900,FList=[InVar],Descr='',Units='umol/(m2 s)')
        
    log.info('NEP-Met Envelopes: All ready')

def LRFBinFilters(cf,ds):
    '''Prepare data into bins for meteorological envelope calculation
    '''
    log.info('Preparing data for meteorological envelope calculation')
    if metutils.cfkeycheck(cf,Base='Params',ThisOne='dependent'):
        InVar = cf['Params']['dependent']
    else:
        InVar = 'Fc'
    
    if (InVar == 'Fc_wpl' or InVar == 'Fc'):
        Fcmg,f = metutils.GetSeriesasMA(ds,InVar)
        Fc1 = ((Fcmg) * (10 ** 6)) / (1000 * 44)
    else:
        Fc1,f = metutils.GetSeriesasMA(ds,InVar)
    
    Ts,f = metutils.GetSeriesasMA(ds,'Ts')
    Sws,f = metutils.GetSeriesasMA(ds,'Sws')
    Fsd,f = metutils.GetSeriesasMA(ds,'Fsd')
    D,f = metutils.GetSeriesasMA(ds,'VPD')
    nRecs = len(Fc1)
    
    #Sws 0-5%
    Sws0 = numpy.ma.masked_where(Sws > 0.05,Fc1)
    Es0 = numpy.ma.masked_where(Fsd > 100,Sws0)
    Sws0Ts6Es0 = numpy.ma.masked_where(Ts > 16,Es0)
    Sws0Ts16Es0 = numpy.ma.masked_where((Ts < 16) | (Ts > 26),Es0)
    Sws0Ts26Es0 = numpy.ma.masked_where((Ts < 26) | (Ts > 36),Es0)
    Sws0Ts36Es0 = numpy.ma.masked_where((Ts < 36) | (Ts > 46),Es0)
    Sws0Ts46Es0 = numpy.ma.masked_where(Ts < 46,Es0)
    metutils.CreateSeries(ds,'Sws0Ts6Es0',Sws0Ts6Es0,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0Ts16Es0',Sws0Ts16Es0,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0Ts26Es0',Sws0Ts26Es0,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0Ts36Es0',Sws0Ts36Es0,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0Ts46Es0',Sws0Ts46Es0,FList=[InVar],Descr='',Units='umol/(m2 s)')
    
    Es100 = numpy.ma.masked_where((Fsd < 100) | (Fsd > 200),Sws0)
    Sws0Ts6Es100 = numpy.ma.masked_where(Ts > 16,Es100)
    Sws0Ts16Es100 = numpy.ma.masked_where((Ts < 16) | (Ts > 26),Es100)
    Sws0Ts26Es100 = numpy.ma.masked_where((Ts < 26) | (Ts > 36),Es100)
    Sws0Ts36Es100 = numpy.ma.masked_where((Ts < 36) | (Ts > 46),Es100)
    Sws0Ts46Es100 = numpy.ma.masked_where(Ts < 46,Es100)
    metutils.CreateSeries(ds,'Sws0Ts6Es100',Sws0Ts6Es100,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0Ts16Es100',Sws0Ts16Es100,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0Ts26Es100',Sws0Ts26Es100,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0Ts36Es100',Sws0Ts36Es100,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0Ts46Es100',Sws0Ts46Es100,FList=[InVar],Descr='',Units='umol/(m2 s)')
    
    Es200 = numpy.ma.masked_where((Fsd < 200) | (Fsd > 300),Sws0)
    Sws0Ts6Es200 = numpy.ma.masked_where(Ts > 16,Es200)
    Sws0Ts16Es200 = numpy.ma.masked_where((Ts < 16) | (Ts > 26),Es200)
    Sws0Ts26Es200 = numpy.ma.masked_where((Ts < 26) | (Ts > 36),Es200)
    Sws0Ts36Es200 = numpy.ma.masked_where((Ts < 36) | (Ts > 46),Es200)
    Sws0Ts46Es200 = numpy.ma.masked_where(Ts < 46,Es200)
    metutils.CreateSeries(ds,'Sws0Ts6Es200',Sws0Ts6Es200,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0Ts16Es200',Sws0Ts16Es200,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0Ts26Es200',Sws0Ts26Es200,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0Ts36Es200',Sws0Ts36Es200,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0Ts46Es200',Sws0Ts46Es200,FList=[InVar],Descr='',Units='umol/(m2 s)')
    
    Es300 = numpy.ma.masked_where((Fsd < 300) | (Fsd > 400),Sws0)
    Sws0Ts6Es300 = numpy.ma.masked_where(Ts > 16,Es300)
    Sws0Ts16Es300 = numpy.ma.masked_where((Ts < 16) | (Ts > 26),Es300)
    Sws0Ts26Es300 = numpy.ma.masked_where((Ts < 26) | (Ts > 36),Es300)
    Sws0Ts36Es300 = numpy.ma.masked_where((Ts < 36) | (Ts > 46),Es300)
    Sws0Ts46Es300 = numpy.ma.masked_where(Ts < 46,Es300)
    metutils.CreateSeries(ds,'Sws0Ts6Es300',Sws0Ts6Es300,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0Ts16Es300',Sws0Ts16Es300,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0Ts26Es300',Sws0Ts26Es300,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0Ts36Es300',Sws0Ts36Es300,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0Ts46Es300',Sws0Ts46Es300,FList=[InVar],Descr='',Units='umol/(m2 s)')
    
    Es400 = numpy.ma.masked_where((Fsd < 400) | (Fsd > 500),Sws0)
    Sws0Ts6Es400 = numpy.ma.masked_where(Ts > 16,Es400)
    Sws0Ts16Es400 = numpy.ma.masked_where((Ts < 16) | (Ts > 26),Es400)
    Sws0Ts26Es400 = numpy.ma.masked_where((Ts < 26) | (Ts > 36),Es400)
    Sws0Ts36Es400 = numpy.ma.masked_where((Ts < 36) | (Ts > 46),Es400)
    Sws0Ts46Es400 = numpy.ma.masked_where(Ts < 46,Es400)
    metutils.CreateSeries(ds,'Sws0Ts6Es400',Sws0Ts6Es400,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0Ts16Es400',Sws0Ts16Es400,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0Ts26Es400',Sws0Ts26Es400,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0Ts36Es400',Sws0Ts36Es400,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0Ts46Es400',Sws0Ts46Es400,FList=[InVar],Descr='',Units='umol/(m2 s)')
    
    Es500 = numpy.ma.masked_where((Fsd < 500) | (Fsd > 600),Sws0)
    Sws0Ts6Es500 = numpy.ma.masked_where(Ts > 16,Es500)
    Sws0Ts16Es500 = numpy.ma.masked_where((Ts < 16) | (Ts > 26),Es500)
    Sws0Ts26Es500 = numpy.ma.masked_where((Ts < 26) | (Ts > 36),Es500)
    Sws0Ts36Es500 = numpy.ma.masked_where((Ts < 36) | (Ts > 46),Es500)
    Sws0Ts46Es500 = numpy.ma.masked_where(Ts < 46,Es500)
    metutils.CreateSeries(ds,'Sws0Ts6Es500',Sws0Ts6Es500,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0Ts16Es500',Sws0Ts16Es500,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0Ts26Es500',Sws0Ts26Es500,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0Ts36Es500',Sws0Ts36Es500,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0Ts46Es500',Sws0Ts46Es500,FList=[InVar],Descr='',Units='umol/(m2 s)')
    
    Es600 = numpy.ma.masked_where((Fsd < 600) | (Fsd > 700),Sws0)
    Sws0Ts6Es600 = numpy.ma.masked_where(Ts > 16,Es600)
    Sws0Ts16Es600 = numpy.ma.masked_where((Ts < 16) | (Ts > 26),Es600)
    Sws0Ts26Es600 = numpy.ma.masked_where((Ts < 26) | (Ts > 36),Es600)
    Sws0Ts36Es600 = numpy.ma.masked_where((Ts < 36) | (Ts > 46),Es600)
    Sws0Ts46Es600 = numpy.ma.masked_where(Ts < 46,Es600)
    metutils.CreateSeries(ds,'Sws0Ts6Es600',Sws0Ts6Es600,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0Ts16Es600',Sws0Ts16Es600,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0Ts26Es600',Sws0Ts26Es600,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0Ts36Es600',Sws0Ts36Es600,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0Ts46Es600',Sws0Ts46Es600,FList=[InVar],Descr='',Units='umol/(m2 s)')
    
    Es700 = numpy.ma.masked_where((Fsd < 700) | (Fsd > 800),Sws0)
    Sws0Ts6Es700 = numpy.ma.masked_where(Ts > 16,Es700)
    Sws0Ts16Es700 = numpy.ma.masked_where((Ts < 16) | (Ts > 26),Es700)
    Sws0Ts26Es700 = numpy.ma.masked_where((Ts < 26) | (Ts > 36),Es700)
    Sws0Ts36Es700 = numpy.ma.masked_where((Ts < 36) | (Ts > 46),Es700)
    Sws0Ts46Es700 = numpy.ma.masked_where(Ts < 46,Es700)
    metutils.CreateSeries(ds,'Sws0Ts6Es700',Sws0Ts6Es700,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0Ts16Es700',Sws0Ts16Es700,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0Ts26Es700',Sws0Ts26Es700,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0Ts36Es700',Sws0Ts36Es700,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0Ts46Es700',Sws0Ts46Es700,FList=[InVar],Descr='',Units='umol/(m2 s)')
    
    Es800 = numpy.ma.masked_where((Fsd < 800) | (Fsd > 900),Sws0)
    Sws0Ts6Es800 = numpy.ma.masked_where(Ts > 16,Es800)
    Sws0Ts16Es800 = numpy.ma.masked_where((Ts < 16) | (Ts > 26),Es800)
    Sws0Ts26Es800 = numpy.ma.masked_where((Ts < 26) | (Ts > 36),Es800)
    Sws0Ts36Es800 = numpy.ma.masked_where((Ts < 36) | (Ts > 46),Es800)
    Sws0Ts46Es800 = numpy.ma.masked_where(Ts < 46,Es800)
    metutils.CreateSeries(ds,'Sws0Ts6Es800',Sws0Ts6Es800,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0Ts16Es800',Sws0Ts16Es800,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0Ts26Es800',Sws0Ts26Es800,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0Ts36Es800',Sws0Ts36Es800,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0Ts46Es800',Sws0Ts46Es800,FList=[InVar],Descr='',Units='umol/(m2 s)')
    
    Es900 = numpy.ma.masked_where(Fsd < 900,Sws0)
    Sws0Ts6Es900 = numpy.ma.masked_where(Ts > 16,Es900)
    Sws0Ts16Es900 = numpy.ma.masked_where((Ts < 16) | (Ts > 26),Es900)
    Sws0Ts26Es900 = numpy.ma.masked_where((Ts < 26) | (Ts > 36),Es900)
    Sws0Ts36Es900 = numpy.ma.masked_where((Ts < 36) | (Ts > 46),Es900)
    Sws0Ts46Es900 = numpy.ma.masked_where(Ts < 46,Es900)
    metutils.CreateSeries(ds,'Sws0Ts6Es900',Sws0Ts6Es900,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0Ts16Es900',Sws0Ts16Es900,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0Ts26Es900',Sws0Ts26Es900,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0Ts36Es900',Sws0Ts36Es900,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0Ts46Es900',Sws0Ts46Es900,FList=[InVar],Descr='',Units='umol/(m2 s)')
        
    #Sws 5-20%
    Sws5 = numpy.ma.masked_where(Sws < 0.05,Fc1)
    Es0 = numpy.ma.masked_where(Fsd > 100,Sws5)
    Sws5Ts6Es0 = numpy.ma.masked_where(Ts > 16,Es0)
    Sws5Ts16Es0 = numpy.ma.masked_where((Ts < 16) | (Ts > 26),Es0)
    Sws5Ts26Es0 = numpy.ma.masked_where((Ts < 26) | (Ts > 36),Es0)
    Sws5Ts36Es0 = numpy.ma.masked_where((Ts < 36) | (Ts > 46),Es0)
    Sws5Ts46Es0 = numpy.ma.masked_where(Ts < 46,Es0)
    metutils.CreateSeries(ds,'Sws5Ts6Es0',Sws5Ts6Es0,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5Ts16Es0',Sws5Ts16Es0,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5Ts26Es0',Sws5Ts26Es0,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5Ts36Es0',Sws5Ts36Es0,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5Ts46Es0',Sws5Ts46Es0,FList=[InVar],Descr='',Units='umol/(m2 s)')
    
    Es100 = numpy.ma.masked_where((Fsd < 100) | (Fsd > 200),Sws5)
    Sws5Ts6Es100 = numpy.ma.masked_where(Ts > 16,Es100)
    Sws5Ts16Es100 = numpy.ma.masked_where((Ts < 16) | (Ts > 26),Es100)
    Sws5Ts26Es100 = numpy.ma.masked_where((Ts < 26) | (Ts > 36),Es100)
    Sws5Ts36Es100 = numpy.ma.masked_where((Ts < 36) | (Ts > 46),Es100)
    Sws5Ts46Es100 = numpy.ma.masked_where(Ts < 46,Es100)
    metutils.CreateSeries(ds,'Sws5Ts6Es100',Sws5Ts6Es100,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5Ts16Es100',Sws5Ts16Es100,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5Ts26Es100',Sws5Ts26Es100,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5Ts36Es100',Sws5Ts36Es100,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5Ts46Es100',Sws5Ts46Es100,FList=[InVar],Descr='',Units='umol/(m2 s)')
    
    Es200 = numpy.ma.masked_where((Fsd < 200) | (Fsd > 300),Sws5)
    Sws5Ts6Es200 = numpy.ma.masked_where(Ts > 16,Es200)
    Sws5Ts16Es200 = numpy.ma.masked_where((Ts < 16) | (Ts > 26),Es200)
    Sws5Ts26Es200 = numpy.ma.masked_where((Ts < 26) | (Ts > 36),Es200)
    Sws5Ts36Es200 = numpy.ma.masked_where((Ts < 36) | (Ts > 46),Es200)
    Sws5Ts46Es200 = numpy.ma.masked_where(Ts < 46,Es200)
    metutils.CreateSeries(ds,'Sws5Ts6Es200',Sws5Ts6Es200,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5Ts16Es200',Sws5Ts16Es200,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5Ts26Es200',Sws5Ts26Es200,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5Ts36Es200',Sws5Ts36Es200,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5Ts46Es200',Sws5Ts46Es200,FList=[InVar],Descr='',Units='umol/(m2 s)')
    
    Es300 = numpy.ma.masked_where((Fsd < 300) | (Fsd > 400),Sws5)
    Sws5Ts6Es300 = numpy.ma.masked_where(Ts > 16,Es300)
    Sws5Ts16Es300 = numpy.ma.masked_where((Ts < 16) | (Ts > 26),Es300)
    Sws5Ts26Es300 = numpy.ma.masked_where((Ts < 26) | (Ts > 36),Es300)
    Sws5Ts36Es300 = numpy.ma.masked_where((Ts < 36) | (Ts > 46),Es300)
    Sws5Ts46Es300 = numpy.ma.masked_where(Ts < 46,Es300)
    metutils.CreateSeries(ds,'Sws5Ts6Es300',Sws5Ts6Es300,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5Ts16Es300',Sws5Ts16Es300,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5Ts26Es300',Sws5Ts26Es300,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5Ts36Es300',Sws5Ts36Es300,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5Ts46Es300',Sws5Ts46Es300,FList=[InVar],Descr='',Units='umol/(m2 s)')
    
    Es400 = numpy.ma.masked_where((Fsd < 400) | (Fsd > 500),Sws5)
    Sws5Ts6Es400 = numpy.ma.masked_where(Ts > 16,Es400)
    Sws5Ts16Es400 = numpy.ma.masked_where((Ts < 16) | (Ts > 26),Es400)
    Sws5Ts26Es400 = numpy.ma.masked_where((Ts < 26) | (Ts > 36),Es400)
    Sws5Ts36Es400 = numpy.ma.masked_where((Ts < 36) | (Ts > 46),Es400)
    Sws5Ts46Es400 = numpy.ma.masked_where(Ts < 46,Es400)
    metutils.CreateSeries(ds,'Sws5Ts6Es400',Sws5Ts6Es400,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5Ts16Es400',Sws5Ts16Es400,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5Ts26Es400',Sws5Ts26Es400,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5Ts36Es400',Sws5Ts36Es400,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5Ts46Es400',Sws5Ts46Es400,FList=[InVar],Descr='',Units='umol/(m2 s)')
    
    Es500 = numpy.ma.masked_where((Fsd < 500) | (Fsd > 600),Sws5)
    Sws5Ts6Es500 = numpy.ma.masked_where(Ts > 16,Es500)
    Sws5Ts16Es500 = numpy.ma.masked_where((Ts < 16) | (Ts > 26),Es500)
    Sws5Ts26Es500 = numpy.ma.masked_where((Ts < 26) | (Ts > 36),Es500)
    Sws5Ts36Es500 = numpy.ma.masked_where((Ts < 36) | (Ts > 46),Es500)
    Sws5Ts46Es500 = numpy.ma.masked_where(Ts < 46,Es500)
    metutils.CreateSeries(ds,'Sws5Ts6Es500',Sws5Ts6Es500,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5Ts16Es500',Sws5Ts16Es500,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5Ts26Es500',Sws5Ts26Es500,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5Ts36Es500',Sws5Ts36Es500,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5Ts46Es500',Sws5Ts46Es500,FList=[InVar],Descr='',Units='umol/(m2 s)')
    
    Es600 = numpy.ma.masked_where((Fsd < 600) | (Fsd > 700),Sws5)
    Sws5Ts6Es600 = numpy.ma.masked_where(Ts > 16,Es600)
    Sws5Ts16Es600 = numpy.ma.masked_where((Ts < 16) | (Ts > 26),Es600)
    Sws5Ts26Es600 = numpy.ma.masked_where((Ts < 26) | (Ts > 36),Es600)
    Sws5Ts36Es600 = numpy.ma.masked_where((Ts < 36) | (Ts > 46),Es600)
    Sws5Ts46Es600 = numpy.ma.masked_where(Ts < 46,Es600)
    metutils.CreateSeries(ds,'Sws5Ts6Es600',Sws5Ts6Es600,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5Ts16Es600',Sws5Ts16Es600,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5Ts26Es600',Sws5Ts26Es600,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5Ts36Es600',Sws5Ts36Es600,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5Ts46Es600',Sws5Ts46Es600,FList=[InVar],Descr='',Units='umol/(m2 s)')
    
    Es700 = numpy.ma.masked_where((Fsd < 700) | (Fsd > 800),Sws5)
    Sws5Ts6Es700 = numpy.ma.masked_where(Ts > 16,Es700)
    Sws5Ts16Es700 = numpy.ma.masked_where((Ts < 16) | (Ts > 26),Es700)
    Sws5Ts26Es700 = numpy.ma.masked_where((Ts < 26) | (Ts > 36),Es700)
    Sws5Ts36Es700 = numpy.ma.masked_where((Ts < 36) | (Ts > 46),Es700)
    Sws5Ts46Es700 = numpy.ma.masked_where(Ts < 46,Es700)
    metutils.CreateSeries(ds,'Sws5Ts6Es700',Sws5Ts6Es700,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5Ts16Es700',Sws5Ts16Es700,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5Ts26Es700',Sws5Ts26Es700,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5Ts36Es700',Sws5Ts36Es700,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5Ts46Es700',Sws5Ts46Es700,FList=[InVar],Descr='',Units='umol/(m2 s)')
    
    Es800 = numpy.ma.masked_where((Fsd < 800) | (Fsd > 900),Sws5)
    Sws5Ts6Es800 = numpy.ma.masked_where(Ts > 16,Es800)
    Sws5Ts16Es800 = numpy.ma.masked_where((Ts < 16) | (Ts > 26),Es800)
    Sws5Ts26Es800 = numpy.ma.masked_where((Ts < 26) | (Ts > 36),Es800)
    Sws5Ts36Es800 = numpy.ma.masked_where((Ts < 36) | (Ts > 46),Es800)
    Sws5Ts46Es800 = numpy.ma.masked_where(Ts < 46,Es800)
    metutils.CreateSeries(ds,'Sws5Ts6Es800',Sws5Ts6Es800,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5Ts16Es800',Sws5Ts16Es800,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5Ts26Es800',Sws5Ts26Es800,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5Ts36Es800',Sws5Ts36Es800,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5Ts46Es800',Sws5Ts46Es800,FList=[InVar],Descr='',Units='umol/(m2 s)')
    
    Es900 = numpy.ma.masked_where(Fsd < 900,Sws5)
    Sws5Ts6Es900 = numpy.ma.masked_where(Ts > 16,Es900)
    Sws5Ts16Es900 = numpy.ma.masked_where((Ts < 16) | (Ts > 26),Es900)
    Sws5Ts26Es900 = numpy.ma.masked_where((Ts < 26) | (Ts > 36),Es900)
    Sws5Ts36Es900 = numpy.ma.masked_where((Ts < 36) | (Ts > 46),Es900)
    Sws5Ts46Es900 = numpy.ma.masked_where(Ts < 46,Es900)
    metutils.CreateSeries(ds,'Sws5Ts6Es900',Sws5Ts6Es900,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5Ts16Es900',Sws5Ts16Es900,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5Ts26Es900',Sws5Ts26Es900,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5Ts36Es900',Sws5Ts36Es900,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5Ts46Es900',Sws5Ts46Es900,FList=[InVar],Descr='',Units='umol/(m2 s)')
        
        
    log.info('NEP-Met Envelopes: All ready')

def LRFDBinFilters(cf,ds):
    '''Prepare data into bins for meteorological envelope calculation
    '''
    log.info('Preparing data for correct meteorological envelope calculation')
    if metutils.cfkeycheck(cf,Base='Params',ThisOne='dependent'):
        InVar = cf['Params']['dependent']
    else:
        InVar = 'Fc'
    
    if (InVar == 'Fc_wpl' or InVar == 'Fc'):
        Fcmg,f = metutils.GetSeriesasMA(ds,InVar)
        Fc1 = ((Fcmg) * (10 ** 6)) / (1000 * 44)
    else:
        Fc1,f = metutils.GetSeriesasMA(ds,InVar)
    
    Ts,f = metutils.GetSeriesasMA(ds,'Ts')
    Sws,f = metutils.GetSeriesasMA(ds,'Sws')
    Fsd,f = metutils.GetSeriesasMA(ds,'Fsd')
    D,f = metutils.GetSeriesasMA(ds,'VPD')
    nRecs = len(Fc1)
    
    #Sws 0-5%
    Sws0 = numpy.ma.masked_where(Sws > 0.05,Fc1)
    Es0 = numpy.ma.masked_where(Fsd > 100,Sws0)
    Sws0D0Es0 = numpy.ma.masked_where(D > 0.65,Es0)
    Sws0D65Es0 = numpy.ma.masked_where((D < 0.65) | (D > 1.3),Es0)
    Sws0D130Es0 = numpy.ma.masked_where((D < 1.3) | (D > 1.95),Es0)
    Sws0D195Es0 = numpy.ma.masked_where((D < 1.95) | (D > 2.60),Es0)
    Sws0D260Es0 = numpy.ma.masked_where((D < 2.60) | (D > 3.25),Es0)
    Sws0D325Es0 = numpy.ma.masked_where((D < 3.25) | (D > 3.90),Es0)
    Sws0D390Es0 = numpy.ma.masked_where((D < 3.90) | (D > 4.55),Es0)
    Sws0D455Es0 = numpy.ma.masked_where((D < 4.55) | (D > 5.20),Es0)
    Sws0D520Es0 = numpy.ma.masked_where((D < 5.20) | (D > 5.85),Es0)
    Sws0D585Es0 = numpy.ma.masked_where(D < 5.85,Es0)
    
    metutils.CreateSeries(ds,'Sws0D0Es0',Sws0D0Es0,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0D65Es0',Sws0D65Es0,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0D130Es0',Sws0D130Es0,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0D195Es0',Sws0D195Es0,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0D260Es0',Sws0D260Es0,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0D325Es0',Sws0D325Es0,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0D390Es0',Sws0D390Es0,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0D455Es0',Sws0D455Es0,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0D520Es0',Sws0D520Es0,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0D585Es0',Sws0D585Es0,FList=[InVar],Descr='',Units='umol/(m2 s)')
        
    Es100 = numpy.ma.masked_where((Fsd < 100) | (Fsd > 200),Sws0)
    Sws0D0Es100 = numpy.ma.masked_where(D > 0.65,Es100)
    Sws0D65Es100 = numpy.ma.masked_where((D < 0.65) | (D > 1.3),Es100)
    Sws0D130Es100 = numpy.ma.masked_where((D < 1.3) | (D > 1.95),Es100)
    Sws0D195Es100 = numpy.ma.masked_where((D < 1.95) | (D > 2.60),Es100)
    Sws0D260Es100 = numpy.ma.masked_where((D < 2.60) | (D > 3.25),Es100)
    Sws0D325Es100 = numpy.ma.masked_where((D < 3.25) | (D > 3.90),Es100)
    Sws0D390Es100 = numpy.ma.masked_where((D < 3.90) | (D > 4.55),Es100)
    Sws0D455Es100 = numpy.ma.masked_where((D < 4.55) | (D > 5.20),Es100)
    Sws0D520Es100 = numpy.ma.masked_where((D < 5.20) | (D > 5.85),Es100)
    Sws0D585Es100 = numpy.ma.masked_where(D < 5.85,Es100)
    
    metutils.CreateSeries(ds,'Sws0D0Es100',Sws0D0Es100,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0D65Es100',Sws0D65Es100,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0D130Es100',Sws0D130Es100,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0D195Es100',Sws0D195Es100,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0D260Es100',Sws0D260Es100,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0D325Es100',Sws0D325Es100,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0D390Es100',Sws0D390Es100,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0D455Es100',Sws0D455Es100,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0D520Es100',Sws0D520Es100,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0D585Es100',Sws0D585Es100,FList=[InVar],Descr='',Units='umol/(m2 s)')
        
    Es200 = numpy.ma.masked_where((Fsd < 200) | (Fsd > 300),Sws0)
    Sws0D0Es200 = numpy.ma.masked_where(D > 0.65,Es200)
    Sws0D65Es200 = numpy.ma.masked_where((D < 0.65) | (D > 1.3),Es200)
    Sws0D130Es200 = numpy.ma.masked_where((D < 1.3) | (D > 1.95),Es200)
    Sws0D195Es200 = numpy.ma.masked_where((D < 1.95) | (D > 2.60),Es200)
    Sws0D260Es200 = numpy.ma.masked_where((D < 2.60) | (D > 3.25),Es200)
    Sws0D325Es200 = numpy.ma.masked_where((D < 3.25) | (D > 3.90),Es200)
    Sws0D390Es200 = numpy.ma.masked_where((D < 3.90) | (D > 4.55),Es200)
    Sws0D455Es200 = numpy.ma.masked_where((D < 4.55) | (D > 5.20),Es200)
    Sws0D520Es200 = numpy.ma.masked_where((D < 5.20) | (D > 5.85),Es200)
    Sws0D585Es200 = numpy.ma.masked_where(D < 5.85,Es200)
    
    metutils.CreateSeries(ds,'Sws0D0Es200',Sws0D0Es200,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0D65Es200',Sws0D65Es200,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0D130Es200',Sws0D130Es200,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0D195Es200',Sws0D195Es200,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0D260Es200',Sws0D260Es200,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0D325Es200',Sws0D325Es200,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0D390Es200',Sws0D390Es200,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0D455Es200',Sws0D455Es200,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0D520Es200',Sws0D520Es200,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0D585Es200',Sws0D585Es200,FList=[InVar],Descr='',Units='umol/(m2 s)')
        
    Es300 = numpy.ma.masked_where((Fsd < 300) | (Fsd > 400),Sws0)
    Sws0D0Es300 = numpy.ma.masked_where(D > 0.65,Es300)
    Sws0D65Es300 = numpy.ma.masked_where((D < 0.65) | (D > 1.3),Es300)
    Sws0D130Es300 = numpy.ma.masked_where((D < 1.3) | (D > 1.95),Es300)
    Sws0D195Es300 = numpy.ma.masked_where((D < 1.95) | (D > 2.60),Es300)
    Sws0D260Es300 = numpy.ma.masked_where((D < 2.60) | (D > 3.25),Es300)
    Sws0D325Es300 = numpy.ma.masked_where((D < 3.25) | (D > 3.90),Es300)
    Sws0D390Es300 = numpy.ma.masked_where((D < 3.90) | (D > 4.55),Es300)
    Sws0D455Es300 = numpy.ma.masked_where((D < 4.55) | (D > 5.20),Es300)
    Sws0D520Es300 = numpy.ma.masked_where((D < 5.20) | (D > 5.85),Es300)
    Sws0D585Es300 = numpy.ma.masked_where(D < 5.85,Es300)
    
    metutils.CreateSeries(ds,'Sws0D0Es300',Sws0D0Es300,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0D65Es300',Sws0D65Es300,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0D130Es300',Sws0D130Es300,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0D195Es300',Sws0D195Es300,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0D260Es300',Sws0D260Es300,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0D325Es300',Sws0D325Es300,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0D390Es300',Sws0D390Es300,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0D455Es300',Sws0D455Es300,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0D520Es300',Sws0D520Es300,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0D585Es300',Sws0D585Es300,FList=[InVar],Descr='',Units='umol/(m2 s)')
        
    Es400 = numpy.ma.masked_where((Fsd < 400) | (Fsd > 500),Sws0)
    Sws0D0Es400 = numpy.ma.masked_where(D > 0.65,Es400)
    Sws0D65Es400 = numpy.ma.masked_where((D < 0.65) | (D > 1.3),Es400)
    Sws0D130Es400 = numpy.ma.masked_where((D < 1.3) | (D > 1.95),Es400)
    Sws0D195Es400 = numpy.ma.masked_where((D < 1.95) | (D > 2.60),Es400)
    Sws0D260Es400 = numpy.ma.masked_where((D < 2.60) | (D > 3.25),Es400)
    Sws0D325Es400 = numpy.ma.masked_where((D < 3.25) | (D > 3.90),Es400)
    Sws0D390Es400 = numpy.ma.masked_where((D < 3.90) | (D > 4.55),Es400)
    Sws0D455Es400 = numpy.ma.masked_where((D < 4.55) | (D > 5.20),Es400)
    Sws0D520Es400 = numpy.ma.masked_where((D < 5.20) | (D > 5.85),Es400)
    Sws0D585Es400 = numpy.ma.masked_where(D < 5.85,Es400)
    
    metutils.CreateSeries(ds,'Sws0D0Es400',Sws0D0Es400,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0D65Es400',Sws0D65Es400,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0D130Es400',Sws0D130Es400,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0D195Es400',Sws0D195Es400,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0D260Es400',Sws0D260Es400,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0D325Es400',Sws0D325Es400,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0D390Es400',Sws0D390Es400,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0D455Es400',Sws0D455Es400,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0D520Es400',Sws0D520Es400,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0D585Es400',Sws0D585Es400,FList=[InVar],Descr='',Units='umol/(m2 s)')
        
    Es500 = numpy.ma.masked_where((Fsd < 500) | (Fsd > 600),Sws0)
    Sws0D0Es500 = numpy.ma.masked_where(D > 0.65,Es500)
    Sws0D65Es500 = numpy.ma.masked_where((D < 0.65) | (D > 1.3),Es500)
    Sws0D130Es500 = numpy.ma.masked_where((D < 1.3) | (D > 1.95),Es500)
    Sws0D195Es500 = numpy.ma.masked_where((D < 1.95) | (D > 2.60),Es500)
    Sws0D260Es500 = numpy.ma.masked_where((D < 2.60) | (D > 3.25),Es500)
    Sws0D325Es500 = numpy.ma.masked_where((D < 3.25) | (D > 3.90),Es500)
    Sws0D390Es500 = numpy.ma.masked_where((D < 3.90) | (D > 4.55),Es500)
    Sws0D455Es500 = numpy.ma.masked_where((D < 4.55) | (D > 5.20),Es500)
    Sws0D520Es500 = numpy.ma.masked_where((D < 5.20) | (D > 5.85),Es500)
    Sws0D585Es500 = numpy.ma.masked_where(D < 5.85,Es500)
    
    metutils.CreateSeries(ds,'Sws0D0Es500',Sws0D0Es500,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0D65Es500',Sws0D65Es500,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0D130Es500',Sws0D130Es500,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0D195Es500',Sws0D195Es500,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0D260Es500',Sws0D260Es500,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0D325Es500',Sws0D325Es500,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0D390Es500',Sws0D390Es500,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0D455Es500',Sws0D455Es500,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0D520Es500',Sws0D520Es500,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0D585Es500',Sws0D585Es500,FList=[InVar],Descr='',Units='umol/(m2 s)')
        
    Es600 = numpy.ma.masked_where((Fsd < 600) | (Fsd > 700),Sws0)
    Sws0D0Es600 = numpy.ma.masked_where(D > 0.65,Es600)
    Sws0D65Es600 = numpy.ma.masked_where((D < 0.65) | (D > 1.3),Es600)
    Sws0D130Es600 = numpy.ma.masked_where((D < 1.3) | (D > 1.95),Es600)
    Sws0D195Es600 = numpy.ma.masked_where((D < 1.95) | (D > 2.60),Es600)
    Sws0D260Es600 = numpy.ma.masked_where((D < 2.60) | (D > 3.25),Es600)
    Sws0D325Es600 = numpy.ma.masked_where((D < 3.25) | (D > 3.90),Es600)
    Sws0D390Es600 = numpy.ma.masked_where((D < 3.90) | (D > 4.55),Es600)
    Sws0D455Es600 = numpy.ma.masked_where((D < 4.55) | (D > 5.20),Es600)
    Sws0D520Es600 = numpy.ma.masked_where((D < 5.20) | (D > 5.85),Es600)
    Sws0D585Es600 = numpy.ma.masked_where(D < 5.85,Es600)
    
    metutils.CreateSeries(ds,'Sws0D0Es600',Sws0D0Es600,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0D65Es600',Sws0D65Es600,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0D130Es600',Sws0D130Es600,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0D195Es600',Sws0D195Es600,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0D260Es600',Sws0D260Es600,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0D325Es600',Sws0D325Es600,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0D390Es600',Sws0D390Es600,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0D455Es600',Sws0D455Es600,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0D520Es600',Sws0D520Es600,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0D585Es600',Sws0D585Es600,FList=[InVar],Descr='',Units='umol/(m2 s)')
        
    Es700 = numpy.ma.masked_where((Fsd < 700) | (Fsd > 800),Sws0)
    Sws0D0Es700 = numpy.ma.masked_where(D > 0.65,Es700)
    Sws0D65Es700 = numpy.ma.masked_where((D < 0.65) | (D > 1.3),Es700)
    Sws0D130Es700 = numpy.ma.masked_where((D < 1.3) | (D > 1.95),Es700)
    Sws0D195Es700 = numpy.ma.masked_where((D < 1.95) | (D > 2.60),Es700)
    Sws0D260Es700 = numpy.ma.masked_where((D < 2.60) | (D > 3.25),Es700)
    Sws0D325Es700 = numpy.ma.masked_where((D < 3.25) | (D > 3.90),Es700)
    Sws0D390Es700 = numpy.ma.masked_where((D < 3.90) | (D > 4.55),Es700)
    Sws0D455Es700 = numpy.ma.masked_where((D < 4.55) | (D > 5.20),Es700)
    Sws0D520Es700 = numpy.ma.masked_where((D < 5.20) | (D > 5.85),Es700)
    Sws0D585Es700 = numpy.ma.masked_where(D < 5.85,Es700)
    
    metutils.CreateSeries(ds,'Sws0D0Es700',Sws0D0Es700,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0D65Es700',Sws0D65Es700,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0D130Es700',Sws0D130Es700,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0D195Es700',Sws0D195Es700,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0D260Es700',Sws0D260Es700,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0D325Es700',Sws0D325Es700,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0D390Es700',Sws0D390Es700,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0D455Es700',Sws0D455Es700,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0D520Es700',Sws0D520Es700,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0D585Es700',Sws0D585Es700,FList=[InVar],Descr='',Units='umol/(m2 s)')
        
    Es800 = numpy.ma.masked_where((Fsd < 800) | (Fsd > 900),Sws0)
    Sws0D0Es800 = numpy.ma.masked_where(D > 0.65,Es800)
    Sws0D65Es800 = numpy.ma.masked_where((D < 0.65) | (D > 1.3),Es800)
    Sws0D130Es800 = numpy.ma.masked_where((D < 1.3) | (D > 1.95),Es800)
    Sws0D195Es800 = numpy.ma.masked_where((D < 1.95) | (D > 2.60),Es800)
    Sws0D260Es800 = numpy.ma.masked_where((D < 2.60) | (D > 3.25),Es800)
    Sws0D325Es800 = numpy.ma.masked_where((D < 3.25) | (D > 3.90),Es800)
    Sws0D390Es800 = numpy.ma.masked_where((D < 3.90) | (D > 4.55),Es800)
    Sws0D455Es800 = numpy.ma.masked_where((D < 4.55) | (D > 5.20),Es800)
    Sws0D520Es800 = numpy.ma.masked_where((D < 5.20) | (D > 5.85),Es800)
    Sws0D585Es800 = numpy.ma.masked_where(D < 5.85,Es800)
    
    metutils.CreateSeries(ds,'Sws0D0Es800',Sws0D0Es800,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0D65Es800',Sws0D65Es800,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0D130Es800',Sws0D130Es800,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0D195Es800',Sws0D195Es800,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0D260Es800',Sws0D260Es800,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0D325Es800',Sws0D325Es800,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0D390Es800',Sws0D390Es800,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0D455Es800',Sws0D455Es800,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0D520Es800',Sws0D520Es800,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0D585Es800',Sws0D585Es800,FList=[InVar],Descr='',Units='umol/(m2 s)')
        
    Es900 = numpy.ma.masked_where(Fsd < 900,Sws0)
    Sws0D0Es900 = numpy.ma.masked_where(D > 0.65,Es900)
    Sws0D65Es900 = numpy.ma.masked_where((D < 0.65) | (D > 1.3),Es900)
    Sws0D130Es900 = numpy.ma.masked_where((D < 1.3) | (D > 1.95),Es900)
    Sws0D195Es900 = numpy.ma.masked_where((D < 1.95) | (D > 2.60),Es900)
    Sws0D260Es900 = numpy.ma.masked_where((D < 2.60) | (D > 3.25),Es900)
    Sws0D325Es900 = numpy.ma.masked_where((D < 3.25) | (D > 3.90),Es900)
    Sws0D390Es900 = numpy.ma.masked_where((D < 3.90) | (D > 4.55),Es900)
    Sws0D455Es900 = numpy.ma.masked_where((D < 4.55) | (D > 5.20),Es900)
    Sws0D520Es900 = numpy.ma.masked_where((D < 5.20) | (D > 5.85),Es900)
    Sws0D585Es900 = numpy.ma.masked_where(D < 5.85,Es900)
    
    metutils.CreateSeries(ds,'Sws0D0Es900',Sws0D0Es900,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0D65Es900',Sws0D65Es900,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0D130Es900',Sws0D130Es900,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0D195Es900',Sws0D195Es900,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0D260Es900',Sws0D260Es900,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0D325Es900',Sws0D325Es900,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0D390Es900',Sws0D390Es900,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0D455Es900',Sws0D455Es900,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0D520Es900',Sws0D520Es900,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws0D585Es900',Sws0D585Es900,FList=[InVar],Descr='',Units='umol/(m2 s)')
        
    #Sws 5-15%
    Sws5 = numpy.ma.masked_where(Sws < 0.05,Fc1)
    Es0 = numpy.ma.masked_where(Fsd > 100,Sws5)
    Sws5D0Es0 = numpy.ma.masked_where(D > 0.65,Es0)
    Sws5D65Es0 = numpy.ma.masked_where((D < 0.65) | (D > 1.3),Es0)
    Sws5D130Es0 = numpy.ma.masked_where((D < 1.3) | (D > 1.95),Es0)
    Sws5D195Es0 = numpy.ma.masked_where((D < 1.95) | (D > 2.60),Es0)
    Sws5D260Es0 = numpy.ma.masked_where((D < 2.60) | (D > 3.25),Es0)
    Sws5D325Es0 = numpy.ma.masked_where((D < 3.25) | (D > 3.90),Es0)
    Sws5D390Es0 = numpy.ma.masked_where((D < 3.90) | (D > 4.55),Es0)
    Sws5D455Es0 = numpy.ma.masked_where((D < 4.55) | (D > 5.20),Es0)
    Sws5D520Es0 = numpy.ma.masked_where((D < 5.20) | (D > 5.85),Es0)
    Sws5D585Es0 = numpy.ma.masked_where(D < 5.85,Es0)
    
    metutils.CreateSeries(ds,'Sws5D0Es0',Sws5D0Es0,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D65Es0',Sws5D65Es0,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D130Es0',Sws5D130Es0,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D195Es0',Sws5D195Es0,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D260Es0',Sws5D260Es0,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D325Es0',Sws5D325Es0,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D390Es0',Sws5D390Es0,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D455Es0',Sws5D455Es0,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D520Es0',Sws5D520Es0,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D585Es0',Sws5D585Es0,FList=[InVar],Descr='',Units='umol/(m2 s)')
        
    Es100 = numpy.ma.masked_where((Fsd < 100) | (Fsd > 200),Sws5)
    Sws5D0Es100 = numpy.ma.masked_where(D > 0.65,Es100)
    Sws5D65Es100 = numpy.ma.masked_where((D < 0.65) | (D > 1.3),Es100)
    Sws5D130Es100 = numpy.ma.masked_where((D < 1.3) | (D > 1.95),Es100)
    Sws5D195Es100 = numpy.ma.masked_where((D < 1.95) | (D > 2.60),Es100)
    Sws5D260Es100 = numpy.ma.masked_where((D < 2.60) | (D > 3.25),Es100)
    Sws5D325Es100 = numpy.ma.masked_where((D < 3.25) | (D > 3.90),Es100)
    Sws5D390Es100 = numpy.ma.masked_where((D < 3.90) | (D > 4.55),Es100)
    Sws5D455Es100 = numpy.ma.masked_where((D < 4.55) | (D > 5.20),Es100)
    Sws5D520Es100 = numpy.ma.masked_where((D < 5.20) | (D > 5.85),Es100)
    Sws5D585Es100 = numpy.ma.masked_where(D < 5.85,Es100)
    
    metutils.CreateSeries(ds,'Sws5D0Es100',Sws5D0Es100,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D65Es100',Sws5D65Es100,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D130Es100',Sws5D130Es100,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D195Es100',Sws5D195Es100,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D260Es100',Sws5D260Es100,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D325Es100',Sws5D325Es100,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D390Es100',Sws5D390Es100,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D455Es100',Sws5D455Es100,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D520Es100',Sws5D520Es100,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D585Es100',Sws5D585Es100,FList=[InVar],Descr='',Units='umol/(m2 s)')
        
    Es200 = numpy.ma.masked_where((Fsd < 200) | (Fsd > 300),Sws5)
    Sws5D0Es200 = numpy.ma.masked_where(D > 0.65,Es200)
    Sws5D65Es200 = numpy.ma.masked_where((D < 0.65) | (D > 1.3),Es200)
    Sws5D130Es200 = numpy.ma.masked_where((D < 1.3) | (D > 1.95),Es200)
    Sws5D195Es200 = numpy.ma.masked_where((D < 1.95) | (D > 2.60),Es200)
    Sws5D260Es200 = numpy.ma.masked_where((D < 2.60) | (D > 3.25),Es200)
    Sws5D325Es200 = numpy.ma.masked_where((D < 3.25) | (D > 3.90),Es200)
    Sws5D390Es200 = numpy.ma.masked_where((D < 3.90) | (D > 4.55),Es200)
    Sws5D455Es200 = numpy.ma.masked_where((D < 4.55) | (D > 5.20),Es200)
    Sws5D520Es200 = numpy.ma.masked_where((D < 5.20) | (D > 5.85),Es200)
    Sws5D585Es200 = numpy.ma.masked_where(D < 5.85,Es200)
    
    metutils.CreateSeries(ds,'Sws5D0Es200',Sws5D0Es200,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D65Es200',Sws5D65Es200,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D130Es200',Sws5D130Es200,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D195Es200',Sws5D195Es200,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D260Es200',Sws5D260Es200,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D325Es200',Sws5D325Es200,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D390Es200',Sws5D390Es200,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D455Es200',Sws5D455Es200,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D520Es200',Sws5D520Es200,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D585Es200',Sws5D585Es200,FList=[InVar],Descr='',Units='umol/(m2 s)')
        
    Es300 = numpy.ma.masked_where((Fsd < 300) | (Fsd > 400),Sws5)
    Sws5D0Es300 = numpy.ma.masked_where(D > 0.65,Es300)
    Sws5D65Es300 = numpy.ma.masked_where((D < 0.65) | (D > 1.3),Es300)
    Sws5D130Es300 = numpy.ma.masked_where((D < 1.3) | (D > 1.95),Es300)
    Sws5D195Es300 = numpy.ma.masked_where((D < 1.95) | (D > 2.60),Es300)
    Sws5D260Es300 = numpy.ma.masked_where((D < 2.60) | (D > 3.25),Es300)
    Sws5D325Es300 = numpy.ma.masked_where((D < 3.25) | (D > 3.90),Es300)
    Sws5D390Es300 = numpy.ma.masked_where((D < 3.90) | (D > 4.55),Es300)
    Sws5D455Es300 = numpy.ma.masked_where((D < 4.55) | (D > 5.20),Es300)
    Sws5D520Es300 = numpy.ma.masked_where((D < 5.20) | (D > 5.85),Es300)
    Sws5D585Es300 = numpy.ma.masked_where(D < 5.85,Es300)
    
    metutils.CreateSeries(ds,'Sws5D0Es300',Sws5D0Es300,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D65Es300',Sws5D65Es300,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D130Es300',Sws5D130Es300,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D195Es300',Sws5D195Es300,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D260Es300',Sws5D260Es300,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D325Es300',Sws5D325Es300,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D390Es300',Sws5D390Es300,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D455Es300',Sws5D455Es300,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D520Es300',Sws5D520Es300,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D585Es300',Sws5D585Es300,FList=[InVar],Descr='',Units='umol/(m2 s)')
        
    Es400 = numpy.ma.masked_where((Fsd < 400) | (Fsd > 500),Sws5)
    Sws5D0Es400 = numpy.ma.masked_where(D > 0.65,Es400)
    Sws5D65Es400 = numpy.ma.masked_where((D < 0.65) | (D > 1.3),Es400)
    Sws5D130Es400 = numpy.ma.masked_where((D < 1.3) | (D > 1.95),Es400)
    Sws5D195Es400 = numpy.ma.masked_where((D < 1.95) | (D > 2.60),Es400)
    Sws5D260Es400 = numpy.ma.masked_where((D < 2.60) | (D > 3.25),Es400)
    Sws5D325Es400 = numpy.ma.masked_where((D < 3.25) | (D > 3.90),Es400)
    Sws5D390Es400 = numpy.ma.masked_where((D < 3.90) | (D > 4.55),Es400)
    Sws5D455Es400 = numpy.ma.masked_where((D < 4.55) | (D > 5.20),Es400)
    Sws5D520Es400 = numpy.ma.masked_where((D < 5.20) | (D > 5.85),Es400)
    Sws5D585Es400 = numpy.ma.masked_where(D < 5.85,Es400)
    
    metutils.CreateSeries(ds,'Sws5D0Es400',Sws5D0Es400,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D65Es400',Sws5D65Es400,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D130Es400',Sws5D130Es400,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D195Es400',Sws5D195Es400,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D260Es400',Sws5D260Es400,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D325Es400',Sws5D325Es400,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D390Es400',Sws5D390Es400,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D455Es400',Sws5D455Es400,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D520Es400',Sws5D520Es400,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D585Es400',Sws5D585Es400,FList=[InVar],Descr='',Units='umol/(m2 s)')
        
    Es500 = numpy.ma.masked_where((Fsd < 500) | (Fsd > 600),Sws5)
    Sws5D0Es500 = numpy.ma.masked_where(D > 0.65,Es500)
    Sws5D65Es500 = numpy.ma.masked_where((D < 0.65) | (D > 1.3),Es500)
    Sws5D130Es500 = numpy.ma.masked_where((D < 1.3) | (D > 1.95),Es500)
    Sws5D195Es500 = numpy.ma.masked_where((D < 1.95) | (D > 2.60),Es500)
    Sws5D260Es500 = numpy.ma.masked_where((D < 2.60) | (D > 3.25),Es500)
    Sws5D325Es500 = numpy.ma.masked_where((D < 3.25) | (D > 3.90),Es500)
    Sws5D390Es500 = numpy.ma.masked_where((D < 3.90) | (D > 4.55),Es500)
    Sws5D455Es500 = numpy.ma.masked_where((D < 4.55) | (D > 5.20),Es500)
    Sws5D520Es500 = numpy.ma.masked_where((D < 5.20) | (D > 5.85),Es500)
    Sws5D585Es500 = numpy.ma.masked_where(D < 5.85,Es500)
    
    metutils.CreateSeries(ds,'Sws5D0Es500',Sws5D0Es500,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D65Es500',Sws5D65Es500,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D130Es500',Sws5D130Es500,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D195Es500',Sws5D195Es500,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D260Es500',Sws5D260Es500,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D325Es500',Sws5D325Es500,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D390Es500',Sws5D390Es500,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D455Es500',Sws5D455Es500,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D520Es500',Sws5D520Es500,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D585Es500',Sws5D585Es500,FList=[InVar],Descr='',Units='umol/(m2 s)')
        
    Es600 = numpy.ma.masked_where((Fsd < 600) | (Fsd > 700),Sws5)
    Sws5D0Es600 = numpy.ma.masked_where(D > 0.65,Es600)
    Sws5D65Es600 = numpy.ma.masked_where((D < 0.65) | (D > 1.3),Es600)
    Sws5D130Es600 = numpy.ma.masked_where((D < 1.3) | (D > 1.95),Es600)
    Sws5D195Es600 = numpy.ma.masked_where((D < 1.95) | (D > 2.60),Es600)
    Sws5D260Es600 = numpy.ma.masked_where((D < 2.60) | (D > 3.25),Es600)
    Sws5D325Es600 = numpy.ma.masked_where((D < 3.25) | (D > 3.90),Es600)
    Sws5D390Es600 = numpy.ma.masked_where((D < 3.90) | (D > 4.55),Es600)
    Sws5D455Es600 = numpy.ma.masked_where((D < 4.55) | (D > 5.20),Es600)
    Sws5D520Es600 = numpy.ma.masked_where((D < 5.20) | (D > 5.85),Es600)
    Sws5D585Es600 = numpy.ma.masked_where(D < 5.85,Es600)
    
    metutils.CreateSeries(ds,'Sws5D0Es600',Sws5D0Es600,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D65Es600',Sws5D65Es600,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D130Es600',Sws5D130Es600,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D195Es600',Sws5D195Es600,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D260Es600',Sws5D260Es600,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D325Es600',Sws5D325Es600,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D390Es600',Sws5D390Es600,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D455Es600',Sws5D455Es600,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D520Es600',Sws5D520Es600,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D585Es600',Sws5D585Es600,FList=[InVar],Descr='',Units='umol/(m2 s)')
        
    Es700 = numpy.ma.masked_where((Fsd < 700) | (Fsd > 800),Sws5)
    Sws5D0Es700 = numpy.ma.masked_where(D > 0.65,Es700)
    Sws5D65Es700 = numpy.ma.masked_where((D < 0.65) | (D > 1.3),Es700)
    Sws5D130Es700 = numpy.ma.masked_where((D < 1.3) | (D > 1.95),Es700)
    Sws5D195Es700 = numpy.ma.masked_where((D < 1.95) | (D > 2.60),Es700)
    Sws5D260Es700 = numpy.ma.masked_where((D < 2.60) | (D > 3.25),Es700)
    Sws5D325Es700 = numpy.ma.masked_where((D < 3.25) | (D > 3.90),Es700)
    Sws5D390Es700 = numpy.ma.masked_where((D < 3.90) | (D > 4.55),Es700)
    Sws5D455Es700 = numpy.ma.masked_where((D < 4.55) | (D > 5.20),Es700)
    Sws5D520Es700 = numpy.ma.masked_where((D < 5.20) | (D > 5.85),Es700)
    Sws5D585Es700 = numpy.ma.masked_where(D < 5.85,Es700)
    
    metutils.CreateSeries(ds,'Sws5D0Es700',Sws5D0Es700,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D65Es700',Sws5D65Es700,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D130Es700',Sws5D130Es700,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D195Es700',Sws5D195Es700,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D260Es700',Sws5D260Es700,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D325Es700',Sws5D325Es700,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D390Es700',Sws5D390Es700,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D455Es700',Sws5D455Es700,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D520Es700',Sws5D520Es700,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D585Es700',Sws5D585Es700,FList=[InVar],Descr='',Units='umol/(m2 s)')
        
    Es800 = numpy.ma.masked_where((Fsd < 800) | (Fsd > 900),Sws5)
    Sws5D0Es800 = numpy.ma.masked_where(D > 0.65,Es800)
    Sws5D65Es800 = numpy.ma.masked_where((D < 0.65) | (D > 1.3),Es800)
    Sws5D130Es800 = numpy.ma.masked_where((D < 1.3) | (D > 1.95),Es800)
    Sws5D195Es800 = numpy.ma.masked_where((D < 1.95) | (D > 2.60),Es800)
    Sws5D260Es800 = numpy.ma.masked_where((D < 2.60) | (D > 3.25),Es800)
    Sws5D325Es800 = numpy.ma.masked_where((D < 3.25) | (D > 3.90),Es800)
    Sws5D390Es800 = numpy.ma.masked_where((D < 3.90) | (D > 4.55),Es800)
    Sws5D455Es800 = numpy.ma.masked_where((D < 4.55) | (D > 5.20),Es800)
    Sws5D520Es800 = numpy.ma.masked_where((D < 5.20) | (D > 5.85),Es800)
    Sws5D585Es800 = numpy.ma.masked_where(D < 5.85,Es800)
    
    metutils.CreateSeries(ds,'Sws5D0Es800',Sws5D0Es800,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D65Es800',Sws5D65Es800,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D130Es800',Sws5D130Es800,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D195Es800',Sws5D195Es800,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D260Es800',Sws5D260Es800,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D325Es800',Sws5D325Es800,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D390Es800',Sws5D390Es800,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D455Es800',Sws5D455Es800,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D520Es800',Sws5D520Es800,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D585Es800',Sws5D585Es800,FList=[InVar],Descr='',Units='umol/(m2 s)')
        
    Es900 = numpy.ma.masked_where(Fsd < 900,Sws5)
    Sws5D0Es900 = numpy.ma.masked_where(D > 0.65,Es900)
    Sws5D65Es900 = numpy.ma.masked_where((D < 0.65) | (D > 1.3),Es900)
    Sws5D130Es900 = numpy.ma.masked_where((D < 1.3) | (D > 1.95),Es900)
    Sws5D195Es900 = numpy.ma.masked_where((D < 1.95) | (D > 2.60),Es900)
    Sws5D260Es900 = numpy.ma.masked_where((D < 2.60) | (D > 3.25),Es900)
    Sws5D325Es900 = numpy.ma.masked_where((D < 3.25) | (D > 3.90),Es900)
    Sws5D390Es900 = numpy.ma.masked_where((D < 3.90) | (D > 4.55),Es900)
    Sws5D455Es900 = numpy.ma.masked_where((D < 4.55) | (D > 5.20),Es900)
    Sws5D520Es900 = numpy.ma.masked_where((D < 5.20) | (D > 5.85),Es900)
    Sws5D585Es900 = numpy.ma.masked_where(D < 5.85,Es900)
    
    metutils.CreateSeries(ds,'Sws5D0Es900',Sws5D0Es900,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D65Es900',Sws5D65Es900,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D130Es900',Sws5D130Es900,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D195Es900',Sws5D195Es900,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D260Es900',Sws5D260Es900,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D325Es900',Sws5D325Es900,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D390Es900',Sws5D390Es900,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D455Es900',Sws5D455Es900,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D520Es900',Sws5D520Es900,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5D585Es900',Sws5D585Es900,FList=[InVar],Descr='',Units='umol/(m2 s)')
        
    log.info('NEP-Met Envelopes: All ready')

def TsBinFilters(cf,ds):
    '''Prepare data into bins for meteorological envelope calculation
    '''
    log.info('Preparing data for meteorological envelope calculation')
    if metutils.cfkeycheck(cf,Base='Params',ThisOne='dependent'):
        InVar = cf['Params']['dependent']
    else:
        InVar = 'Fc'
    
    if (InVar == 'Fc_wpl' or InVar == 'Fc'):
        Fcmg,f = metutils.GetSeriesasMA(ds,InVar)
        Fc1 = ((0-Fcmg) * (10 ** 6)) / (1000 * 44)
    elif (InVar == 'GPP'):
        Fc1,f = metutils.GetSeriesasMA(ds,InVar)
    
    Ts,f = metutils.GetSeriesasMA(ds,'Ts')
    Sws,f = metutils.GetSeriesasMA(ds,'Sws')
    Fsd,f = metutils.GetSeriesasMA(ds,'Fsd')
    D,f = metutils.GetSeriesasMA(ds,'VPD')
    nRecs = len(Fc1)
    
    Ts5 = numpy.ma.masked_where(Ts > 10,Fc1)
    Sws2Ts5 = numpy.ma.masked_where(Sws > 0.05,Ts5)
    Sws5Ts5 = numpy.ma.masked_where((Sws < 0.05) | (Sws > 0.08),Ts5)
    Sws8Ts5 = numpy.ma.masked_where((Sws < 0.08) | (Sws > 0.11),Ts5)
    Sws11Ts5 = numpy.ma.masked_where((Sws < 0.11) | (Sws > 0.14),Ts5)
    Sws14Ts5 = numpy.ma.masked_where((Sws < 0.14) | (Sws > 0.17),Ts5)
    Sws17Ts5 = numpy.ma.masked_where((Sws < 0.17) | (Sws > 0.20),Ts5)
    Sws20Ts5 = numpy.ma.masked_where((Sws < 0.20) | (Sws > 0.23),Ts5)
    Sws23Ts5 = numpy.ma.masked_where((Sws < 0.23) | (Sws > 0.26),Ts5)
    Sws26Ts5 = numpy.ma.masked_where((Sws < 0.26) | (Sws > 0.29),Ts5)
    Sws29Ts5 = numpy.ma.masked_where(Sws < 0.29,Ts5)
    
    metutils.CreateSeries(ds,'Sws2Ts5',Sws2Ts5,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5Ts5',Sws5Ts5,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws8Ts5',Sws8Ts5,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws11Ts5',Sws11Ts5,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws14Ts5',Sws14Ts5,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws17Ts5',Sws17Ts5,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws20Ts5',Sws20Ts5,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws23Ts5',Sws23Ts5,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws26Ts5',Sws26Ts5,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws29Ts5',Sws29Ts5,FList=[InVar],Descr='',Units='umol/(m2 s)')
        
    Ts10 = numpy.ma.masked_where((Ts < 10) | (Ts > 15),Fc1)
    Sws2Ts10 = numpy.ma.masked_where(Sws > 0.05,Ts10)
    Sws5Ts10 = numpy.ma.masked_where((Sws < 0.05) | (Sws > 0.08),Ts10)
    Sws8Ts10 = numpy.ma.masked_where((Sws < 0.08) | (Sws > 0.11),Ts10)
    Sws11Ts10 = numpy.ma.masked_where((Sws < 0.11) | (Sws > 0.14),Ts10)
    Sws14Ts10 = numpy.ma.masked_where((Sws < 0.14) | (Sws > 0.17),Ts10)
    Sws17Ts10 = numpy.ma.masked_where((Sws < 0.17) | (Sws > 0.20),Ts10)
    Sws20Ts10 = numpy.ma.masked_where((Sws < 0.20) | (Sws > 0.23),Ts10)
    Sws23Ts10 = numpy.ma.masked_where((Sws < 0.23) | (Sws > 0.26),Ts10)
    Sws26Ts10 = numpy.ma.masked_where((Sws < 0.26) | (Sws > 0.29),Ts10)
    Sws29Ts10 = numpy.ma.masked_where(Sws < 0.29,Ts10)
    
    metutils.CreateSeries(ds,'Sws2Ts10',Sws2Ts10,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5Ts10',Sws5Ts10,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws8Ts10',Sws8Ts10,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws11Ts10',Sws11Ts10,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws14Ts10',Sws14Ts10,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws17Ts10',Sws17Ts10,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws20Ts10',Sws20Ts10,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws23Ts10',Sws23Ts10,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws26Ts10',Sws26Ts10,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws29Ts10',Sws29Ts10,FList=[InVar],Descr='',Units='umol/(m2 s)')
        
    Ts15 = numpy.ma.masked_where((Ts < 15) | (Ts > 20),Fc1)
    Sws2Ts15 = numpy.ma.masked_where(Sws > 0.05,Ts15)
    Sws5Ts15 = numpy.ma.masked_where((Sws < 0.05) | (Sws > 0.08),Ts15)
    Sws8Ts15 = numpy.ma.masked_where((Sws < 0.08) | (Sws > 0.11),Ts15)
    Sws11Ts15 = numpy.ma.masked_where((Sws < 0.11) | (Sws > 0.14),Ts15)
    Sws14Ts15 = numpy.ma.masked_where((Sws < 0.14) | (Sws > 0.17),Ts15)
    Sws17Ts15 = numpy.ma.masked_where((Sws < 0.17) | (Sws > 0.20),Ts15)
    Sws20Ts15 = numpy.ma.masked_where((Sws < 0.20) | (Sws > 0.23),Ts15)
    Sws23Ts15 = numpy.ma.masked_where((Sws < 0.23) | (Sws > 0.26),Ts15)
    Sws26Ts15 = numpy.ma.masked_where((Sws < 0.26) | (Sws > 0.29),Ts15)
    Sws29Ts15 = numpy.ma.masked_where(Sws < 0.29,Ts15)
    
    metutils.CreateSeries(ds,'Sws2Ts15',Sws2Ts15,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5Ts15',Sws5Ts15,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws8Ts15',Sws8Ts15,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws11Ts15',Sws11Ts15,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws14Ts15',Sws14Ts15,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws17Ts15',Sws17Ts15,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws20Ts15',Sws20Ts15,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws23Ts15',Sws23Ts15,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws26Ts15',Sws26Ts15,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws29Ts15',Sws29Ts15,FList=[InVar],Descr='',Units='umol/(m2 s)')
        
    Ts20 = numpy.ma.masked_where((Ts < 20) | (Ts > 25),Fc1)
    Sws2Ts20 = numpy.ma.masked_where(Sws > 0.05,Ts20)
    Sws5Ts20 = numpy.ma.masked_where((Sws < 0.05) | (Sws > 0.08),Ts20)
    Sws8Ts20 = numpy.ma.masked_where((Sws < 0.08) | (Sws > 0.11),Ts20)
    Sws11Ts20 = numpy.ma.masked_where((Sws < 0.11) | (Sws > 0.14),Ts20)
    Sws14Ts20 = numpy.ma.masked_where((Sws < 0.14) | (Sws > 0.17),Ts20)
    Sws17Ts20 = numpy.ma.masked_where((Sws < 0.17) | (Sws > 0.20),Ts20)
    Sws20Ts20 = numpy.ma.masked_where((Sws < 0.20) | (Sws > 0.23),Ts20)
    Sws23Ts20 = numpy.ma.masked_where((Sws < 0.23) | (Sws > 0.26),Ts20)
    Sws26Ts20 = numpy.ma.masked_where((Sws < 0.26) | (Sws > 0.29),Ts20)
    Sws29Ts20 = numpy.ma.masked_where(Sws < 0.29,Ts20)
    
    metutils.CreateSeries(ds,'Sws2Ts20',Sws2Ts20,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5Ts20',Sws5Ts20,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws8Ts20',Sws8Ts20,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws11Ts20',Sws11Ts20,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws14Ts20',Sws14Ts20,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws17Ts20',Sws17Ts20,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws20Ts20',Sws20Ts20,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws23Ts20',Sws23Ts20,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws26Ts20',Sws26Ts20,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws29Ts20',Sws29Ts20,FList=[InVar],Descr='',Units='umol/(m2 s)')
        
    Ts25 = numpy.ma.masked_where((Ts < 25) | (Ts > 30),Fc1)
    Sws2Ts25 = numpy.ma.masked_where(Sws > 0.05,Ts25)
    Sws5Ts25 = numpy.ma.masked_where((Sws < 0.05) | (Sws > 0.08),Ts25)
    Sws8Ts25 = numpy.ma.masked_where((Sws < 0.08) | (Sws > 0.11),Ts25)
    Sws11Ts25 = numpy.ma.masked_where((Sws < 0.11) | (Sws > 0.14),Ts25)
    Sws14Ts25 = numpy.ma.masked_where((Sws < 0.14) | (Sws > 0.17),Ts25)
    Sws17Ts25 = numpy.ma.masked_where((Sws < 0.17) | (Sws > 0.20),Ts25)
    Sws20Ts25 = numpy.ma.masked_where((Sws < 0.20) | (Sws > 0.23),Ts25)
    Sws23Ts25 = numpy.ma.masked_where((Sws < 0.23) | (Sws > 0.26),Ts25)
    Sws26Ts25 = numpy.ma.masked_where((Sws < 0.26) | (Sws > 0.29),Ts25)
    Sws29Ts25 = numpy.ma.masked_where(Sws < 0.29,Ts25)
    
    metutils.CreateSeries(ds,'Sws2Ts25',Sws2Ts25,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5Ts25',Sws5Ts25,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws8Ts25',Sws8Ts25,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws11Ts25',Sws11Ts25,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws14Ts25',Sws14Ts25,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws17Ts25',Sws17Ts25,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws20Ts25',Sws20Ts25,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws23Ts25',Sws23Ts25,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws26Ts25',Sws26Ts25,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws29Ts25',Sws29Ts25,FList=[InVar],Descr='',Units='umol/(m2 s)')
        
    Ts30 = numpy.ma.masked_where((Ts < 30) | (Ts > 35),Fc1)
    Sws2Ts30 = numpy.ma.masked_where(Sws > 0.05,Ts30)
    Sws5Ts30 = numpy.ma.masked_where((Sws < 0.05) | (Sws > 0.08),Ts30)
    Sws8Ts30 = numpy.ma.masked_where((Sws < 0.08) | (Sws > 0.11),Ts30)
    Sws11Ts30 = numpy.ma.masked_where((Sws < 0.11) | (Sws > 0.14),Ts30)
    Sws14Ts30 = numpy.ma.masked_where((Sws < 0.14) | (Sws > 0.17),Ts30)
    Sws17Ts30 = numpy.ma.masked_where((Sws < 0.17) | (Sws > 0.20),Ts30)
    Sws20Ts30 = numpy.ma.masked_where((Sws < 0.20) | (Sws > 0.23),Ts30)
    Sws23Ts30 = numpy.ma.masked_where((Sws < 0.23) | (Sws > 0.26),Ts30)
    Sws26Ts30 = numpy.ma.masked_where((Sws < 0.26) | (Sws > 0.29),Ts30)
    Sws29Ts30 = numpy.ma.masked_where(Sws < 0.29,Ts30)
    
    metutils.CreateSeries(ds,'Sws2Ts30',Sws2Ts30,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5Ts30',Sws5Ts30,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws8Ts30',Sws8Ts30,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws11Ts30',Sws11Ts30,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws14Ts30',Sws14Ts30,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws17Ts30',Sws17Ts30,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws20Ts30',Sws20Ts30,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws23Ts30',Sws23Ts30,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws26Ts30',Sws26Ts30,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws29Ts30',Sws29Ts30,FList=[InVar],Descr='',Units='umol/(m2 s)')
        
    Ts35 = numpy.ma.masked_where((Ts < 35) | (Ts > 40),Fc1)
    Sws2Ts35 = numpy.ma.masked_where(Sws > 0.05,Ts35)
    Sws5Ts35 = numpy.ma.masked_where((Sws < 0.05) | (Sws > 0.08),Ts35)
    Sws8Ts35 = numpy.ma.masked_where((Sws < 0.08) | (Sws > 0.11),Ts35)
    Sws11Ts35 = numpy.ma.masked_where((Sws < 0.11) | (Sws > 0.14),Ts35)
    Sws14Ts35 = numpy.ma.masked_where((Sws < 0.14) | (Sws > 0.17),Ts35)
    Sws17Ts35 = numpy.ma.masked_where((Sws < 0.17) | (Sws > 0.20),Ts35)
    Sws20Ts35 = numpy.ma.masked_where((Sws < 0.20) | (Sws > 0.23),Ts35)
    Sws23Ts35 = numpy.ma.masked_where((Sws < 0.23) | (Sws > 0.26),Ts35)
    Sws26Ts35 = numpy.ma.masked_where((Sws < 0.26) | (Sws > 0.29),Ts35)
    Sws29Ts35 = numpy.ma.masked_where(Sws < 0.29,Ts35)
    
    metutils.CreateSeries(ds,'Sws2Ts35',Sws2Ts35,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5Ts35',Sws5Ts35,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws8Ts35',Sws8Ts35,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws11Ts35',Sws11Ts35,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws14Ts35',Sws14Ts35,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws17Ts35',Sws17Ts35,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws20Ts35',Sws20Ts35,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws23Ts35',Sws23Ts35,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws26Ts35',Sws26Ts35,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws29Ts35',Sws29Ts35,FList=[InVar],Descr='',Units='umol/(m2 s)')
        
    Ts40 = numpy.ma.masked_where((Ts < 40) | (Ts > 45),Fc1)
    Sws2Ts40 = numpy.ma.masked_where(Sws > 0.05,Ts40)
    Sws5Ts40 = numpy.ma.masked_where((Sws < 0.05) | (Sws > 0.08),Ts40)
    Sws8Ts40 = numpy.ma.masked_where((Sws < 0.08) | (Sws > 0.11),Ts40)
    Sws11Ts40 = numpy.ma.masked_where((Sws < 0.11) | (Sws > 0.14),Ts40)
    Sws14Ts40 = numpy.ma.masked_where((Sws < 0.14) | (Sws > 0.17),Ts40)
    Sws17Ts40 = numpy.ma.masked_where((Sws < 0.17) | (Sws > 0.20),Ts40)
    Sws20Ts40 = numpy.ma.masked_where((Sws < 0.20) | (Sws > 0.23),Ts40)
    Sws23Ts40 = numpy.ma.masked_where((Sws < 0.23) | (Sws > 0.26),Ts40)
    Sws26Ts40 = numpy.ma.masked_where((Sws < 0.26) | (Sws > 0.29),Ts40)
    Sws29Ts40 = numpy.ma.masked_where(Sws < 0.29,Ts40)
    
    metutils.CreateSeries(ds,'Sws2Ts40',Sws2Ts40,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5Ts40',Sws5Ts40,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws8Ts40',Sws8Ts40,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws11Ts40',Sws11Ts40,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws14Ts40',Sws14Ts40,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws17Ts40',Sws17Ts40,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws20Ts40',Sws20Ts40,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws23Ts40',Sws23Ts40,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws26Ts40',Sws26Ts40,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws29Ts40',Sws29Ts40,FList=[InVar],Descr='',Units='umol/(m2 s)')
        
    Ts45 = numpy.ma.masked_where(Ts < 45,Fc1)
    Sws2Ts45 = numpy.ma.masked_where(Sws > 0.05,Ts45)
    Sws5Ts45 = numpy.ma.masked_where((Sws < 0.05) | (Sws > 0.08),Ts45)
    Sws8Ts45 = numpy.ma.masked_where((Sws < 0.08) | (Sws > 0.11),Ts45)
    Sws11Ts45 = numpy.ma.masked_where((Sws < 0.11) | (Sws > 0.14),Ts45)
    Sws14Ts45 = numpy.ma.masked_where((Sws < 0.14) | (Sws > 0.17),Ts45)
    Sws17Ts45 = numpy.ma.masked_where((Sws < 0.17) | (Sws > 0.20),Ts45)
    Sws20Ts45 = numpy.ma.masked_where((Sws < 0.20) | (Sws > 0.23),Ts45)
    Sws23Ts45 = numpy.ma.masked_where((Sws < 0.23) | (Sws > 0.26),Ts45)
    Sws26Ts45 = numpy.ma.masked_where((Sws < 0.26) | (Sws > 0.29),Ts45)
    Sws29Ts45 = numpy.ma.masked_where(Sws < 0.29,Ts45)
    
    metutils.CreateSeries(ds,'Sws2Ts45',Sws2Ts45,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws5Ts45',Sws5Ts45,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws8Ts45',Sws8Ts45,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws11Ts45',Sws11Ts45,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws14Ts45',Sws14Ts45,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws17Ts45',Sws17Ts45,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws20Ts45',Sws20Ts45,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws23Ts45',Sws23Ts45,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws26Ts45',Sws26Ts45,FList=[InVar],Descr='',Units='umol/(m2 s)')
    metutils.CreateSeries(ds,'Sws29Ts45',Sws29Ts45,FList=[InVar],Descr='',Units='umol/(m2 s)')
        
    log.info('Met Envelopes: All ready')
