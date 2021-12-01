#    qcck.py  Provides L2 qc checks
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
import copy
import constants as c
import datetime
import numpy
import time
import qcts
import qcutils
import logging

log = logging.getLogger('qc.ck')

def CoordinateFluxGaps(cf,ds,Fc_in='Fc',Fe_in='Fe',Fh_in='Fh'):
    if qcutils.cfkeycheck(cf,Base='FunctionArgs',ThisOne='gapsvars'):
        vars = ast.literal_eval(cf['FunctionArgs']['gapsvars'])
        Fc_in = vars[0]
        Fe_in = vars[1]
        Fh_in = vars[2]
    Fc,flagC,a = qcutils.GetSeriesasMA(ds,Fc_in)
    Fe,flagE,a = qcutils.GetSeriesasMA(ds,Fe_in)
    Fh,flagH,a = qcutils.GetSeriesasMA(ds,Fh_in)
    index = numpy.ma.where((numpy.mod(flagC,10)!=0) | (numpy.mod(flagE,10)!=0) | (numpy.mod(flagH,10)!=0))[0]
    rC_i = numpy.ma.where((numpy.mod(flagC,10)==0) & ((numpy.mod(flagE,10)!=0) | (numpy.mod(flagH,10)!=0)))[0]
    rE_i = numpy.ma.where((numpy.mod(flagE,10)==0) & ((numpy.mod(flagC,10)!=0) | (numpy.mod(flagH,10)!=0)))[0]
    rH_i = numpy.ma.where((numpy.mod(flagH,10)==0) & ((numpy.mod(flagC,10)!=0) | (numpy.mod(flagE,10)!=0)))[0]
    ds.series[Fc_in]['Flag'][rC_i] = numpy.int32(19)
    ds.series[Fe_in]['Flag'][rE_i] = numpy.int32(19)
    ds.series[Fh_in]['Flag'][rH_i] = numpy.int32(19)
    flux_series = [Fc_in,Fe_in,Fh_in]
    for ThisOne in flux_series:
        index = numpy.where(ds.series[ThisOne]['Flag'] == 19)[0]
        ds.series[ThisOne]['Data'][index] = numpy.float64(c.missing_value)
    log.info(' Finished gap co-ordination')

def do_7500check(cf,ds):
    '''Rejects data values for series specified in LI75List for times when the Diag_7500
       flag is non-zero.  If the Diag_7500 flag is not present in the data structure passed
       to this routine, it is constructed from the QC flags of the series specified in
       LI75Lisat.  Additional checks are done for AGC_7500 (the LI-7500 AGC value),
       Ah_7500_Sd (standard deviation of absolute humidity) and Cc_7500_Sd (standard
       deviation of CO2 concentration).'''
    log.info(' Doing the 7500 check')
    LI75List = ['Ah_7500_Av','Cc_7500_Av','AhAh','CcCc','UzA','UxA','UyA','UzC','UxC','UyC']
    if 'Diag_7500' not in cf['Variables'].keys():
        ds.series[unicode('Diag_7500')] = {}
        nRecs = numpy.size(ds.series['xlDateTime']['Data'])
        ds.series['Diag_7500']['Flag'] = numpy.zeros(nRecs,dtype=numpy.int32)
        for ThisOne in ['Ah_7500_Av','Cc_7500_Av']:
            if ThisOne in ds.series.keys():
                index = numpy.where(ds.series[ThisOne]['Flag']!=0)[0]
                log.info(' do_7500check: ', ThisOne, ' rejected ',len(index))
                ds.series['Diag_7500']['Flag'] = ds.series['Diag_7500']['Flag'] + ds.series[ThisOne]['Flag']
    index = numpy.where(ds.series['Diag_7500']['Flag']!=0)
    log.info('  7500Check: Diag_7500 ' + str(numpy.size(index)))
    for ThisOne in ['AGC_7500','Ah_7500_Sd','Cc_7500_Sd','AhAh','CcCc']:
        if ThisOne in ds.series.keys():
            index = numpy.where(ds.series[ThisOne]['Flag']!=0)
            log.info('  7500Check: ' + ThisOne + ' ' + str(numpy.size(index)))
            ds.series['Diag_7500']['Flag'] = ds.series['Diag_7500']['Flag'] + ds.series[ThisOne]['Flag']
    index = numpy.where((ds.series['Diag_7500']['Flag']!=0))
    log.info('  7500Check: Total ' + str(numpy.size(index)))
    for ThisOne in LI75List:
        if ThisOne in ds.series.keys():
            ds.series[ThisOne]['Data'][index] = numpy.float64(c.missing_value)
            ds.series[ThisOne]['Flag'][index] = numpy.int32(5)
        else:
            log.error('  qcck.do_7500check: series '+str(ThisOne)+' in LI75List not found in ds.series')

def CoordinateAh7500AndFcGaps(cf,ds,Fcvar='Fc'):
    '''Cleans up Ah_7500_Av based upon Fc gaps to for QA check on Ah_7500_Av v Ah_HMP.'''
    log.info(' Doing the Ah_7500 check')
    if qcutils.cfkeycheck(cf,Base='FunctionArgs',ThisOne='AhcheckFc'):
        Fclist = ast.literal_eval(cf['FunctionArgs']['AhcheckFc'])
        Fcvar = Fclist[0]
    
    # index1  Index of bad Ah_7500_Av observations
    index1 = numpy.where((ds.series['Ah_7500_Av']['Flag']!=0) & (ds.series['Ah_7500_Av']['Flag']!=10))
    
    # index2  Index of bad Fc observations
    index2 = numpy.where((ds.series[Fcvar]['Flag']!=0) & (ds.series[Fcvar]['Flag']!=10))
    
    ds.series['Ah_7500_Av']['Data'][index2] = numpy.float64(c.missing_value)
    ds.series['Ah_7500_Av']['Flag'][index2] = ds.series[Fcvar]['Flag'][index2]
    ds.series['Ah_7500_Av']['Flag'][index1] = ds.series['Ah_7500_Av']['Flag'][index1]

def do_CSATcheck(cf,ds):
    '''Rejects data values for series specified in CSATList for times when the Diag_CSAT
       flag is non-zero.  If the Diag_CSAT flag is not present in the data structure passed
       to this routine, it is constructed from the QC flags of the series specified in
       CSATList.'''
    log.info(' Doing the CSAT check')
    if 'Wd_CSAT_Compass' in ds.series.keys():
        Wd = 'Wd_CSAT_Compass'
    else:
        Wd = 'Wd_CSAT'
    CSATList = ['Ux','Uy','Uz','Ws_CSAT',Wd,'Tv_CSAT',
                'UzT','UxT','UyT','UzA','UxA','UyA','UzC','UxC','UyC',
                'UxUz','UyUz','UxUy','UxUx','UyUy']
    if 'Diag_CSAT' not in cf['Variables'].keys():
        ds.series['Diag_CSAT']= {}
        nRecs = numpy.size(ds.series['xlDateTime']['Data'])
        ds.series['Diag_CSAT']['Flag'] = numpy.zeros(nRecs,dtype=numpy.int32)
        for ThisOne in ['Ux','Uy','Uz','Tv_CSAT']:
            if ThisOne in ds.series.keys():
                index = numpy.where(ds.series[ThisOne]['Flag']!=0)[0]
                log.info(' do_CSATcheck: ', ThisOne, ' rejected ',len(index))
                ds.series['Diag_CSAT']['Flag'] = ds.series['Diag_CSAT']['Flag'] + ds.series[ThisOne]['Flag']
    index = numpy.where(ds.series['Diag_CSAT']['Flag']!=0)
    log.info('  CSATCheck: Diag_CSAT ' + str(numpy.size(index)))
    for ThisOne in CSATList:
        if ThisOne in ds.series.keys():
            ds.series[ThisOne]['Data'][index] = numpy.float64(c.missing_value)
            ds.series[ThisOne]['Flag'][index] = numpy.int32(4)
        else:
            log.error('  qcck.do_CSATcheck: series '+str(ThisOne)+' in CSATList not found in ds.series')

def do_diurnalcheck(cf,ds,section='',series='',code=5):
    if 'DiurnalCheck' not in cf[section][series].keys(): return
    if 'NumSd' not in cf[section][series]['DiurnalCheck'].keys(): return
    if code > 80:
        if 'Level' not in cf[section][series]['DiurnalCheck'].keys(): return
        if ds.globalattributes['nc_level'] != str(cf[section][series]['DiurnalCheck']['Level']): return
    dt = numpy.float64(ds.globalattributes['time_step'])
    n = numpy.int32((60./dt) + 0.5)             #Number of timesteps per hour
    nInts = numpy.int32((1440.0/dt)+0.5)        #Number of timesteps per day
    Av = numpy.array([c.missing_value]*nInts,dtype=numpy.float64)
    Sd = numpy.array([c.missing_value]*nInts,dtype=numpy.float64)
    NSd = numpy.array(eval(cf[section][series]['DiurnalCheck']['NumSd']),dtype=numpy.float64)
    for m in range(1,13):
        mindex = numpy.where(ds.series['Month']['Data']==m)[0]
        if len(mindex)!=0:
            lHdh = ds.series['Hdh']['Data'][mindex]
            l2ds = ds.series[series]['Data'][mindex]
            for i in range(nInts):
                li = numpy.where(abs(lHdh-(numpy.float64(i)/numpy.float64(n))<c.eps)&(l2ds!=numpy.float64(c.missing_value)))
                if numpy.size(li)!=0:
                    Av[i] = numpy.mean(l2ds[li])
                    Sd[i] = numpy.std(l2ds[li])
                else:
                    Av[i] = numpy.float64(c.missing_value)
                    Sd[i] = numpy.float64(c.missing_value)
            Lwr = Av - NSd[m-1]*Sd
            Upr = Av + NSd[m-1]*Sd
            hindex = numpy.array(n*lHdh,numpy.int32)
            index = numpy.where(((l2ds!=numpy.float64(c.missing_value))&(l2ds<Lwr[hindex]))|
                                ((l2ds!=numpy.float64(c.missing_value))&(l2ds>Upr[hindex])))[0] + mindex[0]
            ds.series[series]['Data'][index] = numpy.float64(c.missing_value)
            ds.series[series]['Flag'][index] = numpy.int32(code)
            ds.series[series]['Attr']['diurnalcheck_numsd'] = cf[section][series]['DiurnalCheck']['NumSd']

def do_excludedates(cf,ds,section='',series='',code=6):
    if 'ExcludeDates' not in cf[section][series].keys(): return
    #if ds.globalattributes['nc_level'] != 'L2': return
    ldt = ds.series['DateTime']['Data']
    ExcludeList = cf[section][series]['ExcludeDates'].keys()
    NumExclude = len(ExcludeList)
    for i in range(NumExclude):
        ExcludeDateList = ast.literal_eval(cf[section][series]['ExcludeDates'][str(i)])
        try:
            si = ldt.index(datetime.datetime.strptime(ExcludeDateList[0],'%Y-%m-%d %H:%M'))
        except ValueError:
            si = 0
        try:
            ei = ldt.index(datetime.datetime.strptime(ExcludeDateList[1],'%Y-%m-%d %H:%M')) + 1
        except ValueError:
            ei = -1
        ds.series[series]['Data'][si:ei] = numpy.float64(c.missing_value)
        ds.series[series]['Flag'][si:ei] = numpy.int32(code)
        ds.series[series]['Attr']['ExcludeDates_'+str(i)] = cf[section][series]['ExcludeDates'][str(i)]

def do_excludehours(cf,ds,section='',series='',code=7):
    if 'ExcludeHours' not in cf[section][series].keys(): return
    if ds.globalattributes['nc_level'] != 'L2': return
    ldt = ds.series['DateTime']['Data']
    ExcludeList = cf[section][series]['ExcludeHours'].keys()
    NumExclude = len(ExcludeList)
    for i in range(NumExclude):
        ExcludeHourList = ast.literal_eval(cf[section][series]['ExcludeHours'][str(i)])
        try:
            si = ldt.index(datetime.datetime.strptime(ExcludeHourList[0],'%Y-%m-%d %H:%M'))
        except ValueError:
            si = 0
        try:
            ei = ldt.index(datetime.datetime.strptime(ExcludeHourList[1],'%Y-%m-%d %H:%M')) + 1
        except ValueError:
            ei = -1
        for j in range(len(ExcludeHourList[2])):
            ExHr = datetime.datetime.strptime(ExcludeHourList[2][j],'%H:%M').hour
            ExMn = datetime.datetime.strptime(ExcludeHourList[2][j],'%H:%M').minute
            index = numpy.where((ds.series['Hour']['Data'][si:ei]==ExHr)&
                                (ds.series['Minute']['Data'][si:ei]==ExMn))[0] + si
            ds.series[series]['Data'][index] = numpy.float64(c.missing_value)
            ds.series[series]['Flag'][index] = numpy.int32(code)
            ds.series[series]['Attr']['ExcludeHours_'+str(i)] = cf[section][series]['ExcludeHours'][str(i)]

def do_linear(cf,ds):
    log.info(' Applying linear corrections ...')
    level = ds.globalattributes['nc_level']
    for ThisOne in cf['Variables'].keys():
        if qcutils.haskey(cf,ThisOne,'Linear'):
            qcts.ApplyLinear(cf,ds,ThisOne)
        if qcutils.haskey(cf,ThisOne,'Drift'):
            qcts.ApplyLinearDrift(cf,ds,ThisOne)
        if qcutils.haskey(cf,ThisOne,'LocalDrift'):
            qcts.ApplyLinearDriftLocal(cf,ds,ThisOne)

def do_rangecheck(cf,ds,section='',series='',code=2):
    '''Applies a range check to data series listed in the control file.  Data values that
       are less than the lower limit or greater than the upper limit are replaced with
       c.missing_value and the corresponding QC flag element is set to 2.'''
    if 'RangeCheck' not in cf[section][series].keys(): return
    if code > 80:
        if 'Level' not in cf[section][series]['RangeCheck'].keys(): return
        #log.info('   nc_level: '+ds.globalattributes['nc_level']+'; cf_level: '+str(cf[section][series]['Level']))
        if ds.globalattributes['nc_level'] != str(cf[section][series]['RangeCheck']['Level']): return
    if 'Lower' in cf[section][series]['RangeCheck'].keys():
        lwr = numpy.array(eval(cf[section][series]['RangeCheck']['Lower']))
        valid_lower = numpy.min(lwr)
        lwr = lwr[ds.series['Month']['Data']-1]
        index = numpy.where((abs(ds.series[series]['Data']-numpy.float64(c.missing_value))>c.eps)&
                                (ds.series[series]['Data']<lwr))
        ds.series[series]['Data'][index] = numpy.float64(c.missing_value)
        ds.series[series]['Flag'][index] = numpy.int32(code)
        ds.series[series]['Attr']['rangecheck_lower'] = cf[section][series]['RangeCheck']['Lower']
    if 'Upper' in cf[section][series]['RangeCheck'].keys():
        upr = numpy.array(eval(cf[section][series]['RangeCheck']['Upper']))
        valid_upper = numpy.min(upr)
        upr = upr[ds.series['Month']['Data']-1]
        index = numpy.where((abs(ds.series[series]['Data']-numpy.float64(c.missing_value))>c.eps)&
                                (ds.series[series]['Data']>upr))
        ds.series[series]['Data'][index] = numpy.float64(c.missing_value)
        ds.series[series]['Flag'][index] = numpy.int32(code)
        ds.series[series]['Attr']['rangecheck_upper'] = cf[section][series]['RangeCheck']['Upper']
        ds.series[series]['Attr']['valid_range'] = str(valid_lower)+','+str(valid_upper)

def do_qcchecks(cf,ds):
    level = ds.globalattributes['nc_level']
    log.info(' Doing the QC checks at level '+str(level))
    for series in ds.series.keys():
            do_qcchecks_oneseries(cf,ds,series=series)

def do_qcchecks_oneseries(cf,ds,series=''):
    section = qcutils.get_cfsection(cf,series=series,mode='quiet')
    if len(section)==0: return
    level = ds.globalattributes['nc_level']
    if level == 'L2':
        range_code        = 2
        diurnal_code      = 3
        excludedates_code = 6
        excludehours_code = 7
    if level == 'L3':
        excludedates_code = 6
        excludehours_code = 7
        range_code        = 16
        diurnal_code      = 17
    if level == 'L4':
        excludedates_code = 6
        excludehours_code = 7
        range_code        = 82
        diurnal_code      = 83
    if level == 'L5':
        excludedates_code = 6
        excludehours_code = 7
        range_code        = 84
        diurnal_code      = 85
    if level == 'L6':
        excludedates_code = 6
        excludehours_code = 7
        range_code        = 86
        diurnal_code      = 87
    # do the range check
    do_rangecheck(cf,ds,section=section,series=series,code=range_code)
    # do the diurnal check
    do_diurnalcheck(cf,ds,section=section,series=series,code=diurnal_code)
    # do exclude dates
    do_excludedates(cf,ds,section=section,series=series,code=excludedates_code)
    # do exclude hours
    do_excludehours(cf,ds,section=section,series=series,code=excludehours_code)
    try:
        if 'do_qcchecks' not in ds.globalattributes['Functions']:
            ds.globalattributes['Functions'] = ds.globalattributes['Functions']+',do_qcchecks'
    except:
        ds.globalattributes['Functions'] = 'do_qcchecks'
