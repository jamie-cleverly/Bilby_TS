#    qcls.py  Queue up processes for each level
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

import sys
import logging
import ast
import constants as c
import copy
import numpy
import qcck
import qcio
import qcts
import qcutils
import time
import xlrd
import meteorologicalfunctions as mf
import cfg

log = logging.getLogger('qc.ls')

def l1qc(cf,ds):
    ds.globalattributes['nc_level'] = 'L1'
    ds.globalattributes['EPDversion'] = sys.version
    ds.globalattributes['QC_version_history'] = cfg.__doc__
    # put the control file name into the global attributes
    ds.globalattributes['controlfile_name'] = cf['controlfile_name']
    # write the percentage of good data as a variable attribute
    qcutils.get_coverage_individual(ds)
    # write the percentage of good data for groups
    qcutils.get_coverage_groups(ds)
    return ds

def l2qc(cf,ds1):
    """
        Perform initial QA/QC on flux data
        Generates L2 from L1 data
        * check parameters specified in control file
        
        Functions performed:
            qcck.do_rangecheck*
            qcck.do_CSATcheck
            qcck.do_7500check
            qcck.do_diurnalcheck*
            qcck.do_excludedates*
            qcck.do_excludehours*
            qcts.albedo
        """
    # make a copy of the L1 data
    ds2 = copy.deepcopy(ds1)
    ds2.globalattributes['nc_level'] = 'L2'
    ds2.globalattributes['EPDversion'] = sys.version
    ds2.globalattributes['QC_version_history'] = cfg.__doc__
    ds2.globalattributes['L2Functions'] = 'do_qccheck(RangeCheck, diurnalcheck, excludedates, excludehours), CSATcheck, 7500check, albedo'
    ds2.globalattributes['Functions'] = 'do_qccheck(RangeCheck, diurnalcheck, excludedates, excludehours), CSATcheck, 7500check, albedo'
    # put the control file name into the global attributes
    ds2.globalattributes['controlfile_name'] = cf['controlfile_name']
    # apply the quality control checks (range, diurnal, exclude dates and exclude hours
    qcck.do_qcchecks(cf,ds2)
    # do the CSAT diagnostic check
    qcck.do_CSATcheck(cf,ds2)
    # do the LI-7500 diagnostic check
    qcck.do_7500check(cf,ds2)
    # constrain albedo estimates to full sun angles
    qcts.albedo(cf,ds2)
    log.info(' Finished the albedo constraints')    # apply linear corrections to the data
    log.info(' Applying linear corrections ...')
    qcck.do_linear(cf,ds2)
    # write series statistics to file
    qcio.get_seriesstats(cf,ds2)
    # write the percentage of good data as a variable attribute
    qcutils.get_coverage_individual(ds2)
    # write the percentage of good data for groups
    qcutils.get_coverage_groups(ds2)
    return ds2

def l3qc(cf,ds2):
    """
        Corrections
        Generates L3 from L2 data
        
        Functions performed:
            qcts.AddMetVars (optional)
            qcts.CorrectSWC (optional*)
            qcck.do_linear (all sites)
            qcutils.GetMergeList + qcts.MergeSeries Ah_EC (optional)x
            qcts.TaFromTv (optional)
            qcutils.GetMergeList + qcts.MergeSeries Ta_EC (optional)x
            qcts.CoordRotation2D (all sites)
            qcts.MassmanApprox (optional*)y
            qcts.Massman (optional*)y
            qcts.CalculateFluxes (used if Massman not optioned)x
            qcts.CalculateFluxesRM (used if Massman optioned)y
            qcts.FhvtoFh (all sites)
            qcts.Fe_WPL (WPL computed on fluxes, as with Campbell algorithm)+x
            qcts.Fc_WPL (WPL computed on fluxes, as with Campbell algorithm)+x
            qcts.Fe_WPLcov (WPL computed on kinematic fluxes (ie, covariances), as with WPL80)+y
            qcts.Fc_WPLcov (WPL computed on kinematic fluxes (ie, covariances), as with WPL80)+y
            qcts.CalculateNetRadiation (optional)
            qcutils.GetMergeList + qcts.MergeSeries Fsd (optional)
            qcutils.GetMergeList + qcts.MergeSeries Fn (optional*)
            qcts.InterpolateOverMissing (optional)
            AverageSeriesByElements (optional)
            qcts.CorrectFgForStorage (all sites)
            qcts.Average3SeriesByElements (optional)
            qcts.CalculateAvailableEnergy (optional)
            qcck.do_qcchecks (all sites)
            qcck.gaps (optional)
            
            *:  requires ancillary measurements for paratmerisation
            +:  each site requires one pair, either Fe_WPL & Fc_WPL (default) or Fe_WPLCov & FcWPLCov
            x:  required together in option set
            y:  required together in option set
        """
    # make a copy of the L2 data
    ds3 = copy.deepcopy(ds2)
    ds3.globalattributes['nc_level'] = 'L3'
    ds3.globalattributes['EPDversion'] = sys.version
    ds3.globalattributes['QC_version_history'] = cfg.__doc__
    # put the control file name into the global attributes
    ds3.globalattributes['controlfile_name'] = cf['controlfile_name']
    
    # calculate NDVI
    if qcutils.cfkeycheck(cf,Base='Functions',ThisOne='NDVI') and cf['Functions']['NDVI'] == 'True':
        try:
            ds3.globalattributes['L3Functions'] = ds3.globalattributes['L3Functions']+', calculateNDVI'
        except:
            ds3.globalattributes['L3Functions'] = 'calculateNDVI'
        
        log.info(' Calculating NDVI from component reflectances ...')
        qcts.CalculateNDVI(cf,ds3)
    
    # bypass soil temperature correction for Sws (when Ts bad)
    if qcutils.cfkeycheck(cf,Base='Functions',ThisOne='BypassSwsTcorr') and cf['Functions']['BypassSwsTcorr'] == 'True':
        try:
            ds3.globalattributes['L3Functions'] = ds3.globalattributes['L3Functions']+', BypassSwsTcorr'
        except:
            ds3.globalattributes['L3Functions'] = 'BypassSwsTcorr'
        
        log.info(' Re-computing Sws without temperature correction ...')
        qcts.BypassTcorr(cf,ds3)
    
    # correct measured soil water content using empirical relationship to collected samples
    if qcutils.cfkeycheck(cf,Base='Functions',ThisOne='CorrectSWC') and cf['Functions']['CorrectSWC'] == 'True':
        try:
            ds3.globalattributes['L3Functions'] = ds3.globalattributes['L3Functions']+', CorrectSWC'
        except:
            ds3.globalattributes['L3Functions'] = 'CorrectSWC'
        
        log.info(' Correcting soil moisture data ...')
        qcts.CorrectSWC(cf,ds3)
    
    # apply linear corrections to the data
    if qcutils.cfkeycheck(cf,Base='Functions',ThisOne='Corrections') and cf['Functions']['Corrections'] == 'True':
        try:
            ds3.globalattributes['L3Functions'] = ds3.globalattributes['L3Functions']+', do_linear'
        except:
            ds3.globalattributes['L3Functions'] = 'do_linear'
        
        log.info(' Applying linear corrections ...')
        qcck.do_linear(cf,ds3)
    
    # determine HMP Ah if not output by datalogger
    if qcutils.cfkeycheck(cf,Base='Functions',ThisOne='CalculateAh') and cf['Functions']['CalculateAh'] == 'True':
        try:
            ds3.globalattributes['L3Functions'] = ds3.globalattributes['L3Functions']+', CalculateAh'
        except:
            ds3.globalattributes['L3Functions'] = 'CalculateAh'
        
        log.info(' Adding HMP Ah to database')
        qcts.CalculateAhHMP(cf,ds3)
    
    # merge the HMP and corrected 7500 data
    if qcutils.cfkeycheck(cf,Base='Functions',ThisOne='MergeSeriesAhTa') and cf['Functions']['MergeSeriesAhTa'] == 'True':
        try:
            ds3.globalattributes['L3Functions'] = ds3.globalattributes['L3Functions']+', MergeSeriesAhTaCc'
        except:
            ds3.globalattributes['L3Functions'] = 'MergeSeriesAhTaCc'
        
        qcts.MergeSeries(cf,ds3,'Ah',[0,10])
        qcts.MergeSeries(cf,ds3,'Cc',[0,10])
        
    # get the air temperature from the CSAT virtual temperature
        try:
            ds3.globalattributes['L3Functions'] = ds3.globalattributes['L3Functions']+', TaFromTv'
        except:
            ds3.globalattributes['L3Functions'] = 'TaFromTv'
        
        qcts.TaFromTv(cf,ds3)
        
    # merge the HMP and corrected CSAT data
        qcts.MergeSeries(cf,ds3,'Ta',[0,10])
    
    # add relevant meteorological values to L3 data
    if (qcutils.cfkeycheck(cf,Base='Functions',ThisOne='Corrections') and cf['Functions']['Corrections'] == 'True') or (qcutils.cfkeycheck(cf,Base='Functions',ThisOne='CalculateMetVars') and cf['Functions']['CalculateMetVars'] == 'True'):
        try:
            ds3.globalattributes['L3Functions'] = ds3.globalattributes['L3Functions']+', CalculateMetVars'
        except:
            ds3.globalattributes['L3Functions'] = 'CalculateMetVars'
            
        log.info(' Adding standard met variables to database')
        qcts.CalculateMeteorologicalVariables(ds3)
    
    # do the 2D coordinate rotation
    if qcutils.cfkeycheck(cf,Base='Functions',ThisOne='Corrections') and cf['Functions']['Corrections'] == 'True':
        try:
            ds3.globalattributes['L3Functions'] = ds3.globalattributes['L3Functions']+', CoordRotation2D'
        except:
            ds3.globalattributes['L3Functions'] = 'CoordRotation2D'
        
        qcts.CoordRotation2D(cf,ds3)
    
    # do the Massman frequency attenuation correction
    if qcutils.cfkeycheck(cf,Base='Functions',ThisOne='Corrections') and cf['Functions']['Corrections'] == 'True':
        try:
            ds3.globalattributes['L3Functions'] = ds3.globalattributes['L3Functions']+', Massman'
        except:
            ds3.globalattributes['L3Functions'] = 'Massman'
        
        qcts.MassmanStandard(cf,ds3)
    
    # calculate the fluxes
    if qcutils.cfkeycheck(cf,Base='Functions',ThisOne='Corrections') and cf['Functions']['Corrections'] == 'True':
        try:
            ds3.globalattributes['L3Functions'] = ds3.globalattributes['L3Functions']+', CalculateFluxes'
        except:
            ds3.globalattributes['L3Functions'] = 'CalculateFluxes'
        
        qcts.CalculateFluxes(cf,ds3)
    
    # approximate wT from virtual wT using wA (ref: Campbell OPECSystem manual)
    if qcutils.cfkeycheck(cf,Base='Functions',ThisOne='Corrections') and cf['Functions']['Corrections'] == 'True':
        try:
            ds3.globalattributes['L3Functions'] = ds3.globalattributes['L3Functions']+', FhvtoFh'
        except:
            ds3.globalattributes['L3Functions'] = 'FhvtoFh'
        
        qcts.FhvtoFh(cf,ds3)
    
    # correct the H2O & CO2 flux due to effects of flux on density measurements
    if qcutils.cfkeycheck(cf,Base='Functions',ThisOne='Corrections') and cf['Functions']['Corrections'] == 'True':
        if qcutils.cfkeycheck(cf,Base='Functions',ThisOne='WPLcov') and cf['Functions']['WPLcov'] == 'True':
            try:
                ds3.globalattributes['L3Functions'] = ds3.globalattributes['L3Functions']+', WPLcov'
            except:
                ds3.globalattributes['L3Functions'] = 'WPLcov'
            
            qcts.do_WPL(cf,ds3,cov='True')
        else:
            try:
                ds3.globalattributes['L3Functions'] = ds3.globalattributes['L3Functions']+', WPL'
            except:
                ds3.globalattributes['L3Functions'] = 'WPL'
            
            qcts.do_WPL(cf,ds3)
    
    # calculate the net radiation from the Kipp and Zonen CNR1
    if qcutils.cfkeycheck(cf,Base='Functions',ThisOne='CalculateNetRadiation') and cf['Functions']['CalculateNetRadiation'] == 'True':
        try:
            ds3.globalattributes['L3Functions'] = ds3.globalattributes['L3Functions']+', CalculateNetRadiation'
        except:
            ds3.globalattributes['L3Functions'] = 'CalculateNetRadiation'
        
        qcts.MergeSeries(cf,ds3,'Fsd',[0,10])
        qcts.CalculateNetRadiation(ds3,'Fn_KZ','Fsd','Fsu','Fld','Flu')
        qcts.MergeSeries(cf,ds3,'Fn',[0,10])
    
    # combine wind speed from the CSAT and the Wind Sentry
    if qcutils.cfkeycheck(cf,Base='Functions',ThisOne='MergeSeriesWS') and cf['Functions']['MergeSeriesWS'] == 'True':
        try:
            ds3.globalattributes['L3Functions'] = ds3.globalattributes['L3Functions']+', MergeSeriesWS'
        except:
            ds3.globalattributes['L3Functions'] = 'MergeSeriesWS'
        
        qcts.MergeSeries(cf,ds3,'Ws',[0,10])
    
    # average the soil temperature data
    if qcutils.cfkeycheck(cf,Base='Functions',ThisOne='Corrections') and cf['Functions']['Corrections'] == 'True':
        if 'SoilAverage' not in ds3.globalattributes['L3Functions']:
            try:
                ds3.globalattributes['L3Functions'] = ds3.globalattributes['L3Functions']+', SoilAverage'
            except:
                ds3.globalattributes['L3Functions'] = 'SoilAverage'
            
        # interpolate over any ramaining gaps up to 3 hours in length
        qcts.AverageSeriesByElementsI(cf,ds3,'Ts')
        qcts.AverageSeriesByElementsI(cf,ds3,'Sws')
        
    # correct the measured soil heat flux for storage in the soil layer above the sensor
    if qcutils.cfkeycheck(cf,Base='Functions',ThisOne='Corrections') and cf['Functions']['Corrections'] == 'True':
        try:
            ds3.globalattributes['L3Functions'] = ds3.globalattributes['L3Functions']+', CorrectFgForStorage'
        except:
            ds3.globalattributes['L3Functions'] = 'CorrectFgForStorage'
        
        if qcutils.cfkeycheck(cf,Base='Functions',ThisOne='IndividualFgCorrection') and cf['Functions']['IndividualFgCorrection'] == 'True':
            qcts.CorrectIndividualFgForStorage(cf,ds3)
            qcts.AverageSeriesByElementsI(cf,ds3,'Fg')
        else:
            qcts.AverageSeriesByElementsI(cf,ds3,'Fg')
            qcts.CorrectGroupFgForStorage(cf,ds3)
    
    # calculate the available energy
    if qcutils.cfkeycheck(cf,Base='Functions',ThisOne='CalculateAvailableEnergy') and cf['Functions']['CalculateAvailableEnergy'] == 'True':
        try:
            ds3.globalattributes['L3Functions'] = ds3.globalattributes['L3Functions']+', CalculateAvailableEnergy'
        except:
            ds3.globalattributes['L3Functions'] = 'CalculateAvailableEnergy'
        
        qcts.CalculateAvailableEnergy(ds3)
    
    if qcutils.cfkeycheck(cf,Base='Functions',ThisOne='DiagnosticMode'):
        if cf['Functions']['DiagnosticMode'] == 'False':
            qcutils.prepOzFluxVars(cf,ds3)
    else:
        qcutils.prepOzFluxVars(cf,ds3)
    
    # calculate specific humidity and saturated specific humidity profile
    if qcutils.cfkeycheck(cf,Base='Functions',ThisOne='qTprofile') and cf['Functions']['qTprofile'] == 'True':
        try:
            ds3.globalattributes['L3Functions'] = ds3.globalattributes['L3Functions']+', qTprofile'
        except:
            ds3.globalattributes['L3Functions'] = 'qTprofile'
        
        qcts.CalculateSpecificHumidityProfile(cf,ds3)
    
    # calculate Penman-Monteith inversion
    if qcutils.cfkeycheck(cf,Base='Functions',ThisOne='PenmanMonteith') and cf['Functions']['PenmanMonteith'] == 'True':
        try:
            ds3.globalattributes['L3Functions'] = ds3.globalattributes['L3Functions']+', PenmanMonteith'
        except:
            ds3.globalattributes['L3Functions'] = 'PenmanMonteith'
        
        qcts.do_PenmanMonteith(cf,ds3)
    
    # calculate bulk Richardson numbers
    if qcutils.cfkeycheck(cf,Base='Functions',ThisOne='bulkRichardson') and cf['Functions']['bulkRichardson'] == 'True':
        try:
            ds3.globalattributes['L3Functions'] = ds3.globalattributes['L3Functions']+', bulkRichardson'
        except:
            ds3.globalattributes['L3Functions'] = 'bulkRichardson'
        
        qcts.do_bulkRichardson(cf,ds3)
    
    # re-apply the quality control checks (range, diurnal and rules)
    if qcutils.cfkeycheck(cf,Base='Functions',ThisOne='Corrections') and cf['Functions']['Corrections'] == 'True':
        ds3.globalattributes['L3Functions'] = ds3.globalattributes['L3Functions']+', do_qccheck(RangeCheck, diurnalcheck, excludedates, excludehours)'
        qcck.do_qcchecks(cf,ds3)
    
    # quality control checks (range, diurnal and rules) without flux post-processing
    if qcutils.cfkeycheck(cf,Base='Functions',ThisOne='QCChecks') and cf['Functions']['QCChecks'] == 'True':
        qcck.do_qcchecks(cf,ds3)
    
    # apply the ustar filter
    if qcutils.cfkeycheck(cf,Base='Functions',ThisOne='ustarFilter') and cf['Functions']['ustarFilter'] == 'True':
        try:
            ds3.globalattributes['L3Functions'] = ds3.globalattributes['L3Functions']+', ustarFilter'
        except:
            ds3.globalattributes['L3Functions'] = 'ustarFilter'
        
        qcts.FilterFcByUstar(cf,ds3)
    
    # coordinate gaps in the three main fluxes
    if qcutils.cfkeycheck(cf,Base='Functions',ThisOne='CoordinateFluxGaps') and cf['Functions']['CoordinateFluxGaps'] == 'True':
        try:
            ds3.globalattributes['L3Functions'] = ds3.globalattributes['L3Functions']+', CoordinateFluxGaps'
        except:
            ds3.globalattributes['L3Functions'] = 'CoordinateFluxGaps'
        
        qcck.CoordinateFluxGaps(cf,ds3)
    
    # coordinate gaps in Ah_7500_Av with Fc
    if qcutils.cfkeycheck(cf,Base='Functions',ThisOne='Corrections') and cf['Functions']['Corrections'] == 'True':
        try:
            ds3.globalattributes['L3Functions'] = ds3.globalattributes['L3Functions']+', CoordinateAh7500AndFcGaps'
        except:
            ds3.globalattributes['L3Functions'] = 'CoordinateAh7500AndFcGaps'
        
        qcck.CoordinateAh7500AndFcGaps(cf,ds3)
    
    # calcluate ET at observation interval
    if qcutils.cfkeycheck(cf,Base='Functions',ThisOne='CalculateET') and cf['Functions']['CalculateET'] == 'True':
        try:
            ds3.globalattributes['L3Functions'] = ds3.globalattributes['L3Functions']+', CalculateET'
        except:
            ds3.globalattributes['L3Functions'] = 'CalculateET'
        
        log.info(' Calculating ET')
        qcts.CalculateET(cf,ds3,'L3')
    
    # run MOST (Buckingham Pi) 2d footprint model (Kljun et al. 2004)
    if qcutils.cfkeycheck(cf,Base='Functions',ThisOne='footprint') and cf['Functions']['footprint'] == 'True':
        try:
            ds3.globalattributes['L3Functions'] = ds3.globalattributes['L3Functions']+', footprint'
        except:
            ds3.globalattributes['L3Functions'] = 'footprint'
        
        qcts.do_footprint_2d(cf,ds3)
    
    if qcutils.cfkeycheck(cf,Base='Functions',ThisOne='Corrections') and cf['Functions']['Corrections'] == 'True':
        qcio.get_seriesstats(cf,ds3)
    
    # convert Fc [mgCO2 m-2 s-1] to Fc_co2 [mgCO2 m-2 s-1], Fc_c [mgC m-2 s-1], NEE [umol m-2 s-1] and NEP = - NEE
    if qcutils.cfkeycheck(cf,Base='Functions',ThisOne='convertFc') and cf['Functions']['convertFc'] == 'True':
        try:
            ds3.globalattributes['L3Functions'] = ds3.globalattributes['L3Functions']+', convertFc'
        except:
            ds3.globalattributes['L3Functions'] = 'convertFc'
        
        qcts.ConvertFc(cf,ds3)
    
    # convert Fc [mgCO2 m-2 s-1] to Fc [umol m-2 s-1]
    if qcutils.cfkeycheck(cf,Base='Functions',ThisOne='JasonFc') and cf['Functions']['JasonFc'] == 'True':
        try:
            ds3.globalattributes['L3Functions'] = ds3.globalattributes['L3Functions']+', convertFc (umol only)'
        except:
            ds3.globalattributes['L3Functions'] = 'convertFc (umol only)'
        
        qcts.ConvertFcJason(cf,ds3)
    
    # write the percentage of good data as a variable attribute
    qcutils.get_coverage_individual(ds3)
    # write the percentage of good data for groups
    qcutils.get_coverage_groups(ds3)
    
    # compute water-use efficiency from flux-gradient similarity (appendix A, Scanlon & Sahu 2008)
    if qcutils.cfkeycheck(cf,Base='Functions',ThisOne='wue') and cf['Functions']['wue'] == 'True':
        try:
            ds3.globalattributes['L3Functions'] = ds3.globalattributes['L3Functions']+', wue'
        except:
            ds3.globalattributes['L3Functions'] = 'wue'
        
        log.info(' Calculating water-use efficiency from flux-gradient similarity')
        qcts.CalculateWUEfromSimilarity(cf,ds3)
    
    # compute climatology for L3 data
    if qcutils.cfkeycheck(cf,Base='Functions',ThisOne='climatology') and cf['Functions']['climatology'] == 'True':
        try:
            ds3.globalattributes['L3Functions'] = ds3.globalattributes['L3Functions']+', climatology'
        except:
            ds3.globalattributes['L3Functions'] = 'climatology'
        
        qcts.do_climatology(cf,ds3)
    
    if qcutils.cfkeycheck(cf,Base='Functions',ThisOne='Sums') and cf['Functions']['Sums'] == 'L3':
        try:
            ds3.globalattributes['L3Functions'] = ds3.globalattributes['L5Functions']+', Sums'
        except:
            ds3.globalattributes['L3Functions'] = 'Sums'
        
        qcts.do_sums(cf,ds3)
    
    try:
        ds3.globalattributes['Functions'] = ds3.globalattributes['Functions'] + ', ' + ds3.globalattributes['L3Functions']
    except:
        ds3.globalattributes['Functions'] = ds3.globalattributes['L3Functions']
    
    return ds3

def l4qc(cf,ds3,InLevel,x):
    ds4 = copy.deepcopy(ds3)
    ds4.globalattributes['nc_level'] = 'L4'
    if (qcutils.cfkeycheck(cf,Base='Functions',ThisOne='L4_offline') and cf['Functions']['L4_offline'] == 'True') and qcutils.cfkeycheck(cf,Base='Functions',ThisOne='L4_keys'):
        try:
            ds4.globalattributes['L4Functions'] = ds4.globalattributes['L4Functions']+', '+cf['Functions']['L4_keys']
        except:
            ds4.globalattributes['L4Functions'] = cf['Functions']['L4_keys']
        x=x+1
    
    # determine HMP Ah if not output by datalogger
    if qcutils.cfkeycheck(cf,Base='Functions',ThisOne='CalculateAh') and cf['Functions']['CalculateAh'] == 'True':
        try:
            ds4.globalattributes['L4Functions'] = ds4.globalattributes['L4Functions']+', CalculateAh'
        except:
            ds4.globalattributes['L4Functions'] = 'CalculateAh'
        
        log.info(' Adding HMP Ah to database')
        qcts.CalculateAhHMP(cf,ds4)
    
    # add relevant meteorological values to L4 data
    if qcutils.cfkeycheck(cf,Base='Functions',ThisOne='CalculateMetVars') and cf['Functions']['CalculateMetVars'] == 'True':
        try:
            ds4.globalattributes['L4Functions'] = ds4.globalattributes['L4Functions']+', CalculateMetVars'
        except:
            ds4.globalattributes['L4Functions'] = 'CalculateMetVars'
        
        log.info(' Adding standard met variables to database')
        qcts.CalculateMeteorologicalVariables(ds4)
        
    # merge CSAT and wind sentry wind speed
    if qcutils.cfkeycheck(cf,Base='Functions',ThisOne='MergeSeriesWS') and cf['Functions']['MergeSeriesWS'] == 'True':
        try:
            ds4.globalattributes['L4Functions'] = ds4.globalattributes['L4Functions']+', MergeSeriesWS'
        except:
            ds4.globalattributes['L4Functions'] = 'MergeSeriesWS'
        
        qcts.MergeSeries(cf,ds4,'Ws',[0,10])
    
    # linear interpolation to fill missing values over gaps of 1 hour
    if qcutils.cfkeycheck(cf,Base='Functions',ThisOne='InterpolateOverMissing') and cf['Functions']['InterpolateOverMissing'] == 'True':
        try:
            ds4.globalattributes['L4Functions'] = ds4.globalattributes['L4Functions']+', InterpolateOverMissing'
        except:
            ds4.globalattributes['L4Functions'] = 'InterpolateOverMissing'
        
        log.info(' Gap filling by linear interpolation to fill missing values over gaps of 1 hour')
        for ThisOne in cf['InterpolateVars'].keys():
            qcts.InterpolateOverMissing(cf,ds4,series=ThisOne,maxlen=2)
        x=x+1
    
    # re-apply the quality control checks (range, diurnal and rules)
    if x > 0:
        log.info(' Doing QC checks on L4 data')
        qcck.do_qcchecks(cf,ds4)
        try:
            ds4.globalattributes['L4Functions'] = ds4.globalattributes['L4Functions']+', do_qccheck(RangeCheck, diurnalcheck, excludedates, excludehours)'
        except:
            ds4.globalattributes['L4Functions'] = 'do_qccheck(RangeCheck, diurnalcheck, excludedates, excludehours)'
    
    # interpolate over any ramaining gaps up to 3 hours in length
    if qcutils.cfkeycheck(cf,Base='Functions',ThisOne='InterpolateOverMissing') and cf['Functions']['InterpolateOverMissing'] == 'True':
        for ThisOne in cf['InterpolateVars'].keys():
            qcts.InterpolateOverMissing(cf,ds4,series=ThisOne,maxlen=6)
    
    ds4.globalattributes['Functions'] = ds4.globalattributes['Functions'] + ', ' + ds4.globalattributes['L4Functions']
    return ds4,x

def l4to6qc(cf,ds3,AttrLevel,InLevel,OutLevel):
    """
        Fill gaps in met data from other sources
        Integrate SOLO-ANN gap filled fluxes performed externally
        Generates L4 from L3 data
        Generates daily sums excel workbook
        
        Variable Series:
            Meteorological (MList): Ah_EC, Cc_7500_Av, ps, Ta_EC, Ws_CSAT, Wd_CSAT
            Radiation (RList): Fld, Flu, Fn, Fsd, Fsu
            Soil water content (SwsList): all variables containing Sws in variable name
            Soil (SList): Fg, Ts, SwsList
            Turbulent fluxes (FList): Fc_wpl, Fe_wpl, Fh, ustar
            Output (OList): MList, RList, SList, FList
        
        Parameters loaded from control file:
            zmd: z-d
            z0: roughness height
        
        Functions performed:
            qcts.AddMetVars
            qcts.ComputeDailySums
            qcts.InterpolateOverMissing (OList for gaps shorter than 3 observations, OList gaps shorter than 7 observations)
            qcts.GapFillFromAlternate (MList, RList)
            qcts.GapFillFromClimatology (Ah_EC, Fn, Fg, ps, Ta_EC, Ws_CSAT, OList)
            qcts.GapFillFromRatios (Fe, Fh, Fc)
            qcts.ReplaceOnDiff (Ws_CSAT, ustar)
            qcts.UstarFromFh
            qcts.ReplaceWhereMissing (Ustar)
            qcck.do_qcchecks
        """
    if AttrLevel == 'False':
        ds3.globalattributes['Functions'] = ''
        AttrLevel = InLevel
    # check to ensure L4 functions are defined in controlfile
    if qcutils.cfkeycheck(cf,Base='Functions'):
        x=0
        y=0
        z=0
    else:
        log.error('FunctionList not found in control file')
        ds3x = copy.deepcopy(ds3)
        ds3x.globalattributes['nc_level'] = 'L3'
        ds3x.globalattributes['L4Functions'] = 'No L4-L6 functions applied'
        return ds3x
    
    # handle meta-data and import L4-L6 from external process
    if InLevel == 'L3':
        ds3x = copy.deepcopy(ds3)
    else:
        infilename = qcio.get_infilename_from_cf(cf,InLevel)
        ds3x = qcio.nc_read_series(infilename)
        
        for ThisOne in ds3.globalattributes.keys():
            if ThisOne not in ds3x.globalattributes.keys():
                ds3x.globalattributes[ThisOne] = ds3.globalattributes[ThisOne]
        
        for ThisOne in ds3.series.keys():
            if ThisOne in ds3x.series.keys():
                for attr in ds3.series[ThisOne]['Attr'].keys():
                    if attr not in ['ancillary_variables','long_name','standard_name','units']:
                        ds3x.series[ThisOne]['Attr'][attr] = ds3.series[ThisOne]['Attr'][attr]
        
        ds3x.globalattributes['nc_level'] = AttrLevel
        ds3x.globalattributes['EPDversion'] = sys.version
        ds3x.globalattributes['QC_version_history'] = cfg.__doc__
        # put the control file name into the global attributes
        ds3x.globalattributes['controlfile_name'] = cf['controlfile_name']
        if OutLevel == 'L6':
            ds3x.globalattributes['xlL6_datemode'] = ds3x.globalattributes['xl_datemode']
            ds3x.globalattributes['xl_datemode'] = ds3.globalattributes['xl_datemode']
            ds3x.globalattributes['xlL6_filename'] = ds3x.globalattributes['xl_filename']
            ds3x.globalattributes['xl_filename'] = ds3.globalattributes['xl_filename']
            ds3x.globalattributes['xlL6_moddatetime'] = ds3x.globalattributes['xl_moddatetime']
            ds3x.globalattributes['xl_moddatetime'] = ds3.globalattributes['xl_moddatetime']
        elif OutLevel == 'L5':
            ds3x.globalattributes['xlL5_datemode'] = ds3x.globalattributes['xl_datemode']
            ds3x.globalattributes['xl_datemode'] = ds3.globalattributes['xl_datemode']
            ds3x.globalattributes['xlL5_filename'] = ds3x.globalattributes['xl_filename']
            ds3x.globalattributes['xl_filename'] = ds3.globalattributes['xl_filename']
            ds3x.globalattributes['xlL5_moddatetime'] = ds3x.globalattributes['xl_moddatetime']
            ds3x.globalattributes['xl_moddatetime'] = ds3.globalattributes['xl_moddatetime']
        elif OutLevel == 'L4':
            ds3x.globalattributes['xlL4_datemode'] = ds3x.globalattributes['xl_datemode']
            ds3x.globalattributes['xl_datemode'] = ds3.globalattributes['xl_datemode']
            ds3x.globalattributes['xlL4_filename'] = ds3x.globalattributes['xl_filename']
            ds3x.globalattributes['xl_filename'] = ds3.globalattributes['xl_filename']
            ds3x.globalattributes['xlL4_moddatetime'] = ds3x.globalattributes['xl_moddatetime']
            ds3x.globalattributes['xl_moddatetime'] = ds3.globalattributes['xl_moddatetime']
        
        qcutils.prepOzFluxVars(cf,ds3x)
        # convert Fc [mgCO2 m-2 s-1] to Fc_co2 [mgCO2 m-2 s-1], Fc_c [mgC m-2 s-1], NEE [umol m-2 s-1] and NEP = - NEE
        if qcutils.cfkeycheck(cf,Base='Functions',ThisOne='convertFc') and cf['Functions']['convertFc'] == 'True':
            try:
                ds3x.globalattributes['L4Functions'] = ds3x.globalattributes['L4Functions']+', convertFc'
            except:
                ds3x.globalattributes['L4Functions'] = 'convertFc'
            if 'Fc_co2' in ds3x.series.keys():
                qcts.ConvertFc(cf,ds3x,Fco2_in='Fc_co2')
            else:
                qcts.ConvertFc(cf,ds3x)
    
    ds4x = copy.deepcopy(ds3x)
    for ThisOne in ['NEE','NEP','Fc','Fc_co2','Fc_c','Fe','Fh']:
        if ThisOne in ds4x.series.keys() and ThisOne in ds3.series.keys():
            ds4x.series[ThisOne] = ds3.series[ThisOne].copy()
    for ThisOne in ['GPP','CE','ER_night','ER_dark','CE_day','CE_NEEmax','ER_bio','PD','ER_n','ER_LRF']:
        if ThisOne in ds4x.series.keys():
            ds4x.series[ThisOne]['Data'] = numpy.ones(len(ds4x.series[ThisOne]['Data']),dtype=numpy.float64) * numpy.float64(c.missing_value)
            ds4x.series[ThisOne]['Flag'] = numpy.ones(len(ds4x.series[ThisOne]['Data']), dtype=numpy.int32)
    if InLevel == 'L4' or AttrLevel == 'L3':
        ds4,x = l4qc(cf,ds4x,InLevel,x)
        qcutils.get_coverage_individual(ds4)
        qcutils.get_coverage_groups(ds4)
        if qcutils.cfkeycheck(cf,Base='Functions',ThisOne='FlagStats') and cf['Functions']['FlagStats'] == 'True':
            qcio.get_seriesstats(cf,ds4)
    if OutLevel == 'L5' or OutLevel == 'L6':
        try:
            ds4y = copy.deepcopy(ds4)
        except:
            ds4y = copy.deepcopy(ds4x)
        for ThisOne in ['NEE','NEP','Fc','Fc_c','Fc_co2','Fc_c','Fe','Fh']:
            var, var_flag, var_attr = qcutils.GetSeriesasMA(ds3x,ThisOne)
            qcutils.CreateSeries(ds4y,ThisOne,var,Flag=var_flag,Attr=var_attr)
            ds4y.series[ThisOne]['Attr']['long_name'] = var_attr['long_name']
        ds5,y = l5qc(cf,ds4y,y)
        qcutils.get_coverage_individual(ds5)
        qcutils.get_coverage_groups(ds5)
        if qcutils.cfkeycheck(cf,Base='Functions',ThisOne='FlagStats') and cf['Functions']['FlagStats'] == 'True':
            qcio.get_seriesstats(cf,ds5)
    if OutLevel == 'L6':
        ds5z = copy.deepcopy(ds5)
        for ThisOne in ['GPP','CE','ER_night','ER_dark','CE_day','CE_NEEmax','ER_bio','PD','ER_n','ER_LRF']:
            if ThisOne in ds3x.series.keys():
                ds5z.series[ThisOne] = ds3x.series[ThisOne].copy()
        ds6,z = l6qc(cf,ds5z,z)
        qcutils.get_coverage_individual(ds6)
        qcutils.get_coverage_groups(ds6)
        if qcutils.cfkeycheck(cf,Base='Functions',ThisOne='FlagStats') and cf['Functions']['FlagStats'] == 'True':
            qcio.get_seriesstats(cf,ds6)
    
    # calculate daily statistics
    if qcutils.cfkeycheck(cf,Base='Functions',ThisOne='Sums'):
        if cf['Functions']['Sums'] == 'L6':
            ds6.globalattributes['Functions'] = ds6.globalattributes['Functions']+', Sums'
            try:
                ds6.globalattributes['L6Functions'] = ds6.globalattributes['L6Functions']+', Sums'
            except:
                ds6.globalattributes['L6Functions'] = 'Sums'
            
            qcts.do_sums(cf,ds6)
        
        elif cf['Functions']['Sums'] == 'L5':
            ds5.globalattributes['Functions'] = ds5.globalattributes['Functions']+', Sums'
            try:
                ds5.globalattributes['L5Functions'] = ds5.globalattributes['L5Functions']+', Sums'
            except:
                ds5.globalattributes['L5Functions'] = 'Sums'
            
            qcts.do_sums(cf,ds5)
        
        elif cf['Functions']['Sums'] == 'L4':
            ds4.globalattributes['Functions'] = ds4.globalattributes['Functions']+', Sums'
            try:
                ds4.globalattributes['L4Functions'] = ds4.globalattributes['L5Functions']+', Sums'
            except:
                ds4.globalattributes['L4Functions'] = 'Sums'
            
            qcts.do_sums(cf,ds4)
        
    
    # compute climatology
    if qcutils.cfkeycheck(cf,Base='Functions',ThisOne='climatology'):
        if cf['Functions']['climatology'] == 'L6':
            ds6.globalattributes['Functions'] = ds6.globalattributes['Functions']+', climatology'
            try:
                ds6.globalattributes['L6Functions'] = ds6.globalattributes['L6Functions']+', climatology'
            except:
                ds6.globalattributes['L6Functions'] = 'climatology'
            
            qcts.do_climatology(cf,ds6)
        
        elif cf['Functions']['climatology'] == 'L5':
            ds5.globalattributes['Functions'] = ds5.globalattributes['Functions']+', climatology'
            try:
                ds5.globalattributes['L5Functions'] = ds5.globalattributes['L5Functions']+', climatology'
            except:
                ds5.globalattributes['L5Functions'] = 'climatology'
            
            qcts.do_climatology(cf,ds5)
        
        elif cf['Functions']['climatology'] == 'L4':
            ds4.globalattributes['Functions'] = ds4.globalattributes['Functions']+', climatology'
            try:
                ds4.globalattributes['L4Functions'] = ds4.globalattributes['L4Functions']+', climatology'
            except:
                ds4.globalattributes['L4Functions'] = 'climatology'
            
            qcts.do_climatology(cf,ds4)
        
    if OutLevel == 'L4' and (InLevel == 'L3' or InLevel == 'L4'):
        if x == 0:
            ds4.globalattributes['Functions'] = ds4.globalattributes['Functions'] + ', No further L4 gapfilling'
            try:
                ds4.globalattributes['L4Functions'] = ds4.globalattributes['L4Functions'] + ', No further L4 gapfilling'
            except:
                ds4.globalattributes['L4Functions'] = 'No further L4 gapfilling'
            
            log.warn('  L4:  no record of gapfilling functions')
        return ds4
    elif OutLevel == 'L5':
        if x == 0:
            if InLevel == 'L3' or InLevel == 'L4':
                ds4.globalattributes['Functions'] = ds4.globalattributes['Functions'] + ', No further L4 gapfilling'
                try:
                    ds4.globalattributes['L4Functions'] = ds4.globalattributes['L4Functions'] + ', No further L4 gapfilling'
                except:
                    ds4.globalattributes['L4Functions'] = 'No further L4 gapfilling'
                log.warn('  L4:  no record of gapfilling functions')
            ds5.globalattributes['Functions'] = ds5.globalattributes['Functions'] + ', No further L4 gapfilling'
            try:
                ds5.globalattributes['L4Functions'] = ds5.globalattributes['L4Functions'] + ', No further L4 gapfilling'
            except:
                ds5.globalattributes['L4Functions'] = 'No further L4 gapfilling'
        if y == 0:
            ds5.globalattributes['Functions'] = ds5.globalattributes['Functions'] + ', No further L5 gapfilling'
            try:
                ds5.globalattributes['L5Functions'] = ds5.globalattributes['L5Functions'] + ', No further L5 gapfilling'
            except:
                ds5.globalattributes['L5Functions'] = 'No further L5 gapfilling'
            
            log.warn('  L5:  no record of gapfilling functions')
        return ds4,ds5
    elif OutLevel == 'L6':
        if x == 0:
            if InLevel == 'L3' or InLevel == 'L4':
                ds4.globalattributes['Functions'] = ds4.globalattributes['Functions'] + ', No further L4 gapfilling'
                try:
                    ds4.globalattributes['L4Functions'] = ds4.globalattributes['L4Functions'] + ', No further L4 gapfilling'
                except:
                    ds4.globalattributes['L4Functions'] = 'No further L4 gapfilling'
                log.warn('  L4:  no record of gapfilling functions')
            if InLevel == 'L3' or InLevel == 'L4' or InLevel == 'L5':
                ds5.globalattributes['Functions'] = ds5.globalattributes['Functions'] + ', No further L4 gapfilling'
                try:
                    ds5.globalattributes['L4Functions'] = ds5.globalattributes['L4Functions'] + ', No further L4 gapfilling'
                except:
                    ds5.globalattributes['L4Functions'] = 'No further L4 gapfilling'
                log.warn('  L4:  no record of gapfilling functions')
            ds6.globalattributes['Functions'] = ds6.globalattributes['Functions'] + ', No further L4 gapfilling'
            try:
                ds6.globalattributes['L4Functions'] = ds6.globalattributes['L4Functions'] + ', No further L4 gapfilling'
            except:
                ds6.globalattributes['L4Functions'] = 'No further L4 gapfilling'
        
        if y == 0:
            if InLevel == 'L3' or InLevel == 'L4' or InLevel == 'L5':
                ds5.globalattributes['Functions'] = ds5.globalattributes['Functions'] + ', No further L5 gapfilling'
                try:
                    ds5.globalattributes['L5Functions'] = ds5.globalattributes['L5Functions'] + ', No further L5 gapfilling'
                except:
                    ds5.globalattributes['L5Functions'] = 'No further L5 gapfilling'
                log.warn('  L5:  no record of gapfilling functions')
            ds6.globalattributes['Functions'] = ds6.globalattributes['Functions'] + ', No further L5 gapfilling'
            try:
                ds6.globalattributes['L5Functions'] = ds6.globalattributes['L5Functions'] + ', No further L5 gapfilling'
            except:
                ds6.globalattributes['L5Functions'] = 'No further L5 gapfilling'
        if z == 0:
            ds6.globalattributes['Functions'] = ds6.globalattributes['Functions'] + ', No further L6 partitioning'
            try:
                ds6.globalattributes['L6Functions'] = ds5.globalattributes['L6Functions'] + ', No further L6 partitioning'
            except:
                ds6.globalattributes['L6Functions'] = 'No further L6 partitioning'
            log.warn('  L6:  no record of gapfilling functions')
        return ds4,ds5,ds6
    

def l5qc(cf,ds4,y):
    ds5 = copy.deepcopy(ds4)
    ds5.globalattributes['nc_level'] = 'L5'
    if (qcutils.cfkeycheck(cf,Base='Functions',ThisOne='L5_offline') and cf['Functions']['L5_offline'] == 'True') and qcutils.cfkeycheck(cf,Base='Functions',ThisOne='L5_keys'):
        try:
            ds5.globalattributes['L5Functions'] = ds5.globalattributes['L5Functions']+', '+cf['Functions']['L5_keys']
        except:
            ds5.globalattributes['L5Functions'] = cf['Functions']['L5_keys']
        
        y=y+1
    
    # calculate u* from Fh and corrected wind speed
    if qcutils.cfkeycheck(cf,Base='Functions',ThisOne='UstarFromFh') and cf['Functions']['UstarFromFh'] == 'True':
        try:
            ds5.globalattributes['L5Functions'] = ds4.globalattributes['L5Functions']+', UstarFromFh'
        except:
            ds4.globalattributes['L5Functions'] = 'UstarFromFh'
        
        qcts.UstarFromFh(cf,ds5)
        y=y+1
    
    # calcluate ET at observation interval
    if qcutils.cfkeycheck(cf,Base='Functions',ThisOne='CalculateET') and cf['Functions']['CalculateET'] == 'True':
        try:
            ds5.globalattributes['L5Functions'] = ds5.globalattributes['L5Functions']+', CalculateET'
        except:
            ds5.globalattributes['L5Functions'] = 'CalculateET'
        
        log.info(' Calculating ET')
        qcts.CalculateET(cf,ds5,'L5')
    
    # calculate rst, rc and Gst, Gc from Penman-Monteith inversion
    if qcutils.cfkeycheck(cf,Base='Functions',ThisOne='PenmanMonteith') and cf['Functions']['PenmanMonteith'] == 'True':
        try:
            ds5.globalattributes['L5Functions'] = ds5.globalattributes['L5Functions']+', PenmanMonteith'
        except:
            ds5.globalattributes['L5Functions'] = 'PenmanMonteith'
        
        qcts.do_PenmanMonteith(cf,ds5)
    
    # re-calculate the available energy from L5 (gapfilled) fluxes
        try:
            ds5.globalattributes['L5Functions'] = ds5.globalattributes['L5Functions']+', CalculateAvailableEnergy'
        except:
            ds.globalattributes['L5Functions'] = 'CalculateAvailableEnergy'
        
        qcts.CalculateAvailableEnergy(ds5)
    
    # re-apply the quality control checks (range, diurnal and rules)
    if y > 0:
        log.info(' Doing QC checks on L5 data')
        qcck.do_qcchecks(cf,ds5)
        try:
            ds5.globalattributes['L5Functions'] = ds5.globalattributes['L5Functions']+', do_qccheck(RangeCheck, diurnalcheck, excludedates, excludehours)'
        except:
            ds5.globalattributes['L5Functions'] = 'do_qccheck(RangeCheck, diurnalcheck, excludedates, excludehours)'
    
    try:
        ds5.globalattributes['Functions'] = ds5.globalattributes['Functions'] + ', ' + ds5.globalattributes['L5Functions']
    except:
        ds5.globalattributes['Functions'] = ds5.globalattributes['Functions']
    
    return ds5,y

def l6qc(cf,ds5,z):
    ds6 = copy.deepcopy(ds5)
    ds6.globalattributes['nc_level'] = 'L6'
    if (qcutils.cfkeycheck(cf,Base='Functions',ThisOne='L6_offline') and cf['Functions']['L6_offline'] == 'True') and qcutils.cfkeycheck(cf,Base='Functions',ThisOne='L6_keys'):
        try:
            ds6.globalattributes['L6Functions'] = ds6.globalattributes['L6Functions']+', '+cf['Functions']['L6_keys']
        except:
            ds6.globalattributes['L6Functions'] = cf['Functions']['L6_keys']
        z=z+1
    
    # run MOST (Buckingham Pi) 2d footprint model (Kljun et al. 2004)
    if qcutils.cfkeycheck(cf,Base='Functions',ThisOne='footprint') and cf['Functions']['footprint'] == 'True':
        try:
            ds6.globalattributes['L6Functions'] = ds6.globalattributes['L6Functions']+', footprint'
        except:
            ds6.globalattributes['L6Functions'] = 'footprint'
        
        qcts.do_footprint_2d(cf,ds6,datalevel='L6')
    
    try:
        ds6.globalattributes['Functions'] = ds6.globalattributes['Functions'] + ', ' + ds6.globalattributes['L6Functions']
    except:
        ds6.globalattributes['Functions'] = ''
    
    # re-apply the quality control checks (range, diurnal and rules)
    if z > 0:
        log.info(' Doing QC checks on L6 data')
        qcck.do_qcchecks(cf,ds6)
        try:
            ds6.globalattributes['L6Functions'] = ds6.globalattributes['L6Functions']+', do_qccheck(RangeCheck, diurnalcheck, excludedates, excludehours)'
        except:
            ds6.globalattributes['L6Functions'] = 'do_qccheck(RangeCheck, diurnalcheck, excludedates, excludehours)'
    
    return ds6,z
