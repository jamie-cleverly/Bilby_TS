#    pls.py
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
    Partition v0.3 24 July 2013;

    Version History:
    <<v0.0  20 July 2011>>
    <<v0.1: 21 May 2012, stand-alone gui/app set to OzFluxQCv1.8.2.b1 standards>>
    <<v0.2: 25 June 2012, stand-alone gui/app set to OzFluxQCv2.0 standards>>
    <<v0.3: 24 July 2013, stand-alone gui/app set to OzFluxQCPlusv2.5 standards>>
"""

import sys
import logging
import ast
import constants as c
import copy
import numpy
import pio
import pts
import putils
import time
import xlrd
import logging

log = logging.getLogger('partition.ls')

def l3partition(cf,ds2):
    '''Processing OzFlux_Level3 data for partitioning'''
    # make a copy of the OzFlux_Level3 data
    ds3 = copy.deepcopy(ds2)
    ds3.globalattributes['Level'] = 'L5'
    # compute daily statistics
    if putils.cfkeycheck(cf,Base='Params',ThisOne='firstMonth'):
        M1st = int(cf['Params']['firstMonth'])
    else:
        M1st = 1
    if putils.cfkeycheck(cf,Base='Params',ThisOne='secondMonth'):
        M2nd = int(cf['Params']['secondMonth'])
    else:
        M2nd = 12
    # prep nighttime ER observations
    pts.ER_nightL3(cf,ds3,M1st,M2nd,'Fc')
    pts.LightResponseCurves(ds3)
    log.info('L3 Partitioning: All done')
    return ds3
    

def l6partition(cf,ds5):
    '''Processing OzFlux_Level4 data to partition daytime ER, CE and GPP'''
    # make a copy of the OzFlux_Level4 data
    ds6 = copy.deepcopy(ds5)
    ds6.globalattributes['Level'] = 'L6'
    # prep nighttime ER observations
    pts.ER_nightL6(cf,ds6,'ER_night_gapfilled')
    if 'AliceSpringsMulga' in ds6.globalattributes['site_name']:
        pts.DayCEGPP_ASM(cf,ds6,'Fc')
    elif 'TiTreeEast' in ds6.globalattributes['site_name']:
        pts.DayERdark_TTE(cf,ds6,'Fc')
        pts.DayPD_CE_GPP_TTE(cf,ds6,'NEE')
    else:
        log.error(' Site designation undefined in ds6.globalattributes:xl_filename')
    
    log.info('L6 Partitioning: All done')
    return ds6
    
