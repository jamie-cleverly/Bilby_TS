#    metls.py
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
    Partition v0.2 25 June 2012;

    Version History:
    <<v0.0  20 July 2011>>
    <<v0.1: 21 May 2012, stand-alone gui/app set to OzFluxQCv1.8.2.b1 standards>>
    <<v0.2: 25 June 2012, stand-alone gui/app set to OzFluxQCv2.0 standards>>
"""

import sys
import logging
import ast
import constants as c
import copy
import numpy
import metts
import metutils
import time
import xlrd
import logging

log = logging.getLogger('envelope.ls')

def l3envelope(cf,ds2,Params):
    '''Processing OzFlux_Level3 data for enveloping'''
    # make a copy of the OzFlux_Level3 data
    ds3 = copy.deepcopy(ds2)
    # prep nighttime Re observations
    if Params == 'Es':
        metts.EsBinFilters(cf,ds3)
    elif Params == 'Ts':
        metts.TsBinFilters(cf,ds3)
    elif Params == 'D':
        metts.DBinFilters(cf,ds3)
    elif Params == 'LRF':
        metts.LRFBinFilters(cf,ds3)
    elif Params == 'LRFD':
        metts.LRFDBinFilters(cf,ds3)
    log.info('Envelopes: All done')
    return ds3
