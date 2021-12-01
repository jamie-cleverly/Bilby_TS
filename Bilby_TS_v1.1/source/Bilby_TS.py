#    Bilby_TS.py  Starts primary GUI for Bilby TS Time Series analyst
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

import sys, os

sys.path.append(os.path.abspath('scripts'))

import ast
import copy
import datetime
import logging
import numpy
import time
import Tkinter
import cfg
import qcio
import qcls
import qcplot
import qcts
import qcutils
import pdb

class qcgui(Tkinter.Frame):
    """
        QC Data Main GUI
        Used to access read, save, and data processing (qcls) prodecures
        
        Columns: Data levels:
            1:  L1 Raw Data (read excel into NetCDF)
            2:  L2 QA/QC (general QA/QC algorithms, site independent)
            3:  L3 Corrections (Flux data corrections, site dependent based on ancillary measurements available and technical issues)
            4:  L4 Gap Filling (Used for fill met data gaps and ingesting SOLO-ANN Gap Filled fluxes from external processes)
        
        Rows:  function access
            1:  Ingest excel dataset into NetCDF files
            2:  Process data from previous level and generate NetCDF file(s) at current level
            3-6:  Show Timestamp range of dataset and accept date range for graphical plots
            7:  Export excel dataset from NetCDF file
        """
    def __init__(self, master=None):
        Tkinter.Frame.__init__(self, master)
        self.grid()
        self.createWidgets()

    def createWidgets(self):
        self.process1Label = Tkinter.Label(self,text='L1: raw data')
        self.process1Label.grid(row=0,column=1,columnspan=2)
        self.process1Label = Tkinter.Label(self,text='L2: QA/QC')
        self.process1Label.grid(row=0,column=3,columnspan=2)
        self.process2Label = Tkinter.Label(self,text='L3: Corrections')
        self.process2Label.grid(row=0,column=5,columnspan=2)
        self.process3Label = Tkinter.Label(self,text='L4-L6: Gapfilled & C partitioning')
        self.process3Label.grid(row=0,column=7,columnspan=2)
        
        self.fileloadLabel = Tkinter.Label(self,text='Xcel ->')
        self.fileloadLabel.grid(row=1,column=0,columnspan=1)
        self.doxl2nc1Button = Tkinter.Button (self, text="Load L1 Data", command=self.do_xl2ncL1 )
        self.doxl2nc1Button.grid(row=1,column=1,columnspan=2)
        self.doxl2nc2Button = Tkinter.Button (self, text="Load L3 Data", command=self.do_xl2ncL3 )
        self.doxl2nc2Button.grid(row=1,column=5,columnspan=2)
        self.doxl2nc2Button = Tkinter.Button (self, text="Load L4-L6 Data", command=self.do_xl2ncL4 )
        self.doxl2nc2Button.grid(row=1,column=7,columnspan=2)
        
        self.initiateLabel = Tkinter.Label(self,text='Process Data')
        self.initiateLabel.grid(row=2,column=0,columnspan=1)
        self.doL1Button = Tkinter.Button (self, text="L2 Processing", command=self.do_l2qc )
        self.doL1Button.grid(row=2,column=3,columnspan=2)
        self.doL2Button = Tkinter.Button (self, text="L3 Processing", command=self.do_l3qc )
        self.doL2Button.grid(row=2,column=5,columnspan=2)
        self.doL3Button = Tkinter.Button (self, text="L4-L6 Processing", command=self.do_l4to6qc )
        self.doL3Button.grid(row=2,column=7,columnspan=2)
        
        self.initiateLabel = Tkinter.Label(self,text='Re/Load NetCDF Data')
        self.initiateLabel.grid(row=3,column=0,columnspan=1)
        self.doL1Button = Tkinter.Button (self, text="Re/Load L2 Data", command=self.load_l2nc )
        self.doL1Button.grid(row=3,column=3,columnspan=2)
        self.doL2Button = Tkinter.Button (self, text="Re/Load L3 Data", command=self.load_l3nc )
        self.doL2Button.grid(row=3,column=5,columnspan=2)
        self.doL3Button = Tkinter.Button (self, text="Re/Load L4-L6 Data", command=self.load_l4nc )
        self.doL3Button.grid(row=3,column=7,columnspan=2)
        
        self.filestartLabel = Tkinter.Label(self,text='File start date')
        self.filestartLabel.grid(row=4,column=3,columnspan=2)
        self.fileendLabel = Tkinter.Label(self,text='File end date')
        self.fileendLabel.grid(row=4,column=5,columnspan=2)
        
        self.clearcfButton = Tkinter.Button (self, text="Join NetCDF files", command=self.do_ncconcat )
        self.clearcfButton.grid(row=5,column=0,columnspan=1)
        self.filestartValue = Tkinter.Label(self,text='No file loaded ...')
        self.filestartValue.grid(row=5,column=3,columnspan=2)
        self.fileendValue = Tkinter.Label(self,text='No file loaded ...')
        self.fileendValue.grid(row=5,column=5,columnspan=2)
        
        self.closeplotwindowsButton = Tkinter.Button (self, text="Close plot windows", command=self.do_closeplotwindows )
        self.closeplotwindowsButton.grid(row=6,column=0,columnspan=1)
        self.plotstartLabel = Tkinter.Label(self, text='Start date (YYYY-MM-DD)')
        self.plotstartLabel.grid(row=6,column=3,columnspan=2)
        self.plotstartEntry = Tkinter.Entry(self)
        self.plotstartEntry.grid(row=6,column=5,columnspan=2)
        
        self.plotendLabel = Tkinter.Label(self, text='  End date (YYYY-MM-DD)')
        self.plotendLabel.grid(row=7,column=3,columnspan=2)
        self.plotendEntry = Tkinter.Entry(self)
        self.plotendEntry.grid(row=7,column=5,columnspan=2)
        
        self.plotL1L2Button = Tkinter.Button (self, text="Plot L1 & L2 Data", command=self.do_plotL1L2 )
        self.plotL1L2Button.grid(row=8,column=3,columnspan=2)
        self.plotL3L3Button = Tkinter.Button (self, text="Plot L3 Data", command=self.do_plotL3L3 )
        self.plotL3L3Button.grid(row=8,column=5,columnspan=2)
        
        self.filesave2Label = Tkinter.Label(self,text='-> Xcel')
        self.filesave2Label.grid(row=9,column=0,columnspan=1)
        self.savexL2Button = Tkinter.Button (self, text="Save L2 Data", command=self.do_savexL2 )
        self.savexL2Button.grid(row=9,column=3,columnspan=2)
        self.savexL3Button = Tkinter.Button (self, text="Save L3 Data", command=self.do_savexL3 )
        self.savexL3Button.grid(row=9,column=5,columnspan=2)
        self.savexL4Button = Tkinter.Button (self, text="Save L4-L6 Data", command=self.do_savexL4 )
        self.savexL4Button.grid(row=9,column=7,columnspan=2)
        
        self.quitButton = Tkinter.Button (self, text="Quit", command=self.do_quit )
        self.quitButton.grid(row=10,column=0,columnspan=1)
        self.progress = Tkinter.Label(self, text='Waiting for input ...')
        self.progress.grid(row=10,column=1,columnspan=7)
    
    def do_closeplotwindows(self):
        """
            Close plot windows
            """
        import matplotlib
        self.do_progress(text='Closing plot windows ...')             # tell the user what we're doing
        log.info(' Closing plot windows ...')
        fig_numbers = [n.num for n in matplotlib._pylab_helpers.Gcf.get_all_fig_managers()]
        log.info('  Closing plot windows: '+str(fig_numbers))
        for n in fig_numbers:
            matplotlib.pyplot.close(n)
        self.do_progress(text='Waiting for input ...')             # tell the user what we're doing
        log.info(' Waiting for input ...')

    def do_l2qc(self):
        """
            Call qcls.l2qc function
            Performs L2 QA/QC processing on raw data
            Outputs L2 netCDF file to ncData folder
            
            ControlFiles:
                L2_year.txt
                or
                L2.txt
            
            ControlFile contents (see ControlFile/Templates/L2.txt for example):
                [General]:
                    Enter list of functions to be performed
                [Files]:
                    L1 input file name and path
                    L2 output file name and path
                [Variables]:
                    Variable names and parameters for:
                        Range check to set upper and lower rejection limits
                        Diurnal check to reject observations by time of day that
                            are outside specified standard deviation limits
                        Timestamps for excluded dates
                        Timestamps for excluded hours
                [Plots]:
                    Variable lists for plot generation
            """
        self.do_progress(text='Load L2 Control File ...')
        log.info('Beginning L2 process')
        self.cf = qcio.load_controlfile(path='controlfiles')
        if len(self.cf)==0: self.do_progress(text='Waiting for input ...'); return
        if qcutils.cfkeycheck(self.cf,Base='General',ThisOne='InputLevel'):
            InLevel = cf['General']['InputLevel']
        else:
            InLevel = 'L1'
        infilename = qcio.get_infilename_from_cf(self.cf,InLevel)
        if len(infilename)==0: self.do_progress(text='An error occurred, check the console ...'); return
        self.do_progress(text='Doing L2 QC ...')
        self.ds1 = qcio.nc_read_series(infilename)
        if len(self.ds1.series.keys())==0: self.do_progress(text='An error occurred, check the console ...'); del self.ds1; return
        self.update_startenddate(str(self.ds1.series['DateTime']['Data'][0]),
                                 str(self.ds1.series['DateTime']['Data'][-1]))
        self.ds2 = qcls.l2qc(self.cf,self.ds1)
        log.info(' Finished L2 QC process')
        self.do_progress(text='Finished L2 QC process')
        self.do_progress(text='Saving L2 QC ...')                     # put up the progress message
        outfilename = qcio.get_outfilename_from_cf(self.cf,'L2')
        if len(outfilename)==0: self.do_progress(text='An error occurred, check the console ...'); return
        ncFile = qcio.nc_open_write(outfilename)
        log.info(' Writing netCDF file')
        qcio.nc_write_series(ncFile,self.ds2)                   # save the L2 data
        self.do_progress(text='Finished saving L2 QC data')              # tell the user we are done
        log.info(' Finished saving L2 QC data')
        print '\a'

    def do_l3qc(self):
        """
            Call qcls.l3qc_sitename function
            Performs L3 Corrections and QA/QC processing on L2 data
            Outputs L3 netCDF file to ncData folder
            Outputs L3 netCDF file to OzFlux folder
            
            Available corrections:
            * corrections requiring ancillary measurements or samples
              marked with an asterisk
                Linear correction
                    fixed slope
                    linearly shifting slope
                Conversion of virtual temperature to actual temperature
                2D Coordinate rotation
                Massman correction for frequency attenuation*
                Webb, Pearman and Leuning correction for flux effects on density
                    measurements
                Conversion of virtual heat flux to actual heat flux
                Correction of soil moisture content to empirical calibration
                    curve*
                Addition of soil heat storage to ground ground heat flux*
            
            ControlFiles:
                L3_year.txt
                or
                L3a.txt
            
            ControlFile contents (see ControlFile/Templates/L3.txt for example):
                [General]:
                    Python control parameters
                [Files]:
                    L2 input file name and path
                    L3 output file name and ncData folder path
                    L3 OzFlux output file name and OzFlux folder path
                [Massman] (where available):
                    Constants used in frequency attenuation correction
                        zmd: instrument height (z) less zero-plane displacement
                            height (d), m
                        z0: aerodynamic roughness length, m
                        angle: angle from CSAT mounting point between CSAT and
                            IRGA mid-path, degrees
                        CSATarm: distance from CSAT mounting point to CSAT
                            mid-path, m
                        IRGAarm: distance from CSAT mounting point to IRGA
                            mid-path, m
                [Soil]:
                    Constants used in correcting Fg for storage and in empirical
                    corrections of soil water content 
                        FgDepth: Heat flux plate depth, m
                        BulkDensity: Soil bulk density, kg/m3
                        OrganicContent: Soil organic content, fraction
                        SwsDefault
                        Constants for empirical corrections using log(sensor)
                            and exp(sensor) functions (SWC_a0, SWC_a1, SWC_b0,
                            SWC_b1, SWC_t, TDR_a0, TDR_a1, TDR_b0, TDR_b1,
                            TDR_t)
                        Variable and attributes lists (empSWCin, empSWCout,
                            empTDRin, empTDRout, linTDRin, SWCattr, TDRattr)
                [Output]:
                    Variable subset list for OzFlux output file
                [Variables]:
                    Variable names and parameters for:
                        Range check to set upper and lower rejection limits
                        Diurnal check to reject observations by time of day that
                            are outside specified standard deviation limits
                        Timestamps, slope, and offset for Linear correction
                [Plots]:
                    Variable lists for plot generation
            """
        self.do_progress(text='Load L3 Control File ...')
        log.info('Beginning L3 process')
        self.cf = qcio.load_controlfile(path='controlfiles')
        if len(self.cf)==0: self.do_progress(text='Waiting for input ...'); return
        if qcutils.cfkeycheck(self.cf,Base='General',ThisOne='InputLevel'):
            InLevel = self.cf['General']['InputLevel']
        else:
            InLevel = 'L2'
        infilename = qcio.get_infilename_from_cf(self.cf,InLevel)
        if len(infilename)==0: self.do_progress(text='An error occurred, check the console ...'); return
        self.ds2 = qcio.nc_read_series(infilename)
        if len(self.ds2.series.keys())==0: self.do_progress(text='An error occurred, check the console ...'); del self.ds2; return
        self.update_startenddate(str(self.ds2.series['DateTime']['Data'][0]),
                                 str(self.ds2.series['DateTime']['Data'][-1]))
        self.do_progress(text='Doing L3 QC & Corrections ...')
        self.ds3 = qcls.l3qc(self.cf,self.ds2)
        self.do_progress(text='Finished L3')
        try:
            txtstr = ' Finished L3: Standard processing for site: '
            txtstr = txtstr+self.ds3.globalattributes['site_name'].replace(' ','')
            log.info(txtstr)
            
            self.do_progress(text='Saving L3 QC & Corrected NetCDF data ...')                     # put up the progress message
            outfilename = qcio.get_outfilename_from_cf(self.cf,'L3')
            if len(outfilename)==0: self.do_progress(text='An error occurred, check the console ...'); return
            ncFile = qcio.nc_open_write(outfilename)
            log.info(' Writing netCDF file')
            qcio.nc_write_series(ncFile,self.ds3)                   # save the L3 data
            if qcutils.cfkeycheck(self.cf,Base='Output',ThisOne='nc'):
                outfilename2 = qcio.get_outfilename_from_cf(self.cf,'L3_Corrected')
                self.ds3x = copy.deepcopy(self.ds3)
                self.ds3y = qcio.convert_L3Corrected(self.ds3x)
                self.ds3y.globalattributes['Level'] = 'L3_Corrected'
                if len(outfilename2)==0: self.do_progress(text='An error occurred, check the console ...'); return
                ncFile = qcio.nc_open_write(outfilename2)
                outputlist = qcio.get_outputlist_from_cf(self.cf,'nc')
                log.info(' Writing netCDF file')
                qcio.nc_write_series(ncFile,self.ds3y,outputlist=outputlist)
            self.do_progress(text='Finished saving L3 QC & Corrected NetCDF data')              # tell the user we are done
            log.info(' Finished saving L3 QC & Corrected NetCDF data')
            print '\a'
        except:
            txtstr = ' Interupted L3: Error prevented completion for site: '
            txtstr = txtstr+self.ds2.globalattributes['site_name'].replace(' ','')
            log.info(txtstr)
            print '\a'
    
    def do_l4to6qc(self):
        """
            Call qcls.l4qc_gapfill function
            Performs L4 gap filling on L3 met data
            or
            Ingests L4 gap filled fluxes performed in external SOLO-ANN and c
                omputes daily sums
            Outputs L4 netCDF file to ncData folder
            Outputs L4 netCDF file to OzFlux folder
            
            ControlFiles:
                L4_year.txt
                or
                L4b.txt
            
            ControlFile contents (see ControlFile/Templates/L4.txt and
            ControlFile/Templates/L4b.txt for examples):
                [General]:
                    Python control parameters (SOLO)
                    Site characteristics parameters (Gap filling)
                [Files]:
                    L3 input file name and path (Gap filling)
                    L4 input file name and path (SOLO)
                    L4 output file name and ncData folder path (both)
                    L4 OzFlux output file name and OzFlux folder path
                [Variables]:
                    Variable subset list for OzFlux output file (where
                        available)
            """
        log.info('Beginning L4, L5 & L6 processes')
        self.cf = qcio.load_controlfile(path='controlfiles')
        if len(self.cf)==0: self.do_progress(text='Waiting for input ...'); return
        if qcutils.cfkeycheck(self.cf,Base='General',ThisOne='InputLevel'):
            InLevel = self.cf['General']['InputLevel']
        else:
            InLevel = 'L3'
        if qcutils.cfkeycheck(self.cf,Base='General',ThisOne='AttributesLevel'):
            AttrLevel = self.cf['General']['AttributesLevel']
        else:
            AttrLevel = InLevel
        if qcutils.cfkeycheck(self.cf,Base='General',ThisOne='OutputLevel'):
            OutLevel = self.cf['General']['OutputLevel']
        else:
            OutLevel = 'L6'
        try:
            infilename = qcio.get_infilename_from_cf(self.cf,AttrLevel)
        except:
            infilename = qcio.get_infilename_from_cf(self.cf,InLevel)
        if len(infilename)==0: self.do_progress(text='An error occurred, check the console ...'); return
        if not qcutils.file_exists(infilename): self.do_progress(text='An error occurred, check the console ...'); return
        self.ds3 = qcio.nc_read_series(infilename)
        if len(self.ds3.series.keys())==0: self.do_progress(text='An error occurred, check the console ...'); del self.ds3; return
        self.ds3.globalattributes['controlfile_name'] = self.cf['controlfile_name']
        self.update_startenddate(str(self.ds3.series['DateTime']['Data'][0]),
                                 str(self.ds3.series['DateTime']['Data'][-1]))
        sitename = self.ds3.globalattributes['site_name']
        self.do_progress(text='Doing L4 QC: '+sitename+' ...')
        if OutLevel == 'L6':
            self.ds4,self.ds5,self.ds6 = qcls.l4to6qc(self.cf,self.ds3,AttrLevel,InLevel,OutLevel)
        elif OutLevel == 'L5':
            self.ds4,self.ds5 = qcls.l4to6qc(self.cf,self.ds3,AttrLevel,InLevel,OutLevel)
        elif OutLevel == 'L4':
            self.ds4 = qcls.l4to6qc(self.cf,self.ds3,AttrLevel,InLevel,OutLevel)
        else:
            log.error('L4-L6 gapfill:  Invalid output level '+str(OutLevel)+' specified')
        self.do_progress(text='Finished L4-L6: '+sitename)
        log.info(' Finished L4-L6: '+sitename)
        self.do_progress(text='Saving L4-L6 Gap Filled NetCDF data ...')                     # put up the progress message
        outfilename = qcio.get_outfilename_from_cf(self.cf,OutLevel)
        if len(outfilename)==0: self.do_progress(text='An error occurred, check the console ...'); return
        if OutLevel == 'L4':
            ncFile = qcio.nc_open_write(outfilename)
            log.info(' Writing netCDF file')
            qcio.nc_write_series(ncFile,self.ds4)                   # save the L3 data
            self.do_progress(text='Finished saving L4 gap filled NetCDF data')              # tell the user we are done
            log.info(' Finished saving L4 gap filled NetCDF data')
        elif OutLevel == 'L5':
            ncFile = qcio.nc_open_write(outfilename)
            log.info(' Writing netCDF file')
            qcio.nc_write_series(ncFile,self.ds5)                   # save the L3 data
            self.do_progress(text='Finished saving L5 gap filled NetCDF data')              # tell the user we are done
            log.info(' Finished saving L5 gap filled NetCDF data')
        elif OutLevel == 'L6':
            ncFile = qcio.nc_open_write(outfilename)
            log.info(' Writing netCDF file')
            qcio.nc_write_series(ncFile,self.ds6)                   # save the L3 data
            self.do_progress(text='Finished saving L6 gap filled NetCDF data')              # tell the user we are done
            log.info(' Finished saving L6 gap filled NetCDF data')
        if qcutils.cfkeycheck(self.cf,Base='Output',ThisOne='nc4'):
            outfilename2 = qcio.get_outfilename_from_cf(self.cf,'L4_MetFilled')
            self.ds4x = copy.deepcopy(self.ds4)
            self.ds4x.globalattributes['Level'] = 'L4_MetFilled'
            if len(outfilename2)==0: self.do_progress(text='An error occurred, check the console ...'); return
            ncFile = qcio.nc_open_write(outfilename2)
            outputlist = qcio.get_outputlist_from_cf(self.cf,'nc4')
            log.info(' Writing netCDF file')
            qcio.nc_write_series(ncFile,self.ds4x,outputlist=outputlist)
        if qcutils.cfkeycheck(self.cf,Base='Output',ThisOne='nc5'):
            outfilename3 = qcio.get_outfilename_from_cf(self.cf,'L5_FluxFilled')
            self.ds5x = copy.deepcopy(self.ds5)
            self.ds5x.globalattributes['Level'] = 'L5_FluxFilled'
            if len(outfilename3)==0: self.do_progress(text='An error occurred, check the console ...'); return
            ncFile = qcio.nc_open_write(outfilename3)
            outputlist = qcio.get_outputlist_from_cf(self.cf,'nc5')
            log.info(' Writing netCDF file')
            qcio.nc_write_series(ncFile,self.ds5x,outputlist=outputlist)
        if qcutils.cfkeycheck(self.cf,Base='Output',ThisOne='nc6'):
            outfilename4 = qcio.get_outfilename_from_cf(self.cf,'L6_CPartitioned')
            self.ds6x = copy.deepcopy(self.ds6)
            self.ds6x.globalattributes['Level'] = 'L6_CPartitioned'
            if len(outfilename4)==0: self.do_progress(text='An error occurred, check the console ...'); return
            ncFile = qcio.nc_open_write(outfilename4)
            outputlist = qcio.get_outputlist_from_cf(self.cf,'nc6')
            log.info(' Writing netCDF file')
            qcio.nc_write_series(ncFile,self.ds6x,outputlist=outputlist)
        print '\a'

    def do_ncconcat(self):
        """
            Join NetCDF files
            """
        self.do_progress(text='Concatenating NetCDF files ...')             # tell the user what we're doing
        log.info(' Concatenating NetCDF files ...')
        self.cf = qcio.load_controlfile(path='controlfiles')
        if len(self.cf)==0: self.do_progress(text='Waiting for input ...'); return
        qcio.nc_concatenate(self.cf)
        self.do_progress(text='NetCDF files joined')                         # tell the user what we're doing
        log.info(' NetCDF files joined')
        print '\a'

    def do_plotL1L2(self):
        """
            Plot L1 (raw) and L2 (QA/QC) data in blue and red, respectively
            
            Control File for do_l2qc function used.
            If L2 Control File not loaded, requires control file selection.
            """
        if 'ds1' not in dir(self) or 'ds2' not in dir(self):
            self.cf = qcio.load_controlfile(path='controlfiles')
            if len(self.cf)==0: self.do_progress(text='Waiting for input ...'); return
            l1filename = qcio.get_infilename_from_cf(self.cf,'L1')
            if len(l1filename)==0: return
            self.ds1 = qcio.nc_read_series(l1filename)
            if len(self.ds1.series.keys())==0: self.do_progress(text='An error occurred, check the console ...'); del self.ds1; return
            l2filename = qcio.get_outfilename_from_cf(self.cf,'L2')
            self.ds2 = qcio.nc_read_series(l2filename)
            if len(self.ds2.series.keys())==0: self.do_progress(text='An error occurred, check the console ...'); del self.ds2; return
            self.update_startenddate(str(self.ds1.series['DateTime']['Data'][0]),
                                     str(self.ds1.series['DateTime']['Data'][-1]))
        self.do_progress(text='Plotting L1 & L2 QC ...')
        cfname = self.ds2.globalattributes['controlfile_name']
        self.cf = qcio.get_controlfilecontents(cfname)
        for nFig in self.cf['Plots'].keys():
            si = qcutils.GetDateIndex(self.ds1.series['DateTime']['Data'],self.plotstartEntry.get(),
                                      ts=self.ds1.globalattributes['time_step'],default=0,match='exact')
            ei = qcutils.GetDateIndex(self.ds1.series['DateTime']['Data'],self.plotendEntry.get(),
                                      ts=self.ds1.globalattributes['time_step'],default=-1,match='exact')
            plt_cf = self.cf['Plots'][str(nFig)]
            if 'Type' in plt_cf.keys():
                if str(plt_cf['Type']).lower() =='xy':
                    self.do_progress(text='Plotting L1 and L2 XY ...')
                    qcplot.plotxy(self.cf,nFig,plt_cf,self.ds1,self.ds2,si,ei)
                else:
                    self.do_progress(text='Plotting L1 and L2 QC ...')
                    qcplot.plottimeseries(self.cf,nFig,self.ds1,self.ds2,si,ei)
            else:
                self.do_progress(text='Plotting L1 and L2 QC ...')
                qcplot.plottimeseries(self.cf,nFig,self.ds1,self.ds2,si,ei)
        self.do_progress(text='Finished plotting L1 and L2')
        log.info(' Finished plotting L1 and L2, check the GUI')
        print '\a'

    def do_plotL3L3(self):
        """
            Plot L3 (QA/QC and Corrected) data
            
            Control File for do_l3qc function used.
            If L3 Control File not loaded, requires control file selection.
            """
        if 'ds3' not in dir(self):
            self.cf = qcio.load_controlfile(path='controlfiles')
            if len(self.cf)==0: self.do_progress(text='Waiting for input ...'); return
            l3filename = qcio.get_outfilename_from_cf(self.cf,'L3')
            self.ds3 = qcio.nc_read_series(l3filename)
            if len(self.ds3.series.keys())==0: self.do_progress(text='An error occurred, check the console ...'); del self.ds3; return
            self.update_startenddate(str(self.ds3.series['DateTime']['Data'][0]),
                                     str(self.ds3.series['DateTime']['Data'][-1]))
        self.do_progress(text='Plotting L3 QC ...')
        cfname = self.ds3.globalattributes['controlfile_name']
        self.cf = qcio.get_controlfilecontents(cfname)
        for nFig in self.cf['Plots'].keys():
            si = qcutils.GetDateIndex(self.ds3.series['DateTime']['Data'],self.plotstartEntry.get(),
                                      ts=self.ds3.globalattributes['time_step'],default=0,match='exact')
            ei = qcutils.GetDateIndex(self.ds3.series['DateTime']['Data'],self.plotendEntry.get(),
                                      ts=self.ds3.globalattributes['time_step'],default=-1,match='exact')
            plt_cf = self.cf['Plots'][str(nFig)]
            if 'Type' in plt_cf.keys():
                if str(plt_cf['Type']).lower() =='xy':
                    self.do_progress(text='Plotting L3 XY ...')
                    qcplot.plotxy(self.cf,nFig,plt_cf,self.ds3,self.ds3,si,ei)
                else:
                    self.do_progress(text='Plotting L3 QC ...')
                    SeriesList = ast.literal_eval(plt_cf['Variables'])
                    qcplot.plottimeseries(self.cf,nFig,self.ds3,self.ds3,si,ei)
            else:
                self.do_progress(text='Plotting L3 QC ...')
                qcplot.plottimeseries(self.cf,nFig,self.ds3,self.ds3,si,ei)
        self.do_progress(text='Finished plotting L3')
        log.info(' Finished plotting L3, check the GUI')
        print '\a'

    def do_progress(self,text):
        """
            Update progress message in QC Data GUI
            """
        self.progress.destroy()
        self.progress = Tkinter.Label(self, text=text)
        self.progress.grid(row=10,column=1,columnspan=7)
        self.update()

    def do_quit(self):
        """
            Close plot windows and quit QC Data GUI
            """
        import matplotlib
        self.do_progress(text='Closing plot windows ...')             # tell the user what we're doing
        log.info(' Closing plot windows ...')
        matplotlib.pyplot.close('all')
        self.do_progress(text='Quitting ...')                         # tell the user what we're doing
        log.info(' Quitting ...')
        self.quit()

    def do_reloadcf(self):
        """
            Re-load a controlfile (replace current cf or load new cf without processing instructions)
            """
        self.do_progress(text='Loading/Re-loading controlfile ...')             # tell the user what we're doing
        log.info(' Loading/Re-loading controlfile ...')
        self.cf = qcio.load_controlfile(path='controlfiles')
        self.do_progress(text='Controlfile loaded/re-loaded')                         # tell the user what we're doing
        log.info(' Controlfile loaded/re-loaded')

    def do_savexL2(self):
        """
            Call interior_nc2xl.autonc2xl function
            Exports excel data from NetCDF file
            
            Outputs L2 Excel file containing Data and Flag worksheets
            ControlFile:
                nc2xl_year
            ControlFile contents (see ControlFile/Templates/nc2xl.txt for example):
                [General] (optional):
                    InputLevel
                [Files]:
                    Input netCDF file name and path
                    Output excel file name and path
                [Output]:
                    Variables included in output file
            """
        check = 'True'
        self.do_progress(text='Exporting L2 NetCDF -> Xcel ...')                     # put up the progress message
        log.info('Saving L2 to Excel')
        try:
            self.cf
            if len(self.cf)==0:
                self.cf = qcio.load_controlfile(path='controlfiles')
                if len(self.cf)==0: self.do_progress(text='Waiting for input ...'); return
                check = 'False'
        except:
            self.cf = qcio.load_controlfile(path='controlfiles')
            if len(self.cf)==0: self.do_progress(text='Waiting for input ...'); return
            check = 'False'
        
        if (qcutils.cfkeycheck(self.cf,Base='Output',ThisOne='DefaultXl') and self.cf['Output']['DefaultXl'] == 'False') and check == 'True':
            self.cf = qcio.load_controlfile(path='controlfiles')
            if len(self.cf)==0: self.do_progress(text='Waiting for input ...'); return
        
        if qcutils.cfkeycheck(self.cf,Base='General',ThisOne='InputLevel') and qcutils.cfkeycheck(self.cf,Base='General',ThisOne='OutputLevel'):
            InLevel = self.cf['General']['InputLevel']
            OutLevel = self.cf['General']['OutputLevel']
        else:
            InLevel = 'L2'
            OutLevel = 'L2'
        
        if 'in_file_path' in self.cf['Files'][OutLevel]:
            infilename = qcio.get_infilename_from_cf(self.cf,OutLevel)
        else:
            infilename = qcio.get_infilename_from_cf(self.cf,InLevel)
        if len(infilename)==0: self.do_progress(text='An error occurred, check the console ...'); return
        outfilename = qcio.get_outfilename_from_cf(self.cf,OutLevel)
        if len(outfilename)==0: self.do_progress(text='An error occurred, check the console ...'); return
        outputlist = qcio.get_outputlist_from_cf(self.cf,'xl')
        if qcutils.cfkeycheck(self.cf,Base='Output',ThisOne='Sort') and self.cf['Output']['Sort'] == 'False':
            sort = False
        else:
            sort = True
        qcio.nc_2xls(self.cf,infilename,outfilename,outputlist=outputlist,sortoption=sort)
        self.do_progress(text='Finished L2 Data Export')              # tell the user we are done
        log.info(' Finished saving L2 data')
        print '\a'

    def do_savexL3(self):
        """
            Call interior_nc2xl.autonc2xl function
            Exports excel data from NetCDF file
            
            Outputs L3 Excel file containing Data and Flag worksheets
            ControlFile:
                nc2xl_year
                L3d_nc2xl_OzFlux_Lvl3.txt
            ControlFile contents (see ControlFile/Templates/nc2xl.txt for
            example):
                [General] (optional):
                    InputLevel
                [Files]:
                    Input netCDF file name and path
                    Output excel file name and path
                [Output]:
                    Variables included in output file
            """
        self.do_progress(text='Exporting L3 NetCDF -> Xcel ...')                     # put up the progress message
        log.info('Saving L3 to Excel')
        # pdb.set_trace()
        check = 'True'
        try:
            self.cf
            if len(self.cf)==0:
                self.cf = qcio.load_controlfile(path='controlfiles')
                if len(self.cf)==0: self.do_progress(text='Waiting for input ...'); return
                check = 'False'
        except:
            self.cf = qcio.load_controlfile(path='controlfiles')
            if len(self.cf)==0: self.do_progress(text='Waiting for input ...'); return
            check = 'False'
        
        if (qcutils.cfkeycheck(self.cf,Base='Output',ThisOne='DefaultXl') and self.cf['Output']['DefaultXl'] == 'False') and check == 'True':
            self.cf = qcio.load_controlfile(path='controlfiles')
            if len(self.cf)==0: self.do_progress(text='Waiting for input ...'); return
        
        if qcutils.cfkeycheck(self.cf,Base='General',ThisOne='InputLevel') and qcutils.cfkeycheck(self.cf,Base='General',ThisOne='OutputLevel'):
            InLevel = self.cf['General']['InputLevel']
            OutLevel = self.cf['General']['OutputLevel']
        else:
            InLevel = 'L3'
            OutLevel = 'L3'
        
        if 'in_file_path' in self.cf['Files'][OutLevel]:
            infilename = qcio.get_infilename_from_cf(self.cf,OutLevel)
        else:
            infilename = qcio.get_infilename_from_cf(self.cf,InLevel)
        if len(infilename)==0: self.do_progress(text='An error occurred, check the console ...'); return
        outfilename = qcio.get_outfilename_from_cf(self.cf,OutLevel)
        if len(outfilename)==0: self.do_progress(text='An error occurred, check the console ...'); return
        outputlist = qcio.get_outputlist_from_cf(self.cf,'xl')
        if qcutils.cfkeycheck(self.cf,Base='Output',ThisOne='Sort') and self.cf['Output']['Sort'] == 'False':
            sort = False
        else:
            sort = True
        qcio.nc_2xls(self.cf,infilename,outfilename,outputlist=outputlist,sortoption=sort)
        self.do_progress(text='Finished L3 Data Export')              # tell the user we are done
        log.info(' Finished saving L3 data')
        print '\a'

    def do_savexL4(self):
        """
            Call interior_nc2xl.autonc2xl function
            Exports excel data from NetCDF file
            
            Outputs L4 Excel file containing Data and Flag worksheets
            ControlFile:
                nc2xl_year
                L4c_nc2xl_OzFlux_Lvl4.txt
            ControlFile contents (see ControlFile/Templates/nc2xl.txt for
            example):
                [General] (optional):
                    InputLevel
                [Files]:
                    Input netCDF file name and path
                    Output excel file name and path
                [Output]:
                    Variables included in output file
            """
        self.do_progress(text='Exporting L4 NetCDF -> Xcel ...')                     # put up the progress message
        log.info('Saving L4, L5 or L6 to Excel')
        check = 'True'
        try:
            self.cf
            if len(self.cf)==0:
                self.cf = qcio.loadc_ontrolfile('controlfiles')
                if len(self.cf)==0: self.do_progress(text='Waiting for input ...'); return
                check = 'False'
        except:
            self.cf = qcio.load_controlfile(path='controlfiles')
            if len(self.cf)==0: self.do_progress(text='Waiting for input ...'); return
            check = 'False'
        
        if (qcutils.cfkeycheck(self.cf,Base='Output',ThisOne='DefaultXl') and self.cf['Output']['DefaultXl'] == 'False') and check == 'True':
            self.cf = qcio.load_controlfile(path='controlfiles')
            if len(self.cf)==0: self.do_progress(text='Waiting for input ...'); return
        
        if qcutils.cfkeycheck(self.cf,Base='General',ThisOne='InputLevel') and qcutils.cfkeycheck(self.cf,Base='General',ThisOne='OutputLevel'):
            InLevel = self.cf['General']['InputLevel']
            OutLevel = self.cf['General']['OutputLevel']
        else:
            InLevel = 'L4'
            OutLevel = 'L4'
        
        if 'in_file_path' in self.cf['Files'][OutLevel]:
            infilename = qcio.get_infilename_from_cf(self.cf,OutLevel)
        else:
            infilename = qcio.get_infilename_from_cf(self.cf,InLevel)
        if len(infilename)==0: self.do_progress(text='An error occurred, check the console ...'); return
        outfilename = qcio.get_outfilename_from_cf(self.cf,OutLevel)
        if len(outfilename)==0: self.do_progress(text='An error occurred, check the console ...'); return
        outputlist = qcio.get_outputlist_from_cf(self.cf,'xl')
        if qcutils.cfkeycheck(self.cf,Base='Output',ThisOne='Sort') and self.cf['Output']['Sort'] == 'False':
            sort = False
        else:
            sort = True
        qcio.nc_2xls(self.cf,infilename,outfilename,outputlist=outputlist,sortoption=sort)
        outtext = 'Finished ' + str(OutLevel) + ' Data Export'
        self.do_progress(text=outtext)              # tell the user we are done
        log.info(' '+outtext)
        print '\a'

    def do_xl2ncL1(self):
        """
        Calls do_xl2nc with in_level set to L1
            Level 1:
                Read L1 Excel workbook
                Generate flags for missing observations
                Output L1 netCDF file to ncData folder
                Control file: L1.txt
        """
        log.info('Reading L1 from Excel')
        self.in_level = 'L1'
        self.do_xl2nc()
        infilename = qcio.get_outfilename_from_cf(self.cf,'L1')
        if len(infilename)==0: self.do_progress(text='An error occurred, check the console ...'); return
        self.ds0 = qcio.nc_read_series(infilename)
        self.ds1 = qcls.l1qc(self.cf,self.ds0)
        self.update_startenddate(str(self.ds1.series['DateTime']['Data'][0]),
                                 str(self.ds1.series['DateTime']['Data'][-1]))
        ncFile = qcio.nc_open_write(infilename)
        qcio.nc_write_series(ncFile,self.ds1)                   # save the L2 data
        self.do_progress(text='Finished saving L1 raw data/metadata')              # tell the user we are done
        log.info(' Finished saving L1 raw data/metadata')
        print '\a'

    def do_xl2ncL3(self):
        """
        Calls do_xl2nc with in_level set to L3
            Level 3:
                Read L3 Excel workbook
                Ingest flags generated in L3
                Output L3 netCDF file to ncData folder
                Control file: L3.txt
        """
        log.info('Reading L3 from Excel')
        self.in_level = 'L3'
        self.do_xl2nc()
        print '\a'

    def do_xl2ncL4(self):
        """
        Calls do_xl2nc with in_level set to L4
            Level 4:
                Read L4 Excel workbook
                Ingest flags generated in L3 or L4
                Output L4 netCDF file to ncData folder
                Control file: L4.txt
        """
        log.info('Reading L4, L5 or L6 from Excel')
        self.in_level = 'L4'
        self.do_xl2nc()
        print '\a'

    def do_xl2nc(self):
        """
        Calls qcio.xl2nc
        """
        self.do_progress(text='Loading control file ...')
        self.cf = qcio.load_controlfile(path='controlfiles')
        if len(self.cf)==0: self.do_progress(text='Waiting for input ...'); return
        if qcutils.cfkeycheck(self.cf,Base='General',ThisOne='InLevel'):
            self.in_level = self.cf['General']['InLevel']
        self.do_progress(text='Reading Excel file & writing to netCDF')
        rcode = qcio.xl2nc(self.cf,self.in_level)
        if rcode==0:
            self.do_progress(text='Finished writing to netCDF ...')
            log.info(' Finished writing to netCDF ...')
        else:
            self.do_progress(text='An error occurred, check the console ...')

    def load_l1nc(self):
        log.info('Loading L1 ncData to data structure')
        self.do_progress(text='Load L2 Control File ...')
        self.cf = qcio.load_controlfile(path='controlfiles')
        if len(self.cf)==0: self.do_progress(text='Waiting for input ...'); return
        self.do_progress(text='Re/Loading L1 NetCDF Data ...')
        outfilename = qcio.get_outfilename_from_cf(self.cf,'L1')
        if len(outfilename)==0: self.do_progress(text='An error occurred, check the console ...'); return
        self.ds1 = qcio.nc_read_series(outfilename)
        self.update_startenddate(str(self.ds1.series['DateTime']['Data'][0]),
                                 str(self.ds1.series['DateTime']['Data'][-1]))
        log.info(' Finished Re/Loading L1 NetCDF Data')
        self.do_progress(text='Finished Re/Loading L1 NetCDF Data')
        print '\a'
    
    def load_l2nc(self):
        log.info('Loading L2 ncData to data structure')
        self.do_progress(text='Load L2 Control File ...')
        self.cf = qcio.load_controlfile(path='controlfiles')
        if len(self.cf)==0: self.do_progress(text='Waiting for input ...'); return
        self.do_progress(text='Re/Loading L2 NetCDF Data ...')
        infilename = qcio.get_infilename_from_cf(self.cf,'L1')
        if len(infilename)==0: self.do_progress(text='An error occurred, check the console ...'); return
        outfilename = qcio.get_outfilename_from_cf(self.cf,'L2')
        if len(outfilename)==0: self.do_progress(text='An error occurred, check the console ...'); return
        self.ds1 = qcio.nc_read_series(infilename)
        self.ds2 = qcio.nc_read_series(outfilename)
        self.update_startenddate(str(self.ds2.series['DateTime']['Data'][0]),
                                 str(self.ds2.series['DateTime']['Data'][-1]))
        log.info(' Finished Re/Loading L2 NetCDF Data')
        self.do_progress(text='Finished Re/Loading L2 NetCDF Data')
        print '\a'
    
    def load_l3nc(self):
        log.info('Loading L3 ncData to data structure')
        self.do_progress(text='Load L3 Control File ...')
        self.cf = qcio.load_controlfile(path='controlfiles')
        if len(self.cf)==0: self.do_progress(text='Waiting for input ...'); return
        self.do_progress(text='Re/Loading L3 NetCDF Data ...')
        infilename = qcio.get_infilename_from_cf(self.cf,'L2')
        self.ds2 = qcio.nc_read_series(infilename)
        try:
            infilename = qcio.get_infilename_from_cf(self.cf,'L3',fail=True)
        except:
            infilename = qcio.get_outfilename_from_cf(self.cf,'L3')
        self.ds3 = qcio.nc_read_series(infilename)
        self.update_startenddate(str(self.ds3.series['DateTime']['Data'][0]),
                                 str(self.ds3.series['DateTime']['Data'][-1]))
        log.info(' Finished Re/Loading L3 NetCDF Data')
        self.do_progress(text='Finished Re/Loading L3 NetCDF Data')
        print '\a'
    
    def load_l4nc(self):
        log.info('Loading L4 ncData to data structure')
        self.do_progress(text='Load L4 Control File ...')
        self.cf = qcio.load_controlfile(path='controlfiles')
        if len(self.cf)==0: self.do_progress(text='Waiting for input ...'); return
        self.do_progress(text='Re/Loading L4 NetCDF Data ...')
        try:
            infilename = qcio.get_infilename_from_cf(self.cf,'L3')
        except:
            infilename = qcio.get_outfilename_from_cf(self.cf,'L3')
        self.ds3 = qcio.nc_read_series(infilename)
        try:
            infilename = qcio.get_infilename_from_cf(self.cf,'L4')
        except:
            infilename = qcio.get_outfilename_from_cf(self.cf,'L4')
        self.ds4 = qcio.nc_read_series(infilename)
        self.update_startenddate(str(self.ds4.series['DateTime']['Data'][0]),
                                 str(self.ds4.series['DateTime']['Data'][-1]))
        log.info(' Finished Re/Loading L4 NetCDF Data')
        self.do_progress(text='Finished Re/Loading L4 NetCDF Data')
        print '\a'
    
    def update_startenddate(self,startstr,endstr):
        """
            Read start and end timestamps from data and report in QC Data GUI
            """
        self.filestartValue.destroy()
        self.fileendValue.destroy()
        self.filestartValue = Tkinter.Label(self,text=startstr)
        self.filestartValue.grid(row=5,column=3,columnspan=2)
        self.fileendValue = Tkinter.Label(self,text=endstr)
        self.fileendValue.grid(row=5,column=5,columnspan=2)
        self.update()


if __name__ == "__main__":
    log = qcutils.startlog('qc','logfiles/qc.log')
    qcGUI = qcgui()
    main_title = cfg.version_name+' '+cfg.version_number+' Main GUI'
    qcGUI.master.title(main_title)
    qcGUI.mainloop()
    qcGUI.master.destroy()

    log.info('QC: All done')
