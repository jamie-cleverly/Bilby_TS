#    Partition.py
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
    Partition Data Main GUI
    Used to partition NEE into GPP and CE (CE = ER + abiotic decomposition)
    
    Nighttime ER:  determined from exponential relationship between temperature (nocturnal mean) and Fc (nocturnal sum) on nights without missing observations; binned in soil moisture classes
    Daytime ER:  determined from light response curve between Fsi and Fc; binned in temperature classes
    
    """

import sys

sys.path.append('scripts')

import ast
import copy
import datetime
import logging
import numpy
import time
import Tkinter
import pio
import pls
import putils
import pts


class qcgui(Tkinter.Frame):
    def __init__(self, master=None):
        Tkinter.Frame.__init__(self, master)
        self.grid()
        self.createWidgets()

    def createWidgets(self):
        self.process2Label = Tkinter.Label(self,text='L3: ER_night')
        self.process2Label.grid(row=0,column=1,columnspan=1)
        self.process3Label = Tkinter.Label(self,text='L6: ER & GPP')
        self.process3Label.grid(row=0,column=2,columnspan=1)
        self.process3Label = Tkinter.Label(self,text='L6: ER + AD & GPP (gaps)')
        self.process3Label.grid(row=0,column=3,columnspan=1)
        
        self.fileloadLabel = Tkinter.Label(self,text='Xcel ->')
        self.fileloadLabel.grid(row=1,column=0,columnspan=1)
        self.doxl2nc1Button = Tkinter.Button (self, text="Load Corrected Data", command=self.do_xl2ncCall )
        self.doxl2nc1Button.grid(row=1,column=1,columnspan=1)
        self.doxl2nc2Button = Tkinter.Button (self, text="Load Gapfilled Data", command=self.do_xl2ncCall )
        self.doxl2nc2Button.grid(row=1,column=2,columnspan=1)
        #self.doxl2nc2Button = Tkinter.Button (self, text="Conditional Correlation 10Hz", command=self.do_cc )
        #self.doxl2nc2Button.grid(row=1,column=3,columnspan=1)
        
        self.initiateLabel = Tkinter.Label(self,text='Process Data')
        self.initiateLabel.grid(row=2,column=0,columnspan=1)
        self.doL3Button = Tkinter.Button (self, text="ER_night", command=self.do_l3_er_night )
        self.doL3Button.grid(row=2,column=1,columnspan=1)
        self.doL6aButton = Tkinter.Button (self, text="Day CE & GPP", command=self.do_l6_partition )
        self.doL6aButton.grid(row=2,column=2,columnspan=1)
        
        self.filesave2Label = Tkinter.Label(self,text='-> Xcel')
        self.filesave2Label.grid(row=3,column=0,columnspan=1)
        self.savexL3Button = Tkinter.Button (self, text="Save L3 Day Data", command=self.do_savexL3 )
        self.savexL3Button.grid(row=3,column=1,columnspan=1)
        self.savexL6aButton = Tkinter.Button (self, text="Save L6 Data", command=self.do_savexL6 )
        self.savexL6aButton.grid(row=3,column=2,columnspan=1)
        
        self.quitButton = Tkinter.Button (self, text="Quit", command=self.do_quit )
        self.quitButton.grid(row=4,column=0,columnspan=1)
        self.progress = Tkinter.Label(self, text='Waiting for input ...')
        self.progress.grid(row=4,column=1,columnspan=2)

    def do_cc(self):
        self.do_progress(text='Performing conditional correlation')
        self.cf = pio.loadcontrolfile('../controlfiles')
        InLevel = 'L1'
        self.ds = pio.xl_read_series_hf(self.cf,InLevel)
        pts.conditional_correlation(self.cf,self.ds)
        self.do_progress(text='Finished conditional correlation')
        print '\a'
    
    def do_l3_er_night(self):
        self.cf = pio.loadcontrolfile('../controlfiles')
        if len(self.cf)==0:
            self.do_progress(text='Waiting for input ...')
            return
        if putils.cfkeycheck(self.cf,Base='General',ThisOne='InputLevel'):
            InLevel = self.cf['General']['InputLevel']
        else:
            InLevel = 'L3'
        self.ds2 = pio.nc_read_series(self.cf,InLevel)
        self.do_progress(text='Doing partitioning prep...')
        self.ds3 = pls.l3partition(self.cf,self.ds2)
        self.do_progress(text='Finished partitioning prep')
        self.do_progress(text='Saving L3 Partitioning NetCDF data ...')                     # put up the progress message
        if putils.cfkeycheck(self.cf,Base='General',ThisOne='OutputLevel'):
            OutLevel = self.cf['General']['OutputLevel']
        else:
            OutLevel = 'Partitioning'
        pio.nc_write_series(self.cf,self.ds3,OutLevel)                   # save the L3 data
        self.do_progress(text='Finished L3 ER_night Partitioning')              # tell the user we are done
        log.info(' Finished saving L3 Partitioning NetCDF data')
        print '\a'
    
    def do_l6_partition(self):
        self.cf = pio.loadcontrolfile('../controlfiles')
        if len(self.cf)==0:
            self.do_progress(text='Waiting for input ...')
            return
        if putils.cfkeycheck(self.cf,Base='General',ThisOne='InputLevel'):
            InLevel = self.cf['General']['InputLevel']
        else:
            InLevel = 'L5'
        self.ds5 = pio.nc_read_series(self.cf,InLevel)
        self.do_progress(text='Partitioning daytime ER_dark and GPP...')
        self.ds6 = pls.l6partition(self.cf,self.ds5)
        self.do_progress(text='Finished partitioning daytime ER_dark and GPP')
        self.do_progress(text='Saving L4 Partitioning NetCDF data ...')                     # put up the progress message
        if putils.cfkeycheck(self.cf,Base='General',ThisOne='OutputLevel'):
            OutLevel = self.cf['General']['OutputLevel']
        else:
            OutLevel = 'Partitioning'
        pio.nc_write_series(self.cf,self.ds6,OutLevel)                   # save the L3 data
        self.do_progress(text='Finished L6 Partitioning')              # tell the user we are done
        log.info(' Finished saving L6 Partitioning NetCDF data')
        print '\a'
    
    def do_progress(self,text):
        self.progress.destroy()
        self.progress = Tkinter.Label(self, text=text)
        self.progress.grid(row=4,column=1,columnspan=3)
        self.update()
        print '\a'
    
    def do_quit(self):
        self.do_progress(text='Quitting ...')                         # tell the user what we're doing
        log.info(' Quitting ...')
        self.quit()
    
    def do_savexL3(self):
        self.do_progress(text='Exporting L3 Diurnal NetCDF -> Xcel ...')                     # put up the progress message
        if (putils.cfkeycheck(self.cf,'Output','DefaultXl') and self.cf['Output']['DefaultXl'] == 'False'):
            self.cf = pio.loadcontrolfile('../controlfiles')
            if len(self.cf)==0:
                self.do_progress(text='Waiting for input ...')
                return
        if putils.cfkeycheck(self.cf,Base='General',ThisOne='InputLevel') and putils.cfkeycheck(self.cf,Base='General',ThisOne='OutputLevel'):
            InLevel = self.cf['General']['InputLevel']
            OutLevel = self.cf['General']['OutputLevel']
        else:
            InLevel = 'L3'
            OutLevel = 'Partitioning'
        
        pio.autonc2xl(self.cf,InLevel,OutLevel)
        self.do_progress(text='Finished L3 Data Export')              # tell the user we are done
        log.info(' Finished saving L3 data')
        print '\a'
    
    def do_savexL6(self):
        self.do_progress(text='Exporting L4 Diurnal NetCDF -> Xcel ...')                     # put up the progress message
        if (putils.cfkeycheck(self.cf,'Output','DefaultXl') and self.cf['Output']['DefaultXl'] == 'False'):
            self.cf = pio.loadcontrolfile('../controlfiles')
            if len(self.cf)==0:
                self.do_progress(text='Waiting for input ...')
                return
        if putils.cfkeycheck(self.cf,Base='General',ThisOne='InputLevel') and putils.cfkeycheck(self.cf,Base='General',ThisOne='OutputLevel'):
            InLevel = self.cf['General']['InputLevel']
            OutLevel = self.cf['General']['OutputLevel']
        else:
            InLevel = 'L3'
            OutLevel = 'Partitioning'
        
        pio.autonc2xl(self.cf,InLevel,OutLevel)
        self.do_progress(text='Finished L4 Data Export')              # tell the user we are done
        log.info(' Finished saving L4 data')
        print '\a'
    
    def do_xl2ncCall(self):
        self.do_progress(text='Load xl2nc Control File ...')
        self.cf = pio.loadcontrolfile('../controlfiles')
        if len(self.cf)==0:
            self.do_progress(text='Waiting for input ...')
            return
        self.do_progress(text='Importing Xcel file -> NetCDF v4 ...')
        if putils.cfkeycheck(self.cf,'General','InLevel') and putils.cfkeycheck(self.cf,'General','OutLevel'):
            InLevel = self.cf['General']['InLevel']
            OutLevel = self.cf['General']['OutLevel']
        else:
            InLevel = 'L3'
            OutLevel = 'L3'
        pio.autoxl2nc(self.cf,InLevel,OutLevel)
        self.do_progress(text='Finished Data Ingest')
        log.info(' Finished Data Ingest')
        print '\a'


if __name__ == "__main__":
    log = putils.startlog('partition','../logfiles/partition.log')
    qcGUI = qcgui()
    qcGUI.master.title("Carbon Partitioning Main GUI")
    qcGUI.mainloop()
    qcGUI.master.destroy()

    print 'QC: All done'
