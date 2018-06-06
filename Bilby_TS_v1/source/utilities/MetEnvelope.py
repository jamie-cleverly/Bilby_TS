#    MetEnvelope.py
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
    Used to NEE into GPP and Re
    
    Nighttime Re:  determined from exponential relationship between temperature (nocturnal mean) and Fc (nocturnal sum) on nights without missing observations; binned in soil moisture classes
    Daytime Re:  determined from light response curve between Fsi and Fc; binned in temperature classes
    
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
import metio
import metls
import metutils


class qcgui(Tkinter.Frame):
    def __init__(self, master=None):
        Tkinter.Frame.__init__(self, master)
        self.grid()
        self.createWidgets()

    def createWidgets(self):
        self.process2Label = Tkinter.Label(self,text='Ts x Sws Blocks')
        self.process2Label.grid(row=0,column=2,columnspan=1)
        self.process2Label = Tkinter.Label(self,text='D x Sws Blocks')
        self.process2Label.grid(row=0,column=3,columnspan=1)
        self.process2Label = Tkinter.Label(self,text='Fsd x Sws Blocks')
        self.process2Label.grid(row=0,column=4,columnspan=1)
        self.process2Label = Tkinter.Label(self,text='Fsd x Sws x Ts Blocks')
        self.process2Label.grid(row=0,column=5,columnspan=1)
        self.process2Label = Tkinter.Label(self,text='Fsd x Sws x VPD Blocks')
        self.process2Label.grid(row=0,column=6,columnspan=1)
        
        self.fileloadLabel = Tkinter.Label(self,text='Xcel ->')
        self.fileloadLabel.grid(row=1,column=0,columnspan=1)
        self.doxl2nc1Button = Tkinter.Button (self, text="Load Corrected Data", command=self.do_xl2ncCall )
        self.doxl2nc1Button.grid(row=1,column=1,columnspan=1)
        
        self.initiateLabel = Tkinter.Label(self,text='Process Data')
        self.initiateLabel.grid(row=2,column=0,columnspan=1)
        self.doL2Button = Tkinter.Button (self, text="Apply Envelope Filters", command=self.do_qcTs )
        self.doL2Button.grid(row=2,column=2,columnspan=1)
        self.doL2Button = Tkinter.Button (self, text="Apply Envelope Filters", command=self.do_qcD )
        self.doL2Button.grid(row=2,column=3,columnspan=1)
        self.doL2Button = Tkinter.Button (self, text="Apply Envelope Filters", command=self.do_qcEs )
        self.doL2Button.grid(row=2,column=4,columnspan=1)
        self.doL2Button = Tkinter.Button (self, text="Apply Envelope Filters", command=self.do_qcLRF )
        self.doL2Button.grid(row=2,column=5,columnspan=1)
        self.doL2Button = Tkinter.Button (self, text="Apply Envelope Filters", command=self.do_qcLRFD )
        self.doL2Button.grid(row=2,column=6,columnspan=1)
        
        self.filesave2Label = Tkinter.Label(self,text='-> Xcel')
        self.filesave2Label.grid(row=3,column=0,columnspan=1)
        self.savexLButton = Tkinter.Button (self, text="Save Filtered Data", command=self.do_savexL )
        self.savexLButton.grid(row=3,column=2,columnspan=1)
        self.savexL4Button = Tkinter.Button (self, text="Save Filtered Data", command=self.do_savexL )
        self.savexL4Button.grid(row=3,column=3,columnspan=1)
        self.savexL5Button = Tkinter.Button (self, text="Save Filtered Data", command=self.do_savexL )
        self.savexL5Button.grid(row=3,column=4,columnspan=1)
        self.savexL5Button = Tkinter.Button (self, text="Save Filtered Data", command=self.do_savexL )
        self.savexL5Button.grid(row=3,column=5,columnspan=1)
        self.savexL5Button = Tkinter.Button (self, text="Save Filtered Data", command=self.do_savexL )
        self.savexL5Button.grid(row=3,column=6,columnspan=1)
        
        self.quitButton = Tkinter.Button (self, text="Quit", command=self.do_quit )
        self.quitButton.grid(row=4,column=0,columnspan=1)
        self.progress = Tkinter.Label(self, text='Waiting for input ...')
        self.progress.grid(row=4,column=1,columnspan=6)

    def do_qcEs(self):
        self.cf = metio.loadcontrolfile('../controlfiles')
        if len(self.cf)==0:
            self.do_progress(text='Waiting for input ...')
            return
        if metutils.cfkeycheck(self.cf,Base='General',ThisOne='InputLevel'):
            InLevel = self.cf['General']['InputLevel']
        else:
            InLevel = 'L3'
        self.ds2 = metio.nc_read_series(self.cf,InLevel)
        self.do_progress(text='Doing met envelope prep...')
        self.ds3 = metls.l3envelope(self.cf,self.ds2,Params='Es')
        self.do_progress(text='Finished met envelope prep')
        self.do_progress(text='Saving met envelope NetCDF data ...')                     # put up the progress message
        if metutils.cfkeycheck(self.cf,Base='General',ThisOne='OutputLevel'):
            OutLevel = self.cf['General']['OutputLevel']
        else:
            OutLevel = 'MetEnvelope'
        metio.nc_write_series(self.cf,self.ds3,OutLevel)                   # save the L3 data
        self.do_progress(text='Finished Blocking')              # tell the user we are done
        log.info(' Finished saving NetCDF data')
    
    def do_qcLRF(self):
        self.cf = metio.loadcontrolfile('../controlfiles')
        if len(self.cf)==0:
            self.do_progress(text='Waiting for input ...')
            return
        if metutils.cfkeycheck(self.cf,Base='General',ThisOne='InputLevel'):
            InLevel = self.cf['General']['InputLevel']
        else:
            InLevel = 'L3'
        self.ds2 = metio.nc_read_series(self.cf,InLevel)
        self.do_progress(text='Doing met envelope prep...')
        self.ds3 = metls.l3envelope(self.cf,self.ds2,Params='LRF')
        self.do_progress(text='Finished met envelope prep')
        self.do_progress(text='Saving met envelope NetCDF data ...')                     # put up the progress message
        if metutils.cfkeycheck(self.cf,Base='General',ThisOne='OutputLevel'):
            OutLevel = self.cf['General']['OutputLevel']
        else:
            OutLevel = 'MetEnvelope'
        metio.nc_write_series(self.cf,self.ds3,OutLevel)                   # save the L3 data
        self.do_progress(text='Finished Blocking')              # tell the user we are done
        log.info(' Finished saving NetCDF data')
    
    def do_qcD(self):
        self.cf = metio.loadcontrolfile('../controlfiles')
        if len(self.cf)==0:
            self.do_progress(text='Waiting for input ...')
            return
        if metutils.cfkeycheck(self.cf,Base='General',ThisOne='InputLevel'):
            InLevel = self.cf['General']['InputLevel']
        else:
            InLevel = 'L3'
        self.ds2 = metio.nc_read_series(self.cf,InLevel)
        self.do_progress(text='Doing met envelope prep...')
        self.ds3 = metls.l3envelope(self.cf,self.ds2,Params='D')
        self.do_progress(text='Finished met envelope prep')
        self.do_progress(text='Saving met envelope NetCDF data ...')                     # put up the progress message
        if metutils.cfkeycheck(self.cf,Base='General',ThisOne='OutputLevel'):
            OutLevel = self.cf['General']['OutputLevel']
        else:
            OutLevel = 'MetEnvelope'
        metio.nc_write_series(self.cf,self.ds3,OutLevel)                   # save the L3 data
        self.do_progress(text='Finished Blocking')              # tell the user we are done
        log.info(' Finished saving NetCDF data')
    
    def do_qcTs(self):
        self.cf = metio.loadcontrolfile('../controlfiles')
        if len(self.cf)==0:
            self.do_progress(text='Waiting for input ...')
            return
        if metutils.cfkeycheck(self.cf,Base='General',ThisOne='InputLevel'):
            InLevel = self.cf['General']['InputLevel']
        else:
            InLevel = 'L3'
        self.ds2 = metio.nc_read_series(self.cf,InLevel)
        self.do_progress(text='Doing met envelope prep...')
        self.ds3 = metls.l3envelope(self.cf,self.ds2,Params='Ts')
        self.do_progress(text='Finished met envelope prep')
        self.do_progress(text='Saving met envelope NetCDF data ...')                     # put up the progress message
        if metutils.cfkeycheck(self.cf,Base='General',ThisOne='OutputLevel'):
            OutLevel = self.cf['General']['OutputLevel']
        else:
            OutLevel = 'MetEnvelope'
        metio.nc_write_series(self.cf,self.ds3,OutLevel)                   # save the L3 data
        self.do_progress(text='Finished Blocking')              # tell the user we are done
        log.info(' Finished saving NetCDF data')
    
    def do_qcLRFD(self):
        self.cf = metio.loadcontrolfile('../controlfiles')
        if len(self.cf)==0:
            self.do_progress(text='Waiting for input ...')
            return
        if metutils.cfkeycheck(self.cf,Base='General',ThisOne='InputLevel'):
            InLevel = self.cf['General']['InputLevel']
        else:
            InLevel = 'L3'
        self.ds2 = metio.nc_read_series(self.cf,InLevel)
        self.do_progress(text='Doing met envelope prep...')
        self.ds3 = metls.l3envelope(self.cf,self.ds2,Params='LRFD')
        self.do_progress(text='Finished met envelope prep')
        self.do_progress(text='Saving met envelope NetCDF data ...')                     # put up the progress message
        if metutils.cfkeycheck(self.cf,Base='General',ThisOne='OutputLevel'):
            OutLevel = self.cf['General']['OutputLevel']
        else:
            OutLevel = 'MetEnvelope'
        metio.nc_write_series(self.cf,self.ds3,OutLevel)                   # save the L3 data
        self.do_progress(text='Finished Blocking')              # tell the user we are done
        log.info(' Finished saving NetCDF data')
    
    def do_progress(self,text):
        self.progress.destroy()
        self.progress = Tkinter.Label(self, text=text)
        self.progress.grid(row=4,column=1,columnspan=3)
        self.update()
    
    def do_quit(self):
        self.do_progress(text='Quitting ...')                         # tell the user what we're doing
        log.info(' Quitting ...')
        self.quit()
    
    def do_savexL(self):
        check = 'True'
        self.do_progress(text='Exporting Diurnal NetCDF -> Xcel ...')                     # put up the progress message
        #if metutils.cfkeycheck(self.cf,'Output','noDefaultXl'):
        #    self.cf = metio.loadcontrolfile('../controlfiles')
        #    if len(self.cf)==0:
        #        self.do_progress(text='Waiting for input ...')
        #        return
        try:
            self.cf
            if len(self.cf)==0:
                self.cf = metio.loadcontrolfile('../controlfiles')
                if len(self.cf)==0: self.do_progress(text='Waiting for input ...'); return
                check = 'False'
        except:
            self.cf = metio.loadcontrolfile('../controlfiles')
            if len(self.cf)==0: self.do_progress(text='Waiting for input ...'); return
            check = 'False'
        
        if (metutils.cfkeycheck(self.cf,Base='Output',ThisOne='DefaultXl') and self.cf['Output']['DefaultXl'] == 'False') and check == 'True':
            self.cf = metio.loadcontrolfile('../controlfiles')
            if len(self.cf)==0: self.do_progress(text='Waiting for input ...'); return
        
        
        if metutils.cfkeycheck(self.cf,Base='General',ThisOne='OutputLevel'):
            OutLevel = self.cf['General']['OutputLevel']
        else:
            OutLevel = 'Envelope'
        metio.autonc2xl(self.cf,OutLevel)
        self.do_progress(text='Finished Data Export')              # tell the user we are done
        log.info(' Finished saving data')
    
    def do_xl2ncCall(self):
        self.do_progress(text='Load xl2nc Control File ...')
        self.cf = metio.loadcontrolfile('../controlfiles')
        if len(self.cf)==0:
            self.do_progress(text='Waiting for input ...')
            return
        self.do_progress(text='Importing Xcel file -> NetCDF v4 ...')
        if metutils.cfkeycheck(self.cf,'General','InLevel') and metutils.cfkeycheck(self.cf,'General','OutLevel'):
            InLevel = self.cf['General']['InLevel']
            OutLevel = self.cf['General']['OutLevel']
        else:
            InLevel = 'L3'
            OutLevel = 'L3'
        metio.autoxl2nc(self.cf,InLevel,OutLevel)
        self.do_progress(text='Finished Data Ingest')
        log.info(' Finished Data Ingest')


if __name__ == "__main__":
    log = metutils.startlog('envelope','../logfiles/envelope.log')
    qcGUI = qcgui()
    qcGUI.master.title("Meteorology-Carbon envelope main GUI")
    qcGUI.mainloop()
    qcGUI.master.destroy()

    print 'QC: All done'
