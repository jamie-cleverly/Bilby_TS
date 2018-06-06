#    pio.py
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

from configobj import ConfigObj
import ast
import datetime
import numpy
import os
import sys
import time
import constants as c
import Tkinter, tkFileDialog
import xlrd
import xlwt
import netCDF4
import logging
import pts
import putils

log = logging.getLogger('partition.io')

class DataStructure(object):
    def __init__(self):
        self.series = {}
        self.globalattributes = {}

def autonc2xl(cf,InLevel,OutLevel):
    # get the variables
    ds = nc_read_series(cf,InLevel)
    # write the variables to the excel file
    if putils.cfkeycheck(cf,Base='Output',ThisOne='Format') and cf['Output']['Format'] == 'Subset':
        xl_write_series(cf,ds,OutLevel)
    else:
        xl_write_library(cf,ds,OutLevel)

def autoxl2nc(cf,InLevel,OutLevel):
    # get the data series from the Excel file
    ds = xl_read_series(cf,InLevel)
    ds.globalattributes['Level'] = str(InLevel)
    # get the year, month, day, hour, minute and second from the xl date/time
    pts.get_yearmonthdayhourminutesecond(cf,ds)
    # get the quality control flags
    if InLevel == 'L1':
        pts.get_qcflag(ds)
    # get the flags from gap filled 'L3' or 'L4' Excel file
    if InLevel != 'L1':
        VariablesInFile = ds.series.keys()
        for ThisOne in ['xlDateTime','Gap','Year','Month','Day','Hour','Minute','Second','Hdh','Ddd']:
            if ThisOne in VariablesInFile:
                VariablesInFile.remove(ThisOne)
        ds1 = xl_read_flags(cf,ds,InLevel,VariablesInFile)
        if InLevel == 'L4':
            for ThisOne in ['Fc_gapfilled','Fe_gapfilled','Fh_gapfilled','Fc','Fe','Fh']:
                if ThisOne in ds.series.keys():
                    ds1.series[ThisOne]['Flag'] = numpy.int32(ds.series['Gap']['Data'])
        ds = ds1
    # do any functions to create new series
    pts.do_functions(cf,ds)
    # get the netCDF attributes from the control file
    pts.do_attributes(cf,ds)
    # write the data to the netCDF file
    nc_write_series(cf,ds,OutLevel)

def get_controlfilecontents(ControlFileName):
    log.info(' Processing the control file ')
    if len(ControlFileName)!=0:
        cf = ConfigObj(ControlFileName)
        cf['ControlFileName'] = ControlFileName
    else:
        cf = ConfigObj()
    return cf

def get_controlfilename(ControlFilePath):
    log.info(' Choosing the control file ')
    root = Tkinter.Tk(); root.withdraw()
    ControlFileName = tkFileDialog.askopenfilename(initialdir=ControlFilePath)
    root.destroy()
#    if len(ControlFileName)==0:
#        sys.exit()
    return ControlFileName

def get_datetime(ds):
    ''' Creates a series of Python datetime objects from the year, month,
    day, hour, minute and second series stored in the netCDF file.'''
    log.info(' Getting the date and time series')
    nRecs = len(ds.series['Year']['Data'])
    ds.series[unicode('DateTime')] = {}
    ds.series['DateTime']['Data'] = []
    for i in range(nRecs):
        ds.series['DateTime']['Data'].append(datetime.datetime(int(ds.series['Year']['Data'][i]),
                                                       int(ds.series['Month']['Data'][i]),
                                                       int(ds.series['Day']['Data'][i]),
                                                       int(ds.series['Hour']['Data'][i]),
                                                       int(ds.series['Minute']['Data'][i]),
                                                       int(ds.series['Second']['Data'][i])))
    ds.series['DateTime']['Flag'] = numpy.zeros(nRecs)
    ds.series['DateTime']['Attr'] = {}
    ds.series['DateTime']['Attr']['long_name'] = 'Date-time object'
    ds.series['DateTime']['Attr']['units'] = 'None'

def get_infilename_from_cf(cf,InLevel,fail=False):
    filename = ""
    if fail == True:
        filename = cf['Files'][InLevel]['in_file_path']+cf['Files'][InLevel]['in_filename']
        return str(filename)
    else:
        if "Files" in cf.keys():
            if "in_file_path" in cf['Files'][InLevel].keys():
                if "in_filename" in cf['Files'][InLevel].keys():
                    filename = cf['Files'][InLevel]['in_file_path']+cf['Files'][InLevel]['in_filename']
                else:
                    log.error("get_infilename_from_cf: 'in_filename' key not found in 'Files' section of control file")
            else:
                log.error("get_infilename_from_cf: 'file_path' key not found in 'Files' section of control file")
        else:
            log.error("get_infilename_from_cf: 'Files' section not found in control file")
        return str(filename)

def get_keyvalue_from_cf(section,key):
    try:
        value = section[key]
    except:
        log.error('get_keyvalue_from_cf: '+str(key)+' not found in '+str(section.name)+' section of control file')
        value = ''
    return value

def get_ncdtype(Series):
    sd = Series.dtype.name
    dt = 'f'
    if sd=='float64': dt = 'd'
    if sd=='int32': dt = 'i'
    if sd=='int64': dt = 'l'
    return dt

def get_ncfilename(path='.',title='Choose a netCDF file to open'):
    '''Get a netCDF file name'''
    root = Tkinter.Tk(); root.withdraw()
    ncFileName = tkFileDialog.askopenfilename(parent=root,initialdir=path,title=title)
    root.destroy()
    if len(ncFileName)==0:
        sys.exit()
    return ncFileName

def get_saveasfilename(path='.',title='Save file as'):
    '''Get a file name for saving'''
    root = Tkinter.Tk(); root.withdraw()
    SaveAsFileName = tkFileDialog.asksaveasfilename(parent=root,initialdir=path,title=title)
    root.destroy()
    if len(SaveAsFileName)==0:
        sys.exit()
    return SaveAsFileName

def get_xlfilename(path='.'):
    '''Get an Excel file name'''
    root = Tkinter.Tk(); root.withdraw()
    xlFileName = tkFileDialog.askopenfilename(parent=root,initialdir=path,title='Choose an Excel file to open')
    root.destroy
    return str(xlFileName)

def loadcontrolfile(ControlFilePath):
    ControlFileName = get_controlfilename(ControlFilePath)
    cf = get_controlfilecontents(ControlFileName)
    return cf

def nc_read_series(cf,level):
    ''' Read a netCDF file and put the data and meta-data into a DataStructure'''
    ncFullName = cf['Files'][level]['ncFilePath']+cf['Files'][level]['ncFileName']
    log.info(' Reading netCDF file '+ncFullName)
    netCDF4.default_encoding = 'latin-1'
    ncFile = netCDF4.Dataset(ncFullName,'r')
    ds = DataStructure()
    gattrlist = ncFile.ncattrs()
    if len(gattrlist)!=0:
        for gattr in gattrlist:
            ds.globalattributes[gattr] = getattr(ncFile,gattr)
            if 'time_step' in ds.globalattributes: c.ts = ds.globalattributes['time_step']
    for ThisOne in ncFile.variables.keys():
        if '_QCFlag' not in ThisOne:
            # create the series in the data structure
            ds.series[unicode(ThisOne)] = {}
            # get the data variable object
            ds.series[ThisOne]['Data'] = ncFile.variables[ThisOne][:]
            # check for a QC flag and if it exists, load it
            if ThisOne+'_QCFlag' in ncFile.variables.keys():
                ds.series[ThisOne]['Flag'] = ncFile.variables[ThisOne+'_QCFlag'][:]
            # get the variable attributes
            vattrlist = ncFile.variables[ThisOne].ncattrs()
            if len(vattrlist)!=0:
                ds.series[ThisOne]['Attr'] = {}
                for vattr in vattrlist:
                    ds.series[ThisOne]['Attr'][vattr] = getattr(ncFile.variables[ThisOne],vattr)
    ncFile.close()
    # get a series of Python datetime objects
    get_datetime(ds)
    return ds

def nc_read_file(ncFullName):
    netCDF4.default_encoding = 'latin-1'
    ncFile = netCDF4.Dataset(ncFullName,'r')
    ds = DataStructure()
    gattrlist = ncFile.ncattrs()
    if len(gattrlist)!=0:
        for gattr in gattrlist:
            ds.globalattributes[gattr] = getattr(ncFile,gattr)
            if 'time_step' in ds.globalattributes: c.ts = ds.globalattributes['time_step']
    for ThisOne in ncFile.variables.keys():
        if '_QCFlag' not in ThisOne:
            # create the series in the data structure
            ds.series[unicode(ThisOne)] = {}
            # get the data variable object
            ds.series[ThisOne]['Data'] = ncFile.variables[ThisOne][:]
            # check for a QC flag and if it exists, load it
            if ThisOne+'_QCFlag' in ncFile.variables.keys():
                ds.series[ThisOne]['Flag'] = ncFile.variables[ThisOne+'_QCFlag'][:]
            # get the variable attributes
            vattrlist = ncFile.variables[ThisOne].ncattrs()
            if len(vattrlist)!=0:
                ds.series[ThisOne]['Attr'] = {}
                for vattr in vattrlist:
                    ds.series[ThisOne]['Attr'][vattr] = getattr(ncFile.variables[ThisOne],vattr)
    ncFile.close()
    # get a series of Python datetime objects
    get_datetime(ds)
    return ds

def nc_write_OzFlux_series(cf,ds,level):
    ncFullName = cf['Files'][level]['ncFilePath']+cf['Files'][level]['ncFileName']
    log.info(' Writing netCDF file '+ncFullName)
    ncFile = netCDF4.Dataset(ncFullName,'w',format='NETCDF3_CLASSIC')
    for ThisOne in ds.globalattributes.keys():
        setattr(ncFile,ThisOne,ds.globalattributes[ThisOne])
    t = time.localtime()
    RunDateTime = str(datetime.datetime(t[0],t[1],t[2],t[3],t[4],t[5]))
    setattr(ncFile,'RunDateTime',RunDateTime)
    nRecs = len(ds.series['xlDateTime']['Data'])
    setattr(ncFile,'NumRecs',str(nRecs))
    setattr(ncFile,'Level',level)
    ncFile.createDimension('Time',nRecs)
    SeriesList = ast.literal_eval(cf['Output']['OFL2'])
    VariableList = ds.series.keys()
    for ThisOne in ds.series.keys():
        if ThisOne not in SeriesList:
            VariableList.remove(ThisOne)
    SeriesList = VariableList
    for ThisOne in ['xlDateTime','Year','Month','Day','Hour','Minute','Second','Hdh']:
        if ThisOne in SeriesList:
            dt = get_ncdtype(ds.series[ThisOne]['Data'])
            ncVar = ncFile.createVariable(ThisOne,dt,('Time',))
            ncVar[:] = ds.series[ThisOne]['Data'].tolist()
            setattr(ncVar,'long_name',ThisOne)
            setattr(ncVar,'units','none')
            SeriesList.remove(ThisOne)
    if 'DateTime' in SeriesList:
        SeriesList.remove('DateTime')
    for ThisOne in SeriesList:
        if 'Data' in ds.series[ThisOne].keys():
            dt = get_ncdtype(ds.series[ThisOne]['Data'])
            ncVar = ncFile.createVariable(ThisOne,dt,('Time',))
            ncVar[:] = ds.series[ThisOne]['Data'].tolist()
        if 'Attr' in ds.series[ThisOne].keys():
            for attr in ds.series[ThisOne]['Attr']:
                setattr(ncVar,attr,ds.series[ThisOne]['Attr'][attr])
        if 'Flag' in ds.series[ThisOne].keys():
            dt = get_ncdtype(ds.series[ThisOne]['Flag'])
            ncVar = ncFile.createVariable(ThisOne+'_QCFlag',dt,('Time',))
            ncVar[:] = ds.series[ThisOne]['Flag'].tolist()
            setattr(ncVar,'long_name','QC flag')
            setattr(ncVar,'units','none')
    ncFile.close()

def nc_write_series(cf,ds,level):
    ncFullName = cf['Files'][level]['ncFilePath']+cf['Files'][level]['ncFileName']
    log.info(' Writing netCDF file '+ncFullName)
    if putils.cfkeycheck(cf,Base='General',ThisOne='netCDFv3') and cf['General']['netCDFv3'] == 'False':
        ncFile = netCDF4.Dataset(ncFullName,'w')
    else:
        ncFile = netCDF4.Dataset(ncFullName,'w',format='NETCDF3_CLASSIC')
    for ThisOne in ds.globalattributes.keys():
        setattr(ncFile,ThisOne,ds.globalattributes[ThisOne])
    t = time.localtime()
    RunDateTime = str(datetime.datetime(t[0],t[1],t[2],t[3],t[4],t[5]))
    setattr(ncFile,'RunDateTime',RunDateTime)
    nRecs = len(ds.series['xlDateTime']['Data'])
    setattr(ncFile,'NumRecs',str(nRecs))
    setattr(ncFile,'Level',level)
    ncFile.createDimension('Time',nRecs)
    SeriesList = ds.series.keys()
    for ThisOne in ['xlDateTime','Year','Month','Day','Hour','Minute','Second','Hdh']:
        if ThisOne in SeriesList:
            dt = get_ncdtype(ds.series[ThisOne]['Data'])
            ncVar = ncFile.createVariable(ThisOne,dt,('Time',))
            ncVar[:] = ds.series[ThisOne]['Data'].tolist()
            setattr(ncVar,'long_name',ThisOne)
            setattr(ncVar,'units','none')
            SeriesList.remove(ThisOne)
    if 'DateTime' in SeriesList:
        SeriesList.remove('DateTime')
    for ThisOne in SeriesList:
        if 'Data' in ds.series[ThisOne].keys():
            dt = get_ncdtype(ds.series[ThisOne]['Data'])
            ncVar = ncFile.createVariable(ThisOne,dt,('Time',))
            ncVar[:] = ds.series[ThisOne]['Data'].tolist()
        if 'Attr' in ds.series[ThisOne].keys():
            for attr in ds.series[ThisOne]['Attr']:
                setattr(ncVar,attr,ds.series[ThisOne]['Attr'][attr])
        if 'Flag' in ds.series[ThisOne].keys():
            dt = get_ncdtype(ds.series[ThisOne]['Flag'])
            ncVar = ncFile.createVariable(ThisOne+'_QCFlag',dt,('Time',))
            ncVar[:] = ds.series[ThisOne]['Flag'].tolist()
            setattr(ncVar,'long_name','QC flag')
            setattr(ncVar,'units','none')
    ncFile.close()

def xl_read_flags(cf,ds,level,VariablesInFile):
    # First data row in Excel worksheets.
    FirstDataRow = int(cf['Files'][level]['xl1stDataRow']) - 1
    # Get the full name of the Excel file from the control file.
    xlFullName = cf['Files'][level]['xlFilePath']+cf['Files'][level]['xlFileName']
    # Get the Excel workbook object.
    if os.path.isfile(xlFullName):
        xlBook = xlrd.open_workbook(xlFullName)
    else:
        log.error(' Excel file '+xlFullName+' not found, choose another')
        xlFullName = get_xlfilename()
        if len(xlFullName)==0:
            return
        xlBook = xlrd.open_workbook(xlFullName)
    ds.globalattributes['xlFullName'] = xlFullName
    
    for ThisOne in VariablesInFile:
        if 'xl' in cf['Variables'][ThisOne].keys():
            log.info(' Getting flags for '+ThisOne+' from spreadsheet')
            ActiveSheet = xlBook.sheet_by_name('Flag')
            LastDataRow = int(ActiveSheet.nrows)
            HeaderRow = ActiveSheet.row_values(int(cf['Files'][level]['xlHeaderRow'])-1)
            if cf['Variables'][ThisOne]['xl']['name'] in HeaderRow:
                xlCol = HeaderRow.index(cf['Variables'][ThisOne]['xl']['name'])
                Values = ActiveSheet.col_values(xlCol)[FirstDataRow:LastDataRow]
                Types = ActiveSheet.col_types(xlCol)[FirstDataRow:LastDataRow]
                ds.series[ThisOne]['Flag'] = numpy.array([-9999]*len(Values),numpy.int32)
                for i in range(len(Values)):
                    if Types[i]==2: #xlType=3 means a date/time value, xlType=2 means a number
                        ds.series[ThisOne]['Flag'][i] = numpy.int32(Values[i])
                    else:
                        log.error('  xl_read_flags: flags for '+ThisOne+' not found in xl file')
    return ds

def xl_read_series(cf,level):
    # Instance the data structure object.
    ds = DataStructure()
    # First data row in Excel worksheets.
    FirstDataRow = int(cf['Files'][level]['xl1stDataRow']) - 1
    # Get the full name of the Excel file from the control file.
    xlFullName = cf['Files'][level]['xlFilePath']+cf['Files'][level]['xlFileName']
    # Get the Excel workbook object.
    if os.path.isfile(xlFullName):
        log.info(' Opening and reading Excel file '+xlFullName)
        xlBook = xlrd.open_workbook(xlFullName)
        log.info(' Opened and read Excel file '+xlFullName)
    else:
        log.error(' Excel file '+xlFullName+' not found, choose another')
        xlFullName = get_xlfilename()
        if len(xlFullName)==0:
            return
        log.info(' Opening and reading Excel file '+xlFullName)
        xlBook = xlrd.open_workbook(xlFullName)
        log.info(' Opened and read Excel file '+xlFullName)
    ds.globalattributes['xlFullName'] = xlFullName
    # Get the Excel file modification date and time, these will be
    # written to the netCDF file to uniquely identify the version
    # of the Excel file used to create this netCDF file.
    s = os.stat(xlFullName)
    t = time.localtime(s.st_mtime)
    ds.globalattributes['xlModDateTime'] = str(datetime.datetime(t[0],t[1],t[2],t[3],t[4],t[5]))
    # Loop over the variables defined in the 'Variables' section of the
    # configuration file.
    for ThisOne in cf['Variables'].keys():
        if 'xl' in cf['Variables'][ThisOne].keys():
            log.info(' Getting data for '+ThisOne+' from spreadsheet')
            ActiveSheet = xlBook.sheet_by_name(cf['Variables'][ThisOne]['xl']['sheet'])
            LastDataRow = int(ActiveSheet.nrows)
            HeaderRow = ActiveSheet.row_values(int(cf['Files'][level]['xlHeaderRow'])-1)
            if cf['Variables'][ThisOne]['xl']['name'] in HeaderRow:
                ds.series[unicode(ThisOne)] = {}
                xlCol = HeaderRow.index(cf['Variables'][ThisOne]['xl']['name'])
                Values = ActiveSheet.col_values(xlCol)[FirstDataRow:LastDataRow]
                Types = ActiveSheet.col_types(xlCol)[FirstDataRow:LastDataRow]
                ds.series[ThisOne]['Data'] = numpy.array([-9999]*len(Values),numpy.float64)
                for i in range(len(Values)):
                    if (Types[i]==3) or (Types[i]==2): #xlType=3 means a date/time value, xlType=2 means a number
                        ds.series[ThisOne]['Data'][i] = numpy.float64(Values[i])
            else:
                log.error('  xl_read_series: series '+ThisOne+' not found in xl file')
    return ds

def xl_read_series_hf(cf,InLevel):
    # Instance the data structure object.
    ds = DataStructure()
    # get the filename
    FileName = get_infilename_from_cf(cf,InLevel)
    if len(FileName)==0:
        log.error(' in_filename not found in control file')
        return ds
    if not os.path.exists(FileName):
        log.error(' Input file '+FileName+' specified in control file not found')
        return ds
    # convert from Excel row number to xlrd row number
    FirstDataRow = numpy.int32(get_keyvalue_from_cf(cf['Files'][InLevel],'in_firstdatarow')) - 1
    HeaderRow = numpy.int32(get_keyvalue_from_cf(cf['Files'][InLevel],'in_headerrow')) - 1
    # get the Excel workbook object.
    log.info(' Opening and reading Excel file '+FileName)
    xlBook = xlrd.open_workbook(FileName)
    log.info(' Opened and read Excel file '+FileName)
    ds.globalattributes['featureType'] = 'timeseries'
    ds.globalattributes['xl_filename'] = FileName
    ds.globalattributes['xl_datemode'] = str(xlBook.datemode)
    xlsheet_names = [x.lower() for x in xlBook.sheet_names()]
    # Get the Excel file modification date and time, these will be
    # written to the netCDF file to uniquely identify the version
    # of the Excel file used to create this netCDF file.
    s = os.stat(FileName)
    t = time.localtime(s.st_mtime)
    ds.globalattributes['xl_moddatetime'] = str(datetime.datetime(t[0],t[1],t[2],t[3],t[4],t[5]))
    # Loop over the variables defined in the 'Variables' section of the
    # configuration file.
    for ThisOne in cf['Variables'].keys():
        if 'xl' in cf['Variables'][ThisOne].keys():
            if 'sheet' in cf['Variables'][ThisOne]['xl'].keys():
                xlsheet_name = cf['Variables'][ThisOne]['xl']['sheet']
                if xlsheet_name.lower() in xlsheet_names:
                    log.info(' Getting data for '+ThisOne+' from spreadsheet')
                    xlsheet_index = xlsheet_names.index(xlsheet_name.lower())
                    ActiveSheet = xlBook.sheet_by_index(xlsheet_index)
                    HeaderList = [x.lower() for x in ActiveSheet.row_values(HeaderRow)]
                    if cf['Variables'][ThisOne]['xl']['name'].lower() in HeaderList:
                        LastDataRow = numpy.int32(ActiveSheet.nrows)
                        ds.series[unicode(ThisOne)] = {}
                        xlCol = HeaderList.index(cf['Variables'][ThisOne]['xl']['name'].lower())
                        Values = ActiveSheet.col_values(xlCol)[FirstDataRow:LastDataRow]
                        Types = ActiveSheet.col_types(xlCol)[FirstDataRow:LastDataRow]
                        ds.series[ThisOne]['Data'] = numpy.ones(len(Values),dtype=numpy.float64)*numpy.float64(c.missing_value)
                        ds.series[ThisOne]['Flag'] = numpy.ones(len(Values),dtype=numpy.int32)
                        # we could use "where" and get rid of this for loop
                        for i in range(len(Values)):
                            if (Types[i]==3) or (Types[i]==2): #xlType=3 means a date/time value, xlType=2 means a number
                                ds.series[ThisOne]['Data'][i] = numpy.float64(Values[i])
                                if ds.series[ThisOne]['Data'][i] != numpy.float64(c.missing_value):
                                    ds.series[ThisOne]['Flag'][i] = numpy.int32(0)
                    else:
                        log.error('  xl_read_series: series '+ThisOne+' not found in xl file')
                else:
                    log.error('  xl_read_series: sheet '+xlsheet_name+' not found in xl file')
            else:
                log.error('  xl_read_series: key "sheet" not found in control file entry for '+ThisOne)
        else:
            log.error('  xl_read_series: key "xl" not found in control file entry for '+ThisOne)
    ds.globalattributes['nc_nrecs'] = str(len(ds.series['xlDateTime']['Data']))
    return ds

def xl_write_library(cf,ds,level):
    log.info(' Opening the Excel file ')
    xlfullname = cf['Files'][level]['xlFilePath']+cf['Files'][level]['xlFileName']
    xlfile = xlwt.Workbook()
    if putils.cfkeycheck(cf,Base='General',ThisOne='Platform') and cf['General']['Platform'] == 'Mac':
        xlfile.dates_1904 = True
    # add sheets to the Excel file
    xlAttrSheet = xlfile.add_sheet('Attr')
    xlDataSheet = xlfile.add_sheet('Data')
    xlFlagSheet = xlfile.add_sheet('Flag')
    # write the global attributes
    log.info(' Writing the global attributes to Excel file '+xlfullname)
    xlcol = 0
    xlrow = 0
    xlAttrSheet.write(xlrow,xlcol,'Global attributes')
    xlrow = xlrow + 1
    globalattrlist = ds.globalattributes.keys()
    globalattrlist.sort()
    for ThisOne in [x for x in globalattrlist if 'Flag' not in x]:
        xlAttrSheet.write(xlrow,xlcol,ThisOne)
        xlAttrSheet.write(xlrow,xlcol+1,str(ds.globalattributes[ThisOne]))
        xlrow = xlrow + 1
    for ThisOne in [x for x in globalattrlist if 'Flag' in x]:
        xlAttrSheet.write(xlrow,xlcol,ThisOne)
        xlAttrSheet.write(xlrow,xlcol+1,str(ds.globalattributes[ThisOne]))
        xlrow = xlrow + 1
    # write the variable attributes
    log.info(' Writing the variable attributes to Excel file '+xlfullname)
    xlrow = xlrow + 1
    xlAttrSheet.write(xlrow,xlcol,'Variable attributes')
    xlrow = xlrow + 1
    xlcol_varname = 0
    xlcol_attrname = 1
    xlcol_attrvalue = 2
    variablelist = ds.series.keys()
    variablelist.sort()
    if 'DateTime' in variablelist:
        variablelist.remove('DateTime')
    for ThisOne in variablelist:
        xlAttrSheet.write(xlrow,xlcol_varname,ThisOne)
        attributelist = ds.series[ThisOne]['Attr'].keys()
        attributelist.sort()
        for Attr in attributelist:
            xlAttrSheet.write(xlrow,xlcol_attrname,Attr)
            xlAttrSheet.write(xlrow,xlcol_attrvalue,ds.series[ThisOne]['Attr'][Attr])
            xlrow = xlrow + 1
    # write the Excel date/time to the data and the QC flags as the first column
    log.info(' Writing the datetime to Excel file '+xlfullname)
    nRecs = len(ds.series['xlDateTime']['Data'])
    d_xf = xlwt.easyxf(num_format_str='dd/mm/yyyy hh:mm')
    xlDataSheet.write(2,xlcol,'xlDateTime')
    for j in range(nRecs):
        xlDataSheet.write(j+3,xlcol,ds.series['xlDateTime']['Data'][j],d_xf)
        xlFlagSheet.write(j+3,xlcol,ds.series['xlDateTime']['Data'][j],d_xf)
    # remove xlDateTime from the list of variables to be written to the Excel file
    variablelist.remove('xlDateTime')
    # now start looping over the other variables in the xl file
    xlcol = xlcol + 1
    # loop over variables to be output to xl file
    for ThisOne in variablelist:
        # put up a progress message
        log.info(' Writing '+ThisOne+' into column '+str(xlcol)+' of the Excel file')
        # write the units and the variable name to the header rows in the xl file
        attrlist = ds.series[ThisOne]['Attr'].keys()
        if 'long_name' in attrlist:
            longname = ds.series[ThisOne]['Attr']['long_name']
        elif 'Description' in attrlist:
            longname = ds.series[ThisOne]['Attr']['Description']
        else:
            longname = None
        if 'units' in attrlist:
            units = ds.series[ThisOne]['Attr']['units']
        elif 'Units' in attrlist:
            units = ds.series[ThisOne]['Attr']['Units']
        else:
            units = None
        xlDataSheet.write(0,xlcol,longname)
        xlDataSheet.write(1,xlcol,units)
        xlDataSheet.write(2,xlcol,ThisOne)
        # loop over the values in the variable series (array writes don't seem to work)
        for j in range(nRecs):
            xlDataSheet.write(j+3,xlcol,float(ds.series[ThisOne]['Data'][j]))
        # check to see if this variable has a quality control flag
        if 'Flag' in ds.series[ThisOne].keys():
            # write the QC flag name to the xk file
            xlFlagSheet.write(2,xlcol,ThisOne)
            # specify the format of the QC flag (integer)
            d_xf = xlwt.easyxf(num_format_str='0')
            # loop over QV flag values and write to xl file
            for j in range(nRecs):
                xlFlagSheet.write(j+3,xlcol,int(ds.series[ThisOne]['Flag'][j]),d_xf)
        # increment the column pointer
        xlcol = xlcol + 1
    
    log.info(' Saving the Excel file '+xlfullname)
    xlfile.save(xlfullname)

def xl_write_series(cf,ds,level):
    log.info(' Opening the Excel file ')
    nRecs = len(ds.series['xlDateTime']['Data'])
    xlFileName = cf['Files'][level]['xlFilePath']+cf['Files'][level]['xlFileName']
    xlFile = xlwt.Workbook()
    if cf['General']['Platform'] == 'Mac':
        xlFile.dates_1904 = True
    xlDataSheet = xlFile.add_sheet('Data')
    xlFlagSheet = xlFile.add_sheet('Flag')
    xlCol = 0
    VariablesInFile = ds.series.keys()
    VariablesToOutput = ast.literal_eval(cf['Output'][level])
    # write the xl date/time value to the first column of the worksheets
    d_xf = xlwt.easyxf(num_format_str='dd/mm/yyyy hh:mm')
    for j in range(nRecs):
        xlDataSheet.write(j+10,xlCol,ds.series['xlDateTime']['Data'][j],d_xf)
    xlFlagSheet.write(9,xlCol,'TIMESTAMP')
    xlDataSheet.write(9,xlCol,'TIMESTAMP')
    if not putils.cfkeycheck(cf,Base='Output',ThisOne='FlagDefs') or (putils.cfkeycheck(cf,Base='Output',ThisOne='FlagDefs') and cf['Output']['FlagDefs'] == 'True'):
        xlDataSheet.write(0,xlCol,'Site:')
        xlDataSheet.write(0,xlCol+1,ds.globalattributes['site'])
        xlDataSheet.write(2,xlCol,'Institution:')
        xlDataSheet.write(2,xlCol+1,ds.globalattributes['institution'])
        xlDataSheet.write(1,xlCol,'Latitude:')
        xlDataSheet.write(1,xlCol+1,ds.globalattributes['latitude'])
        xlDataSheet.write(1,xlCol+2,'Longitude:')
        xlDataSheet.write(1,xlCol+3,ds.globalattributes['longitude'])
        xlDataSheet.write(3,xlCol,'Contact:')
        xlDataSheet.write(3,xlCol+1,ds.globalattributes['contact'])
        xlFlagSheet.write(0,xlCol,'0:')
        xlFlagSheet.write(0,xlCol+1,ds.globalattributes['Flag000'])
        xlFlagSheet.write(0,xlCol+2,'1:')
        xlFlagSheet.write(0,xlCol+3,ds.globalattributes['Flag001'])
        xlFlagSheet.write(0,xlCol+4,'2:')
        xlFlagSheet.write(0,xlCol+5,ds.globalattributes['Flag002'])
        xlFlagSheet.write(0,xlCol+6,'3:')
        xlFlagSheet.write(0,xlCol+7,ds.globalattributes['Flag003'])
        xlFlagSheet.write(0,xlCol+8,'4:')
        xlFlagSheet.write(0,xlCol+9,ds.globalattributes['Flag004'])
        xlFlagSheet.write(0,xlCol+10,'5:')
        xlFlagSheet.write(0,xlCol+11,ds.globalattributes['Flag005'])
        xlFlagSheet.write(0,xlCol+12,'6:')
        xlFlagSheet.write(0,xlCol+13,ds.globalattributes['Flag006'])
        xlFlagSheet.write(0,xlCol+14,'7:')
        xlFlagSheet.write(0,xlCol+15,ds.globalattributes['Flag007'])
        xlFlagSheet.write(0,xlCol+16,'8:')
        xlFlagSheet.write(0,xlCol+17,ds.globalattributes['Flag008'])
        xlFlagSheet.write(0,xlCol+18,'9:')
        xlFlagSheet.write(0,xlCol+19,ds.globalattributes['Flag009'])
        xlFlagSheet.write(1,xlCol,'10:')
        xlFlagSheet.write(1,xlCol+1,ds.globalattributes['Flag010'])
        xlFlagSheet.write(1,xlCol+2,'11:')
        xlFlagSheet.write(1,xlCol+3,ds.globalattributes['Flag011'])
        xlFlagSheet.write(1,xlCol+4,'12:')
        xlFlagSheet.write(1,xlCol+5,ds.globalattributes['Flag012'])
        xlFlagSheet.write(1,xlCol+6,'13:')
        xlFlagSheet.write(1,xlCol+7,ds.globalattributes['Flag013'])
        xlFlagSheet.write(1,xlCol+8,'14:')
        xlFlagSheet.write(1,xlCol+9,ds.globalattributes['Flag014'])
        xlFlagSheet.write(1,xlCol+10,'15:')
        xlFlagSheet.write(1,xlCol+11,ds.globalattributes['Flag015'])
        xlFlagSheet.write(1,xlCol+12,'16:')
        xlFlagSheet.write(1,xlCol+13,ds.globalattributes['Flag016'])
        xlFlagSheet.write(1,xlCol+14,'17:')
        xlFlagSheet.write(1,xlCol+15,ds.globalattributes['Flag017'])
        xlFlagSheet.write(1,xlCol+16,'18:')
        xlFlagSheet.write(1,xlCol+17,ds.globalattributes['Flag018'])
        xlFlagSheet.write(1,xlCol+18,'19:')
        xlFlagSheet.write(1,xlCol+19,ds.globalattributes['Flag019'])
        xlFlagSheet.write(1,xlCol+20,'21:')
        xlFlagSheet.write(1,xlCol+21,ds.globalattributes['Flag021'])
        xlFlagSheet.write(1,xlCol+22,'22:')
        xlFlagSheet.write(1,xlCol+23,ds.globalattributes['Flag022'])
        xlFlagSheet.write(2,xlCol,'30:')
        xlFlagSheet.write(2,xlCol+1,ds.globalattributes['Flag030'])
        xlFlagSheet.write(2,xlCol+2,'31:')
        xlFlagSheet.write(2,xlCol+3,ds.globalattributes['Flag031'])
        xlFlagSheet.write(2,xlCol+4,'40:')
        xlFlagSheet.write(2,xlCol+5,ds.globalattributes['Flag040'])
        xlFlagSheet.write(2,xlCol+6,'50:')
        xlFlagSheet.write(2,xlCol+7,ds.globalattributes['Flag050'])
        xlFlagSheet.write(2,xlCol+8,'60:')
        xlFlagSheet.write(2,xlCol+9,ds.globalattributes['Flag060'])
        xlFlagSheet.write(2,xlCol+10,'70:')
        xlFlagSheet.write(2,xlCol+11,ds.globalattributes['Flag070'])
        xlFlagSheet.write(2,xlCol+12,'80:')
        xlFlagSheet.write(2,xlCol+13,ds.globalattributes['Flag080'])
        xlFlagSheet.write(2,xlCol+14,'81:')
        xlFlagSheet.write(2,xlCol+15,ds.globalattributes['Flag081'])
        xlFlagSheet.write(2,xlCol+16,'82:')
        xlFlagSheet.write(2,xlCol+17,ds.globalattributes['Flag082'])
        xlFlagSheet.write(2,xlCol+18,'83:')
        xlFlagSheet.write(2,xlCol+19,ds.globalattributes['Flag083'])
        xlFlagSheet.write(2,xlCol+20,'84:')
        xlFlagSheet.write(2,xlCol+21,ds.globalattributes['Flag084'])
        xlFlagSheet.write(2,xlCol+22,'85:')
        xlFlagSheet.write(2,xlCol+23,ds.globalattributes['Flag085'])
        xlFlagSheet.write(2,xlCol+24,'86:')
        xlFlagSheet.write(2,xlCol+25,ds.globalattributes['Flag086'])
        xlFlagSheet.write(2,xlCol+26,'87:')
        xlFlagSheet.write(2,xlCol+27,ds.globalattributes['Flag087'])
        xlFlagSheet.write(3,xlCol,'90:')
        xlFlagSheet.write(3,xlCol+1,ds.globalattributes['Flag090'])
        xlFlagSheet.write(3,xlCol+2,'100:')
        xlFlagSheet.write(3,xlCol+3,ds.globalattributes['Flag100'])
        xlFlagSheet.write(3,xlCol+4,'110:')
        xlFlagSheet.write(3,xlCol+5,ds.globalattributes['Flag110'])
        xlFlagSheet.write(3,xlCol+6,'120:')
        xlFlagSheet.write(3,xlCol+7,ds.globalattributes['Flag120'])
        xlFlagSheet.write(3,xlCol+8,'130:')
        xlFlagSheet.write(3,xlCol+9,ds.globalattributes['Flag130'])
        xlFlagSheet.write(3,xlCol+10,'140:')
        xlFlagSheet.write(3,xlCol+11,ds.globalattributes['Flag140'])
        xlFlagSheet.write(3,xlCol+12,'150:')
        xlFlagSheet.write(3,xlCol+13,ds.globalattributes['Flag150'])
        xlFlagSheet.write(3,xlCol+14,'151:')
        xlFlagSheet.write(3,xlCol+15,ds.globalattributes['Flag151'])
        xlFlagSheet.write(4,xlCol,'161:')
        xlFlagSheet.write(4,xlCol+1,ds.globalattributes['Flag161'])
        xlFlagSheet.write(4,xlCol+2,'162:')
        xlFlagSheet.write(4,xlCol+3,ds.globalattributes['Flag162'])
        xlFlagSheet.write(5,xlCol,'171:')
        xlFlagSheet.write(5,xlCol+1,ds.globalattributes['Flag171'])
        xlFlagSheet.write(5,xlCol+2,'173:')
        xlFlagSheet.write(5,xlCol+3,ds.globalattributes['Flag173'])
        xlFlagSheet.write(5,xlCol+4,'174:')
        xlFlagSheet.write(5,xlCol+5,ds.globalattributes['Flag174'])
        xlFlagSheet.write(5,xlCol+6,'180:')
        xlFlagSheet.write(5,xlCol+7,ds.globalattributes['Flag180'])
        xlFlagSheet.write(5,xlCol+8,'190:')
        xlFlagSheet.write(5,xlCol+9,ds.globalattributes['Flag190'])
        xlFlagSheet.write(5,xlCol+10,'191:')
        xlFlagSheet.write(5,xlCol+11,ds.globalattributes['Flag191'])
        xlFlagSheet.write(5,xlCol+12,'200:')
        xlFlagSheet.write(5,xlCol+13,ds.globalattributes['Flag200'])
        xlFlagSheet.write(5,xlCol+14,'201:')
        xlFlagSheet.write(5,xlCol+15,ds.globalattributes['Flag201'])
        xlFlagSheet.write(6,xlCol,'211:')
        xlFlagSheet.write(6,xlCol+1,ds.globalattributes['Flag211'])
    
    #d_xf = xlwt.easyxf('font: height 160',num_format_str='dd/mm/yyyy hh:mm')
    d_xf = xlwt.easyxf(num_format_str='dd/mm/yyyy hh:mm')
    for j in range(nRecs):
        xlFlagSheet.write(j+10,xlCol,ds.series['xlDateTime']['Data'][j],d_xf)
    # remove the date and time variables from the list to output
    for ThisOne in ['xlDateTime','Year','Month','Day','Hour','Minute','Second','Hdh']:
        if ThisOne in VariablesInFile:
            VariablesInFile.remove(ThisOne)
    # now start looping over the other variables in the xl file
    xlCol = xlCol + 1
    # list of variables in the data structure
    #VariablesInFile.sort(key=str.lower)
    # list of variables to write out (specified in the control file)
    #VariablesToOutput.sort(key=str.lower)
    # loop over variables to be output to xl file
    for ThisOne in VariablesToOutput:
        if ThisOne in VariablesInFile:
            # put up a progress message
            log.info(' Writing '+ThisOne+' into column '+str(xlCol)+' of the Excel file')
            # specify the style of the output
            #d_xf = xlwt.easyxf('font: height 160')
            d_xf = xlwt.easyxf()
            # write the units and the variable name to the header rows in the xl file
            longname = ds.series[ThisOne]['Attr']['long_name']
            units = ds.series[ThisOne]['Attr']['units']
            #xlDataSheet.write(8,xlCol,Units,d_xf)
            #xlDataSheet.write(9,xlCol,ThisOne,d_xf)
            xlDataSheet.write(7,xlCol,longname)
            xlDataSheet.write(8,xlCol,units)
            xlDataSheet.write(9,xlCol,ThisOne)
            # loop over the values in the variable series (array writes don't seem to work)
            for j in range(nRecs):
                #xlDataSheet.write(j+10,xlCol,float(ds.series[ThisOne]['Data'][j]),d_xf)
                xlDataSheet.write(j+10,xlCol,float(ds.series[ThisOne]['Data'][j]))
            # check to see if this variable has a quality control flag
            if 'Flag' in ds.series[ThisOne].keys():
                # write the QC flag name to the xk file
                #xlFlagSheet.write(9,xlCol,ThisOne,d_xf)
                xlFlagSheet.write(9,xlCol,ThisOne)
                # specify the format of the QC flag (integer)
                #d_xf = xlwt.easyxf('font: height 160',num_format_str='0')
                d_xf = xlwt.easyxf(num_format_str='0')
                # loop over QV flag values and write to xl file
                for j in range(nRecs):
                    xlFlagSheet.write(j+10,xlCol,ds.series[ThisOne]['Flag'][j],d_xf)
            # increment the column pointer
            xlCol = xlCol + 1
    # tell the user what we are doing
    log.info(' Saving the Excel file '+xlFileName)
    # save the xl file
    xlFile.save(xlFileName)

