#    qcio.py  Handles input/output functions, data ingest and file saving
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

from configobj import ConfigObj
import ast
import cfg
import copy
import constants as c
import csv
import datetime
import dateutil
import meteorologicalfunctions as mf
import numpy
import os
import platform
import sys
import time
import Tkinter, tkFileDialog
import xlrd
import xlwt
import xlsxwriter
import netCDF4
import logging
import qcts
import qcutils
import pdb

log = logging.getLogger('qc.io')

class DataStructure(object):
    def __init__(self):
        self.series = {}
        self.globalattributes = {}
        self.dimensions = {}
        self.mergeserieslist = []
        self.averageserieslist = []
        self.soloserieslist = []
        self.climatologyserieslist = []

def convert_v27tov28():
    """ Convert V2.7 (1D) netCDF files to V2.8 (3D). """
    # get the file names
    ncV27name = get_filename_dialog(path="../Sites")
    ncV28name = ncV27name.replace(".nc","_V28.nc")
    # read the V2.7 file
    ds = nc_read_series(ncV27name)
    # add the "time_zone" global attribute if it is not present
    if "time_zone" not in ds.globalattributes.keys():
        for gattr in ["site_name","SiteName"]:
            if gattr in ds.globalattributes.keys():
                time_zone = qcutils.get_timezone(ds.globalattributes[gattr])
        ds.globalattributes["time_zone"] = time_zone
    # add the "missing_value" attribute if it is not present
    for ThisOne in ds.series.keys():
        if "missing_value" not in ds.series[ThisOne]["Attr"].keys():
            ds.series[ThisOne]["Attr"]["missing_value"] = c.missing_value
    # write the V2.8 file
    ncFile = nc_open_write(ncV28name,nctype='NETCDF4')
    nc_write_series(ncFile, ds)

def convert_L3Corrected(ds):
    if ds.globalattributes['site_name'] == "AliceSpringsMulga":
        for ThisOne in ['Ts','dTs','Sws','Fg_bs','Fg_ms','Fg_mu','Sws_bs','Sws_ms','Sws_mu','Ts_bs','Ts_ms','Ts_mu','dTs_bs','dTs_ms','dTs_mu','svwc_s_baresoil','svwc_s_understory','svwc_s_mulga']:
            index = numpy.where(ds.series[ThisOne]['Flag']>10)[0]
            if len(index) > 0:
                ds.series[ThisOne]['Data'][index] = numpy.float64(c.missing_value)
                ds.series[ThisOne]['Flag'][index] = numpy.int32(1)
                ds.series['Fg']['Data'][index] = numpy.float64(c.missing_value)
                ds.series['Fg']['Flag'][index] = numpy.int32(1)
                ds.series['Fa']['Data'][index] = numpy.float64(c.missing_value)
                ds.series['Fa']['Flag'][index] = numpy.int32(1)
    
    elif ds.globalattributes['site_name'] == "TiTreeEast":
        for ThisOne in ['Ts','dTs','Sws','Fg_spinifex','Fg_mulga','Fg_grass','Sws_spinifex_0cm','Sws_mulga_0cm','Sws_grass_0cm','Ts_spinifex_0cm','Ts_mulga_0cm','dTs_mulga','dTs_spinifex']:
            index = numpy.where(ds.series[ThisOne]['Flag']>10)[0]
            if len(index) > 0:
                ds.series[ThisOne]['Data'][index] = numpy.float64(c.missing_value)
                ds.series[ThisOne]['Flag'][index] = numpy.int32(1)
                ds.series['Fg']['Data'][index] = numpy.float64(c.missing_value)
                ds.series['Fg']['Flag'][index] = numpy.int32(1)
                ds.series['Fa']['Data'][index] = numpy.float64(c.missing_value)
                ds.series['Fa']['Flag'][index] = numpy.int32(1)
    
    # write the percentage of good data as a variable attribute
    qcutils.get_coverage_individual(ds)
    # write the percentage of good data for groups
    qcutils.get_coverage_groups(ds)
    return ds


def copy_datastructure(cf,ds_in):
    '''
    Return a copy of a data structure based on the following rules:
     1) if the netCDF file at the "copy_to" level does not exist
        then copy the existing data structure at the "input" level
        to create a new data structure at the "output" level.
    '''
    # assumptions that need to be checked are:
    #  - the start datetime of the two sets of data are the same
    #  - the end datetime of the L3 data is the same or after the
    #    end datetime of the the L4 data
    #    - if the end datetimes are the same then we are just re-processing something
    #    - if the end datetime for the L3 data is after the end date of the L4 data
    #      then more data has been added to this year and the user wants to gap fill
    #      the new data
    # modificatons to be made:
    #  - check the modification datetime of the L3 and L4 files:
    #     - if the L3 file is newer than the L4 file the disregard the "UseExistingOutFile" setting
    # get the output (L4) file name
    ct_filename = cf['Files']['file_path']+cf['Files']['out_filename']
    # if the L4 file does not exist then create the L4 data structure as a copy
    # of the L3 data structure
    if not os.path.exists(ct_filename):
        ds_out = copy.deepcopy(ds_in)
    # if the L4 file does exist ...
    if os.path.exists(ct_filename):
        # check to see if the user wants to use it
        if cf['Options']['UseExistingOutFile']!='Yes':
            # if the user doesn't want to use the existing L4 data then create
            # the L4 data structure as a copy of the L3 data structure
            ds_out = copy.deepcopy(ds_in)
        else:
            # the user wants to use the data from an existing L4 file
            # get the netCDF file name at the "input" level
            outfilename = get_outfilename_from_cf(cf,OutLevel)
            # read the netCDF file at the "input" level
            ds_file = nc_read_series(outfilename)
            dt_file = ds_file.series['DateTime']['Data']
            sd_file = str(dt_file[0])
            ed_file = str(dt_file[-1])
            # create a copy of the data
            ds_out = copy.deepcopy(ds_in)
            dt_out = ds_out.series['DateTime']['Data']
            ts = ds_out.globalattributes['time_step']
            sd_out = str(dt_out[0])
            ed_out = str(dt_out[-1])
            # get the start and end indices based on the start and end dates
            si = qcutils.GetDateIndex(dt_out,sd_file,ts=ts,default=0,match='exact')
            ei = qcutils.GetDateIndex(dt_out,ed_file,ts=ts,default=-1,match='exact')
            # now replace parts of ds_out with the data read from file
            for ThisOne in ds_file.series.keys():
                # check to see if the L4 series exists in the L3 data
                if ThisOne in ds_out.series.keys():
                    # ds_out is the copy of the L3 data, now fill it with the L4 data read from file
                    ds_out.series[ThisOne]['Data'][si:ei+1] = ds_file.series[ThisOne]['Data']
                    ds_out.series[ThisOne]['Flag'][si:ei+1] = ds_file.series[ThisOne]['Flag']
                else:
                    # if it doesn't, create the series and put the data into it
                    ds_out.series[ThisOne] = {}
                    ds_out.series[ThisOne] = ds_file.series[ThisOne].copy()
                    # check to see if we have to append data to make the copy of the L4 data now
                    # in the L3 data structure the same length as the existing L3 data
                    nRecs_file = numpy.int32(ds_file.globalattributes['nc_nrecs'])
                    nRecs_out = numpy.int32(ds_out.globalattributes['nc_nrecs'])
                    if nRecs_file < nRecs_out:
                        # there is more data at L3 than at L4
                        # append missing data to make the series the same length
                        nRecs_append = nRecs_out - nRecs_file
                        data = numpy.array([c.missing_value]*nRecs_append,dtype=numpy.float64)
                        flag = numpy.ones(nRecs_append,dtype=numpy.int32)
                        ds_out.series[ThisOne]['Data'] = numpy.concatenate((ds_out.series[ThisOne]['Data'],data))
                        ds_out.series[ThisOne]['Flag'] = numpy.concatenate((ds_out.series[ThisOne]['Flag'],flag))
                    elif nRecs_file > nRecs_out:
                        # tell the user something is wrong
                        log.error('copy_datastructure: L3 file contains less data than L4 file')
                        # return an empty dictionary
                        ds_out = {}
                    else:
                        # nRecs_file and nRecs_out are equal so we do not need to do anything
                        pass
    return ds_out

def nc_2xls(cf,ncfilename,xlfilename,outputlist=None,sortoption=True):
    # read the netCDF file
    ds = nc_read_series(ncfilename)
    # write the variables to the Excel 97/2003 file
    xlfilename= xlfilename.replace('.nc','.xls')
    if "nc_nrecs" in ds.globalattributes.keys():
        nRecs = numpy.int32(ds.globalattributes["nc_nrecs"])
    else:
        variablelist = ds.series.keys()
        nRecs = len(ds.series[variablelist[0]]["Data"])
    # write the variables to the Excel file
    # pdb.set_trace()
    if sortoption == False:
        if nRecs<65535:
            xl_write_series_unsorted(cf,ds,xlfilename,outputlist=outputlist)
        else:
            xlsx_write_series_unsorted(cf,ds,xlfilename,outputlist=outputlist)
    else:
        if nRecs<65535:
            xl_write_series(ds,xlfilename,outputlist=outputlist)
        else:
            xlsx_write_series(ds,xlfilename,outputlist=outputlist)

def read_eddypro_full(csvname):
    ds = DataStructure()
    csvfile = open(csvname,'rb')
    csvreader = csv.reader(csvfile)
    n = 0
    adatetime = []
    us_data_list = []
    us_flag_list = []
    Fh_data_list = []
    Fh_flag_list = []
    Fe_data_list = []
    Fe_flag_list = []
    Fc_data_list = []
    Fc_flag_list = []
    for row in csvreader:
        if n==0:
            header=row
        elif n==1:
            varlist=row
            us_data_col = varlist.index('u*')
            us_flag_col = varlist.index('qc_Tau')
            Fh_data_col = varlist.index('H')
            Fh_flag_col = varlist.index('qc_H')
            Fe_data_col = varlist.index('LE')
            Fe_flag_col = varlist.index('qc_LE')
            Fc_data_col = varlist.index('co2_flux')
            Fc_flag_col = varlist.index('qc_co2_flux')
        elif n==2:
            unitlist=row
        else:
            adatetime.append(datetime.datetime.strptime(row[1]+' '+row[2],'%Y-%m-%d %H:%M'))
            us_data_list.append(numpy.float64(row[us_data_col]))
            us_flag_list.append(numpy.float64(row[us_flag_col]))
            Fh_data_list.append(numpy.float64(row[Fh_data_col]))
            Fh_flag_list.append(numpy.float64(row[Fh_flag_col]))
            Fe_data_list.append(numpy.float64(row[Fe_data_col]))
            Fe_flag_list.append(numpy.float64(row[Fe_flag_col]))
            Fc_data_list.append(numpy.float64(row[Fc_data_col]))
            Fc_flag_list.append(numpy.float64(row[Fc_flag_col]))
        n = n + 1
    nRecs = len(adatetime)
    adatetime = qcutils.RoundDateTime(adatetime,dt=30)
    ds.series['DateTime'] = {}
    ds.series['DateTime']['Data'] = adatetime
    ds.series['ustar'] = {}
    ds.series['ustar']['Data'] = numpy.array(us_data_list,dtype=numpy.float64)
    ds.series['ustar']['Flag'] = numpy.array(us_flag_list,dtype=numpy.int32)
    ds.series['Fh'] = {}
    ds.series['Fh']['Data'] = numpy.array(Fh_data_list,dtype=numpy.float64)
    ds.series['Fh']['Flag'] = numpy.array(Fh_flag_list,dtype=numpy.int32)
    ds.series['Fe'] = {}
    ds.series['Fe']['Data'] = numpy.array(Fe_data_list,dtype=numpy.float64)
    ds.series['Fe']['Flag'] = numpy.array(Fe_flag_list,dtype=numpy.int32)
    ds.series['Fc'] = {}
    ds.series['Fc']['Data'] = numpy.array(Fc_data_list,dtype=numpy.float64)
    ds.series['Fc']['Flag'] = numpy.array(Fc_flag_list,dtype=numpy.int32)
    ds.globalattributes["nc_nrecs"] = nRecs
    return ds
    
def xl2nc(cf,InLevel):
    # get the data series from the Excel file
    ds = xl_read_series(cf,InLevel)
    if len(ds.series.keys())==0: return 1
    if InLevel != 'L1':
        VariablesInFile = ds.series.keys()
        for ThisOne in ['xlDateTime','Gap','Year','Month','Day','Hour','Minute','Second','Hdh','Ddd']:
            if ThisOne in VariablesInFile:
                VariablesInFile.remove(ThisOne)
        ds1 = xl_read_flags(cf,ds,InLevel,VariablesInFile)
        if 'Gap' in ds.series.keys():
            for ThisOne in ['NEE','NEP','Fc_gapfilled','Fco2_gapfilled','Fe_gapfilled','Fh_gapfilled','Fc','Fc_c','Fc_co2','Fe','Fh']:
                if ThisOne in ds.series.keys():
                    ds1.series[ThisOne]['Flag'] = numpy.int32(ds.series['Gap']['Data'])
        ds = ds1
    # get the netCDF attributes from the control file
    qcts.do_attributes(cf,ds)
    # get a series of Python datetime objects from the Excel datetime
    qcutils.get_datetimefromxldate(ds)
    # get series of UTC datetime
    qcutils.get_UTCfromlocaltime(ds)
    if ((qcutils.cfkeycheck(cf,Base='General',ThisOne='FixTimeGaps')) and (cf['General']['FixTimeGaps'] == 'False')):
        ftype = 'skip'
    else:
        ftype = 'gaps'
    #check for gaps in the Excel datetime series
    has_gaps = qcutils.CheckTimeStep(ds,fix=ftype)
    # write the processing level to a global attribute
    ds.globalattributes['nc_level'] = str(InLevel)
    # get the start and end date from the datetime series unless they were
    # given in the control file
    if 'start_date' not in ds.globalattributes.keys():
        ds.globalattributes['start_date'] = str(ds.series['DateTime']['Data'][0])
    if 'end_date' not in ds.globalattributes.keys():
        ds.globalattributes['end_date'] = str(ds.series['DateTime']['Data'][-1])
    # get the year, month, day, hour, minute and second from the xl date/time
    qcutils.get_ymdhmsfromxldate(ds)
    # do any functions to create new series
    qcts.do_functions(cf,ds)
    # create a series of synthetic downwelling shortwave radiation
    qcts.get_synthetic_fsd(ds)
    # write the data to the netCDF file
    outfilename = get_outfilename_from_cf(cf,InLevel)
    ncFile = nc_open_write(outfilename)
    nc_write_series(ncFile,ds)
    return 0

def fn_write_csv(cf):
    # get the file names
    ncFileName = get_infilename_from_cf(cf)
    csvFileName = get_outfilename_from_cf(cf)
    # open the csv file
    csvfile = open(csvFileName,'wb')
    writer = csv.writer(csvfile)
    # read the netCDF file
    ds = nc_read_series(ncFileName)
    # Tumbarumba doesn't have RH in the netCDF files
    if "RH" not in ds.series.keys():
        Ah,f,a = qcutils.GetSeriesasMA(ds,'Ah')
        Ta,f,a = qcutils.GetSeriesasMA(ds,'Ta')
        RH = mf.RHfromabsolutehumidity(Ah, Ta)
        attr = qcutils.MakeAttributeDictionary(long_name='Relative humidity',units='%',standard_name='relative_humidity')
        qcutils.CreateSeries(ds,"RH",RH,FList=['Ta','Ah'],Attr=attr)
    ts = numpy.int32(ds.globalattributes["time_step"])
    ts_delta = datetime.timedelta(minutes=ts)
    # get the datetime series
    dt = ds.series["DateTime"]["Data"]
    # check the start datetime of the series and adjust if necessary
    start_datetime = dateutil.parser.parse(str(cf["General"]["start_datetime"]))
    if dt[0]<start_datetime:
        # requested start_datetime is after the start of the file
        log.info(" Truncating start of file")
        si = qcutils.GetDateIndex(dt,str(start_datetime),ts=ts,match="exact")
        for thisone in ds.series.keys():
            ds.series[thisone]["Data"] = ds.series[thisone]["Data"][si:]
            ds.series[thisone]["Flag"] = ds.series[thisone]["Flag"][si:]
        ds.globalattributes["nc_nrecs"] = str(len(ds.series["DateTime"]["Data"]))
    elif dt[0]>start_datetime:
        # requested start_datetime is before the start of the file
        log.info(" Padding start of file")
        dt_patched = [ldt for ldt in qcutils.perdelta(start_datetime, dt[0], ts_delta)]
        data_patched = numpy.ones(len(dt_patched),dtype=numpy.float64)*numpy.float64(c.missing_value)
        flag_patched = numpy.ones(len(dt_patched),dtype=numpy.int32)
        # list of series in the data structure
        series_list = ds.series.keys()
        # ds.series["DateTime"]["Data"] is a list not a numpy array so we must treat it differently
        ds.series["DateTime"]["Data"] = dt_patched+ds.series["DateTime"]["Data"]
        ds.series["DateTime"]["Flag"] = numpy.concatenate((flag_patched,ds.series["DateTime"]["Flag"]))
        series_list.remove("DateTime")
        for thisone in series_list:
            ds.series[thisone]["Data"] = numpy.concatenate((data_patched,ds.series[thisone]["Data"]))
            ds.series[thisone]["Flag"] = numpy.concatenate((flag_patched,ds.series[thisone]["Flag"]))
        ds.globalattributes["nc_nrecs"] = str(len(ds.series["DateTime"]["Data"]))
        # refresh the year, month, day etc arrays now that we have padded the datetime series
        qcutils.get_ymdhms_from_datetime(ds)
    # now check the end datetime of the file
    end_datetime = dateutil.parser.parse(str(cf["General"]["end_datetime"]))
    if dt[-1]>end_datetime:
        # requested end_datetime is before the end of the file
        msg = " Truncating end of file "+dt[-1].strftime("%Y-%m-%d %H:%M")+" "+end_datetime.strftime("%Y-%m-%d %H:%M")
        log.info(msg)
        ei = qcutils.GetDateIndex(dt,str(end_datetime),ts=ts,match="exact")
        for thisone in ds.series.keys():
            ds.series[thisone]["Data"] = ds.series[thisone]["Data"][:ei+1]
            ds.series[thisone]["Flag"] = ds.series[thisone]["Flag"][:ei+1]
        ds.globalattributes["nc_nrecs"] = str(len(ds.series["DateTime"]["Data"]))
    elif dt[-1]<end_datetime:
        # requested start_datetime is before the start of the file
        msg = " Padding end of file "+dt[-1].strftime("%Y-%m-%d %H:%M")+" "+end_datetime.strftime("%Y-%m-%d %H:%M")
        log.info(msg)
        dt_patched = [ldt for ldt in qcutils.perdelta(dt[-1]+ts_delta, end_datetime+ts_delta, ts_delta)]
        data_patched = numpy.ones(len(dt_patched),dtype=numpy.float64)*numpy.float64(c.missing_value)
        flag_patched = numpy.ones(len(dt_patched),dtype=numpy.int32)
        # list of series in the data structure
        series_list = ds.series.keys()
        # ds.series["DateTime"]["Data"] is a list not a numpy array so we must treat it differently
        ds.series["DateTime"]["Data"] = ds.series["DateTime"]["Data"]+dt_patched
        ds.series["DateTime"]["Flag"] = numpy.concatenate((ds.series["DateTime"]["Flag"],flag_patched))
        series_list.remove("DateTime")
        for thisone in series_list:
            ds.series[thisone]["Data"] = numpy.concatenate((ds.series[thisone]["Data"],data_patched))
            ds.series[thisone]["Flag"] = numpy.concatenate((ds.series[thisone]["Flag"],flag_patched))
        ds.globalattributes["nc_nrecs"] = str(len(ds.series["DateTime"]["Data"]))
        # refresh the year, month, day etc arrays now that we have padded the datetime series
        qcutils.get_ymdhms_from_datetime(ds)
    if ts==30:
        nRecs_year = 17520
        nRecs_leapyear = 17568
    elif ts==60:
        nRecs_year = 8760
        nRecs_leapyear = 8784
    else:
        log.error(" Unrecognised time step ("+str(ts)+")")
        return
    if (numpy.int32(ds.globalattributes["nc_nrecs"])!=nRecs_year) & (numpy.int32(ds.globalattributes["nc_nrecs"])!=nRecs_leapyear):
        log.error(" Number of records in file does not equal "+str(nRecs_year)+" or "+str(nRecs_leapyear))
        log.error(len(ds.series["DateTime"]["Data"]),ds.series["DateTime"]["Data"][0],ds.series["DateTime"]["Data"][-1])
        return
    # get the date and time data
    Day,flag,attr = qcutils.GetSeries(ds,'Day')
    Month,flag,attr = qcutils.GetSeries(ds,'Month')
    Year,flag,attr = qcutils.GetSeries(ds,'Year')
    Hour,flag,attr = qcutils.GetSeries(ds,'Hour')
    Minute,flag,attr = qcutils.GetSeries(ds,'Minute')
    # get the data
    data = {}
    series_list = cf["Variables"].keys()
    for series in series_list:
        ncname = cf["Variables"][series]["ncname"]
        if ncname not in ds.series.keys():
            log.error("Series "+ncname+" not in netCDF file, skipping ...")
            series_list.remove(series)
            continue
        data[series] = ds.series[ncname]
        fmt = cf["Variables"][series]["format"]
        if "." in fmt:
            numdec = len(fmt) - (fmt.index(".") + 1)
            strfmt = "{0:."+str(numdec)+"f}"
        else:
            strfmt = "{0:d}"
        data[series]["fmt"] = strfmt
    #adjust units if required
    for series in series_list:
        if series=="FC" and data[series]["Attr"]["units"]=='mg/m2/s':
            data[series]["Data"] = mf.Fc_umolpm2psfrommgpm2ps(data[series]["Data"])
            data[series]["Attr"]["units"] = "umol/m2/s"
        if series=="CO2" and data[series]["Attr"]["units"]=='mg/m3':
            CO2 = data["CO2"]["Data"]
            TA = data["TA"]["Data"]
            PA = data["PA"]["Data"]
            data[series]["Data"] = mf.co2_ppmfrommgpm3(CO2,TA,PA)
            data[series]["Attr"]["units"] = "umol/mol"
        if series=="H2O" and data[series]["Attr"]["units"]=='g/m3':
            H2O = data["H2O"]["Data"]
            TA = data["TA"]["Data"]
            PA = data["PA"]["Data"]
            data[series]["Data"] = mf.h2o_mmolpmolfromgpm3(H2O,TA,PA)
            data[series]["Attr"]["units"] = "mmol/mol"
        if series=="RH" and data[series]["Attr"]["units"] in ["fraction","frac"]:
            data[series]["Data"] = numpy.float64(100)*data[series]["Data"]
            data[series]["Attr"]["units"] = "%"
    # write the general information to csv file
    for item in cf["General"]:
        writer.writerow([item,str(cf['General'][item])])
    # write the variable names to the csv file
    row_list = ['DateTime','Year','Month','Day','HHMM']
    for item in series_list:
        row_list.append(item)
    writer.writerow(row_list)
    # write the units line to the csv file
    units_list = ["-","-","-","-","-"]
    for item in series_list:
        units_list.append(data[item]["Attr"]["units"])
    writer.writerow(units_list)
    # now write the data
    for i in range(len(Year)):
        # get the datetime string
        dtstr = '%02d/%02d/%d %02d:%02d'%(Day[i],Month[i],Year[i],Hour[i],Minute[i])
        hrmn = '%02d%02d'%(Hour[i],Minute[i])
        dttup = datetime.datetime(Year[i],Month[i],Day[i],Hour[i],Minute[i]).timetuple()
        doy = numpy.float64(dttup.tm_yday) + numpy.float64(dttup.tm_hour)/24 + numpy.float64(dttup.tm_min)/1440
        data_list = [dtstr,'%d'%(Year[i]),'%02d'%(Month[i]),'%02d'%(Day[i]),hrmn]
        for series in series_list:
            strfmt = data[series]["fmt"]
            if "d" in strfmt:
                data_list.append(strfmt.format(numpy.int32(round(data[series]["Data"][i]))))
            else:
                data_list.append(strfmt.format(data[series]["Data"][i]))
        writer.writerow(data_list)
    # close the csv file
    csvfile.close()
    return

def get_controlfilecontents(ControlFileName,mode="verbose"):
    if mode!="quiet": log.info(' Processing the control file '+ControlFileName)
    if len(ControlFileName)!=0:
        cf = ConfigObj(ControlFileName)
        cf['controlfile_name'] = ControlFileName
    else:
        cf = ConfigObj()
    return cf

def get_controlfilename(path=''):
    log.info(' Choosing the control file ')
    root = Tkinter.Tk(); root.withdraw()
    name = tkFileDialog.askopenfilename(initialdir=path)
    root.destroy()
    return name

def get_ncdtype(Series):
    sd = Series.dtype.name
    dt = 'f8'
    if sd=='int32' or sd== 'int64': dt = 'i4'
    return dt

def get_filename_dialog(path='.',title='Choose a file'):
    '''
    Put up a file open dialog.
    USEAGE:
     fname = qcio.get_filename_dialog(path=<path_to_file>,title=<tile>)
    INPUT:
     path  - the path to the file location, optional
     title - the title for the file open dialog, optional
    RETURNS:
     fname - the full file name including path, string
    '''
    root = Tkinter.Tk(); root.withdraw()
    FileName = tkFileDialog.askopenfilename(parent=root,initialdir=path,title=title)
    root.destroy()
    return str(FileName)

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

def get_outfilename_from_cf(cf,OutLevel):
    try:
        filename = cf['Files'][OutLevel]['out_file_path']+cf['Files'][OutLevel]['out_filename']
    except:
        log.error('get_outfilename_from_cf: Error getting out_filename from control file')
        filename = ''
    return str(filename)

def get_keyvalue_from_cf(section,key):
    try:
        value = section[key]
    except:
        log.error('get_keyvalue_from_cf: '+str(key)+' not found in '+str(section.name)+' section of control file')
        value = ''
    return value

def get_outputlist_from_cf(cf,filetype):
    # pdb.set_trace()
    try:
        outputlist = ast.literal_eval(cf['Output'][filetype])
    except:
        #log.info('get_outputlist_from_cf: Unable to get output list from Output section in control file')
        outputlist = None
    return outputlist

def get_seriesstats(cf,ds):
    # open an Excel file for the flag statistics
    level = ds.globalattributes['nc_level']
    prereplace = cf['Files'][level]['out_filename']
    xl_filename = cf['Files'][level]['flag_file_path']+prereplace.replace('.nc','_FlagStats.xls')
    log.info(' Writing flag stats to Excel file '+xl_filename)
    xlFile = xlwt.Workbook()
    if platform.system()=="darwin": xlFile.dates_1904 = True
    xlFlagSheet = xlFile.add_sheet('Flag')
    # get the flag statistics
    xlRow = 0
    xlCol = 1
    xlFlagSheet.write(xlRow,xlCol,ds.globalattributes['Flag000'])
    xlFlagSheet.write(xlRow,xlCol+1,ds.globalattributes['Flag001'])
    xlFlagSheet.write(xlRow,xlCol+2,ds.globalattributes['Flag002'])
    xlFlagSheet.write(xlRow,xlCol+3,ds.globalattributes['Flag003'])
    xlFlagSheet.write(xlRow,xlCol+4,ds.globalattributes['Flag004'])
    xlFlagSheet.write(xlRow,xlCol+5,ds.globalattributes['Flag005'])
    xlFlagSheet.write(xlRow,xlCol+6,ds.globalattributes['Flag006'])
    xlFlagSheet.write(xlRow,xlCol+7,ds.globalattributes['Flag007'])
    xlFlagSheet.write(xlRow,xlCol+8,ds.globalattributes['Flag008'])
    xlFlagSheet.write(xlRow,xlCol+9,ds.globalattributes['Flag009'])
    xlFlagSheet.write(xlRow,xlCol+10,ds.globalattributes['Flag010'])
    xlFlagSheet.write(xlRow,xlCol+11,ds.globalattributes['Flag011'])
    xlFlagSheet.write(xlRow,xlCol+12,ds.globalattributes['Flag012'])
    xlFlagSheet.write(xlRow,xlCol+13,ds.globalattributes['Flag013'])
    xlFlagSheet.write(xlRow,xlCol+14,ds.globalattributes['Flag014'])
    xlFlagSheet.write(xlRow,xlCol+15,ds.globalattributes['Flag015'])
    xlFlagSheet.write(xlRow,xlCol+16,ds.globalattributes['Flag016'])
    xlFlagSheet.write(xlRow,xlCol+17,ds.globalattributes['Flag017'])
    xlFlagSheet.write(xlRow,xlCol+18,ds.globalattributes['Flag018'])
    xlFlagSheet.write(xlRow,xlCol+19,ds.globalattributes['Flag019'])
    xlFlagSheet.write(xlRow,xlCol+20,ds.globalattributes['Flag021'])
    xlFlagSheet.write(xlRow,xlCol+21,ds.globalattributes['Flag022'])
    xlFlagSheet.write(xlRow,xlCol+22,ds.globalattributes['Flag030'])
    xlFlagSheet.write(xlRow,xlCol+23,ds.globalattributes['Flag031'])
    xlFlagSheet.write(xlRow,xlCol+24,ds.globalattributes['Flag040'])
    xlFlagSheet.write(xlRow,xlCol+25,ds.globalattributes['Flag050'])
    xlFlagSheet.write(xlRow,xlCol+26,ds.globalattributes['Flag060'])
    xlFlagSheet.write(xlRow,xlCol+27,ds.globalattributes['Flag070'])
    xlFlagSheet.write(xlRow,xlCol+28,ds.globalattributes['Flag080'])
    xlFlagSheet.write(xlRow,xlCol+29,ds.globalattributes['Flag081'])
    xlFlagSheet.write(xlRow,xlCol+30,ds.globalattributes['Flag082'])
    xlFlagSheet.write(xlRow,xlCol+31,ds.globalattributes['Flag083'])
    xlFlagSheet.write(xlRow,xlCol+32,ds.globalattributes['Flag084'])
    xlFlagSheet.write(xlRow,xlCol+33,ds.globalattributes['Flag085'])
    xlFlagSheet.write(xlRow,xlCol+34,ds.globalattributes['Flag086'])
    xlFlagSheet.write(xlRow,xlCol+35,ds.globalattributes['Flag087'])
    xlFlagSheet.write(xlRow,xlCol+36,ds.globalattributes['Flag090'])
    xlFlagSheet.write(xlRow,xlCol+37,ds.globalattributes['Flag100'])
    xlFlagSheet.write(xlRow,xlCol+38,ds.globalattributes['Flag110'])
    xlFlagSheet.write(xlRow,xlCol+39,ds.globalattributes['Flag120'])
    xlFlagSheet.write(xlRow,xlCol+40,ds.globalattributes['Flag130'])
    xlFlagSheet.write(xlRow,xlCol+41,ds.globalattributes['Flag140'])
    xlFlagSheet.write(xlRow,xlCol+42,ds.globalattributes['Flag150'])
    xlFlagSheet.write(xlRow,xlCol+43,ds.globalattributes['Flag151'])
    xlFlagSheet.write(xlRow,xlCol+44,ds.globalattributes['Flag161'])
    xlFlagSheet.write(xlRow,xlCol+45,ds.globalattributes['Flag162'])
    xlFlagSheet.write(xlRow,xlCol+46,ds.globalattributes['Flag171'])
    xlFlagSheet.write(xlRow,xlCol+47,ds.globalattributes['Flag172'])
    xlFlagSheet.write(xlRow,xlCol+48,ds.globalattributes['Flag173'])
    xlFlagSheet.write(xlRow,xlCol+49,ds.globalattributes['Flag174'])
    xlFlagSheet.write(xlRow,xlCol+50,ds.globalattributes['Flag180'])
    xlFlagSheet.write(xlRow,xlCol+51,ds.globalattributes['Flag190'])
    xlFlagSheet.write(xlRow,xlCol+52,ds.globalattributes['Flag191'])
    xlFlagSheet.write(xlRow,xlCol+53,ds.globalattributes['Flag200'])
    xlFlagSheet.write(xlRow,xlCol+54,ds.globalattributes['Flag201'])
    xlFlagSheet.write(xlRow,xlCol+55,ds.globalattributes['Flag211'])
    bins0 = numpy.array([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,21,22,30,31,40,50,60,70,80,81,82,83,84,85,86,87,90,100,110,120,130,140,150,151,161,162,171,172,173,174,180,190,191,200,201,211,212])
    bins = bins0 - 0.5
    xlRow = xlRow + 1
    xlCol = 1
    for Value in bins[:len(bins)-1]:
        xlFlagSheet.write(xlRow,xlCol,numpy.int32(Value+0.5))
        xlCol = xlCol + 1
    xlRow = xlRow + 1
    xlCol = 0
    dsVarNames = ds.series.keys()
    dsVarNames.sort(key=unicode.lower)
    for ThisOne in dsVarNames:
        data,flag,attr = qcutils.GetSeries(ds, ThisOne)
        hist, bin_edges = numpy.histogram(flag, bins=bins)
        xlFlagSheet.write(xlRow,xlCol,ThisOne)
        xlCol = xlCol + 1
        for Value in hist:
            xlFlagSheet.write(xlRow,xlCol,numpy.float64(Value))
            xlCol = xlCol + 1
        xlCol = 0
        xlRow = xlRow + 1
    xlFile.save(xl_filename)

def load_controlfile(path=''):
    ''' 
    Returns a control file object.
    USAGE: cf = load_controlfile([path=<some_path_to_a_controlfile>])
    The "path" keyword is optional.
    '''
    name = get_controlfilename(path=path)
    cf = get_controlfilecontents(name)
    return cf

def nc_concatenate(cf):
    # initialise logicals
    TimeGap = False
    # get an instance of the data structure
    ds = DataStructure()
    # get the input file list
    InFile_list = cf['Files']['In'].keys()
    # read in the first file
    ncFileName = cf['Files']['In'][InFile_list[0]]
    log.info('nc_concatenate: Reading data from '+ncFileName)
    ds_n = nc_read_series(ncFileName)
    if len(ds_n.series.keys())==0:
        log.error('nc_concatenate: An error occurred reading netCDF file: '+ncFileName)
        return
    # fill the global attributes
    for ThisOne in ds_n.globalattributes.keys():
        ds.globalattributes[ThisOne] = ds_n.globalattributes[ThisOne]
    # fill the variables
    for ThisOne in ds_n.series.keys():
        ds.series[ThisOne] = {}
        ds.series[ThisOne]['Data'] = ds_n.series[ThisOne]['Data']
        if 'Flag' in ds_n.series[ThisOne].keys():
            ds.series[ThisOne]['Flag'] = ds_n.series[ThisOne]['Flag']
        ds.series[ThisOne]['Attr'] = {}
        for attr in ds_n.series[ThisOne]['Attr'].keys():
            ds.series[ThisOne]['Attr'][attr] = ds_n.series[ThisOne]['Attr'][attr]
    ts = numpy.int32(ds.globalattributes['time_step'])
    # loop over the remaining files given in the control file
    for n in InFile_list[1:]:
        ncFileName = cf['Files']['In'][InFile_list[numpy.int32(n)]]
        log.info('nc_concatenate: Reading data from '+ncFileName)
        #print 'ncconcat: reading data from '+ncFileName
        ds_n = nc_read_series(ncFileName)
        if len(ds.series.keys())==0:
            log.error('nc_concatenate: An error occurred reading the netCDF file: '+ncFileName)
            return
        dt_n = ds_n.series['DateTime']['Data']
        dt = ds.series['DateTime']['Data']
        nRecs_n = len(ds_n.series['xlDateTime']['Data'])
        nRecs = len(ds.series['xlDateTime']['Data'])
        #print ds.series['DateTime']['Data'][-1],ds_n.series['DateTime']['Data'][-1]
        #print dt[-1],dt[-1]+datetime.timedelta(minutes=ts),dt_n[0]
        #if dt_n[0]<dt[-1]+datetime.timedelta(minutes=ts):
        #    log.info('nc_concatenate: Overlapping times detected in consecutive files')
        #    si = qcutils.GetDateIndex(dt_n,str(dt[-1]),ts=ts)+1
        #    ei = -1
        #if dt_n[0]==dt[-1]+datetime.timedelta(minutes=ts):
        #    log.info('nc_concatenate: Start and end times OK in consecutive files')
        #    si = 0; ei = -1
        #if dt_n[0]>dt[-1]+datetime.timedelta(minutes=ts):
        #    log.info('nc_concatenate: Gap between start and end times in consecutive files')
        #    si = 0; ei = -1
        #    TimeGap = True
        # loop over the data series in the concatenated file
        si = 0; ei = -1
        for ThisOne in ds.series.keys():
            # does this series exist in the file being added to the concatenated file
            if ThisOne in ds_n.series.keys():
                # if so, then append this series to the concatenated series
                ds.series[ThisOne]['Data'] = numpy.append(ds.series[ThisOne]['Data'],ds_n.series[ThisOne]['Data'][si:ei])
                ds.series[ThisOne]['Flag'] = numpy.append(ds.series[ThisOne]['Flag'],ds_n.series[ThisOne]['Flag'][si:ei])
            else:
                # if not, then create a dummy series and concatenate that
                ds_n.series[ThisOne] = {}
                ds_n.series[ThisOne]['Data'] = numpy.array([c.missing_value]*nRecs_n,dtype=numpy.float64)
                ds_n.series[ThisOne]['Flag'] = numpy.array([1]*nRecs_n,dtype=numpy.int32)
                ds.series[ThisOne]['Data'] = numpy.append(ds.series[ThisOne]['Data'],ds_n.series[ThisOne]['Data'][si:ei])
                ds.series[ThisOne]['Flag'] = numpy.append(ds.series[ThisOne]['Flag'],ds_n.series[ThisOne]['Flag'][si:ei])
        # and now loop over the series in the file being concatenated
        for ThisOne in ds_n.series.keys():
            # does this series exist in the concatenated data
            if ThisOne not in ds.series.keys():
                # if not then add it
                ds.series[ThisOne] = {}
                ds.series[ThisOne]['Data'] = numpy.array([c.missing_value]*nRecs,dtype=numpy.float64)
                ds.series[ThisOne]['Flag'] = numpy.array([1]*nRecs,dtype=numpy.int32)
                ds.series[ThisOne]['Data'] = numpy.append(ds.series[ThisOne]['Data'],ds_n.series[ThisOne]['Data'][si:ei])
                ds.series[ThisOne]['Flag'] = numpy.append(ds.series[ThisOne]['Flag'],ds_n.series[ThisOne]['Flag'][si:ei])
                ds.series[ThisOne]['Attr'] = {}
                for attr in ds_n.series[ThisOne]['Attr'].keys():
                    ds.series[ThisOne]['Attr'][attr] = ds_n.series[ThisOne]['Attr'][attr]
        # now sort out any time gaps
        if TimeGap:
            qcutils.FixTimeGaps(ds)
            TimeGap = False
    
    ds.globalattributes['nc_nrecs'] = str(len(ds.series['xlDateTime']['Data']))
    
    # write the netCDF file
    ncFileName = get_keyvalue_from_cf(cf['Files']['Out'],'ncFileName')
    log.info('nc_concatenate: Writing data to '+ncFileName)
    ncFile = nc_open_write(ncFileName)
    nc_write_series(ncFile,ds)

def nc_read_series(ncFullName):
    ''' Read a netCDF file and put the data and meta-data into a DataStructure'''
    log.info(' Reading netCDF file '+ncFullName)
    netCDF4.default_encoding = 'latin-1'
    ds = DataStructure()
    # check to see if the requested file exists, return empty ds if it doesn't
    if not qcutils.file_exists(ncFullName,mode="quiet"):
        log.error(' netCDF file '+ncFullName+' not found')
        raise Exception("nc_read_series: file not found")
    # file probably exists, so let's read it
    ncFile = netCDF4.Dataset(ncFullName,'r')
    # now deal with the global attributes
    gattrlist = ncFile.ncattrs()
    if len(gattrlist)!=0:
        for gattr in gattrlist:
            ds.globalattributes[gattr] = getattr(ncFile,gattr)
            if "time_step" in ds.globalattributes: c.ts = ds.globalattributes["time_step"]
    # get a list of the variables in the netCDF file (not their QC flags)
    varlist = [x for x in ncFile.variables.keys() if "_QCFlag" not in x]
    for ThisOne in varlist:
        # skip variables that do not have time as a dimension
        dimlist = [x.lower() for x in ncFile.variables[ThisOne].dimensions]
        if "time" not in dimlist: continue
        # create the series in the data structure
        ds.series[unicode(ThisOne)] = {}
        # get the data and the QC flag
        data,flag,attr = nc_read_var(ncFile,ThisOne)
        ds.series[ThisOne]["Data"] = data
        ds.series[ThisOne]["Flag"] = flag
        ds.series[ThisOne]["Attr"] = attr
    ncFile.close()
    # make sure all values of -9999 have non-zero QC flag
    qcutils.CheckQCFlags(ds)
    # get a series of Python datetime objects
    qcutils.get_datetimefromymdhms(ds)
    # get series of UTC datetime
    qcutils.get_UTCfromlocaltime(ds)
    # tell the user when the data starts and ends
    ldt = ds.series["DateTime"]["Data"]
    msg = " Got data from "+ldt[0].strftime("%Y-%m-%d %H:%M")+" to "+ldt[-1].strftime("%Y-%m-%d %H:%M")
    log.info(msg)
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

def nc_read_var(ncFile,ThisOne):
    """ Reads a variable from a netCDF file and returns the data, the QC flag and the variable
        attribute dictionary.
    """
    # check the number of dimensions
    nDims = len(ncFile.variables[ThisOne].shape)
    if nDims not in [1,3]:
        msg = "nc_read_var: unrecognised number of dimensions ("+str(nDims)
        msg = msg+") for netCDF variable "+ ThisOne
        raise Exception(msg)
    if nDims==1:
        # single dimension
        data = ncFile.variables[ThisOne][:]
        # netCDF4 returns a masked array if the "missing_variable" attribute has been set
        # for the variable, here we trap this and force the array in ds.series to be ndarray
        if numpy.ma.isMA(data): data,dummy = qcutils.MAtoSeries(data)
        # check for a QC flag
        if ThisOne+'_QCFlag' in ncFile.variables.keys():
            # load it from the netCDF file
            flag = ncFile.variables[ThisOne+'_QCFlag'][:]
        else:
            # create an empty flag series if it does not exist
            nRecs = numpy.size(data)
            flag = numpy.zeros(nRecs,dtype=numpy.int32)
    elif nDims==3:
        # 3 dimensions
        data = ncFile.variables[ThisOne][:,0,0]
        # netCDF4 returns a masked array if the "missing_variable" attribute has been set
        # for the variable, here we trap this and force the array in ds.series to be ndarray
        if numpy.ma.isMA(data): data,dummy = qcutils.MAtoSeries(data)
        # check for a QC flag
        if ThisOne+'_QCFlag' in ncFile.variables.keys():
            # load it from the netCDF file
            flag = ncFile.variables[ThisOne+'_QCFlag'][:,0,0]
        else:
            # create an empty flag series if it does not exist
            nRecs = numpy.size(data)
            flag = numpy.zeros(nRecs,dtype=numpy.int32)
    # get the variable attributes
    vattrlist = ncFile.variables[ThisOne].ncattrs()
    attr = {}
    if len(vattrlist)!=0:
        for vattr in vattrlist:
            attr[vattr] = getattr(ncFile.variables[ThisOne],vattr)
    return data,flag,attr

def nc_open_write(ncFullName,nctype='NETCDF4'):
    log.info(' Opening netCDF file '+ncFullName+' for writing')
    try:
        ncFile = netCDF4.Dataset(ncFullName,'w',format=nctype)
    except:
        log.error(' Unable to open netCDF file '+ncFullName+' for writing')
        ncFile = ''
    return ncFile

def nc_write_series(ncFile,ds,outputlist=None,ndims=3):
    ldt = ds.series["DateTime"]["Data"]
    ds.globalattributes['QC_version'] = str(cfg.version_name)+' '+str(cfg.version_number)
    for ThisOne in ds.globalattributes.keys():
        setattr(ncFile,ThisOne,ds.globalattributes[ThisOne])
    t = time.localtime()
    rundatetime = str(datetime.datetime(t[0],t[1],t[2],t[3],t[4],t[5]))
    setattr(ncFile,'nc_rundatetime',rundatetime)
    # we specify the size of the Time dimension because netCDF4 is slow to write files
    # when the Time dimension is unlimited
    nRecs = numpy.int32(ds.globalattributes['nc_nrecs'])
    ncFile.createDimension("time",nRecs)
    if ndims==3:
        ncFile.createDimension("latitude",1)
        ncFile.createDimension("longitude",1)
        dims = ("time","latitude","longitude")
    else:
        dims = ("time",)
    if outputlist==None:
        outputlist = ds.series.keys()
    else:
        for ThisOne in outputlist:
            if ThisOne not in ds.series.keys():
                log.warn(' Requested series '+ThisOne+' not found in data structure')
                outputlist.remove(ThisOne)
        if len(outputlist)==0: outputlist = ds.series.keys()
    # can't write an array of Python datetime objects to a netCDF file
    # actually, this could be written as characters
    for ThisOne in ["DateTime","DateTime_UTC"]:
        if ThisOne in outputlist: outputlist.remove(ThisOne)
    # write the time variable
    if "time" not in outputlist:
        nc_time = netCDF4.date2num(ldt,"days since 1800-01-01 00:00:00.0",calendar="gregorian")
        ncVar = ncFile.createVariable("time","d",("time",))
        ncVar[:] = nc_time
        setattr(ncVar,"long_name","time")
        setattr(ncVar,"standard_name","time")
        setattr(ncVar,"units","days since 1800-01-01 00:00:00.0")
        setattr(ncVar,"calendar","gregorian")
    # now write the latitude and longitude variables
    if ndims==3:
        if "latitude" not in outputlist:
            ncVar = ncFile.createVariable("latitude","d",("latitude",))
            ncVar[:] = numpy.float64(ds.globalattributes["latitude"])
            setattr(ncVar,'long_name','latitude')
            setattr(ncVar,'standard_name','latitude')
            setattr(ncVar,'units','degrees north')
        if "longitude" not in outputlist:
            ncVar = ncFile.createVariable("longitude","d",("longitude",))
            ncVar[:] = numpy.float64(ds.globalattributes["longitude"])
            setattr(ncVar,'long_name','longitude')
            setattr(ncVar,'standard_name','longitude')
            setattr(ncVar,'units','degrees east')
    # now make sure the date and time series are in outputlist
    datetimelist = ['xlDateTime','xlDateTime_UTC','Year','Month','Day','Hour','Minute','Second','Hdh']
    # and write them to the netCDF file
    for ThisOne in sorted(datetimelist):
        if ThisOne in ds.series.keys(): nc_write_var(ncFile,ds,ThisOne,dims)
        if ThisOne in outputlist: outputlist.remove(ThisOne)
    # write everything else to the netCDF file
    for ThisOne in sorted(outputlist):
        #pdb.set_trace()
        nc_write_var(ncFile,ds,ThisOne,dims)
    # write the coordinate reference system (crs) variable
    if "crs" not in outputlist:
        ncVar = ncFile.createVariable("crs","i",())
        setattr(ncVar,"grid_mapping_name","latitude_longitude")
        setattr(ncVar,"long_name","WGS 1984 datum")
        setattr(ncVar,"longitude_of_prime_meridian","0.0")
        setattr(ncVar,"semi_major_axis","6378137.0")
        setattr(ncVar,"inverse_flattening","298.257223563")
    ncFile.close()

def nc_write_var(ncFile,ds,ThisOne,dim):
    """
    PURPOSE:
     Function to write data from a series in the data structure to a netCDF variable.
    USAGE:
     nc_write_var(ncFile,ds,ThisOne,("time","latitude","longitude"))
      where ncFile is a netCDF file object
            ds is the data structure
            ThisOne is the label of a series in ds
            ("time","latitude","longitude") is the dimension tuple
    AUTHOR: PRI
    DATE: August 2014
    """
    # get the data type of the series in ds
    dt = get_ncdtype(ds.series[ThisOne]['Data'])
    # create the netCDF variable
    if ThisOne not in ncFile.variables.keys():
        ncVar = ncFile.createVariable(ThisOne,dt,dim)
    else:
        log.warn('   write variable to netCDF file:  '+str(ThisOne)+' duplicated in controlfile')
        return
    # different writes to the variable depending on whether it is 1D or 3D
    if len(dim)==1: ncVar[:] = ds.series[ThisOne]['Data'].tolist()
    if len(dim)==3: ncVar[:,0,0] = ds.series[ThisOne]['Data'].tolist()
    # write the attributes
    for attr in ds.series[ThisOne]['Attr']:
        setattr(ncVar,attr,ds.series[ThisOne]['Attr'][attr])
    # get the data type of the QC flag
    dt = get_ncdtype(ds.series[ThisOne]['Flag'])
    # create the variable
    ncVar = ncFile.createVariable(ThisOne+'_QCFlag',dt,dim)
    # write 1D or 3D
    if len(dim)==1: ncVar[:] = ds.series[ThisOne]['Flag'].tolist()
    if len(dim)==3: ncVar[:,0,0] = ds.series[ThisOne]['Flag'].tolist()
    # set the attributes
    setattr(ncVar,'long_name',ThisOne+'QC flag')
    setattr(ncVar,'units','none')

def xl_read_flags(cf,ds,level,VariablesInFile):
    # First data row in Excel worksheets.
    FirstDataRow = numpy.int32(cf['Files'][level]['in_firstdatarow']) - 1
    # Get the full name of the Excel file from the control file.
    xlFullName = get_infilename_from_cf(cf,level)
    # Get the Excel workbook object.
    if os.path.isfile(xlFullName):
        xlBook = xlrd.open_workbook(xlFullName)
    else:
        log.error(' Excel file '+xlFullName+' not found, choose another')
        xlFullName = get_filename_dialog(path='.',title='Choose an Excel file')
        if len(xlFullName)==0:
            return
        xlBook = xlrd.open_workbook(xlFullName)
    
    for ThisOne in VariablesInFile:
        if 'xl' in cf['Variables'][ThisOne].keys():
            log.info(' Getting flags for '+ThisOne+' from spreadsheet')
            ActiveSheet = xlBook.sheet_by_name('Flag')
            LastDataRow = numpy.int32(ActiveSheet.nrows)
            HeaderRow = ActiveSheet.row_values(numpy.int32(cf['Files'][level]['in_headerrow'])-1)
            if cf['Variables'][ThisOne]['xl']['name'] in HeaderRow:
                xlCol = HeaderRow.index(cf['Variables'][ThisOne]['xl']['name'])
                Values = ActiveSheet.col_values(xlCol)[FirstDataRow:LastDataRow]
                ds.series[ThisOne]['Flag'] = numpy.array([c.missing_value]*len(Values),numpy.int32)
                for i in range(len(Values)):
                    ds.series[ThisOne]['Flag'][i] = numpy.int32(Values[i])
    
    return ds

def xl_read_series(cf,InLevel):
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

def xl_write_ACCESSStats(ds):
    if "access" not in dir(ds): return
    # open an Excel file for the fit statistics
    cfname = ds.globalattributes["controlfile_name"]
    cf = get_controlfilecontents(cfname,mode="quiet")
    out_filename = get_outfilename_from_cf(cf)
    # get the Excel file name
    xl_filename = out_filename.replace('.nc','_ACCESSStats.xls')
    log.info(' Writing ACCESS fit statistics to Excel file '+xl_filename)
    # open the Excel file
    xlfile = xlwt.Workbook()
    # list of outputs to write to the Excel file
    date_list = ["startdate","enddate"]
    # loop over the series that have been gap filled using ACCESS data
    d_xf = xlwt.easyxf(num_format_str='dd/mm/yyyy hh:mm')
    for label in ds.access.keys():
        # get the list of values to output with the start and end dates removed
        output_list = ds.access[label]["results"].keys()
        for item in date_list:
            if item in output_list: output_list.remove(item)
        # add a sheet with the series label
        xlResultsSheet = xlfile.add_sheet(label)
        xlRow = 9
        xlCol = 0
        for dt in date_list:
            xlResultsSheet.write(xlRow,xlCol,dt)
            for item in ds.access[label]["results"][dt]:
                xlRow = xlRow + 1
                xlResultsSheet.write(xlRow,xlCol,item,d_xf)
            xlRow = 9
            xlCol = xlCol + 1
        for output in output_list:
            xlResultsSheet.write(xlRow,xlCol,output)
            for item in ds.access[label]["results"][output]:
                xlRow = xlRow + 1
                # xlwt under Anaconda seems to only allow float64!
                xlResultsSheet.write(xlRow,xlCol,numpy.float64(item))
            xlRow = 9
            xlCol = xlCol + 1
    xlfile.save(xl_filename)

def xl_write_SOLOStats(ds):
    if "solo" not in dir(ds): return
    # open an Excel file for the fit statistics
    cfname = ds.globalattributes["controlfile_name"]
    cf = get_controlfilecontents(cfname)
    out_filename = get_outfilename_from_cf(cf)
    # get the Excel file name
    xl_filename = out_filename.replace('.nc','_SOLOStats.xls')
    log.info(' Writing SOLO fit statistics to Excel file '+xl_filename)
    # open the Excel file
    xlfile = xlwt.Workbook()
    # list of outputs to write to the Excel file
    date_list = ["startdate","enddate"]
    output_list = ["n","r_max","bias","rmse","var_obs","var_mod","m_ols","b_ols"]
    # loop over the series that have been gap filled using ACCESS data
    d_xf = xlwt.easyxf(num_format_str='dd/mm/yyyy hh:mm')
    for label in ds.solo.keys():
        # get the list of values to output with the start and end dates removed
        output_list = ds.solo[label]["results"].keys()
        for item in date_list:
            if item in output_list: output_list.remove(item)
        # add a sheet with the series label
        xlResultsSheet = xlfile.add_sheet(label)
        xlRow = 10
        xlCol = 0
        for dt in date_list:
            xlResultsSheet.write(xlRow,xlCol,dt)
            for item in ds.solo[label]["results"][dt]:
                xlRow = xlRow + 1
                xlResultsSheet.write(xlRow,xlCol,item,d_xf)
            xlRow = 10
            xlCol = xlCol + 1
        for output in output_list:
            xlResultsSheet.write(xlRow,xlCol,output)
            for item in ds.solo[label]["results"][output]:
                xlRow = xlRow + 1
                # xlwt under Anaconda seems to only allow float64!
                xlResultsSheet.write(xlRow,xlCol,numpy.float64(item))
            xlRow = 10
            xlCol = xlCol + 1
    xlfile.save(xl_filename)

def xl_write_series(ds, xlfullname, outputlist=None):
    if "nc_nrecs" in ds.globalattributes.keys():
        nRecs = numpy.int32(ds.globalattributes["nc_nrecs"])
    else:
        variablelist = ds.series.keys()
        nRecs = len(ds.series[variablelist[0]]["Data"])
    # open the Excel file
    log.info(' Opening and writing Excel file '+xlfullname)
    xlfile = xlwt.Workbook(encoding="latin-1")
    # set the datemode
    xlfile.dates_1904 = numpy.int32(ds.globalattributes['xl_datemode'])
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
    for ThisOne in sorted([x for x in globalattrlist if 'Flag' not in x]):
        xlAttrSheet.write(xlrow,xlcol,ThisOne)
        try:
            xlAttrSheet.write(xlrow,xlcol+1,ds.globalattributes[ThisOne])
        except:
            xlAttrSheet.write(xlrow,xlcol+1,int(ds.globalattributes[ThisOne]))
        xlrow = xlrow + 1
    for ThisOne in sorted([x for x in globalattrlist if 'Flag' in x]):
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
    if outputlist==None:
        outputlist = variablelist
    else:
        for ThisOne in outputlist:
            if ThisOne not in variablelist:
                log.warn(' Requested series '+ThisOne+' not found in data structure')
                outputlist.remove(ThisOne)
        if len(outputlist)==0:
            outputlist = variablelist
    outputlist.sort()
    for ThisOne in ["DateTime","DateTime_UTC"]:
        if ThisOne in outputlist: outputlist.remove(ThisOne)
    for ThisOne in outputlist:
        xlAttrSheet.write(xlrow,xlcol_varname,ThisOne)
        attributelist = ds.series[ThisOne]['Attr'].keys()
        attributelist.sort()
        for Attr in attributelist:
            xlAttrSheet.write(xlrow,xlcol_attrname,Attr)
            xlAttrSheet.write(xlrow,xlcol_attrvalue,str(ds.series[ThisOne]['Attr'][Attr]))
            xlrow = xlrow + 1
    # write the Excel date/time to the data and the QC flags as the first column
    ldt = ds.series["DateTime"]["Data"]
    # get the datemode of the original spreadsheet
    datemode = numpy.int32(ds.globalattributes['xl_datemode'])
    xlDateTime = qcutils.get_xldate_from_datetime(ldt,datemode=datemode)
    log.info(' Writing the datetime to Excel file '+xlfullname)
    d_xf = xlwt.easyxf(num_format_str='dd/mm/yyyy hh:mm')
    xlDataSheet.write(2,xlcol,'xlDateTime')
    for j in range(nRecs):
        xlDataSheet.write(j+3,xlcol,xlDateTime[j],d_xf)
        xlFlagSheet.write(j+3,xlcol,xlDateTime[j],d_xf)
    # output the xl datetime as UTC if it exists in the file
    #if "xlDateTime_UTC" in ds.series.keys():
        #xlcol = xlcol + 1
        #xlDateTime = ds.series["xlDateTime_UTC"]["Data"]
        #xlDataSheet.write(2,xlcol,"xlDateTime_UTC")
        #for j in range(nRecs):
            #xlDataSheet.write(j+3,xlcol,xlDateTime[j],d_xf)
            #xlFlagSheet.write(j+3,xlcol,xlDateTime[j],d_xf)
    # remove xlDateTime from the list of variables to be written to the Excel file
    if "xlDateTime" in outputlist: outputlist.remove("xlDateTime")
    if "xlDateTime_UTC" in outputlist: outputlist.remove("xlDateTime_UTC")
    # now start looping over the other variables in the xl file
    xlcol = xlcol + 1
    # loop over variables to be output to xl file
    for ThisOne in outputlist:
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
            xlDataSheet.write(j+3,xlcol,numpy.float64(ds.series[ThisOne]['Data'][j]))
        # check to see if this variable has a quality control flag
        if 'Flag' in ds.series[ThisOne].keys():
            # write the QC flag name to the xk file
            xlFlagSheet.write(2,xlcol,ThisOne)
            # specify the format of the QC flag (integer)
            d_xf = xlwt.easyxf(num_format_str='0')
            # loop over QV flag values and write to xl file
            for j in range(nRecs):
                xlFlagSheet.write(j+3,xlcol,numpy.int32(ds.series[ThisOne]['Flag'][j]),d_xf)
        # increment the column pointer
        xlcol = xlcol + 1
    
    log.info(' Saving the Excel file '+xlfullname)
    xlfile.save(xlfullname)

def xlsx_write_series(ds, xlsxfullname, outputlist=None):
    if "nc_nrecs" in ds.globalattributes.keys():
        nRecs = numpy.int32(ds.globalattributes["nc_nrecs"])
    else:
        variablelist = ds.series.keys()
        nRecs = len(ds.series[variablelist[0]]["Data"])
    # open the Excel file
    log.info(' Opening and writing Excel file '+xlsxfullname)
    if int(ds.globalattributes["xl_datemode"])==1:
        xlfile = xlsxwriter.Workbook(xlsxfullname, {'date_1904': True})
    else:
        xlfile = xlsxwriter.Workbook(xlsxfullname, {'date_1904': False})
    # set the datemode
    xlfile.dates_1904 = numpy.int32(ds.globalattributes['xl_datemode'])
    # add sheets to the Excel file
    xlAttrSheet = xlfile.add_worksheet('Attr')
    xlDataSheet = xlfile.add_worksheet('Data')
    xlFlagSheet = xlfile.add_worksheet('Flag')
    # write the global attributes
    log.info(' Writing the global attributes to Excel file '+xlsxfullname)
    xlcol = 0
    xlrow = 0
    xlAttrSheet.write(xlrow,xlcol,'Global attributes')
    xlrow = xlrow + 1
    globalattrlist = ds.globalattributes.keys()
    globalattrlist.sort()
    for ThisOne in sorted([x for x in globalattrlist if 'Flag' not in x]):
        xlAttrSheet.write(xlrow,xlcol,ThisOne)
        xlAttrSheet.write(xlrow,xlcol+1,ds.globalattributes[ThisOne])
        xlrow = xlrow + 1
    for ThisOne in sorted([x for x in globalattrlist if 'Flag' in x]):
        xlAttrSheet.write(xlrow,xlcol,ThisOne)
        xlAttrSheet.write(xlrow,xlcol+1,str(ds.globalattributes[ThisOne]))
        xlrow = xlrow + 1
    # write the variable attributes
    log.info(' Writing the variable attributes to Excel file '+xlsxfullname)
    xlrow = xlrow + 1
    xlAttrSheet.write(xlrow,xlcol,'Variable attributes')
    xlrow = xlrow + 1
    xlcol_varname = 0
    xlcol_attrname = 1
    xlcol_attrvalue = 2
    variablelist = ds.series.keys()
    if outputlist==None:
        outputlist = variablelist
    else:
        for ThisOne in outputlist:
            if ThisOne not in variablelist:
                log.warn(' Requested series '+ThisOne+' not found in data structure')
                outputlist.remove(ThisOne)
        if len(outputlist)==0:
            outputlist = variablelist
    outputlist.sort()
    for ThisOne in ["DateTime","DateTime_UTC"]:
        if ThisOne in outputlist: outputlist.remove(ThisOne)
    for ThisOne in outputlist:
        xlAttrSheet.write(xlrow,xlcol_varname,ThisOne)
        attributelist = ds.series[ThisOne]['Attr'].keys()
        attributelist.sort()
        for Attr in attributelist:
            xlAttrSheet.write(xlrow,xlcol_attrname,Attr)
            xlAttrSheet.write(xlrow,xlcol_attrvalue,str(ds.series[ThisOne]['Attr'][Attr]))
            xlrow = xlrow + 1
    # write the Excel date/time to the data and the QC flags as the first column
    datemode = numpy.int32(ds.globalattributes['xl_datemode'])
    ldt = ds.series["DateTime"]["Data"]
    xlDateTime = qcutils.get_xldate_from_datetime(ldt,datemode=datemode)
    log.info(' Writing the datetime to Excel file '+xlsxfullname)
    dt_format = xlfile.add_format({'num_format': 'dd/mm/yyyy hh:mm'})
    xlDataSheet.write(2,xlcol,'xlDateTime')
    xlFlagSheet.write(2,xlcol,'xlDateTime')
    for j in range(nRecs):
        xlDataSheet.write_datetime(j+3,xlcol,ldt[j],dt_format)
        xlFlagSheet.write_datetime(j+3,xlcol,ldt[j],dt_format)
    # remove xlDateTime from the list of variables to be written to the Excel file
    if "xlDateTime" in outputlist: outputlist.remove("xlDateTime")
    if "xlDateTime_UTC" in outputlist: outputlist.remove("xlDateTime_UTC")
    # now start looping over the other variables in the xl file
    xlcol = xlcol + 1
    # loop over variables to be output to xl file
    for ThisOne in outputlist:
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
            xlDataSheet.write(j+3,xlcol,numpy.float64(ds.series[ThisOne]['Data'][j]))
        # check to see if this variable has a quality control flag
        if 'Flag' in ds.series[ThisOne].keys():
            # write the QC flag name to the Excel file
            xlFlagSheet.write(2,xlcol,ThisOne)
            # specify the format of the QC flag (integer)
            flag_format = xlfile.add_format({'num_format': '0'})
            # loop over QC flag values and write to xl file
            for j in range(nRecs):
                xlFlagSheet.write(j+3,xlcol,numpy.int32(ds.series[ThisOne]['Flag'][j]),flag_format)
        # increment the column pointer
        xlcol = xlcol + 1
    
    log.info(' Saving the Excel file '+xlsxfullname)
    xlfile.close()

def xl_write_series_unsorted(cf, ds, xlfullname, outputlist=None):
    if "nc_nrecs" in ds.globalattributes.keys():
        nRecs = numpy.int32(ds.globalattributes["nc_nrecs"])
    else:
        variablelist = ds.series.keys()
        nRecs = len(ds.series[variablelist[0]]["Data"])
    # open the Excel file
    log.info(' Opening and writing Excel file '+xlfullname)
    xlfile = xlwt.Workbook(encoding="latin-1")
    # set the datemode
    xlfile.dates_1904 = numpy.int32(ds.globalattributes['xl_datemode'])
    # add sheets to the Excel file
    xlDataSheet = xlfile.add_sheet('Data')
    xlFlagSheet = xlfile.add_sheet('Flag')
    xlCol = 0
    variablelist = ds.series.keys()
    if outputlist==None:
        outputlist = variablelist
    else:
        for ThisOne in outputlist:
            if ThisOne not in variablelist:
                log.warn(' Requested series '+ThisOne+' not found in data structure')
                outputlist.remove(ThisOne)
        if len(outputlist)==0:
            outputlist = variablelist
    # write the xl date/time value to the first column of the worksheets
    ldt = ds.series["DateTime"]["Data"]
    # get the datemode of the original spreadsheet
    datemode = numpy.int32(ds.globalattributes['xl_datemode'])
    xlDateTime = qcutils.get_xldate_from_datetime(ldt,datemode=datemode)
    log.info(' Writing the datetime to Excel file '+xlfullname)
    d_xf = xlwt.easyxf(num_format_str='dd/mm/yyyy hh:mm')
    for j in range(nRecs):
        xlDataSheet.write(j+10,xlCol,xlDateTime[j],d_xf)
    xlFlagSheet.write(9,xlCol,'TIMESTAMP')
    xlDataSheet.write(9,xlCol,'TIMESTAMP')
    xlDataSheet.write(0,xlCol+1,'Dataset not licensed for distribution.  For internal use only.')
    xlDataSheet.write(1,xlCol,'Site:')
    try:
        xlDataSheet.write(1,xlCol+1,ds.globalattributes['site'])
    except:
        xlDataSheet.write(1,xlCol+1,ds.globalattributes['site_name'])
    xlDataSheet.write(3,xlCol,'Institution:')
    xlDataSheet.write(3,xlCol+1,ds.globalattributes['institution'])
    xlDataSheet.write(2,xlCol,'Latitude:')
    xlDataSheet.write(2,xlCol+1,ds.globalattributes['latitude'])
    xlDataSheet.write(2,xlCol+2,'Longitude:')
    xlDataSheet.write(2,xlCol+3,ds.globalattributes['longitude'])
    xlDataSheet.write(4,xlCol,'Contact:')
    xlDataSheet.write(4,xlCol+1,ds.globalattributes['contact'])
    xlFlagSheet.write(0,xlCol,'0:')
    if ((qcutils.cfkeycheck(cf,Base='Output',ThisOne='FlagList')) and (cf['Output']['FlagList'] == 'False')):
        xlFlagSheet.write(0,xlCol+1,'FlagList disabled')
    else:
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
        xlFlagSheet.write(j+10,xlCol,xlDateTime[j],d_xf)
    # remove the date and time variables from the list to output
    for ThisOne in ['xlDateTime','Year','Month','Day','Hour','Minute','Second','Hdh']:
        if ThisOne in variablelist:
            variablelist.remove(ThisOne)
    # now start looping over the other variables in the xl file
    xlCol = xlCol + 1
    # list of variables in the data structure
    #VariablesInFile.sort(key=str.lower)
    # list of variables to write out (specified in the control file)
    #VariablesToOutput.sort(key=str.lower)
    # loop over variables to be output to xl file
    for ThisOne in outputlist:
        if ThisOne in variablelist:
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
                xlDataSheet.write(j+10,xlCol,numpy.float64(ds.series[ThisOne]['Data'][j]))
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
    log.info(' Saving the Excel file '+xlfullname)
    # save the xl file
    xlfile.save(xlfullname)

def xlsx_write_series_unsorted(cf, ds, xlsxfullname, outputlist=None):
    if "nc_nrecs" in ds.globalattributes.keys():
        nRecs = numpy.int32(ds.globalattributes["nc_nrecs"])
    else:
        variablelist = ds.series.keys()
        nRecs = len(ds.series[variablelist[0]]["Data"])
    # open the Excel file
    log.info(' Opening and writing Excel file '+xlsxfullname)
    if int(ds.globalattributes["xl_datemode"])==1:
        xlfile = xlsxwriter.Workbook(xlsxfullname, {'date_1904': True})
    else:
        xlfile = xlsxwriter.Workbook(xlsxfullname, {'date_1904': False})
    # set the datemode
    xlfile.dates_1904 = numpy.int32(ds.globalattributes['xl_datemode'])
    # add sheets to the Excel file
    xlDataSheet = xlfile.add_worksheet('Data')
    xlFlagSheet = xlfile.add_worksheet('Flag')
    xlcol = 0
    variablelist = ds.series.keys()
    if outputlist==None:
        outputlist = variablelist
    else:
        for ThisOne in outputlist:
            if ThisOne not in variablelist:
                log.warn(' Requested series '+ThisOne+' not found in data structure')
                outputlist.remove(ThisOne)
        if len(outputlist)==0:
            outputlist = variablelist
    for ThisOne in ["DateTime","DateTime_UTC"]:
        if ThisOne in outputlist: outputlist.remove(ThisOne)
    # write the Excel date/time to the data and the QC flags as the first column
    datemode = numpy.int32(ds.globalattributes['xl_datemode'])
    ldt = ds.series["DateTime"]["Data"]
    xlDateTime = qcutils.get_xldate_from_datetime(ldt,datemode=datemode)
    log.info(' Writing the datetime to Excel file '+xlsxfullname)
    dt_format = xlfile.add_format({'num_format': 'dd/mm/yyyy hh:mm'})
    xlDataSheet.write(9,xlcol,'xlDateTime')
    xlFlagSheet.write(9,xlcol,'xlDateTime')
    for j in range(nRecs):
        xlDataSheet.write_datetime(j+10,xlcol,ldt[j],dt_format)
        xlFlagSheet.write_datetime(j+10,xlcol,ldt[j],dt_format)
    
    xlDataSheet.write(0,xlcol+1,'Dataset not licensed for distribution.  For internal use only.')
    xlDataSheet.write(1,xlcol,'Site:')
    try:
        xlDataSheet.write(1,xlcol+1,ds.globalattributes['site'])
    except:
        xlDataSheet.write(1,xlcol+1,ds.globalattributes['site_name'])
    xlDataSheet.write(3,xlcol,'Institution:')
    xlDataSheet.write(3,xlcol+1,ds.globalattributes['institution'])
    xlDataSheet.write(2,xlcol,'Latitude:')
    xlDataSheet.write(2,xlcol+1,ds.globalattributes['latitude'])
    xlDataSheet.write(2,xlcol+2,'Longitude:')
    xlDataSheet.write(2,xlcol+3,ds.globalattributes['longitude'])
    xlDataSheet.write(4,xlcol,'Contact:')
    xlDataSheet.write(4,xlcol+1,ds.globalattributes['contact'])
    xlFlagSheet.write(0,xlcol,'0:')
    if ((qcutils.cfkeycheck(cf,Base='Output',ThisOne='FlagList')) and (cf['Output']['FlagList'] == 'False')):
        xlFlagSheet.write(0,xlcol+1,'FlagList disabled')
    else:
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
    
    # remove the date and time variables from the list to output
    for ThisOne in ['xlDateTime','Year','Month','Day','Hour','Minute','Second','Hdh','DateTime_UTC','DateTime','xlDateTime_UTC']:
        if ThisOne in variablelist:
            variablelist.remove(ThisOne)
    # now start looping over the other variables in the xl file
    xlcol = xlcol + 1
    # list of variables in the data structure
    #VariablesInFile.sort(key=str.lower)
    # list of variables to write out (specified in the control file)
    #VariablesToOutput.sort(key=str.lower)
    # loop over variables to be output to xl file
    for ThisOne in outputlist:
        if ThisOne in variablelist:
            # put up a progress message
            log.info(' Writing '+ThisOne+' into column '+str(xlcol)+' of the Excel file')
            # write the units and the variable name to the header rows in the xl file
            longname = ds.series[ThisOne]['Attr']['long_name']
            units = ds.series[ThisOne]['Attr']['units']
            #xlDataSheet.write(8,xlcol,Units,d_xf)
            #xlDataSheet.write(9,xlcol,ThisOne,d_xf)
            xlDataSheet.write(7,xlcol,longname)
            xlDataSheet.write(8,xlcol,units)
            xlDataSheet.write(9,xlcol,ThisOne)
            # loop over the values in the variable series (array writes don't seem to work)
            for j in range(nRecs):
                xlDataSheet.write(j+10,xlcol,numpy.float64(ds.series[ThisOne]['Data'][j]))
            # check to see if this variable has a quality control flag
            if 'Flag' in ds.series[ThisOne].keys():
                # write the QC flag name to the Excel file
                xlFlagSheet.write(9,xlcol,ThisOne)
                # specify the format of the QC flag (integer)
                flag_format = xlfile.add_format({'num_format': '0'})
                # loop over QC flag values and write to xl file
                for j in range(nRecs):
                    xlFlagSheet.write(j+10,xlcol,numpy.int32(ds.series[ThisOne]['Flag'][j]),flag_format)
            # increment the column pointer
            xlcol = xlcol + 1
    # tell the user what we are doing
    log.info(' Saving the Excel file '+xlsxfullname)
    # save the xl file
    xlfile.close()
