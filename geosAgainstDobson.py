#!/usr/bin/env python3
import os, h5py, argparse, glob, math,sys
import numpy as np
from lib.dobson_io import readWoudc
import pytz
import matplotlib
matplotlib.use('Agg')
def go ( a ):
    """
    Main program for this thing.
    Input:
            a: argparse which gives command line arguements to do stuff with.
    Output:
            a plot tagged with experiment and control, along with and hdf5 with stats plotted. 
    """

    # *Sigh* This is not great, but work your way through and figure out which data we're reading 
    dobsonFilesCsv = glob.glob(os.path.join(a.dobson_path,'*.csv'))
    ndobsonCsv = len(dobsonFilesCsv)

    controlAnalysisFiles = []
    experimentAnalysisFiles = []
    # A little bit inefficient, but make two passes. First, make a list of files to dmget off dirac, and grab sonde data. 
    # Second loop, actually do stuff.

    dOb = {}
    dOb['lon'] = []
    dOb['lat'] = []  
    dOb['date'] = []
    dOb['time'] = []
    dOb['O3'] = []   
    for s in dobsonFilesCsv:
        lonDobson, latDobson, datetimeDobson, ozoneDobson = readDobson(s, 'woudc')
        for t in datetimeDobson:
            dateDobson = t.astimezone(pytz.UTC).strftime('%Y%m%d')
            timeDobson = t.astimezone(pytz.UTC).strftime('%H%M') 
            # Do stuff for the experimental run (get idx for the experiment and the control along the way) 
            experimentAnalysisFiles.append( getFileName(a.ops, a.experiment, dateDobson, timeDobson) )

            # now for the control
            controlAnalysisFiles.append( getFileName(a.ops, a.control, dateDobson, timeDobson) )
            dOb['lon'].extend(lonDobson)
            dOb['lat'].extend(latDobson)
            dOb['date'].extend(dateDobson)
            dOb['time'].extend(timeDobson)
            dOb['O3'].extend(ozoneDobson)    
    uniqueE = []
    uniqueC = {}
    for e in experimentAnalysisFiles:
        if e not in uniqueE:
            uniqueE.append(e)
    for c in controlAnalysisFiles:
        if c not in uniqueC:
            uniqueC.append(c)
 
    # do dmget
    controlString = " ".join(uniqueC)
    print("dmget on control analysis files "+controlString) 
    os.system('dmget '+ controlString)
    print('done dmget.')
    experimentString = " ".join(uniqueE)
    print("dmget on experiment analysis files "+experimentString) 
    os.system('dmget '+ experimentString)
    print('done dmget...for good!')
    controlOzone = []
    experimentOzone = []
    #Actually do stuff.
    for i,O3 in enumerate(dObs['O3']):
        
        # use sonde latitude to get x,y from analysis.
        idxLon,idxLat =  getIndexFromAnalysis(experimentAnalysisFiles[i], s['lat'], s['lon'])
        print('Reading Experiment Analysis File: {}'.format(experimentAnalysisFiles[i]))
        h5 = h5py.File(experimentAnalysisFiles[i],'r')
        experimentOzone.append(np.asarray(h5['TO3'][0,idxLat,idxLon]))
        h5.close()   
        #same grid, don't need to interpolate that again...

        print('Reading Control Analysis File: {}'.format(controlAnalysisFiles[i]))
        h5 = h5py.File(controlAnalysisFiles[i],'r')
        controlOzone.append(np.asarray(h5['TO3'][0,idxLat,idxLon]))
        h5.close()   
        
        #ss = updateStats(ss, interpolatedSondeOzone, controlOzone, experimentOzone )
    #ss = finishStats( ss )
    #writeH5(a, ss, press_int)
    plt.plot(np.asarray(expermientOzone),'rx')
    plt.plot(np.asarray(controlOzone),'ko')
    plt.savefig('whir.png')

def initProfileStats (nint):
    """
    Initialize dictionary with number of levels we want to interpolate against. (nint)
    Input:
            nint: Number of pressure levels in the grid
    Output:
            o: dictionary of stats
            o['count_sonde'] :        counts (# of sondes used at each interpolated profile gridpoint)
            o['count_both'] :         counts where there are valid analysis profile, and sonde profile points
            o['av_sonde'] :           average of sonde ozone (using count_both)
            o['av_sonde_standalone']: average of sonde ozone (using count_sonde)
            o['av_ana1']:             average of control analysis profile ozone
            o['av_ana2']:             average of experiment analysis profile ozone
            o['std_ana1']:            rms between control analysis and sondes
            o['std_ana2']:            rms between experiment analysis and sondes 
            o['std_sonde']:           std of the sondes
           
    """
    o = {} 
    o['count_sonde'] = np.zeros(nint)
    o['count_both'] = np.zeros(nint)
    o['av_sonde'] = np.zeros(nint)
    o['av_sonde_standalone'] = np.zeros(nint)
    o['av_ana1']  = np.zeros(nint)
    o['av_ana2']  = np.zeros(nint)
    o['std_ana1']   = np.zeros(nint)
    o['std_ana2']   = np.zeros(nint)
    o['std_sonde'] = np.zeros(nint)
    return o

def readDobson ( s, dobsonType ):
    """
    Given a filename, and type of sonde, read the data into a common set of variables. 
    If you're bored, probably should make this more "classy" by using a factory class
    and avoid if statement, but meh... 
    Input:
            s: path to the file to be read in (sonde)
    Output:
            lon                : longitude of the sonde 
            lat                : latitude of the sonde 
            time               : time datetime
            Dobson unit             : Dobson Unit

    """
    lon = []
    lat = []
    if(dobsonType == 'woudc'): 
        d = readWoudc(s)
        alon = float(d['LOCATION']['Longitude'])
        alat = float(d['LOCATION']['Latitude'])
        time = d['OBSERVATIONS']['Time']
        dobsonO3 = d['OBSERVATIONS']['ColumnO3']
        for dd in dobsonO3:
            lon.append(alon)
            lat.append(alat)
    else:
        sys.exit("error don't know what dobson filetype this is!")
   
    return lon, lat, time, dobsonO3

def timeLookup(dateString, timeString):
    """
    Not 100% sure how this guy works, but looks up nearest analysis time for a given sonde.

    Input:
            dateString: date of sonde YYYYMMDD
            timeString: time of sonde hh:mm
    Output:
            YYYY: Year (integer) for analysis lookup
            MM:   Month    "      "     "        "
            DD:   Day      "      "     "        "
            hh:   Hour     "      "     "        " 
            mm:   Minute   "      "     "        "
    
    """
    YYYY = int(dateString[0:4])
    MM = int(dateString[4:6])
    DD = int(dateString[6:8])
    hh = int(timeString.split(':')[0])
    mm = int(timeString.split(':')[1])
    # following Kris' script to look up analysis time.
    timefrac = (hh + mm/60)/24.0
    ihh = math.floor(timefrac*8.0+ 0.5)
    if(ihh > 7): ihh = 7
    hh = '{:02d}'.format(ihh*3)
    hh = int(hh)
    return YYYY, MM, DD, hh, mm 

def getFileName(ops, experiment, dateSonde, timeSonde):
    """   
    Get the sonde ilename given path to users experiments, experiment, datestring, timetring
    Input:
            ops:        path to experiments (i.e. /archive/$USER) 
            experiment: name of GEOS experiment 
            dateSonde:  date string (YYYYMMDD) of the sonde
            timeSonde:  time string (hh:mm) of the sonde
    Output:
            analysisFile: give the analysis file to read in
    """
    YYYY, MM, DD, hh, mm =  timeLookup(dateSonde, timeSonde)
    analysisFile = os.path.join(ops, experiment,'diag','Y{:04d}'.format(YYYY), 'M{:02d}'.format(MM),\
                   experiment+'.inst3_3d_asm_Np.'+dateSonde+'_'+'{:02d}00z.nc4'.format(hh))
    return analysisFile 

def getIndexFromAnalysis(analysisFile, latSonde, lonSonde):
    """
    Look up index on grid for nearest i,j for a given sonde lat/lon
    Input:
            analysisFile: analysis file to read in
            latSonde: latitude of the sonde 
            lonSonde: longitude of the sonde
    Output:
            idxLat: index associated with latitude in GEOS analysis file
            idxLon: index associated with longitude in GEOS analysis file
    """
    h5 = h5py.File(analysisFile,'r')
    lons = h5['lon']
    lats = h5['lat']
    """
    #Kris' way of finding nearest neighbor.
    idx, = np.where(np.asarray(lons) <= float(lonSonde))
    idxLon = idx[max(idx)]
    if idxLon >= nlon: idxLon = 0
    idx, = np.where(np.asarray(lats) <= float(latSonde))
    idxLat = idx[max(idx)]
    """
    # using python and the internets...
    # https://stackoverflow.com/questions/30873844/identifying-the-nearest-grid-point
    idxLon = np.nanargmin((np.asarray(lons)-float(lonSonde))**2.0)
    idxLat = np.nanargmin((np.asarray(lats)-float(latSonde))**2.0)
    print('Index Lon, Index Lat, grid lon, grid lat, sonde lon, sonde lat')
    print(idxLon, idxLat, lons[idxLon], lats[idxLat], lonSonde, latSonde)
    h5.close()
    return idxLon, idxLat


def updateStats(ss, interpolatedSondeOzone, controlOzone, experimentOzone ):
    """
    Compute and update stats for a given profile. 
    These are really running tallies. 
    They aren't the actual statistics until finishStats() has been run.

    Input:
            ss: stats dictionary
            ss['count_sonde'] :        counts (# of sondes used at each interpolated profile gridpoint)
            ss['count_both'] :         counts where there are valid analysis profile, and sonde profile points
            ss['av_sonde'] :           average of sonde ozone (using count_both)
            ss['av_sonde_standalone']: average of sonde ozone (using count_sonde)
            ss['av_ana1']:             average of control analysis profile ozone
            ss['av_ana2']:             average of experiment analysis profile ozone
            ss['std_ana1']:            rms between control analysis and sondes
            ss['std_ana2']:            rms between experiment analysis and sondes 
            ss['std_sonde']:           std of the sondes
            interpolatedSondeOzone: ozone sonde profile interpolated to grid (mPa)
            controlOzone: ozone profile from control analysis interpolated to grid (mPa)
            experimentOzone: ozone profile from experiment analysis interpolated to grid (mPa)
    Output:
            ss: stats dictionary
 
    """

    # The introduction of potentially bad portions of analysis profile leads to having to break out index
    # into sondes for sonde and "both" for rms error stats.
    
    # using nans may not be the greatest idea, as these lines will give you a warning, but will work.
    idxGoodSonde, = np.where( np.isfinite(interpolatedSondeOzone) & (interpolatedSondeOzone > 0.0) )
    idxGoodBoth, = np.where ( np.isfinite(interpolatedSondeOzone) & (interpolatedSondeOzone > 0.0) &\
                        np.isfinite(controlOzone) & np.isfinite(experimentOzone) & (controlOzone > 0.0) & (experimentOzone > 0.0) ) 
        
    ss['count_sonde'][idxGoodSonde] = ss['count_sonde'][idxGoodSonde] + 1.0 
    ss['count_both'][idxGoodBoth] = ss['count_both'][idxGoodBoth] + 1.0

    ss['av_sonde'][idxGoodBoth] = ss['av_sonde'][idxGoodBoth] + interpolatedSondeOzone[idxGoodBoth]
    ss['av_sonde_standalone'][idxGoodSonde] = ss['av_sonde_standalone'][idxGoodSonde] + interpolatedSondeOzone[idxGoodSonde]
    ss['av_ana1'][idxGoodBoth] = ss['av_ana1'][idxGoodBoth] + controlOzone[idxGoodBoth]
    ss['av_ana2'][idxGoodBoth] = ss['av_ana2'][idxGoodBoth] + experimentOzone[idxGoodBoth]

    #here std is really rms error. and we need Both
    ss['std_ana1'][idxGoodBoth] =  ss['std_ana1'][idxGoodBoth] + (controlOzone[idxGoodBoth] -interpolatedSondeOzone[idxGoodBoth])**2
    ss['std_ana2'][idxGoodBoth] =  ss['std_ana2'][idxGoodBoth] + (experimentOzone[idxGoodBoth] -interpolatedSondeOzone[idxGoodBoth])**2
    # here std is std
    ss['std_sonde'][idxGoodSonde] =  ss['std_sonde'][idxGoodSonde] + (interpolatedSondeOzone[idxGoodSonde])**2
    return ss

def finishStats( ss ):
    """
    Take running tally values from stats dictionary and convert to desired statistics.
    Input:
            ss: stats dictionary
            ss['count_sonde'] :        counts (# of sondes used at each interpolated profile gridpoint)
            ss['count_both'] :         counts where there are valid analysis profile, and sonde profile points
            ss['av_sonde'] :           average of sonde ozone (using count_both)
            ss['av_sonde_standalone']: average of sonde ozone (using count_sonde)
            ss['av_ana1']:             average of control analysis profile ozone
            ss['av_ana2']:             average of experiment analysis profile ozone
            ss['std_ana1']:            rms between control analysis and sondes
            ss['std_ana2']:            rms between experiment analysis and sondes 
            ss['std_sonde']:           std of the sondes
            interpolatedSondeOzone: ozone sonde profile interpolated to grid (mPa)
            controlOzone: ozone profile from control analysis interpolated to grid (mPa)
            experimentOzone: ozone profile from experiment analysis interpolated to grid (mPa)
    Output:
            ss: stats dictionary
 
    """
    diff1 = np.zeros(ss['count_sonde'].shape[0])
    diff2 = np.zeros(ss['count_sonde'].shape[0])
    idxGoodSonde, = np.where( ss['count_sonde']>1 )
    idxGoodBoth, = np.where( ss['count_both'] > 1 )

    ss['av_ana1'][idxGoodBoth] = ss['av_ana1'][idxGoodBoth]/ss['count_both'][idxGoodBoth]
    ss['av_ana2'][idxGoodBoth] = ss['av_ana2'][idxGoodBoth]/ss['count_both'][idxGoodBoth]
    ss['av_sonde'][idxGoodBoth] = ss['av_sonde'][idxGoodBoth]/ss['count_both'][idxGoodBoth]

    ss['av_sonde_standalone'][idxGoodBoth] = ss['av_sonde_standalone'][idxGoodBoth]/ss['count_sonde'][idxGoodBoth]

    diff1[idxGoodBoth] = ss['av_ana1'][idxGoodBoth] - ss['av_sonde'][idxGoodBoth]
    diff2[idxGoodBoth] = ss['av_ana2'][idxGoodBoth] - ss['av_sonde'][idxGoodBoth]

    ss['std_ana1'][idxGoodBoth] =  ss['std_ana1'][idxGoodBoth] - ss['count_both'][idxGoodBoth]*(diff1[idxGoodBoth])**2
    ss['std_ana1'][idxGoodBoth] =  np.sqrt(ss['std_ana1'][idxGoodBoth]/(ss['count_both'][idxGoodBoth]-1) ) 

    ss['std_ana2'][idxGoodBoth] =  ss['std_ana2'][idxGoodBoth] - ss['count_both'][idxGoodBoth]*(diff2[idxGoodBoth])**2
    ss['std_ana2'][idxGoodBoth] =  np.sqrt(ss['std_ana2'][idxGoodBoth]/(ss['count_both'][idxGoodBoth]-1) ) 
    ss['std_sonde'][idxGoodSonde] = ss['std_sonde'][idxGoodSonde] - ss['count_sonde'][idxGoodSonde]*(ss['av_sonde_standalone'][idxGoodSonde])**2
    ss['std_sonde'][idxGoodSonde] = np.sqrt( ss['std_sonde'][idxGoodSonde] / (ss['count_sonde'][idxGoodSonde]-1) )
    return ss
def writeH5( a, ss, press_int ):
    """
    Write hdf5 tagged with input experiment and control.
    Input:
            a: command line argparse
            ss: stats dictionary
    Output: 
            hdf5 file: {experiment}_{control}_output_stats.h5 output to wherever the script is running. (cwd)
    """
    with h5py.File(a.experiment +'_'+ a.control +'_output_stats.h5','w') as f:
        for k in list(ss.keys()):
            dset = f.create_dataset(k,data=ss[k])
        dset = f.create_dataset('pressure', data = press_int )
 
if __name__ == "__main__":
    parser = argparse.ArgumentParser( description = 'compare ozone sondes')
    parser.add_argument('--experiment', help = 'experiment', required = True, dest = 'experiment')
    parser.add_argument('--control',help = 'control', required = True, dest='control')
    parser.add_argument('--ops', help = 'Optional arg to specify ops archive.', required = False, dest = 'ops',default="/archive/u/bkarpowi/")
    parser.add_argument('--cname', help="control name.", dest='cname', default='control' )
    parser.add_argument('--ename', help="experiment name.", dest='ename', default='experiment' )
    parser.add_argument('--dobson', help="path to dobson measurements.", dest='dobson_path', required=True)
    a = parser.parse_args()
    go ( a ) 
