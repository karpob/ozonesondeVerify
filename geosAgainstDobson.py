#!/usr/bin/env python3
import os, h5py, argparse, glob, math,sys
import numpy as np
from lib.dobson_io import readWoudc
import pytz
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from datetime import datetime, timedelta
from scipy.stats import pearsonr
def hour_rounder(t):
    # Rounds to nearest hour by adding a timedelta hour if minute >= 30
    return (t.replace(second=0, microsecond=0, minute=0, hour=t.hour)
               +timedelta(hours=t.minute//30))

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
    dOb['datetime'] = []
    dOb['O3'] = []   
    for s in dobsonFilesCsv:
        lonDobson, latDobson, datetimeDobson, ozoneDobson = readDobson(s, 'woudc')
        dOb['lon'].extend(lonDobson)
        dOb['lat'].extend(latDobson)
        dOb['O3'].extend(ozoneDobson)    
        dOb['datetime'].extend(datetimeDobson) 

        for t in datetimeDobson:
            dateDobson = t.astimezone(pytz.UTC).strftime('%Y%m%d')
            timeDobson = t.astimezone(pytz.UTC).strftime('%H:%M') 
            # Do stuff for the experimental run (get idx for the experiment and the control along the way) 
            experimentAnalysisFiles.append( getFileName(a.ops, a.experiment, dateDobson, timeDobson) )
            dOb['date'].extend(dateDobson)
            dOb['time'].extend(timeDobson)
 
            # now for the control
            controlAnalysisFiles.append( getFileName(a.ops, a.control, dateDobson, timeDobson) )
    uniqueE = []
    uniqueC = []
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

    idxLon,idxLat =  getIndexFromAnalysis(experimentAnalysisFiles[0], dOb['lat'][0], dOb['lon'][0])
    print(len(dOb['O3']),len(experimentAnalysisFiles))
    for i,O3 in enumerate(dOb['O3']):
        
        # use sonde latitude to get x,y from analysis.
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
    diffE = np.asarray(dOb['O3']) - np.asarray(experimentOzone)
    diffC = np.asarray(dOb['O3']) - np.asarray(controlOzone)

    print('std experiment, control',np.std(diffE), np.std(diffC)) 
    print('rms experiment, control',np.sqrt(diffE**2).mean(), np.sqrt(diffC**2).mean()) 
    plt.plot(np.asarray(dOb['datetime']),diffE,'rx')
    plt.plot(np.asarray(dOb['datetime']),diffC,'bx')
    plt.xticks(rotation=90)
    #plt.plot(np.asarray(controlOzone),'ko')
    plt.savefig('whir.png')
    plt.close()
    plt.plot(np.asarray(dOb['datetime']),np.asarray(experimentOzone),'rx')
    plt.plot(np.asarray(dOb['datetime']),np.asarray(controlOzone),'bx')
    plt.plot(np.asarray(dOb['datetime']),dOb['O3'],'kx')
    plt.xticks(rotation=90)
    plt.savefig('whir2.png')
    plt.close()
    xmin,xmax = min(dOb['O3']), max(dOb['O3'])
    x= np.linspace(xmin,xmax,100)
    plt.plot(dOb['O3'],np.asarray(experimentOzone),'rx')
    plt.plot(dOb['O3'],np.asarray(controlOzone),'bx')
    plt.plot(x,x,'k')
    plt.savefig('whir3.png')
    r1 = pearsonr(dOb['O3'], np.asarray(experimentOzone))
    r2 = pearsonr(dOb['O3'], np.asarray(controlOzone))
    print(r1,r2)

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
        idx = []
        
        i=0
        for c in list(d['OBSERVATIONS']['StdDevO3'][:]):
            if (float(c) < 0.7):
                idx.append(i)
            i+=1
        idx = np.asarray(idx)
        time = np.asarray(d['OBSERVATIONS']['Time'])[idx].tolist()
        dobsonO3 = np.asarray(d['OBSERVATIONS']['ColumnO3'])[idx].tolist()
        
        for dd in dobsonO3:
            lon.append(alon)
            lat.append(alat)
    else:
        sys.exit("error don't know what dobson filetype this is!")
   
    return lon, lat, time, dobsonO3
def timeLookupHour(dateString, timeString):
    """
    Lookup hourly. 

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
    t = datetime(YYYY, MM, DD, hh,mm,0)
    tt = hour_rounder(t)
    return tt.year, tt.month, tt.day, tt.hour, tt.minute 

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
    YYYY, MM, DD, hh, mm =  timeLookupHour(dateSonde, timeSonde)
    analysisFile = os.path.join(ops, experiment,'diag','Y{:04d}'.format(YYYY), 'M{:02d}'.format(MM),\
                   experiment+'.inst1_2d_asm_Nx.'+dateSonde+'_'+'{:02d}00z.nc4'.format(hh))
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

 
if __name__ == "__main__":
    parser = argparse.ArgumentParser( description = 'compare ozone sondes')
    parser.add_argument('--experiment', help = 'experiment', required = True, dest = 'experiment')
    parser.add_argument('--control',help = 'control', required = True, dest='control')
    parser.add_argument('--ops', help = 'Optional arg to specify ops archive.', required = False, dest = 'ops',default="/discover/nobackup/projects/gmao/obsdev/bkarpowi/")
    parser.add_argument('--cname', help="control name.", dest='cname', default='control' )
    parser.add_argument('--ename', help="experiment name.", dest='ename', default='experiment' )
    parser.add_argument('--dobson', help="path to dobson measurements.", dest='dobson_path', required=True)
    a = parser.parse_args()
    go ( a ) 
