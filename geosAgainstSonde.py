import os, h5py, argparse, glob, math,sys
import numpy as np
from scipy.interpolate import interp1d
from scipy.interpolate import InterpolatedUnivariateSpline
from lib.sonde_io import readWoudc,readShadoz
# note: you'll need matplotlib ;)
from lib.plots import plotSondeAndAnalysisStats


def go ( a ):
    """
    Main program for this thing.
    Input:
            a: argparse which gives command line arguements to do stuff with.
    Output:
            a plot tagged with experiment and control, along with and hdf5 with stats plotted. 
    """

    startLon = convertLongitude360(a.start_lon)
    endLon = convertLongitude360(a.end_lon)

    startLat = float(a.start_lat)
    endLat = float(a.end_lat)

    press_edges, press_int, nint = getPressureGridForOutput()
    
    sondeFiles = glob.glob(os.path.join(a.sonde_path,'*.csv'))
    nsondes = len(sondeFiles)

    if (nsondes == 0):
        sondeType = 'shadoz'
        sondeFiles = glob.glob(os.path.join(a.sonde_path,'*'))
        nsondes = len(sondeFiles)
    else: sondeType = 'woudc'

    sondeFiles.sort()
    # stats dictionary to populate

    ss = initProfileStats(nint)
    fcnt = 0
    for s in sondeFiles:
        lonSonde, latSonde, dateSonde, timeSonde, ozSonde_mPa, pressSonde_hPa, tempSonde_C = readProfile(s,sondeType)
        
        if(dateSonde >= a.start and dateSonde <= a.end and\
           convertLongitude360(lonSonde) >= startLon and convertLongitude360(lonSonde) <= endLon and\
           float(latSonde) >= startLat and float(latSonde) <= endLat):
            fcnt+=1
            print('file count',fcnt)
            print("Sonde read in: {}".format(s) ) 

            # Do stuff for the experimental run (get idx for the experiment and the control along the way) 
            experimentFile = getFileName(a.ops, a.experiment, dateSonde, timeSonde)
            print("Reading {}".format(experimentFile))
            v="ssh dirac 'dmget "+experimentFile+"'"
            os.system(v)
            # comment these two for dmget hackery...
            idxLon,idxLat =  getIndexFromAnalysis(experimentFile, latSonde, lonSonde)
            experimentOzone = getInterpolatedOzoneFromAnalysis(experimentFile, press_int, idxLon, idxLat)

            # now for the control
            controlFile = getFileName(a.ops, a.control, dateSonde, timeSonde)
            print("Reading {}".format(controlFile))
            v="ssh dirac 'dmget "+controlFile+"'"
            os.system(v)
            print('done dirac')

#"""
            #same grid, don't need to interpolate that again...
            controlOzone = getInterpolatedOzoneFromAnalysis(controlFile, press_int, idxLon, idxLat)

            # now for the sonde 
            interpolatedSondeOzone = interpolateSonde(1.0e15, pressSonde_hPa, ozSonde_mPa, press_edges)
            if(a.strict and ( any(controlOzone>1e3) or any(experimentOzone>1e3) ) ):
                # skip this profile in stats. because it has unphysical values for ozone, and we're doing 
                # strict rules
                print("Skipping this sonde because analysis possibly bad:{} {}".format(os.path.basename(controlFile), os.path.basename(experimentFile)))
                continue
            else:
                print("using analysis.")
                
            # s for statistics on profiles
            ss = updateStats(ss, interpolatedSondeOzone, controlOzone, experimentOzone )

    ss = finishStats( ss )
    writeH5(a, ss, press_int)
    # plot only pressure above 10 hPa
    idx = np.where( (ss['count_both'] > 1) & (press_int > 10))

    plotSondeAndAnalysisStats(press_int, ss, idx,  a.control, a.experiment, 'stats_'+a.experiment+'_'+a.control)
#"""
def convertLongitude360(lon):
    """
    convert longitude from -180 to 180 to 0 to 360.
    """
    return float(lon)%360
    
def getPressureGridForOutput():
    """
    Get the pressure grid we want to interpolate to.
    Output:
            press_edges: pressure grid edges in hPa
            press_int: pressures in the middle of the grid hPa
            nint : number of levels
    """

    #global grid stuff
    press_edges = np.asarray([1000., 900., 800., 700., 600., 500., 400., 300., 250.,\
                   200., 150., 100., \
                   90., 80., 70., 60., 50., 40., 30., 25., 20., 15.,  10., \
                   9., 8., 7., 6., 5.])
 
    nint = press_edges.shape[0] - 1
    press_int = np.sqrt(press_edges[0:nint]*press_edges[1::])

    return press_edges, press_int, nint

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

def readProfile ( s, sondeType ):
    """
    Given a filename, and type of sonde, read the data into a common set of variables. 
    If you're bored, probably should make this more "classy" by using a factory class
    and avoid if statement, but meh... 
    Input:
            s: path to the file to be read in (sonde)
    Output:
            lon                : longitude of the sonde (singleton)
            lat                : latitude of the sonde (singleton)
            date               : date string YYYYMMDD  (singleton)
            time               : time string hh:mm (singleton)
            ozPartialPress_mPa : ozone partial pressure in mPa (vector)
            press_hPa          : pressure hPa (vector)
            temp_C             : Temperature in Celsius (vector)
 

    """
    if(sondeType == 'shadoz'): 
        d,names = readShadoz(s)
        lon = d['Longitude (deg)']
        lat = d['Latitude (deg)']
        date = d['Launch Date']
        time = d['Launch Time (UT)']
        ozPartialPress_mPa = d['PROFILE']['O3 mPa']
        press_hPa = d['PROFILE']['Press hPa']
        temp_C = d['PROFILE']['Temp C']
             
    elif(sondeType == 'woudc'): 
        d = readWoudc(s)
        lon = float(d['LOCATION']['Longitude'])
        lat = float(d['LOCATION']['Latitude'])
        date = d['TIMESTAMP']['Date'].replace('-','')
        time = d['TIMESTAMP']['Time']
        ozPartialPress_mPa = d['PROFILE']['O3PartialPressure']
        press_hPa = d['PROFILE']['Pressure']
        temp_C =  d['PROFILE']['Temperature']
    else:
        sys.exit("error don't know what profile this is!")
    return lon, lat, date, time, ozPartialPress_mPa, press_hPa, temp_C

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

    return idxLon, idxLat

def getInterpolatedOzoneFromAnalysis(analysisFile, press_int, idxLon, idxLat):
    """
    Read analysis file, and interpolate ozone to output grid (press_int)
    Input:
            analysisFile: analysis file to read in
            pres_int: pressure levels to interpolate to (hPa)
            idxLat: index associated with latitude in GEOS analysis file
            idxLon: index associated with longitude in GEOS analysis file
    Output:
            interpolatedOzone: ozone partial pressure (mPa) interpolated to press_int grid at idxLon, idxLat  
    """   
    h5 = h5py.File(analysisFile,'r')
    levs = np.asarray(h5['lev'])
    geosO3 = np.asarray(h5['O3'])
    oz = np.asarray(geosO3[0,:,idxLat,idxLon])*604229.0*0.1*np.asarray(levs)
    
    pressureAnalysisFlipped = np.flipud(np.asarray(levs))
    pressureInterpolatedLevelsFlipped = np.flipud(press_int) 
    ozoneAnalysisFlipped = np.flipud(oz)    
    interpolatedOzone = np.zeros(press_int.shape[0]) 
    idxOz, = np.where( (ozoneAnalysisFlipped > 0.0) &  (ozoneAnalysisFlipped <100) )
    fOz = InterpolatedUnivariateSpline(np.log(pressureAnalysisFlipped[idxOz]), ozoneAnalysisFlipped[idxOz], k=3) 
    #fOz = interp1d(np.log(pressureAnalysisFlipped[idxOz]), ozoneAnalysisFlipped[idxOz], kind='cubic' )
    interpolatedOzone = fOz( np.log(pressureInterpolatedLevelsFlipped) )     
    interpolatedOzone = np.flipud(interpolatedOzone)

    return interpolatedOzone

def interpolateSonde(undefinedValue, pressureIn, profileIn, pressEdgesOut):
    """
    Interpolate sonde ozone profile to desired pressure grid(presEdgesOut)
    Input:
            undefinedValue: value which says the value is undefined (nan-like)
            pressureIn: pressure (hPa) of the sonde 
            profileIn: ozone partial pressure (mPa) of the sonde 
            pressEdgesOut: pressure edges used to place pressures on interpolate grid
    Output:
            interpolatedProfile: ozone partial pressure (mPa) interpolated to output pressure grid.  
    """   

    nlevs = pressEdgesOut.shape[0]-1
    interpolatedProfile = np.zeros(nlevs)
    for i in np.arange(0,nlevs):
        idx, = np.where( (pressureIn <= pressEdgesOut[i]) & (pressureIn >= pressEdgesOut[i+1]))

        if(len(idx) == 0): interpolatedProfile[i] = undefinedValue 
        else: interpolatedProfile[i] = np.mean(profileIn[idx])
        
    return interpolatedProfile
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
    idxGoodSonde, = np.where( (interpolatedSondeOzone < 1.0e15) & (interpolatedSondeOzone > 0.0) )
    idxGoodBoth, = np.where ( (interpolatedSondeOzone < 1.0e15) & (interpolatedSondeOzone > 0.0) &\
                        (controlOzone < 1e3) & (experimentOzone <1e3) & (controlOzone > 0.0) & (experimentOzone >0.0) ) 
        
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
    idxGoodSonde, = np.where( ss['count_sonde']>1 )
    idxGoodBoth, = np.where( ss['count_both'] > 1 )

    ss['av_ana1'][idxGoodBoth] = ss['av_ana1'][idxGoodBoth]/ss['count_both'][idxGoodBoth]
    ss['av_ana2'][idxGoodBoth] = ss['av_ana2'][idxGoodBoth]/ss['count_both'][idxGoodBoth]
    ss['av_sonde'][idxGoodBoth] = ss['av_sonde'][idxGoodBoth]/ss['count_both'][idxGoodBoth]

    ss['av_sonde_standalone'][idxGoodBoth] = ss['av_sonde_standalone'][idxGoodBoth]/ss['count_sonde'][idxGoodBoth]

    diff1 = ss['av_ana1'][idxGoodBoth] - ss['av_sonde'][idxGoodBoth]
    diff2 = ss['av_ana2'][idxGoodBoth] - ss['av_sonde'][idxGoodBoth]

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
    parser.add_argument('--start', help = 'start dtg YYYYMMDDhh', required = True, dest = 'start')
    parser.add_argument('--end', help = 'end dtg YYYYMMDDhh', required = True, dest = 'end')
    parser.add_argument('--ops', help = 'Optional arg to specify ops archive.', required = False, dest = 'ops',default="/archive/u/bkarpowi")
    parser.add_argument('--strict', help="reject using any bad analysis levels for ozone.", dest='strict', action='store_false' )
    #Default uses ascension island for SHADOZ
    parser.add_argument('--slat', help = 'southern most latitude', required = False, dest = 'start_lat',default="-10")
    parser.add_argument('--nlat', help = 'northern most latitude', required = False, dest = 'end_lat',default="0")
    parser.add_argument('--wlon', help = 'western most longitude', required = False, dest = 'start_lon',default="-15")
    parser.add_argument('--elon', help = 'eastern most longitude', required = False, dest = 'end_lon',default="-14")
    # tropical pacific
    #parser.add_argument('--slat', help = 'southern most latitude', required = False, dest = 'start_lat',default="-20")
    #parser.add_argument('--nlat', help = 'northern most latitude', required = False, dest = 'end_lat',default="20")
    #parser.add_argument('--wlon', help = 'western most longitude', required = False, dest = 'start_lon',default="60")
    #parser.add_argument('--elon', help = 'eastern most longitude', required = False, dest = 'end_lon',default="-90")

    parser.add_argument('--profiles', help = 'Optional arg to specify profile location.',\
                        required = False, dest = 'sonde_path',default="/archive/u/kwargan/data/SHADOZ/")
    #parser.add_argument('--profiles', help = 'Optional arg to specify profile location.',\
    #                    required = False, dest = 'sonde_path',default="/archive/u/kwargan/data/ozone_sondes/woudc2018/")

    a = parser.parse_args()
    go ( a ) 
