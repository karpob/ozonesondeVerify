import numpy as np
import h5py,os,sys 
from datetime import datetime 
from pyhdf.SD import SD, SDC

def getDatetimeFromMJD(mjdIn):
    """
    Get a useful timestamp from the somewhat useless Modified Julian Day (MJD). 
    *Grumble* Because we needed yet another variation on UNIX time, or Julian day... *Grumble*
    """
    diffBetweenMjdUnix = (datetime(2000,1,1) - datetime(1970,1,1)).total_seconds()
    secondsMjd = 86400.0*mjdIn 
    outArray = []
    for s in secondsMjd:
        outArray.append( datetime.utcfromtimestamp( s + diffBetweenMjdUnix ) )
    
    return np.asarray(outArray)

def h4get(f,dset):
    """
    helper function to pull data from h4 field.
    """
    obj = f.select(dset)
    return obj.get()

def readTolnetH4(fname):
    d = {}
    profileDicts = []
    h4 = SD( fname, SDC.READ )
    #  get available records in hdf (you want this, sadly, because nobody calls o3 mixing ratio the same thing, 
    # and in one case time isn't even called the same thing).
    availableSDS = list(h4.datasets().keys())

    d['ALT'] = np.asarray( h4get(h4,'ALTITUDE') )
    d['Elevation'] = np.asarray( h4get(h4,"ALTITUDE.INSTRUMENT") )

    #annoyingly time is different in one case, so we have ugly "ifs" here...
    if 'DATETIME.START' in availableSDS: 
        d['startTime'] = getDatetimeFromMJD( np.asarray( h4get(h4,"DATETIME.START") ) )
    elif 'DATETIME' in availableSDS:
        d['startTime'] = getDatetimeFromMJD( np.asarray( h4get(h4,"DATETIME") ) ) 
    else:
        sys.exit("Sorry. Time isn't even in this file: {}".format(fname))

    # again for the FTIR Instruments, they only give one time
    if 'DATETIME.STOP' in availableSDS:
        d['endTime'] = getDatetimeFromMJD( np.asarray( h4get(h4,"DATETIME.STOP") ) )
    elif 'DATETIME' in availableSDS: 
        d['endTime'] = getDatetimeFromMJD( np.asarray( h4get(h4,"DATETIME") ) )


    d['dT'] = np.asarray( h4get(h4,"INTEGRATION.TIME") ) 
    d['Latitude']= np.asarray( h4get(h4,"LATITUDE.INSTRUMENT") ) 
    d['Longitude']= np.asarray( h4get(h4,"LONGITUDE.INSTRUMENT") )


    # This gets fun, because everyone decided to call mixing ratio something slightly different.  
    if 'O3.MIXING.RATIO.VOLUME_DERIVED' in availableSDS:
        d['O3MR'] = np.asarray( h4get(h4,"O3.MIXING.RATIO.VOLUME_DERIVED") )
    elif 'O3.MIXING.RATIO.VOLUME_EMISSION' in availableSDS: 
        d['O3MR'] = np.asarray( h4get(h4,"O3.MIXING.RATIO.VOLUME_EMISSION") )
    elif 'O3.MIXING.RATIO.VOLUME_ABSORPTION.SOLAR' in availableSDS:
        d['O3MR'] = np.asarray( h4get(h4,"O3.MIXING.RATIO.VOLUME_ABSORPTION.SOLAR") )
    else:
        sys.exit("Weird file:{} No ozone information.".format(fname))

    # same for uncertianty, although I haven't used it yet.
    if 'O3.MIXING.RATIO.VOLUME_DERIVED_UNCERTAINTY.COMBINED.STANDARD' in availableSDS:
        d['O3MRUncert'] = np.asarray( h4get(h4,"O3.MIXING.RATIO.VOLUME_DERIVED_UNCERTAINTY.COMBINED.STANDARD") )
    elif 'O3.MIXING.RATIO.VOLUME_EMISSION_UNCERTAINTY.COMBINED.STANDARD' in availableSDS:
        d['O3MRUncert'] = np.asarray( h4get(h4,"O3.MIXING.RATIO.VOLUME_EMISSION_UNCERTAINTY.COMBINED.STANDARD") )
    elif 'O3.COLUMN_ABSORPTION.SOLAR_UNCERTAINTY.SYSTEMATIC.STANDARD' in availableSDS: 
        d['O3MRUncert'] = np.asarray( h4get(h4,"O3.COLUMN_ABSORPTION.SOLAR_UNCERTAINTY.SYSTEMATIC.STANDARD") )
    else:
        sys.exit("Weird file:{} No ozone uncertianty information.".format(fname)) 

    # if it's a lidar get number density. 
    
    if 'O3.NUMBER.DENSITY_ABSORPTION.DIFFERENTIAL' in availableSDS:
        d['O3ND'] = np.asarray ( h4get(h4,"O3.NUMBER.DENSITY_ABSORPTION.DIFFERENTIAL") )
    
    d['Press'] = np.asarray( h4get(h4,"PRESSURE_INDEPENDENT") )
    d['Temp'] = np.asarray( h4get(h4,"TEMPERATURE_INDEPENDENT") )
    if(len(d['Press'].shape)==2):    
        if (d['Press'][0,0]<d['Press'][0,1]): flipIt = True
        else:flipIt = False
    else:
        if (d['Press'][0]<d['Press'][1]): flipIt = True
        else:flipIt = False
 

    dims = d['O3MR'].shape
    if(len(dims) > 1):
        nProfiles, nLevels = dims[0],dims[1]
        for i in range(0,nProfiles):
            dd = {}
            dd['startTime'] = d['startTime'][i]
            dd['endTime'] = d['endTime'][i]
            if (flipIt): dd['O3MR'] = d['O3MR'][i,::-1]
            else: dd['O3MR'] = d['O3MR'][i,:]
            if('O3ND' in list( d.keys() ) ): 
                if(flipIt):dd['O3ND'] = d['O3ND'][i,::-1]
                else: dd['O3ND'] = d['O3ND'][i,:]
            if(flipIt):dd['Press'] = d['Press'][i,::-1]
            else:dd['Press'] = d['Press'][i,:]
            if(flipIt): dd['Temp'] = d['Temp'][i,::-1]
            else:dd['Temp'] = d['Temp'][i,:]
            dd['Longitude'] = d['Longitude']
            dd['Latitude'] = d['Latitude']
            dd['Elevation'] = d['Elevation']
            profileDicts.append(dd)
    else:
        dd = {}
        dd['startTime'] = d['startTime']
        dd['endTime'] = d['endTime']
        if(flipIt): dd['O3MR'] = d['O3MR'][::-1]
        else: dd['O3MR'] = d['O3MR'][:]
        if( 'O3ND' in list(d.keys())): 
            if(flipIt): dd['O3ND'] = d['O3ND'][::-1]
            else: dd['O3ND'] = d['O3ND'][:]
        if(flipIt): dd['Press'] = d['Press'][::-1]
        else: dd['Press'] = d['Press']
        if(flipIt): dd['Temp'] = d['Temp'][::-1]
        else: dd['Temp'] = d['Temp']
        dd['Longitude'] = d['Longitude']
        dd['Latitude'] = d['Latitude']
        dd['Elevation'] = d['Elevation']
        profileDicts.append(dd)
   
 

    # because they had to be different..same as ".close()" would be for any other api in the universe.
    h4.end() 
    return profileDicts


def readTolnetH5(fname):
    """
    Pretty much the same as readTolnetH4, but I wrote this first. I started out converting these things
    to hdf5, then reading them, because I really don't like the pyhdf interface. Keeping this here for 
    the day when pyhdf gets orphaned, or the world finally agrees with me that hdf4 is just too old.
    """
    d = {}
    profileDicts = []
    h5 = h5py.File(fname,'r')

    d['ALT'] = np.asarray( h5['ALTITUDE'] )
    d['Elevation'] = np.asarray( h5["ALTITUDE.INSTRUMENT"])
    d['startTime'] = getDatetimeFromMJD(np.asarray( h5["DATETIME.START"]) ) 
    d['endTime'] = getDatetimeFromMJD( np.asarray(h5["DATETIME.STOP"]) )
    d['dT'] = np.asarray( h5["INTEGRATION.TIME"] ) 
    d['Latitude']= np.asarray( h5["LATITUDE.INSTRUMENT"] ) 
    d['Longitude']= np.asarray( h5["LONGITUDE.INSTRUMENT"] )  
    d['O3MR'] = np.asarray( h5["O3.MIXING.RATIO.VOLUME_DERIVED"] )
    d['O3MRUncert'] = np.asarray( h5["O3.MIXING.RATIO.VOLUME_DERIVED_UNCERTAINTY.COMBINED.STANDARD"] ) 
    d['O3ND'] = np.asarray ( h5["O3.NUMBER.DENSITY_ABSORPTION.DIFFERENTIAL"] )
    d['Press'] = np.asarray( h5["PRESSURE_INDEPENDENT"])
    d['Temp'] = np.asarray(h5["TEMPERATURE_INDEPENDENT"])
    dims = d['O3MR'].shape
    if(len(dims) > 1):
        nProfiles, nLevels = dims[0],dims[1]
        for i in range(0,nProfiles):
            dd = {}
            dd['startTime'] = d['startTime'][i]
            dd['endTime'] = d['endTime'][i]
            dd['O3MR'] = d['O3MR'][i,:]
            dd['O3ND'] = d['O3ND'][i,:]
            dd['Press'] = d['Press']
            dd['Temp'] = d['Temp']
            dd['Longitude'] = d['Longitude']
            dd['Latitude'] = d['Latitude']
            dd['Elevation'] = d['Elevation']
            profileDicts.append(dd)
    else:
        dd = {}
        dd['startTime'] = d['startTime']
        dd['endTime'] = d['endTime']
        dd['O3MR'] = d['O3MR'][:]
        dd['O3ND'] = d['O3ND'][:]
        dd['Press'] = d['Press']
        dd['Temp'] = d['Temp']
        dd['Longitude'] = d['Longitude']
        dd['Latitude'] = d['Latitude']
        dd['Elevation'] = d['Elevation']
        profileDicts.append(dd)

     
    return profileDicts

def readTolnetAscii(fname):
    with open(fname) as f: lines = f.readlines()
    
    recordStarts = []
    recordEnds = []
    cnt = 0
    for i,l in enumerate(lines):
        if('#BEGIN' in l):
            recordStarts.append(i)
    eof = i
    recordEnds = np.asarray(recordStarts[1::])-1
    recordStarts = np.asarray(recordStarts)+1
    recordStarts = recordStarts.tolist()
    recordEnds = recordEnds.tolist()
    recordEnds.append(eof)
    profileDicts = []
    
    #iterate through each profile
    for i,p in enumerate(recordStarts):
        profileDict = {}
        #first line after # Begin is the number of header lines
        start = p
        end = recordEnds[i]
    
        ii = start
        nhead = int(lines[ii].split(';')[0]) 
        ii+=1
        for iii in range(ii,ii+nhead):
            if('ALT, O3ND' not in lines[iii]):
                val,key = lines[iii].strip().split(';')
                if('START' in key):
                    profileDict['startTime'] =datetime.strptime(val.replace(" ",""),"%Y-%m-%d,%H:%M:%S")
                elif('END' in key):
                    profileDict['endTime'] = datetime.strptime(val.replace(" ",""),"%Y-%m-%d,%H:%M:%S")
                elif('LONGITUDE' in key):
                    profileDict['Longitude'] = float(val.split(',')[0])
                    profileDict['Latitude'] = float(val.split(',')[1])
                    profileDict['Elevation'] = float(val.split(',')[2])
            else:
                profileVariablesTmp = lines[iii].strip().replace(";","").split(',')
        iii+=1
        profileVariables = []
        #probably should replace this with a zip of some kind.
        for p in profileVariablesTmp: profileVariables.append(p.replace(" ",""))

        for p in profileVariables:
            profileDict[p] = [] 
        while(iii <= end):
            strings = lines[iii].strip().split(',')
            flts = []
            for i,p in enumerate(profileVariables):
                profileDict[p].append(float(strings[i])) 
            iii+=1
        for p in profileVariables:
            #ascii files are in PPB for some reason, make it ppmv consistent with hdf
            if('O3MR') in p: 
                profileDict[p] = 1.0e-3*np.asarray(profileDict[p])
                # make -9.999 back into -999.9
                idx, = np.where( profileDict[p] < -9.0 )
                profileDict[p][idx] = -999.0
 
            else: profileDict[p] = np.asarray(profileDict[p]) 
        profileDicts.append(profileDict)
    return profileDicts

def readTolnet( h5Files, datFiles ):
    bigList = []
    for f in h5Files:
        bigList.extend(readTolnetH4(f))
    for f in datFiles:
        bigList.extend(readTolnetAscii(f)) 

    return bigList
    
if __name__ == "__main__":
    #fs = 'groundbased_lidar.o3_nasa.jpl003_table.mountain.ca_20180821t203649z_20180821t213630z_002.h5' 
    #fs = 'groundbased_lidar.o3_nasa.larc001_langley.research.center.va_20180730t000014z_20180730t210014z_001.h5'
    fs = '../ftp.cpc.ncep.noaa.gov/ndacc/station/rikubets/hdf/ftir/groundbased_ftir.o3_unagoya001_rikubetsu_20180117t031331z_20180727t010244z_002.hdf'  
    #fs = '../TOLnet/hdf/h5/groundbased_lidar.o3_nasa.larc001_langley.research.center.va_20180701t000010z_20180702t000010z_001.h5'
    profileDicts1 = readTolnetH4(fs)
#    profileDicts2 = readTolnetAscii('../TOLnet/zip/larc/TOLNet-O3Lidar_LaRC_20180701_R0.dat')
    print( list(profileDicts1[0].keys() ) )
    for i in range(0,len(profileDicts1[0]['O3MR'][:])):  
        print (profileDicts1[0]['Press'][i],profileDicts1[0]['O3MR'][i]  )

