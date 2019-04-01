import numpy as np
import h5py,os 
from datetime import datetime 

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

def readTolnetH5(fname):
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
    print(d['O3MR'].shape)
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
                profileVariables = lines[iii].strip().replace(";","").split(',')
        iii+=1
        for p in profileVariables:
            profileDict[p] = [] 
        while(iii <= end):
            strings = lines[iii].strip().split(',')
            flts = []
            for i,p in enumerate(profileVariables):
                profileDict[p].append(float(strings[i])) 
            iii+=1
        for p in profileVariables:
            profileDict[p] = np.asarray(profileDict[p]) 
        profileDicts.append(profileDict)

    return profileDicts

def readTolnet( h5Files, datFiles ):
    bigList = []
    for f in h5Files:
        bigList.extend(readTolnetH5(f))
    for f in datFiles:
        bigList.extend(readTolnetAscii(f)) 

    return bigList
    
if __name__ == "__main__":
    fs = os.listdir('hdf/h5')[0]
    #fs = 'groundbased_lidar.o3_nasa.jpl003_table.mountain.ca_20180821t203649z_20180821t213630z_002.h5' 
    fs = 'groundbased_lidar.o3_nasa.larc001_langley.research.center.va_20180730t000014z_20180730t210014z_001.h5'

    profileDicts = readTolnetH5(os.path.join('hdf/h5',fs))
    #profileDicts = readTolnetAscii('UAH/TOLNet-O3Lidar_UAH_20180806_R0.dat')
    for i,d in enumerate(profileDicts):
        print('profile {:d}'.format(i)) 
        print(d['startTime'],d['endTime'])
        print(d['Longitude'],d['Latitude'],d['Elevation'])
        pvars = []
        for k in list(d.keys()):
            if (k not in ['Latitude', 'Longitude', 'Elevation', 'startTime','endTime']):
                pvars.append(k)
                cnt = d[k].shape[0]
        print(pvars)
        print(cnt)
        for i in range(0,cnt):
            line = ''
            for kk in pvars:
               line+= "{:}".format(d[kk][i]) +','
            print(i,line)
    
    """
    fs = os.listdir('hdf/h5')[0]
    #fs = 'groundbased_lidar.o3_nasa.jpl003_table.mountain.ca_20180821t203649z_20180821t213630z_002.h5' 
    fs = 'groundbased_lidar.o3_nasa.larc001_langley.research.center.va_20180730t000014z_20180730t210014z_001.h5'
    d = readTolnetH5(os.path.join('hdf/h5',fs))
    
    #for p in  np.asarray(d['O3MR']):
    #    print (p)
   
    #print(np.asarray(d['startTime'][:]).shape)

    """

