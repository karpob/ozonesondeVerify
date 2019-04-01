import numpy as np
import h5py,os 
from datetime import datetime 
def readTolnetH5(fname):
    h5 = h5py.File(fname,'r')

    d['ALT'] = h5['ALTITUDE']
    d['Elevation'] = h5["ALTITUDE.INSTRUMENT"]
    d['startTime'] = h5["DATETIME.START"]
    d['endTime'] = h5["DATETIME.STOP"]
    d['dT'] = h5["INTEGRATION.TIME"]
    d['Lonitude']= h5["LATITUDE.INSTRUMENT"]
    d['Longitude']= h5["LONGITUDE.INSTRUMENT"] 
    d['O3MR'] = h5["O3.MIXING.RATIO.VOLUME_DERIVED"]
    d['O3MRUncert'] = h5["O3.MIXING.RATIO.VOLUME_DERIVED_UNCERTAINTY.COMBINED.STANDARD"]
    d['O3ND'] = h5["O3.NUMBER.DENSITY_ABSORPTION.DIFFERENTIAL"]
    #d['O3NDResol'] = h5["O3.NUMBER.DENSITY_ABSORPTION.DIFFERENTIAL_RESOLUTION.ALTITUDE.DF.CUTOFF"]
    #d['Precision'] = h5["O3.NUMBER.DENSITY_ABSORPTION.DIFFERENTIAL_RESOLUTION.ALTITUDE.IMPULSE.RESPONSE"]
    #d['ChRange'] = h5["O3.NUMBER.DENSITY_ABSORPTION.DIFFERENTIAL_RESOLUTION.ALTITUDE.DISTANCE.FROM.IMPULSE"]
    #d['O3NDUncert'] = h5["O3.NUMBER.DENSITY_ABSORPTION.DIFFERENTIAL_UNCERTAINTY.COMBINED.STANDARD"]
    d['Press'] = h5["PRESSURE_INDEPENDENT"]
    d['Temp'] = h5["TEMPERATURE_INDEPENDENT"]

    return d

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
if __name__ == "__main__":
    
    profileDicts = readTolnetAscii('UAH/TOLNet-O3Lidar_UAH_20180806_R0.dat')
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
    
    fs = os.listdir('hdf/h5')[0]
    fs = 'groundbased_lidar.o3_nasa.jpl003_table.mountain.ca_20180821t203649z_20180821t213630z_002.h5' 
    d = readTolnetH5(os.path.join('hdf/h5',fs))
    for p in  np.asarray(d['O3MR']):
        print (p)

    



