import numpy as np
import sys
def readWoudc(a):
    """
    Reads in WOUDC ascii data format for ozonesondes.

    Yeah... should comment this madness at some point.
    Let's just say it reads in an ascii file tags meta data and reads profile data in that order
    """
    o = {}
    profileNames = []
    #open file and read in lines.
    with open(a) as f:
        lines = f.readlines()
    # generate dictionary of dictonaries based on # as headers, and immediate lines as datasets.
    endLine = -1
    for i,l in enumerate(lines):
        if '#' in l:            
            if(',,' in l): l = l.replace(',,',',-999.99,') 
            if('#PRELAUNCH' in l): endLine = i
            #Make dictonary of dictoniaries using # names as first metadata header stripping off #
            o[l.strip().replace('#','')] = {}
            # for everything except profile create sub dictionary with headers on next line, and data on following line.
            if '#PROFILE' not in l:
                for ii,var in enumerate(lines[i+1].strip().split(',')):
                    if (ii < len(lines[i+2].strip().split(','))):
                        o[l.strip().replace('#','')][var] = lines[i+2].strip().split(',')[ii]
            # for #PROFILE use the next line for metadata, and initialize arrays for each metadata tag, mark where profile starts
            else:
                for var in lines[i+1].strip().split(','):
                    o[l.strip().replace('#','')][var] = []
                    profileNames.append(var)
                profileLine = i 
    
    #loop through lines that mark the profile, and populate dictionary.
    for l in lines[profileLine+2:endLine]:
        if(',,,' in l): l = l.replace(',,,',',-999.99,-999.99,') 
        if(',,' in l): l = l.replace(',,',',-999.99,')
        if '*' in l: continue
        for i,v in enumerate(l.strip().split(',')):
            if (profileNames[i] == 'LevelCode'):
                if v ==' ' or v=='-999.99': 
                    o['PROFILE'][profileNames[i]].append(int(999))
                else:
                    o['PROFILE'][profileNames[i]].append(int(v))
            elif len(v.strip())>0:
                    o['PROFILE'][profileNames[i]].append(float(v))
            else:
                o['PROFILE'][profileNames[i]].append(float(-999.99))
        if (len(l.strip().split(','))<10):
            for ii in range(len(l.strip().split(',')),10):
                o['PROFILE'][profileNames[ii]].append(float(-999.99))
    # filter profile to remove -999.99
    for p in profileNames: o['PROFILE'][p] = np.asarray(o['PROFILE'][p])
    idx, = np.where(o['PROFILE']['O3PartialPressure'] != -999.99)
    for k in list(o['PROFILE'].keys()):
        o['PROFILE'][k] = o['PROFILE'][k][idx]
    return o

def readShadoz(a):
    """
    Reads in ozone sonde in ascii format from SHADOZ.
    Yeah... should comment this madness at some point.
    Let's just say it reads in an ascii file tags meta data and reads profile data in that order
    """

    o = {}
    with open(a) as f:
        lines = f.readlines()
    for i,l in enumerate(lines):
        if ': ' in l:
            ll  = l.strip()
            sizeLL = len(ll.split(': '))
            # this line is a little wacky to get rid of whitespace at the end of metadata tag
            if(sizeLL>1):
                o[ " ".join(ll.split(': ')[0].replace('   ','').split()) ] = ll.split(': ')[1]
        elif 'Time' in l:
            variablesAndJunk = l.strip().split('  ')
            variables = []
            for v in variablesAndJunk:
                if(v != ''): variables.append(v)
            units = lines[i+1].strip().split()
            o['PROFILE'] = {}
            profileNames = []
            for ii,v in enumerate(variables):
                if v[0] ==' ': vv=v[1::]
                else : vv = v 
                o['PROFILE'][vv+' '+units[ii]] = []
                profileNames.append(vv+' '+units[ii])
            profileStart = i+2
    for l in lines[profileStart::]:
        for i,v in enumerate(l.strip().split()):
            if('Time' in profileNames[i]):
                o['PROFILE'][profileNames[i]].append(int(v))
            else:
                o['PROFILE'][profileNames[i]].append(float(v))
    for p in profileNames:
        o['PROFILE'][p] = np.asarray(o['PROFILE'][p])
    p = profileNames
    idx, = np.where( (o['PROFILE'][p[0]] != 9000.0) &\
                     (o['PROFILE'][p[1]] != 9000.0) &\
                     (o['PROFILE'][p[2]] != 9000.0) &\
                     (o['PROFILE'][p[3]] != 9000.0) &\
                     (o['PROFILE'][p[4]] != 9000.0) &\
                     (o['PROFILE'][p[5]] != 9000.0) &\
                     (o['PROFILE'][p[6]] != 9000.0) )
    for k in profileNames: o['PROFILE'][k] = np.asarray(o['PROFILE'][k])[idx]


    return o, profileNames
    
if __name__ == "__main__":
    #a = '/archive/u/kwargan/data/ozone_sondes/woudc2018/20181214.ecc.z.z32286.fmi-smna.csv'
    a = '/archive/u/kwargan/data/SHADOZ/ascen_20180103T12_V05.1_R.dat'
    #o = readWoudc(a)
    print('lon, lat, date, hour, minute, oz, press, temp')
    #print( list(o.keys()))
    #print( o['TIMESTAMP'])
    """
    oo = o['PROFILE']                 
    print ( 'Pressure', 'O3PartialPressure', 'Temperature', 'WindSpeed', 'WindDirection', 'LevelCode', 'Duration', 'GPHeight', 'RelativeHumidity', 'SampleTemperature')
    for i in range(len(oo['Pressure'])):
         print ( oo['Pressure'][i], oo['O3PartialPressure'][i],oo['Temperature'][i],oo['WindSpeed'][i], oo['WindDirection'][i], oo['LevelCode'][i],oo['Duration'][i],oo['GPHeight'][i],oo['RelativeHumidity'][i], oo['SampleTemperature'][i] )
    """
    
    o,names = readShadoz(a)  
    oo = o['PROFILE']
    print(list(oo.keys()))
    print(list(o.keys()))
    print( o['Launch Time (UT)'])
    print( o['Launch Date'])
    """ 
    lon = o['Longitude (deg)']
    lat = o['Latitude (deg)']
    date = o['Launch Date']
    time = o['Launch Time (UT)']
    oz = o['PROFILE']['O3 mPa']
    press = o['PROFILE']]['Press hPa']
    temp = o['PROFILE']['Temp C']
['Time sec', 'Press hPa', ' Alt km', 'Temp C', 'RH %', 'O3 mPa', 'O3 ppmv', 'O3 du', 'W Dir deg', ' W Spd m/s', 'T Pump C', 'I O3 uA', ' GPS lon deg', ' GPS lat deg', ' GPS alt km']
['NASA/GSFC/SHADOZ Archive', 'SHADOZ Version', 'SHADOZ format data created', 'STATION', 'SHADOZ Principal Investigator', 'Station Principal Investigator(s)', 'Latitude (deg)', 'Longitude (deg)', 'Elevation (m)', 'Launch Date', 'Launch Time (UT)', 'Highest level reached (hPa)', 'Integrated O3 until EOF (DU)', 'Sonde/MLS Climatology(1988-2010)', 'Sonde Instrument, SN', 'Radiosonde, SN', 'KI Solution', 'Applied pump corrections', 'Pump flow rate (sec/100ml)', 'Background current (uA)', 'Missing or bad values', 'PROFILE']

#    for i in range(len(oo['Press hPa'])):
#       print(  oo['Time sec'][i], oo['O3 ppmv'][i])
    """     
