import numpy as np
import os,sys
from datetime import datetime
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

def readShadoz5(a):
    """
    Reads in ozone sonde in ascii format from SHADOZ. For SHADOZ V05.1 format. 
    Yeah... should comment this madness at some point.
    Let's just say it reads in an ascii file tags meta data and reads profile data in that order

    Input: A SHADOZ V05.1 formatted ascii file. (full path)
    Output:
           o - A dictionary with keys of metadata + subdictionary 'PROFILE' which is the actual profile data.
           profileNames - list of the Profile variables in the second dictionary.
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
        elif ('Time' in l) or ('Level' in l and 'Press' in l):
            variablesAndJunk = l.strip().split('  ')
            variables = []
            for v in variablesAndJunk:
                if(v != '') and 'O3 # DN O3 Res' not in v: variables.append(v)
                if 'O3 # DN O3 Res' in v:
                    variables.append('O3 # DN')
                    variables.append('O3 Res')
                 
            units = lines[i+1].strip().split()
            o['PROFILE'] = {}
            profileNames = []
            for ii,v in enumerate(variables):
                if v[0] ==' ': vv=v[1::]
                else : vv = v 
                o['PROFILE'][vv+' '+units[ii]] = []
                profileNames.append(vv+' '+units[ii])
            profileStart = i+2
    # add extra bit to make sure date comes back as number not something a human would use like April.
    for k in list(o.keys()):
        if 'Launch Date' in k:
            lDateKey = k
        elif 'Launch Time' in k:
            lTimeKey = k 

    if( not o[lDateKey].isdigit() ):
        o[lDateKey] = datetime.strptime(o[lDateKey],"%d %B %Y").strftime("%Y%m%d")
    if ('GMT' in o[lTimeKey]):
        o[lTimeKey] = ':'.join(o[lTimeKey].split(':')[0:2])
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

def readShadoz6(a):
    """
    Reads in ozone sonde in ascii format from SHADOZ.
    Yeah... should comment this madness at some point.
    Let's just say it reads in an ascii file tags meta data and reads profile data in that order
    Input: A SHADOZ V06 formatted ascii file. (full path)
    Output:
           o - A dictionary with keys of metadata + subdictionary 'PROFILE' which is the actual profile data.
           profileNames - list of the Profile variables in the second dictionary.
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
        elif ('Time' in l) or ('Level' in l and 'Press' in l):
            variablesAndJunk = l.strip().split(' ')
            variables = []
            for v in variablesAndJunk:
                if(v != '') and 'O3 # DN O3 Res' not in v: variables.append(v)
                if 'O3 # DN O3 Res' in v:
                    variables.append('O3 # DN')
                    variables.append('O3 Res')
                 
            units = lines[i+1].strip().split()
            o['PROFILE'] = {}
            profileNames = []
            for ii,v in enumerate(variables):
                if v[0] ==' ': vv=v[1::]
                else : vv = v 
                vv = vv.replace('O3_mPa','O3')
                vv = vv.replace('O3_ppmv','O3')
                vv = vv.replace('O3_DU','O3')
                o['PROFILE'][vv+' '+units[ii]] = []
                profileNames.append(vv+' '+units[ii])
            profileStart = i+2
    # add extra bit to make sure date comes back as number not something a human would use like April.
    for k in list(o.keys()):
        if 'Launch Date' in k:
            lDateKey = k
        elif 'Launch Time' in k:
            lTimeKey = k 

    if( not o[lDateKey].isdigit() ):
        o[lDateKey] = datetime.strptime(o[lDateKey],"%d %B %Y").strftime("%Y%m%d")
    if ('GMT' in o[lTimeKey]):
        o[lTimeKey] = ':'.join(o[lTimeKey].split(':')[0:2])
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

def readShadozU1(a):
    """
    Reads in ozone sonde in ascii format from SHADOZ.
    Yeah... should comment this madness at some point.
    Let's just say it reads in an ascii file tags meta data and reads profile data in that order

    Input: A SHADOZ V01_uncertainties formatted ascii file. (full path)
    Output:
           o - A dictionary with keys of metadata + subdictionary 'PROFILE' which is the actual profile data.
           profileNames - list of the Profile variables in the second dictionary.
           **caution** temperature returned here is fake! I hard coded something here as to not break the overarching profile reader!
                       V01_uncertainties doesn't have a Temperature!
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
            profileStart = i+1
    # add extra bit to make sure date comes back as number not something a human would use like April.
    for k in list(o.keys()):
        if 'Launch Date' in k:
            lDateKey = k
        elif 'Launch Time' in k:
            lTimeKey = k 

    if( not o[lDateKey].isdigit() ):
        o[lDateKey] = datetime.strptime(o[lDateKey],"%d %B %Y").strftime("%Y%m%d")
    if ('GMT' in o[lTimeKey]):
        o[lTimeKey] = ':'.join(o[lTimeKey].split(':')[0:2])
    profileNames = []
    for k in list(o.keys()):
        if 'Column' in k:
            profileNames.append(o[k].replace('[','').replace(']','').replace('Pressure','Press'))
    # construct PROFILE entry into dictionary.
    o['PROFILE'] = {}
    # Add empty Temp C, because the main reader I wrote starting out wants Temperature! 
    profileNames.append('Temp C')
    for p in profileNames:
        o['PROFILE'][p] = []
    
    for l in lines[profileStart::]:
        for i,v in enumerate(l.strip().split()):
            if('Time' in profileNames[i]):
                o['PROFILE'][profileNames[i]].append(int(v))
            else:
                o['PROFILE'][profileNames[i]].append(float(v))
        o['PROFILE']['Temp C'].append(float(9000.0))
    for p in profileNames:
        o['PROFILE'][p] = np.asarray(o['PROFILE'][p])
    p = profileNames
    idx, = np.where( (o['PROFILE'][p[0]] != 9000.0) &\
                     (o['PROFILE'][p[1]] != 9000.0) &\
                     (o['PROFILE'][p[2]] != 9000.0) &\
                     (o['PROFILE'][p[3]] != 9000.0) )
    for k in profileNames: o['PROFILE'][k] = np.asarray(o['PROFILE'][k])[idx]


    return o, profileNames

    
def readShadoz(a):
    """
    Overarching reader for SHADOZ. Selects which version based of filename (Yuck, but hey -- you figure out a better way!)
    """
    if ('V05.1_R' in os.path.split(a)[-1] ): return readShadoz5(a)
    elif('V06' in os.path.split(a)[-1] ): return readShadoz6(a)
    elif('SHADOZV01_uncertainty' in os.path.split(a)[-1] ): return readShadozU1(a)
    else:
        print('Unknown version of SHADOZ data. Looking for V05.1, V06, or SHADOZV01_uncertainty in filename.')
        sys.exit(1)
if __name__ == "__main__":
    #a = '/archive/u/kwargan/data/ozone_sondes/woudc2018/20181214.ecc.z.z32286.fmi-smna.csv'
    a = '/discover/nobackup/projects/gmao/obsdev/bkarpowi/SouthPole/sp001_2003_06_07_21.l100'
    o,names = readShadoz(a)
    print (o['PROFILE'])
