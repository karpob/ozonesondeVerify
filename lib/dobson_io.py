import numpy as np
import os,sys
from datetime import datetime
from dateutil.parser import parse
def readWoudc(a):
    """
    Reads in WOUDC ascii data format for total column measurements.

    Yeah... should comment this madness at some point.
    Let's just say it reads in an ascii file tags meta data and reads OBSERVATIONS data in that order
    """
    o = {}
    OBSERVATIONSNames = []
    #open file and read in lines.
    with open(a) as f:
        lines = f.readlines()
    # generate dictionary of dictonaries based on # as headers, and immediate lines as datasets.
    endLine = -1
    for i,l in enumerate(lines):
        if '#' in l:            
            if(',,' in l): l = l.replace(',,',',-999.99,')

            if('#DAILY_SUMMARY' in l):
                endLine = i
                break 
            #Make dictonary of dictoniaries using # names as first metadata header stripping off #
            o[l.strip().replace('#','')] = {}
            # for everything except OBSERVATIONS create sub dictionary with headers on next line, and data on following line.
            if '#OBSERVATIONS' not in l:
                for ii,var in enumerate(lines[i+1].strip().split(',')):
                    if (ii < len(lines[i+2].strip().split(','))):
                        o[l.strip().replace('#','')][var] = lines[i+2].strip().split(',')[ii]
            # for #OBSERVATIONS use the next line for metadata, and initialize arrays for each metadata tag, mark where OBSERVATIONS starts
            else:
                for var in lines[i+1].strip().split(','):
                    o[l.strip().replace('#','')][var] = []
                    OBSERVATIONSNames.append(var)
                OBSERVATIONSLine = i 
    
    #loop through lines that mark the OBSERVATIONS, and populate dictionary.
    for l in lines[OBSERVATIONSLine+2:endLine]:
        if(',,,' in l): l = l.replace(',,,',',-999.99,-999.99,') 
        if(',,' in l): l = l.replace(',,',',-999.99,')
        if(',     ,     ,' in l): l = l.replace(',     ,     ,',',-999.99,-999.99,')
        if(',     ,' in l): l = l.replace(',     ,',',-999.99,')
        if '*' in l: continue
        for i,v in enumerate(l.strip().split(',')):
            if ('Time' in OBSERVATIONSNames[i] ):
                timestring = o['TIMESTAMP']['Date']+'T'+v
                print(l)
                print(v)
                if(len(v)>0):
                    print(timestring+o['TIMESTAMP']['UTCOffset'][0:5])
                    aa = parse(timestring+o['TIMESTAMP']['UTCOffset'][0:5])
                    o['OBSERVATIONS'][OBSERVATIONSNames[i]].append(aa)
            else:            
                o['OBSERVATIONS'][OBSERVATIONSNames[i]].append(v)
    return o

if __name__ == "__main__":
    #a = '/archive/u/kwargan/data/ozone_sondes/woudc2018/20181214.ecc.z.z32286.fmi-smna.csv'
    a = '/Users/bkarpowi/20180928.Brewer.MKII.019.MSC.csv'
    o = readWoudc(a)
    print (o['OBSERVATIONS'])
    print (o['LOCATION'])
    print (o['TIMESTAMP'])
    print(o.keys())

