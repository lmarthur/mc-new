import os
from natsort import natsorted
import sys

dir_name = sys.argv[1]
param = sys.argv[2]
allFile = dir_name+"/all_aero_coefs.dat"
first = True

print("Parameter is "+param)

for run in natsorted(os.listdir(dir_name)):
    if "." in run:
        continue
    
    coeffFile = os.path.join(os.getcwd(), dir_name, run, "MN_aero_coefs.dat")
    
    f = open(coeffFile, 'r')
    if first:
        lines = f.readlines()
        header = lines[3]
        coeffs = lines[4]
    else:
        coeffs = f.readlines()[4]
    f.close()

    inputFile = os.path.join(os.getcwd(), dir_name, run, run+".flap")

    f = open(inputFile, 'r')
    lines = f.readlines()

    for linenum, line in enumerate(lines):
        if param in line:
            paramLine = line
            break
    f.close()

    newParamLine = paramLine.split()
    value = newParamLine[0]
    
    if first:
        a = open(allFile, 'w')
        a.write("# Name            "+param+"          "+header[2:])
        a.close()
        first = False
 
    a = open(allFile, 'a')
    a.write(run+"     "+value+"     "+coeffs)
    a.close()
        
