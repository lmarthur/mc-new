import os
from natsort import natsorted
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--dir",help="Directory to store all runs")
parser.add_argument("--param",help="Parameter to iterate over")
parser.add_argument("--start",type=float,help="Starting value to iterate from")
parser.add_argument("--end",type=float,help="Ending value of iteration")
parser.add_argument("--num",type=int,help="Number of steps to take during iteration")
args = vars(parser.parse_args())

dir_name = args["dir"]
param = args["param"]
start = float(args["start"])
end = int(args["end"])
num = int(args["num"])

step = (end - start) / num

for i, run in enumerate(natsorted(os.listdir(dir_name))):
    if "." in run:
        continue
    
    inputFile = os.path.join(os.getcwd(), dir_name, run, run+".flap")
    
    f = open(inputFile, 'r')
    lines = f.readlines()

    for linenum, line in enumerate(lines):
        if param in line:
            paramLine = line
            paramNum = linenum
            break
    f.close()

    newParamLine = paramLine.split()

    print(newParamLine)
   
    newParamLine[0] = '%.2f' % (start + step*i)

    newParamLine = "  " + newParamLine[0] + "     " + ' '.join(newParamLine[1:]) + "\n"

    print(newParamLine)
    
    lines[paramNum] = newParamLine
    w = open(inputFile, 'w')
    w.writelines(lines)
    w.close()

