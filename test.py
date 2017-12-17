import glob, os
import sys

files = glob.glob("output/weather*.dat")
files = sorted(files)
Narg  = len(sys.argv)
first = 1
last  = len(files) 
if (Narg>1): first=int(sys.argv[1])
if (Narg>2):  last=int(sys.argv[Narg-1])
print first,last

for file in files[first-1:last]:
  print file
