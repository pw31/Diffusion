import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, ScalarFormatter, LogLocator
from matplotlib.backends.backend_pdf import PdfPages
plt.rcParams['axes.linewidth'] = 1.5
pp = PdfPages('tdep.pdf')
import glob, os
import sys

colour = ['blue','black','silver','red','darkorange','gold','darkorchid','aqua','cadetblue','darkkhaki','pink','moccasin','cornflowerblue','chartreuse','limegreen','darkgreen','chocolate','darkgoldenrod','darkolivegreen','darkmagenta','aquamarine','coral','burlywood','beige','darkorange','crimson','darkcyan','bisque']
Ncolour = len(colour)

#================== element col density in atmosphere ====================
file = "history3.out"
data = open(file)
header3 = data.readline()
data.close()
keyword3 = np.array(header3.split())
Ncol3 = len(keyword3)
dat3 = np.loadtxt(file,skiprows=2)
print file," no of columns=",Ncol3

indt = np.where(keyword3=='time[s]')[0][0]
t3 = dat3[:,indt]
yr = 365.25*24*3600. 
Myr = 1.E+6*yr
tmin = np.min(t3/yr)
tmax = np.max(t3/yr)
Narg  = len(sys.argv)
#print Narg,sys.argv[1],sys.argv[2]
print tmin,tmax
if (Narg>1): tmin=10**np.float(sys.argv[1])
if (Narg>2): tmax=10**np.float(sys.argv[2])
print tmin,tmax
iii = np.where((t3/yr>tmin) & (t3/yr<tmax))[0]

Nmax = np.max(dat3[:,4:])
print "maximum atmosphere column density [cm-2]",Nmax
col = 0
for i in range(4,Ncol3):
  Natmos = dat3[:,i]
  if (np.max(Natmos)>Nmax*1.E-10):
    lab = keyword3[i]
    plt.plot(np.log10(t3[iii]/yr),np.log10(Natmos[iii]),label=lab,c=colour[col],lw=3)
    col = col+1
plt.title(r'atmosphere',fontsize=26)
plt.xlabel(r'$\log_{10}\ t\ \mathrm{[yr]}$',fontsize=20)
plt.ylabel(r'$\log_{10}\ N\ \mathrm{[cm^{-2}]}$',fontsize=20)
plt.xlim(np.log10(tmin),np.log10(tmax))
plt.ylim(np.log10(Nmax)-10,np.log10(Nmax)+0.5)
plt.tight_layout()
plt.legend(loc='lower center',fontsize=10)
#plt.yscale('log')
plt.savefig(pp,format='pdf')
plt.clf()

#================== molecular col density in atmosphere ====================
file = "history4.out"
data = open(file)
header4 = data.readline()
data.close()
keyword4 = np.array(header4.split())
Ncol4 = len(keyword4)
dat4 = np.loadtxt(file,skiprows=2)
t4 = dat4[:,indt]
iii = np.where((t4/yr>tmin) & (t4/yr<tmax))[0]
print file," no of columns=",Ncol4
Nmax = np.max(dat4[:,4:])
print "maximum atmosphere column density [cm-2]",Nmax
col = 0
for i in range(4,Ncol4):
  Natmos = dat4[:,i]
  if (np.max(Natmos)>Nmax*3.E-3):
    lab = keyword4[i]
    #print lab,iii
    plt.plot(np.log10(t4[iii]/yr),np.log10(Natmos[iii]),label=lab,c=colour[col],lw=3)
    col = col+1
plt.title(r'atmosphere',fontsize=26)
plt.xlabel(r'$\log_{10}\ t\ \mathrm{[yr]}$',fontsize=20)
plt.ylabel(r'$\log_{10}\ N\ \mathrm{[cm^{-2}]}$',fontsize=20)
plt.xlim(np.log10(tmin),np.log10(tmax))
plt.ylim(np.log10(Nmax)-5,np.log10(Nmax)+0.2)
plt.tight_layout()
plt.legend(loc='lower center',fontsize=10)
#plt.yscale('log')
plt.savefig(pp,format='pdf')
plt.clf()

#================== col density in crust ====================
file = "history1.out"
data = open(file)
header1 = data.readline()
data.close()
keyword1 = np.array(header1.split())
Ncol1 = len(keyword1)
dat1 = np.loadtxt(file,skiprows=2)
t1 = dat1[:,indt]
iii = np.where((t1/yr>tmin) & (t1/yr<tmax))[0]
print file," no of columns=",Ncol1
Nmax = np.max(dat1[:,4:])
print "maximum crust column density [cm-2]",Nmax
col = 0
for i in range(4,Ncol1):
  Nsolid = dat1[:,i]
  if (np.max(Nsolid)>Nmax*1.E-5):
    Nsolid[np.where(Nsolid==0)]=1.E-99
    lab = keyword1[i]
    plt.plot(np.log10(t1[iii]/yr),np.log10(Nsolid[iii]),label=lab,c=colour[col],lw=3)
    col = col+1
plt.title(r'crust',fontsize=26)
plt.xlabel(r'$\log_{10}\ t\ \mathrm{[yr]}$',fontsize=20)
plt.ylabel(r'$\log_{10}\ N\ \mathrm{[cm^{-2}]}$',fontsize=20)
plt.xlim(np.log10(tmin),np.log10(tmax))
plt.ylim(np.log10(Nmax)-5,np.log10(Nmax)+0.5)
plt.tight_layout()
plt.legend(loc='lower center',fontsize=10)
#plt.yscale('log')
plt.savefig(pp,format='pdf')
plt.clf()

#================== material crust col density ====================
file = "history2.out"
data = open(file)
header2 = data.readline()
data.close()
keyword2 = np.array(header2.split())
Ncol2 = len(keyword2)
dat2 = np.loadtxt(file,skiprows=2)
t2 = dat2[:,indt]
iii = np.where((t2/yr>tmin) & (t2/yr<tmax))[0]
print file," no of columns=",Ncol2
Nmax = np.max(dat2[:,4:])
col = 0
for i in range(4,Ncol2):
  Nsolid = dat2[:,i]
  if (np.max(Nsolid)>Nmax*1.E-5):  
    Nsolid[np.where(Nsolid==0)]=1.E-99
    lab = keyword2[i]
    plt.plot(np.log10(t2[iii]/yr),np.log10(Nsolid[iii]),label=lab,c=colour[col],lw=3)
    col = col+1
plt.title(r'crust',fontsize=26)
plt.xlabel(r'$\log_{10}\ t\ \mathrm{[yr]}$',fontsize=20)
plt.ylabel(r'$\log_{10}\ N\ \mathrm{[cm^{-2}]}$',fontsize=20)
plt.xlim(np.log10(tmin),np.log10(tmax))
plt.ylim(np.log10(Nmax)-5,np.log10(Nmax)+0.5)
plt.tight_layout()
plt.legend(loc='lower center',fontsize=10)
#plt.yscale('log')
plt.savefig(pp,format='pdf')
plt.clf()

pp.close()
print '... written output to tdep.pdf.'


