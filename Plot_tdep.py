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

file  = "col.out"
print file
data   = open(file)
header = data.readline()
data.close()
keyword = np.array(header.split())
dat = np.loadtxt(file,skiprows=1)

indt = np.where(keyword=='time')[0][0]
t = dat[:,indt]
yr = 365.25*24*3600. 
tmin = np.min(t/yr)
tmax = np.max(t/yr)
Narg  = len(sys.argv)
#print Narg,sys.argv[1],sys.argv[2]
print tmin,tmax
if (Narg>1): tmin=np.float(sys.argv[1])
if (Narg>2): tmax=np.float(sys.argv[2])
print tmin,tmax
iii = np.where((t/yr>tmin) & (t/yr<tmax))[0]

#colo = ['blue','black','silver','red','darkorange','gold','darkorchid','aqua','cadetblue','darkkhaki','pink','moccasin','cornflowerblue','chartreuse','limegreen','darkgreen','chocolate','darkgoldenrod']
#'darkolivegreen','darkmagenta','aquamarine','coral','burlywood',
#'beige','darkorange','crimson','darkcyan','bisque'

#================== col density ====================

colours = ('Navy','DeepSkyBlue','Orange','Orange','Maroon','Maroon','Green','Green','Goldenrod','Purple','Pink','Pink')
style = ('solid','dotted','dotted','dashed','dashed','dashdot','dashdot','dashed','dashdot','solid','dashed','dotted')
material = ('N_TiO2','N_Al2O3','N_MgSiO3','N_Mg2SiO4','N_SiO', 'N_SiO2','N_Fe','N_FeO','N_MgO', 'N_KCl','N_NaCl','N_Na2S')
label = (r'$N_{TiO_2}$',r'$N_{Al_{2}O_{3}}$',r'$N_{MgSiO_{3}}$',r'$N_{MgSiO_{4}}$',r'$N_{SiO}$',r'$N_{SiO_2}$',r'$N_{Fe}$',r'$N_{FeO}$',r'$N_{MgO}$', r'$N_{KCl}$',r'$N_{NaCl}$',r'$N_{Na_2S}$')

for l,m,c,lin in zip(label,material,colours,style):
  fig,ax = plt.subplots()
  indn = np.where(keyword==m)[0]
  if (len(indn)==1):
    n=dat[:,indn[0]]
    lab = l+" $\mathrm{[10^{-5}g\,cm^{-2}]}$" 
    plt.plot(t[iii]/yr,n[iii]/1.E-5,label=l,color=c,linestyle=lin,linewidth=3)
    plt.xlabel(r'$time\ \mathrm{[yrs]}$',fontsize=20)
    plt.ylabel(lab,fontsize=20)
    plt.xlim(tmin,tmax)
    #plt.legend(loc='lower right')
    plt.tight_layout()
    plt.savefig(pp,format='pdf')
    plt.clf()

fig,ax = plt.subplots()
for l,m,c,lin in zip(label,material,colours,style):
  indn = np.where(keyword==m)[0]
  if (len(indn)==1):
    n=dat[:,indn[0]]
    plt.plot(t[iii]/yr,n[iii]/1.E-5,label=l,color=c,linestyle=lin,linewidth=3)
plt.xlabel(r'$time\ \mathrm{[yrs]}$',fontsize=20)
plt.ylabel(r'$\mathrm{col.dens.\ [10^{-5}g\,cm^{-2}]}$',fontsize=20)
plt.xlim(tmin,tmax)
#plt.yscale('log')
plt.xlabel(r'$time [s]$',fontsize=20)
plt.legend(loc='upper right')
plt.tight_layout()
plt.savefig(pp,format='pdf')
plt.clf()

pp.close()
print '... written output to tdep.pdf.'


