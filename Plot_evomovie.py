import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, ScalarFormatter, LogLocator
plt.rcParams['axes.linewidth'] = 1.5
import glob, os
import sys

files = glob.glob("structure*.dat")
files = sorted(files)
Narg  = len(sys.argv)
first = 1
last  = len(files)
step  = 1
if (Narg>1): first=int(sys.argv[1])
if (Narg>2):  last=int(sys.argv[2])
if (Narg>3):  step=int(sys.argv[3])
print first,last,step

bar = 1.E+6                    # 1 bar in dyn/cm2 
day = 3600.*24.
yr  = 365.25*day

data = open('history2.out')
titel = data.readline()
data.close()
solids = titel.split()
NDUST = len(solids)
dat = np.loadtxt('history2.out',skiprows=1)
species = []
Ncond = 0
for i in range(4,NDUST):
  if (np.max(dat[:,i])>0.0):
    name = solids[i]
    if ('[s]' in name): name=name[:len(name)-3]
    species.append(name)
    Ncond = Ncond+1
species = np.array(species)
print Ncond,species
Nmax = np.log10(np.max(dat[:,4:]))+0.1
Nmin = Nmax-7
print "Nmax=",Nmax

nnn = first
while (nnn<=last):
  file = files[nnn-1]
  nnn += step
  print file

  data   = open(file)
  titel  = data.readline()
  dimens = data.readline()
  dimens = np.array(dimens.split())
  NELEM  = int(dimens[0])
  NMOLE  = int(dimens[1])
  NDUST  = int(dimens[2])
  NPOINT = int(dimens[3])
  depth  = data.readline()
  crust_thick = float(depth.split()[0])
  print "crust thickness = ",crust_thick
  crust_Ncond = 1.E-20*np.ones(Ncond)
  for i in range(0,NDUST):
    dat = data.readline().split()
    if (float(dat[1])>0.0):
      ind = np.where(dat[0] == species)[0]
      print ind,dat[0]
      crust_Ncond[ind[0]] = float(dat[1])
  header = data.readline()
  keyword = np.array(header.split())
  data.close()

  dat = np.loadtxt(file,skiprows=4+NDUST)
  NPOINT = len(dat[0:])
  time = np.float(titel.split()[1])
  out = "t ={:10.3e} days".format(time/day)
  print out

  Tg    = dat[:,0]                 # T [K]
  nHtot = dat[:,1]                 # n<H> [cm-3]
  press = dat[:,2]/bar             # p [bar]
  Diff  = dat[:,3]                 # diffusion coefficient cm2/s
  lognH = np.log10(nHtot)          
  lp    = np.log10(press)
  Tsurf = Tg[0]
  pmin  = np.min(lp)
  pmax  = np.max(lp)
  iii   = np.where((lp>pmin) & (lp<pmax))[0]
  colo = ['blue','black','silver','red','darkorange','gold','darkorchid','aqua','cadetblue','darkkhaki','pink','moccasin','cornflowerblue','chartreuse','limegreen','darkgreen','chocolate','darkgoldenrod']
  #'darkolivegreen','darkmagenta','aquamarine','coral','burlywood',
  #'beige','darkorange','crimson','darkcyan','bisque'
  Ncolor = len(colo)
  colo = colo*10
  styl = ['-']*Ncolor + ['--']*Ncolor + [':']*Ncolor + ['-.']*Ncolor*7 
  widt = [2]*Ncolor*10

  #================ the gas phase element abundances ===================
  fig,ax = plt.subplots(1,2,figsize=(13,7),sharex=False)
  #ax[0].subplots_adjust(right=0.25)
  count = 0
  for i in range(5+NELEM+NMOLE+NDUST,5+NELEM+NMOLE+NDUST+NELEM,1):
    elm = keyword[i]
    element = elm[3:]
    yy = dat[:,i]               # log10 eps
    ax[0].plot(lp,yy,c=colo[count],ls=styl[count],lw=widt[count],label=element)
    count = count+1
  ax[0].set_xlabel(r'$\log_{10}\ p\ \mathrm{[bar]}$',fontsize=20)
  ax[0].set_ylabel(r'$\log\,\epsilon_{\rm gas}$',fontsize=20)
  ax[0].set_xlim(pmin,pmax)
  ax[0].set_ylim(-20,0.3)
  ax[0].tick_params(axis='both', labelsize=15)
  ax[0].tick_params('both', length=6, width=1.5, which='major')
  ax[0].tick_params('both', length=3, width=1, which='minor')
  leg = ax[0].legend(loc='lower left',fontsize=15,fancybox=True)
  leg.get_frame().set_alpha(0.7)
  titel2 = "time =%9.2e s" %(time)
  ax[0].set_title(titel2,fontsize=18)
  #================== barplot crust column densities ====================
  xpos = np.arange(Ncond)
  print Ncond,xpos
  print species
  print crust_Ncond
  ax[1].bar(xpos, np.log10(crust_Ncond), align='center', alpha=0.5, color=colo)
  ax[1].set_xlim(-0.6,Ncond-0.4)
  ax[1].set_ylim(Nmin,Nmax)
  ax[1].set_xticks(xpos)
  ax[1].set_xticklabels(species, rotation=70,fontsize=16)
  ax[1].set_ylabel(r'$N_{\rm cond}\ \rm[cm^{-2}]$',fontsize=20)
  titel2 = "Tsurf =%8.2f K,  thickness =%10.3e cm" %(Tsurf,crust_thick)
  ax[1].set_title(titel2,fontsize=18)

  plt.tight_layout()

  print file
  i1 = file.index('structure_')+10
  i2 = file.index('.dat')
  print file[i1:i2]
  no = file[i1:i2]
  outfile = 'frame_'+no+'.png'
  fig.savefig(outfile)
  plt.clf()
  plt.close()
  print '... written output to '+outfile
