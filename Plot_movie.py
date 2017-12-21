import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, ScalarFormatter, LogLocator
plt.rcParams['axes.linewidth'] = 1.5
import glob, os
import sys

files = glob.glob("weather*.dat")
files = sorted(files)
Narg  = len(sys.argv)
first = 1
last  = len(files)
step  = 1
if (Narg>1): first=int(sys.argv[1])
if (Narg>2):  last=int(sys.argv[2])
if (Narg>3):  step=int(sys.argv[3])
print first,last,step

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
  NEPS   = int(dimens[3])
  NNUC   = int(dimens[4])
  NPOINT = int(dimens[5])
  header = data.readline()
  data.close()

  dat = np.loadtxt(file,skiprows=3)
  keyword = np.array(header.split())
  NPOINT = len(dat[0:])

  bar   = 1.E+6                    # 1 bar in dyn/cm2 
  Tg    = dat[:,0]                 # T [K]
  nHtot = dat[:,1]                 # n<H> [cm-3]
  lognH = np.log10(nHtot)          
  press = dat[:,2]/bar             # p [bar]
  lp    = np.log10(press)
  pmin  = np.min(lp)
  pmax  = np.max(lp)
  #pmin  = -5.0
  #pmax  = 1.0
  iii   = np.where((lp>pmin) & (lp<pmax))[0]
  Tmin  = np.min(Tg[iii])
  Tmax  = np.max(Tg[iii])
  Tmin  = Tmin*0.9
  Tmax  = Tmax*1.1
  nHmin = np.min(nHtot[iii])
  nHmax = np.max(nHtot[iii])
  nHmin = nHmin*0.9
  nHmax = nHmax*1.1
  if (nHmax>nHmin*5): 
    nHmin = nHmin/2.0
    nHmax = nHmax*2.0
  #sep = 20
  colo = ['blue','black','silver','red','darkorange','gold','darkorchid','aqua','cadetblue','darkkhaki','pink','moccasin','cornflowerblue','chartreuse','limegreen','darkgreen','chocolate','darkgoldenrod']
  #'darkolivegreen','darkmagenta','aquamarine','coral','burlywood',
  #'beige','darkorange','crimson','darkcyan','bisque'
  Ncolor = len(colo)
  colo = colo*10
  styl = ['-']*Ncolor + ['--']*Ncolor + [':']*Ncolor + ['-.']*Ncolor*7 
  widt = [2]*Ncolor*10

  #================== solid particle densities ===================
  fig,ax = plt.subplots()
  count = 0
  for i in range(5+NELEM+NMOLE+NDUST,5+NELEM+NMOLE+2*NDUST,1):
    solid = keyword[i]
    yy = dat[:,i]               # log10 nsolid/n<H>
    plt.plot(lp[iii],yy[iii],c=colo[count],ls=styl[count],lw=widt[count],label=solid)
    count = count+1

  time = np.float(titel.split()[1])
  day = 3600.*24.
  out = "t ={:10.5f} days".format(time/day)
  plt.title(out,fontsize=15)
  plt.xlabel(r'$\log_{10}\ p\ \mathrm{[bar]}$',fontsize=20)
  plt.ylabel(r'$\mathrm{log}_{10}\ n_\mathrm{solid}/n_\mathrm{\langle H\rangle}$',fontsize=20)
  plt.xlim(pmin,pmax)
  plt.ylim(-16,-3)
  plt.tick_params(axis='both', labelsize=14)
  plt.tick_params('both', length=6, width=1.5, which='major')
  plt.tick_params('both', length=3, width=1, which='minor')
  sz = np.min([11,1+195.0/count])
  leg = plt.legend(loc='lower right',fontsize=11,fancybox=True,
                   handlelength=2.5,prop={'size':sz})
  leg.get_frame().set_alpha(0.7)
  plt.tight_layout()
  ii = file.index('weather_')+8
  print file[ii:ii+5]
  no = file[ii:ii+5]
  outfile = 'nsolid_'+no+'.png'
  fig.savefig(outfile)
  plt.clf()
  plt.close()
  print '... written output to '+outfile

  #================ the gas phase element abundances ===================
  fig,ax = plt.subplots()
  count = 0
  for i in range(5+NELEM+NMOLE+2*NDUST,5+NELEM+NMOLE+2*NDUST+NELEM,1):
    elm = keyword[i]
    element = elm[3:]
    yy = dat[:,i]               # log10 eps
    plt.plot(lp,yy,c=colo[count],ls=styl[count],lw=widt[count],label=element)
    count = count+1
  plt.title(out,fontsize=15)
  plt.xlabel(r'$\log_{10}\ p\ \mathrm{[bar]}$',fontsize=20)
  plt.ylabel(r'$\log\,\epsilon_{\rm gas}$',fontsize=20)
  plt.xlim(pmin,pmax)
  plt.ylim(-20,0.3)
  plt.tick_params(axis='both', labelsize=15)
  plt.tick_params('both', length=6, width=1.5, which='major')
  plt.tick_params('both', length=3, width=1, which='minor')
  sz = np.min([11,1+195.0/count])
  leg = plt.legend(loc='lower right',fontsize=sz,fancybox=True)
  leg.get_frame().set_alpha(0.7)
  plt.tight_layout()
  outfile = 'eps_'+no+'.png'
  fig.savefig(outfile)
  plt.clf()
  plt.close()
  print '... written output to '+outfile

  #================== nucleation rates ... =====================
  fig,ax = plt.subplots(3,figsize=(8,9),sharex=True)
  ymax = -100.0
  count = 3
  for i in range(8+2*NELEM+NMOLE+2*NDUST,8+2*NELEM+NMOLE+2*NDUST+NNUC,1):
    log10_Jstar = dat[:,i]
    ym = np.max(log10_Jstar)
    ax[0].plot(lp,log10_Jstar,lw=4,c=colo[count],label=keyword[i])
    count = count+1
  ax[0].set_title(out,fontsize=15)
  #ax[0].set_xlabel(r'$\log_{10}\ p\ \mathrm{[bar]}$',fontsize=20)
  ax[0].set_ylabel(r'$\log_{10}\ J_{\!\star}\ \mathrm{[cm^{-3}s^{-1}]}$',fontsize=20)
  #ax[0].set_xlim(pmin,pmax)
  ax[0].set_ylim(-19.5,0)
  ax[0].tick_params(axis='both', labelsize=15)
  ax[0].tick_params('both', length=6, width=1.5, which='major')
  ax[0].tick_params('both', length=3, width=1, which='minor')
  leg = ax[0].legend(loc='upper right',fontsize=11,fancybox=True,
                     handlelength=2.5)
  leg.get_frame().set_alpha(0.7)
  #================ ... dust/gas mass ratio ... ===============
  ind = np.where(keyword=='dust/gas')[0][0]
  log10_dust_gas = dat[:,ind]
  ax[1].plot(lp,10**log10_dust_gas,lw=4,label='dust/gas')
  #ax[1].set_xlabel(r'$\log_{10}\ p\ \mathrm{[bar]}$',fontsize=20)
  ax[1].set_ylabel(r'$\mathrm{dust/gas}$',fontsize=20)
  #ax[1].set_xlim(pmin,pmax)
  ax[1].set_ylim(1.1E-9,1.E-2)
  ax[1].set_yscale('log')
  ax[1].tick_params(axis='both', labelsize=15)
  ax[1].tick_params('both', length=6, width=1.5, which='major')
  ax[1].tick_params('both', length=3, width=1, which='minor')
  #===================== ... mean size ========================
  ind = np.where(keyword=='<a>[mic]')[0][0]
  amean = dat[:,ind]
  ax[2].plot(lp,amean,lw=4,label='<a>')
  ax[2].set_xlabel(r'$\log_{10}\ p\ \mathrm{[bar]}$',fontsize=20)
  ax[2].set_ylabel(r'$\langle a\rangle\ \mathrm{[\mu m]}$',fontsize=20)
  ax[2].set_xlim(pmin,pmax)
  ax[2].set_ylim(1.E-3,1.E+2)
  ax[2].set_yscale('log')
  ax[2].tick_params(axis='both', labelsize=15)
  ax[2].tick_params('both', length=6, width=1.5, which='major')
  ax[2].tick_params('both', length=3, width=1, which='minor')
  plt.tight_layout()
  plt.subplots_adjust(wspace=0, hspace=0)
  outfile = 'dust_'+no+'.png'
  fig.savefig(outfile)
  plt.clf()
  plt.close()
  print '... written output to '+outfile
