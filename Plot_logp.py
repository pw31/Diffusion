import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, ScalarFormatter, LogLocator
from matplotlib.backends.backend_pdf import PdfPages
plt.rcParams['axes.linewidth'] = 1.5
pp = PdfPages('weather.pdf')
import glob, os
import sys

files = glob.glob("weather*.dat")
files = sorted(files)
Narg  = len(sys.argv)
last  = len(files)-1 
if (Narg>1):  last=int(sys.argv[1])-1
print Narg,last,len(files)
file  = files[last]
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
time = np.float(titel.split()[1])
day = 3600.*24.
out = "t ={:10.5f} days".format(time/day)

bar   = 1.E+6                    # 1 bar in dyn/cm2 
Tg    = dat[:,0]                 # T [K]
nHtot = dat[:,1]                 # n<H> [cm-3]
press = dat[:,2]/bar             # p [bar]
Diff  = dat[:,3]                 # diffusion coefficient cm2/s
lognH = np.log10(nHtot)          
lp    = np.log10(press)
pmin  = np.min(lp)
pmax  = np.max(lp)
#pmin  = -5.3
#pmax  = +2.1
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
#if (Tmax-Tmin>1500): sep=100
#if (Tmax-Tmin>1000): sep=50
#if (Tmax-Tmin<600): sep=20
#if (Tmax-Tmin<400): sep=10
  colo = ['blue','black','silver','red','darkorange','gold','darkorchid','aqua','cadetblue','darkkhaki','pink','moccasin','cornflowerblue','chartreuse','limegreen','darkgreen','chocolate','darkgoldenrod']
#'darkolivegreen','darkmagenta','aquamarine','coral','burlywood',
#'beige','darkorange','crimson','darkcyan','bisque'
Ncolor = len(colo)
colo = colo*10
styl = ['-']*Ncolor + ['--']*Ncolor + [':']*Ncolor + ['-.']*Ncolor*7 
widt = [2]*Ncolor*10

#================== temperature-pressure structure ====================
fig,ax = plt.subplots()
plt.plot(lp,Tg,lw=4)
plt.title(out,fontsize=15)
plt.xlabel(r'$\log_{10}\ p\ \mathrm{[bar]}$',fontsize=20)
plt.ylabel(r'$T\ \mathrm{[K]}$',fontsize=20)
plt.xlim(pmin,pmax)
plt.ylim(Tmin,Tmax)
plt.tick_params(axis='both', labelsize=15)
plt.tick_params('both', length=6, width=1.5, which='major')
plt.tick_params('both', length=3, width=1, which='minor')
#minorLocator = MultipleLocator(sep)
#ax.xaxis.set_minor_locator(minorLocator)
plt.tight_layout()
plt.savefig(pp,format='pdf')
plt.clf()

#================== density-pressure structure ====================
fig,ax = plt.subplots()
plt.plot(lp,nHtot,lw=4)
plt.xlabel(r'$\log_{10}\ p\ \mathrm{[bar]}$',fontsize=20)
plt.ylabel(r'$n_\mathrm{\langle H\rangle}\ \mathrm{[cm^{-3}]}$',fontsize=20)
plt.xlim(pmin,pmax)
plt.ylim(nHmin,nHmax)
if (nHmax>nHmin*5): plt.yscale('log')
plt.tick_params(axis='both', labelsize=15)
plt.tick_params('both', length=6, width=1.5, which='major')
plt.tick_params('both', length=3, width=1, which='minor')
#minorLocator = MultipleLocator(sep)
#ax.xaxis.set_minor_locator(minorLocator)
plt.tight_layout()
plt.savefig(pp,format='pdf')
plt.clf()

#================== the diffusion coefficient ====================
fig,ax = plt.subplots()
plt.plot(lp,np.log10(Diff),lw=4)
plt.xlabel(r'$\log_{10}\ p\ \mathrm{[bar]}$',fontsize=20)
plt.ylabel(r'$\log_{10}\ D\ \mathrm{[cm^{2}s^{-1}]}$',fontsize=20)
plt.xlim(pmin,pmax)
plt.tick_params(axis='both', labelsize=15)
plt.tick_params('both', length=6, width=1.5, which='major')
plt.tick_params('both', length=3, width=1, which='minor')
plt.tight_layout()
plt.savefig(pp,format='pdf')
plt.clf()

#================== the nucleation rates ====================
fig,ax = plt.subplots()
ymax = -100.0
count = 3
for i in range(8+2*NELEM+NMOLE+2*NDUST,8+2*NELEM+NMOLE+2*NDUST+NNUC,1):
  log10_Jstar = dat[:,i]
  ym = np.max(log10_Jstar)
  plt.plot(lp,log10_Jstar,c=colo[count],lw=4,label=keyword[i])
  count = count+1
  ymax = np.max([ymax,ym])
plt.xlabel(r'$\log_{10}\ p\ \mathrm{[bar]}$',fontsize=20)
plt.ylabel(r'$\log_{10}\ J_{\!\star}\ \mathrm{[cm^{-3}s^{-1}]}$',fontsize=20)
plt.xlim(pmin,pmax)
plt.ylim(ymax-15,ymax+1)
plt.tick_params(axis='both', labelsize=15)
plt.tick_params('both', length=6, width=1.5, which='major')
plt.tick_params('both', length=3, width=1, which='minor')
leg = plt.legend(loc='upper right',fontsize=11,fancybox=True,handlelength=2.5)
leg.get_frame().set_alpha(0.7)
plt.tight_layout()
plt.savefig(pp,format='pdf')
plt.clf()

#================== the dust/gas mass ratio ====================
fig,ax = plt.subplots()
ind = np.where(keyword=='dust/gas')[0][0]
log10_dust_gas = dat[:,ind]
ymax = np.max(log10_dust_gas)
plt.plot(lp,10**log10_dust_gas,lw=4)
plt.xlabel(r'$\log_{10}\ p\ \mathrm{[bar]}$',fontsize=20)
plt.ylabel(r'$\mathrm{dust/gas}$',fontsize=20)
plt.xlim(pmin,pmax)
plt.ylim(10**(ymax-10),10**ymax*3)
plt.yscale('log')
plt.tick_params(axis='both', labelsize=15)
plt.tick_params('both', length=6, width=1.5, which='major')
plt.tick_params('both', length=3, width=1, which='minor')
ax.yaxis.set_minor_locator(LogLocator(subs=[2,3,4,5,6,7,8,9]))
plt.tight_layout()
plt.savefig(pp,format='pdf')
plt.clf()

#================== mean size ====================
fig,ax = plt.subplots()
ind = np.where(keyword=='<a>[mic]')[0][0]
amean = dat[:,ind]
ymax = np.max(amean)
if (ymax>0):
  plt.plot(lp,amean,lw=4)
  plt.xlabel(r'$\log_{10}\ p\ \mathrm{[bar]}$',fontsize=20)
  plt.ylabel(r'$\langle a\rangle\ \mathrm{[\mu m]}$',fontsize=20)
  plt.xlim(pmin,pmax)
  plt.ylim(1.E-3,ymax*3)
  plt.yscale('log')
  plt.tick_params(axis='both', labelsize=15)
  plt.tick_params('both', length=6, width=1.5, which='major')
  plt.tick_params('both', length=3, width=1, which='minor')
  ax.yaxis.set_minor_locator(LogLocator(subs=[2,3,4,5,6,7,8,9]))
  #minorLocator = MultipleLocator(sep)
  #ax.xaxis.set_minor_locator(minorLocator)
  plt.tight_layout()
  plt.savefig(pp,format='pdf')
  plt.clf()

#================ the gas phase element abundances ===================
fig,ax = plt.subplots()
count = 0
ymax = -100.0
for i in range(5+NELEM+NMOLE+2*NDUST,5+NELEM+NMOLE+2*NDUST+NELEM,1):
  elm = keyword[i]
  element = elm[3:]
  yy = dat[:,i]               # log10 eps
  ymax=np.max([ymax,np.max(yy)])            
  plt.plot(lp,yy,c=colo[count],ls=styl[count],lw=widt[count],label=element)
  count = count+1
plt.xlabel(r'$\log_{10}\ p\ \mathrm{[bar]}$',fontsize=20)
plt.ylabel(r'$\log\,\epsilon_{\rm gas}$',fontsize=20)
plt.xlim(pmin,pmax)
plt.ylim(ymax-20,0.3)
plt.tick_params(axis='both', labelsize=15)
plt.tick_params('both', length=6, width=1.5, which='major')
plt.tick_params('both', length=3, width=1, which='minor')
#minorLocator = MultipleLocator(sep)
#ax.xaxis.set_minor_locator(minorLocator)
#minorLocator = MultipleLocator(1)
#ax.yaxis.set_minor_locator(minorLocator)
sz = np.min([11,1+195.0/count])
leg = plt.legend(loc='upper left',fontsize=sz,fancybox=True)
leg.get_frame().set_alpha(0.7)
plt.tight_layout()
plt.savefig(pp,format='pdf')
plt.clf()

#================== solid particle densities ===================
solids = []
ymax = -100.0
fig,ax = plt.subplots()
count = 0
nsolid = np.zeros(NPOINT)
i1 = 5+NELEM+NMOLE
i2 = 5+NELEM+NMOLE+NDUST
ymax = np.max(dat[iii,i1+NDUST:i2+NDUST])
ymin = ymax-8
ymax = ymax+0.3
for i in range(i1,i2,1):
  solid = keyword[i]
  solid = solid[1:]
  ind = np.where(keyword == 'n'+solid)[0]
  if (np.size(ind) == 0): continue
  ind = ind[0]
  yy = dat[:,ind]               # log10 nsolid/n<H>
  nsolid = nsolid + 10**yy
  #print solid,np.max(yy[iii]),ymax,ymin
  if (np.max(yy[iii])>ymin):
    plt.plot(lp[iii],yy[iii],c=colo[count],ls=styl[count],lw=widt[count],label=solid)
  count = count + 1
plt.title('condensates',fontsize=20)
plt.xlabel(r'$\log_{10}\ p\ \mathrm{[bar]}$',fontsize=20)
plt.ylabel(r'$\mathrm{log}_{10}\ n_\mathrm{solid}/n_\mathrm{\langle H\rangle}$',fontsize=20)
#plt.xscale('log')
plt.xlim(pmin,pmax)
plt.ylim(ymin,ymax)
plt.tick_params(axis='both', labelsize=14)
plt.tick_params('both', length=6, width=1.5, which='major')
plt.tick_params('both', length=3, width=1, which='minor')
sz = np.min([11,1+195.0/count])
leg = plt.legend(loc='upper left',fontsize=11,fancybox=True,
                 handlelength=2.5,prop={'size':sz})
leg.get_frame().set_alpha(0.7)
plt.tight_layout()
plt.savefig(pp,format='pdf')
plt.clf()

#================== solid composition ===================
ipos = np.where(nsolid>1.E-99)
fig,ax = plt.subplots()
count = 0
for i in range(5+NELEM+NMOLE,5+NELEM+NMOLE+NDUST,1):
  solid = keyword[i]
  solid = solid[1:]
  ind = np.where(keyword == 'n'+solid)[0]
  if (np.size(ind) == 0): continue
  ind = ind[0]
  yy = 10**dat[:,ind]/nsolid               # log10 nsolid/ntotsolid
  if (np.max(yy[ipos])>1.E-10):
    plt.plot(lp[ipos],yy[ipos],c=colo[count],ls=styl[count],lw=widt[count],label=solid)
  count = count + 1
plt.title('condensates',fontsize=20)
plt.xlabel(r'$\log_{10}\ p\ \mathrm{[bar]}$',fontsize=20)
plt.ylabel(r'$\mathrm{log}_{10}\ b_\mathrm{mix}$',fontsize=20)
plt.yscale('log')
plt.xlim(pmin,pmax)
plt.ylim(10**(-10.0),10**(0.3))
plt.tick_params(axis='both', labelsize=14)
plt.tick_params('both', length=6, width=1.5, which='major')
plt.tick_params('both', length=3, width=1, which='minor')
sz = np.min([11,1+195.0/count])
leg = plt.legend(loc='lower left',fontsize=11,fancybox=True,
                 handlelength=2.5,prop={'size':sz})
leg.get_frame().set_alpha(0.7)
plt.tight_layout()
plt.savefig(pp,format='pdf')
plt.clf()

#================== effective supersaturation ratios ===================
fig,ax = plt.subplots()
count = 0
for i in range(5+NELEM+NMOLE,5+NELEM+NMOLE+NDUST,1):
  solid = keyword[i]
  solid = solid[1:]
  ind = np.where(keyword == 'S'+solid)[0]
  if (np.size(ind) == 0): continue
  ind = ind[0]
  logS = dat[:,ind]              # log10 S
  plt.plot(lp,logS,c=colo[count],ls=styl[count],lw=widt[count],label=solid)
  count = count + 1
plt.title('effective supersaturation ratios',fontsize=20)
plt.xlabel(r'$\log_{10}\ p\ \mathrm{[bar]}$',fontsize=20)
plt.ylabel(r'$\mathrm{log}_{10}\ S$',fontsize=20)
plt.xlim(pmin,pmax)
plt.ylim(-10,10)
plt.tick_params(axis='both', labelsize=14)
plt.tick_params('both', length=6, width=1.5, which='major')
plt.tick_params('both', length=3, width=1, which='minor')
sz = np.min([13,1+195.0/count])
col = 1
if (count>30): 
  sz = np.min([13,1+195.0/count*2])
  col = 2
leg = plt.legend(loc='lower left',fontsize=10,fancybox=True,
                 handlelength=3,prop={'size':sz},ncol=col)
leg.get_frame().set_alpha(0.7)
plt.tight_layout()
plt.savefig(pp,format='pdf')
plt.clf()

#================== some important molecules ====================
fig,ax = plt.subplots()
mols  = ['H2','H','N2','H2O','O2','CO','CO2','CH4','NH3','C2H2','el']
mols  = np.array(mols)
ntot  = 0.0*nHtot
for i in range(4,5+NELEM+NMOLE): # electrons, all atoms, ions and cations
  ntot = ntot + 10**dat[:,i]
lntot = np.log10(ntot)
count = 0
for i in range(4,5+NELEM+NMOLE): 
  mol = keyword[i]
  yy = dat[:,i]-lntot            # log10 nmol/ntot
  crit = -1.5
  ind = np.where(mols == mol)[0]
  if (np.size(ind)>0): crit=-5
  #print i,mol,ind,np.size(ind)
  if (np.max(yy[iii])>crit):
    plt.plot(lp,yy,c=colo[count],ls=styl[count],lw=widt[count],label=mol)
    count = count + 1
plt.title('important molecules',fontsize=20)
plt.xlabel(r'$\log_{10}\ p\ \mathrm{[bar]}$',fontsize=20)
plt.ylabel(r'$\mathrm{log}_{10}\ n_\mathrm{mol}/n_\mathrm{tot}$',fontsize=20)
plt.xlim(pmin,pmax)
plt.ylim(-6.2,0.2)
plt.tick_params(axis='both', labelsize=14)
plt.tick_params('both', length=6, width=1.5, which='major')
plt.tick_params('both', length=3, width=1, which='minor')
leg = plt.legend(loc='lower right',fontsize=11,fancybox=True)
leg.get_frame().set_alpha(0.7)
plt.tight_layout()
plt.savefig(pp,format='pdf')
plt.clf()

#================== where are the elements? ================
ellist = ['H','C','O','N','SI','S','NA','CL','CA','TI','K','AL','MG','FE','LI','F','P','NI','MN','CR','ZN','ZR','RB','CU','B','BR','V','SR','W','el']
allist = [' ',' ',' ',' ','Si',' ','Na','Cl','Ca','Ti',' ','Al','Mg','Fe','Li',' ',' ','Ni','Mn','Cr','Zn','Zr','Rb','Cu',' ','Br',' ','Sr',' ','+']
exlist = [' He ',' Cl CL Ca CA Cr CR Co Cu CU ',' ',' Na NA Ni NI ',' ',' Si SI Sr SR ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' Fe FE ',' ',' ',' ',' ',' ',' ',' ',' ',' Br BR ',' ',' ',' ',' ',' ']
titels = ['hydrogen','carbon','oxygen','nitrogen','silicon','sulphur','sodium','chlorine','calcium','titanium','potassium','aluminum','magnesium','iron','lithium','fluorine','phosphorus','nickel','manganese','chromium','zinc','zirconium','rubidium','copper','boron','bromine','vanadium','strontium','tungston','charge carriers']
limits = [2,5,2.5,6,6,5,6,4,7,8,6,6,6,6,7,6,6,6,6,6,6,6,6,6,6,6,6,6,6,5]   
abund = np.zeros(5+NELEM+NMOLE, dtype=np.int)
for i in range(0,30):
  fig,ax = plt.subplots()
  el = ellist[i]
  al = allist[i]
  ex = exlist[i]
  limit = limits[i]
  titel = titels[i]
  print titel+" ..."
  nmax = np.float(-100)
  nmin = np.float(0)
  mollist = []
  abulist = []
  maxy = 0.0*dat[:,0]
  for mol in range(4,5+NELEM+NMOLE,1):
    molname = keyword[mol]
    ind = str.find(molname,el)
    if (ind < 0): 
      ind = str.find(molname,al)
    if (ind < 0 and el=='el'): 
      ind = str.find(molname,'-')
    if (ind >= 0):
      next1 = molname[ind:ind+2]
      next2 = molname[ind-1:ind+1]
      #print keyword[mol],next1,str.find(ex,next1),len(next1)
      if (len(next1)==1 or str.find(ex,next1)==-1 or molname=='SIS'):
        if (next2!='MN' and next2!='ZN'):
          yy = dat[:,mol]                # log10 nmol [cm-3]
          yy = yy - lognH                # log10 nmol/n<H>
          nmax = np.max([nmax,np.max(yy[iii])])
          maxy = maxy + 10**yy
          if (molname=='el'): nmin = np.min([nmin,np.min(yy[iii])])
          mollist.append(mol)   
          abulist.append(np.mean(yy))
  if (nmax==-100): continue
  indices = np.argsort(abulist)
  count = 0
  maxy = np.log10(maxy)
  nmin = np.min([nmin,np.min(maxy[iii])-limit,nmax-12])
  for ind in reversed(indices):
    mol = mollist[ind]
    abu = abulist[ind]
    molname = keyword[mol]
    yy = dat[:,mol]                # log10 nmol [cm-3]
    yy = yy - lognH                # log10 nmol/n<H>
    if (np.max(yy[iii]-maxy[iii])>-3 or molname=='el'):
      abund[mol] = 1
    if (np.max(yy[iii]-maxy[iii])>-limit or molname=='el'):
      print molname,abu
      plt.plot(lp,yy,c=colo[count],ls=styl[count],lw=widt[count],label=molname)
      count = count + 1
  plt.title(titel,fontsize=20)
  plt.xlabel(r'$\log_{10}\ p\ \mathrm{[bar]}$',fontsize=20)
  plt.ylabel(r'$\mathrm{log}_{10}\ n_\mathrm{mol}/n_\mathrm{\langle H\rangle}$',fontsize=20)
  plt.xlim(pmin,pmax)
  plt.ylim(nmin,nmax+1)
  plt.tick_params(axis='both', labelsize=14)
  plt.tick_params('both', length=6, width=1.5, which='major')
  plt.tick_params('both', length=3, width=1, which='minor')
  if (Tmax/Tmin>10):
    plt.xscale('log')
    #else:  
    #minorLocator = MultipleLocator(sep)
    #ax.xaxis.set_minor_locator(minorLocator)

  #minorLocator = MultipleLocator(1.0)
  #if (nmax-nmin>50): minorLocator = MultipleLocator(2.0)
  #if (nmax-nmin>100): minorLocator = MultipleLocator(5.0)
  #if (nmax-nmin>200): minorLocator = MultipleLocator(10.0)
  #ax.yaxis.set_minor_locator(minorLocator)
  sz = np.min([11,1+195.0/count])
  col = 1
  if (count>30): 
    sz = np.min([9,1+195.0/count*2])
    col = 2
  leg = plt.legend(loc='lower left',fontsize=10,fancybox=True,
                   handlelength=3,prop={'size':sz},ncol=col)
  leg.get_frame().set_alpha(0.7)
  plt.tight_layout()
  plt.savefig(pp,format='pdf')
  plt.clf()

print 'having not plotted these molecules ...'
Nnot = 0
out = ' '
for mol in range(4,5+NELEM+NMOLE,1):
  if (abund[mol]==0): 
    out = out + keyword[mol] + ' '
    Nnot = Nnot+1
print NMOLE," molecules,  ",Nnot," unimportant, namely ..."
print out

pp.close()
print '... written output to weather.pdf.'


