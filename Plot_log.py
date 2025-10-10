import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, ScalarFormatter, LogLocator
from matplotlib.backends.backend_pdf import PdfPages
plt.rcParams['lines.linewidth'] = 3
plt.rcParams['axes.linewidth'] = 1.5
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['xtick.top'] = True
plt.rcParams['ytick.right'] = True 
plt.rcParams['xtick.labelsize'] = plt.rcParams['ytick.labelsize'] = 12
plt.rcParams['xtick.major.size'] = plt.rcParams['ytick.major.size'] = 7
plt.rcParams['xtick.minor.size'] = plt.rcParams['ytick.minor.size'] = 4
plt.rcParams['xtick.major.width'] = plt.rcParams['ytick.major.width'] = 1.6
plt.rcParams['font.size'] = 16
pp = PdfPages('out.pdf')

colo = ['gold','cadetblue','coral','blue','chartreuse','cornflowerblue','darkgray','darkgreen','red','darkorchid','aqua','burlywood','chocolate','darkkhaki','pink','moccasin']

file   = 'out.dat'
data   = open(file)
dimens = data.readline()
dimen  = np.array(dimens.split())
Np     = np.int(dimen[0])
init   = np.int(dimen[1])
zz     = np.array(data.readline().split(),dtype=np.float)
zunit  = 'cm'
tunit  = 's'
yr     = 365.25*24.0*3600.0 
if (np.max(zz)>1.E+5):
  zz    = zz/1.E+5
  zunit = 'km'
  tunit = 'yr'
#--- for analytic solution ---
N = 1
L = np.max(zz)-np.min(zz)
D = 1.0
k = np.pi/(N*L)
w = D*k**2

xmax = 0.0
for i in range(0,999):
  timedat = data.readline()
  print timedat
  if (len(timedat)==0): break
  tt = np.array(timedat.split()[1],dtype=np.float)
  xx = np.array(data.readline().split(),dtype=np.float)
  xmax = np.max([xmax,np.max(xx)])
  #================== x over z plot ====================
  if (tunit=='s'): titel = r'$t\rm = %9.2E$ s' % (tt)
  if (tunit=='yr'): titel = r'$t\rm = %9.2E$ yr' % (tt/yr)
  if (tt/yr>1.E+10): titel = 't-independent'
  plt.plot(zz,xx,lw=3,label=titel,c=colo[i])
  #--- analytic solution ---
  if (init==3):
    x0 = np.exp(-w*tt) * np.sin(k*zz)
    plt.plot(zz,x0,lw=1,c='black')
  if (init==4):  
    if (i==0): t0=tt
    ww = 2*np.sqrt(D*tt)
    AA = np.sqrt(t0/tt)
    x0 = AA*np.exp(-((zz-0.5)/ww)**2)
    plt.plot(zz,x0,lw=1,c='black')

xmin = xmax*1.E-10    
plt.ylim(xmin,xmax*2)
if (zunit=='cm'): plt.xlabel(r'$z\ \mathrm{[cm]}$')
if (zunit=='km'): plt.xlabel(r'$z\ \mathrm{[km]}$')
plt.yscale('log')
plt.ylabel(r'$\mathrm{concentration}$')
plt.legend(loc='best',fontsize=9,fancybox=True)
plt.tight_layout()
plt.savefig(pp,format='pdf')
data.close()

pp.close()
print '... written output to out.pdf.'


