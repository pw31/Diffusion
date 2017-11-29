import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
plt.rcParams['axes.linewidth'] = 1.5
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
  if (tunit=='s'): titel = 't = %5.3f s' % (tt)
  if (tunit=='yr'): titel = 't = %7.4f yr' % (tt/yr)
  if (tt/yr>1.E+10): titel = 't-independent'
  plt.plot(zz,xx,lw=4,label=titel,c=colo[i])
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

plt.ylim(0.0,xmax+0.05)
if (zunit=='cm'): plt.xlabel(r'$z\ \mathrm{[cm]}$',fontsize=22)
if (zunit=='km'): plt.xlabel(r'$z\ \mathrm{[km]}$',fontsize=22)
plt.ylabel(r'$\mathrm{concentration}$',fontsize=22)
plt.legend(loc='lower left',fontsize=10,fancybox=True)
plt.tight_layout()
plt.savefig(pp,format='pdf')
data.close()

pp.close()
print '... written output to out.pdf.'


