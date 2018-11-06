import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
ppp = PdfPages('Tstruc.pdf')
#from matplotlib import rc
#rc('text', usetex=True)   # makes it nice but VERY slow

bar = 1.E+6         
data = np.loadtxt('./Structures/HotRockyPlanets/55cnc-e.txt')
ppdat = (10**data[:,0])/bar
TTdat = data[:,1]
print ppdat
print TTdat

#--------------------------------
p3 = 10.0         # surface
p2 = 0.0          # tropopause (will be derived below)
p1 = 5.E-5        # stratopause
p0 = 1.E-9        # exobase
p4 = 10.0*p3
#--------------------------------
T3 = 2990.        # surface
T2 = 1210.        # tropopause
T1 = 1280.        # stratopause
T0 = 1050.
#--------------------------------
a1 = np.log(p1/p0)/np.sqrt(T1-T0)
a2 = np.log(p3/p1)/(np.sqrt(T1-T2)+np.sqrt(T3-T2))
p2 = p1*np.exp(a2*np.sqrt(T1-T2))
print "a1,a2=",a1,a2


pp = np.zeros(1001, dtype=np.float)
TT = np.zeros(1001, dtype=np.float)
dp = np.log(p4/p0)/1000
nsmooth = 50
for i in range(0,1001):
  lnp = np.log(p0) + dp*i
  pp[i] = np.exp(lnp)
  #--- boxcar filter ---
  Tsum = 0.0
  for j in range(-nsmooth,+nsmooth+1):
    lnp = np.log(p0) + dp*(i+j)
    p = np.exp(lnp)
    if (p<p1):
      Tsum += T0 + (np.log(p/p0)/a1)**2
    elif (p<p3):
      Tsum += T2 + (np.log(p/p2)/a2)**2
    else:  
      Tsum += T3
  TT[i] = Tsum/(2*nsmooth)
  
fig,ax = plt.subplots(figsize=(6,7))
plt.plot(TT,np.log10(pp),color='blue',lw=3,ls='solid',label='model')
plt.scatter(TTdat,np.log10(ppdat),color='black',label='55 Cancri e')
plt.ylim(np.log10(p4),np.log10(p0))
plt.xlabel(r'$T\ \mathrm{[K]}$',fontsize=20)
plt.ylabel(r'$\log_{10}\ p\ \mathrm{[bar]}$',fontsize=20)
plt.tick_params(axis='both', labelsize=15)
plt.tick_params('both', length=6, width=1.5, which='major')
plt.tick_params('both', length=3, width=1, which='minor')
plt.legend(loc='upper right',fontsize=11,fancybox=True)
plt.tight_layout()
#plt.show()
plt.savefig(ppp,format='pdf')

ppp.close()
print '... written output to Tstruc.pdf.'


