import corner
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import astropy.io.ascii as ascii

G = 6.67259e-8 #cm^3 g^-1 s^-2
m1 = 1.148*1.989e33
mjup=1.898e30

acc = fits.getdata('accept.fits')
samples = fits.getdata('samples.fits')
chi = fits.getdata('Chi2vals.fits')

print plt.hist(acc,bins=2,normed=1)[0]/2.
plt.savefig('acc.pdf')
plt.clf()

plt.plot(chi)
plt.xlabel("Step number")
plt.ylabel(r"$\chi^2$")
plt.savefig('chiAll.pdf')
plt.clf()

plt.plot(chi[20000:])
plt.xlabel("Step number after 20,000")
plt.ylabel(r"$\chi^2$")
plt.savefig('chiMinusBurn.pdf')
plt.clf()

samples[:,0]=samples[:,0]-2451725.
samples[:,1]=samples[:,1]-3.5
samples[:,2]=samples[:,2]/mjup
fig,axarr = plt.subplots(3,3)
fig = corner.corner(samples[20000:],fig=fig,truths=[None,0.02472,0.69],labels=[r"$\tau_0 [\mathrm{HJD}-2451725]$",r"$P [\mathrm{days}-3.5]$",r"$M\sin i [M_\mathrm{J}]$"],quantiles=[0.16,0.5,0.84],verbose=1,use_math_text=1,levels=[0.68,0.95,0.997],tick_labelsize=8)
axarr[1,0].tick_params(axis='both',labelsize=8)
axarr[2,1].tick_params(axis='both',labelsize=8)
axarr[2,0].tick_params(axis='both',labelsize=10)
fig.savefig('corner.pdf')
plt.clf()

samples = fits.getdata('samples.fits')
medfit = np.median(samples[20000:],axis=0)
print medfit
T0,P,Msini=medfit

def bfrv(dates):
    Ps = P*86400
    n = (2.*np.pi)/Ps
    a = ((G*m1*Ps**2)/(4.*np.pi**2))**(1./3.)
    vsini = (Msini/m1)*n*a
    model = (vsini * np.sin(2.*np.pi*(dates)/P))
    return model

data = ascii.read('RV.dat')
date = np.array(data['date'])
RV = np.array(data['RV'])*1e2 #m->cm
sig = np.array(data['sig'])*1e2 #m->cm

t=np.linspace(0,P,1000)
m = bfrv(t)

phdate = (date - T0)%P
plt.errorbar(phdate/P,RV*1e-2,yerr=sig*1e-2,fmt='o')
plt.plot(t/P,m*1e-2)
plt.xlabel('Phase')
plt.ylabel(r"$RV [\mathrm{m/s}]$")
plt.savefig('phrv.pdf')
plt.clf()

def npbfrv(dates):
    Ps = P*86400
    n = (2.*np.pi)/Ps
    a = ((G*m1*Ps**2)/(4.*np.pi**2))**(1./3.)
    vsini = (Msini/m1)*n*a
    model = (vsini * np.sin(2.*np.pi*(dates-T0)/P))
    return model

t=np.linspace(np.min(date),np.max(date),500000)

plt.errorbar(date,RV*1e-2,yerr=sig*1e-2,fmt='o')
plt.plot(t,npbfrv(t)*1e-2)
plt.xlim((2.451754e6,P+2.451754e6))
plt.xlabel(r"HJD")
plt.ylabel(r"$RV [\mathrm{m/s}]$")
plt.savefig('rverror.pdf')
plt.clf()
