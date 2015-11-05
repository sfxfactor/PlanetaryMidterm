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

plt.hist(acc)
plt.savefig('acc.pdf')
plt.clf()

plt.plot(chi)
plt.savefig('chiAll.pdf')
plt.clf()

plt.plot(chi[20000:])
plt.savefig('chiMinusBurn.pdf')
plt.clf()

uSamples=samples
uSamples[:,0]=samples[:,0]-2451725.
uSamples[:,1]=samples[:,1]-3.5246
uSamples[:,2]=samples[:,2]/mjup
fig = corner.corner(uSamples[20000:],labels=[r"$\tau_0~[\mathrm{HJD}]$",r"$P~[\mathrm{days}]$",r"$M\sin i~[M_\mathrm{Jup}]$"],title_fmt=".6f",quantiles=[0.16,0.5,0.84],show_titles=1,verbose=1,use_math_text=1,levels=[0.68,0.95,0.997])
fig.savefig('corner.pdf')
plt.clf()

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
plt.errorbar(phdate,RV,yerr=sig,fmt='o')
plt.plot(t,m)
plt.savefig('phrf.pdf')
plt.clf()
