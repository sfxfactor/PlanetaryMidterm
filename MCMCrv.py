import numpy as np
import astropy.io.ascii as ascii
from astropy.io import fits
G = 6.67259e-8 #cm^3 g^-1 s^-2

#get data
dataTitle = 'RV'
data = ascii.read(dataTitle+'.dat')
date = np.array(data['date'])
RV = np.array(data['RV'])*1e2 #m->cm
sig = np.array(data['sig'])*1e2 #m->cm
m1 = 1.148*1.989e33
niter=100001

def calcX2(T0,P,Msini):#,v0):
    Ps=P*86400.
    n = (2.*np.pi)/Ps
    a = ((G*m1*Ps**2)/(4.*np.pi**2))**(1./3.)
    vsini = (Msini/m1)*n*a
    model = vsini * np.sin(2.*np.pi*(date-T0)/P) #+v0
    X2 = sum(((RV-model)/sig)**2)
    return X2


#initial model
T0 = 2451725.4
P = 3.524
Msini =1.*1.898e30
#v0 = 0.0
sT = 0.01
sP = 0.0006
sM = 4e28
#sv = 1.

Chi2 = calcX2(T0,P,Msini) #,v0)

accept = np.array([])
samples = np.zeros((niter,3))
Chi2vals = Chi2

T0n = T0
Pn = P
Msinin=Msini
change=[]

for i in range(0,niter):
    if (i%10000 == 0):
        print i/1000,'%'
        #generate new model
    n=np.random.random_integers(1,high=3)
    #n=3 #used for fine tuning widths to get acceptance fraction 20-40%
    change=np.append(change,n)
    if n==1:
        T0n=np.random.normal(T0,sT)
        Pn=P
        Msinin=Msini
        #v0n = v0
    elif n==2:
        T0n=T0
        Pn=np.random.normal(P,sP)
        Msinin=Msini
        #v0n = v0
    elif n==3:
        T0n=T0
        Pn=P
        Msinin=np.random.normal(Msini,sM)
        #v0n = v0
    #else:
    #    T0n=T0
    #    Pn=P
    #    Msinin=Msini
    #    v0n = np.random.normal(v0,sv)

    Chi2n=calcX2(T0n,Pn,Msinin) #,v0n)
    
    #random number
    n=np.random.random()
    if n<np.exp(0.5*(Chi2-Chi2n)):
        #accept
        accept=np.append(accept,1)
        T0=T0n
        P=Pn
        Msini=Msinin
        Chi2=Chi2n
        #v0=v0n
        samples[i]=np.array([T0,P,Msini]) #,v0])
        Chi2vals=np.append(Chi2vals,Chi2)
    else:    
        #reject
        accept=np.append(accept,0)
        samples[i]=np.array([T0,P,Msini]) #,v0])
        Chi2vals=np.append(Chi2vals,Chi2)

fits.writeto('samples.fits',samples)
fits.writeto('Chi2vals.fits',Chi2vals)
fits.writeto('accept.fits',accept)
print 'done!'

