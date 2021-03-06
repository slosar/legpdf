#!/usr/bin/env python

from optparse import OptionParser
from numpy.polynomial.legendre import Legendre as legser
import LegendrePDF
import random
import numpy as np
import matplotlib.pyplot as plt
import math
from statsmodels.compat.python import range
import numpy as np
from scipy.misc import comb
from scipy.interpolate import interp1d
from scipy.optimize import minimize
from scipy.special import erf

def main():
    o=getOptions()
    if o.pdffile:
        global pdfi, pdft
        pdfx,pdfy=np.loadtxt(o.pdffile,unpack=True)
        cpdfy=np.cumsum(pdfy)/len(pdfy)
        cpdfy=[0]+list(cpdfy)+[1.0]        
        pdfy=[0]+list(pdfy)+[0.0]
        pdfx=[0]+list(pdfx)+[1.0]
        pdft=interp1d(pdfx,pdfy)
        pdfi=interp1d(cpdfy,pdfx)
    hist,resM,resC,resL,resD=getBootStraps(o)
    ana=analyzeBootStraps(o,[resM,resC,resL,resD])
    if o.show:
        plotResults(o,hist,ana)    

def gauss_noise_prefactor(n):
    if n==0:
        return 0;
    if (n%2==1):
        return 0
    return math.gamma((n+1)/2.)*2**(n/2.)/np.sqrt(np.pi)


def cum2mc(kappa):
    '''convert non-central moments to cumulants
    recursive formula produces as many cumulants as moments

    References
    ----------
    Kenneth Lange: Numerical Analysis for Statisticians, page 40
    (http://books.google.ca/books?id=gm7kwttyRT0C&pg=PA40&lpg=PA40&dq=convert+cumulants+to+moments&source=web&ots=qyIaY6oaWH&sig=cShTDWl-YrWAzV7NlcMTRQV6y0A&hl=en&sa=X&oi=book_result&resnum=1&ct=result)


    '''
    mc = [1,0.0] #_kappa[0]]  #insert 0-moment and mean
    kappa0 = kappa[0]
    kappa = [1] + list(kappa)
    for nn,m in enumerate(kappa[2:]):
        n = nn+2
        mc.append(0)
        for k in range(n-1):
            mc[n] += comb(n-1,k,exact=1) * kappa[n-k]*mc[k]

    mc[1] = kappa0 # insert mean as first moments by convention
    return np.array(mc[1:])

def cum2mnc(kappa):
    mnc=[1]
    for k in range(len(kappa)):
        a1=comb(k,np.arange(k+1))
        a2=np.array(mnc)
        a3=kappa[k::-1]
        mnc.append( (a1*a2*a3).sum())
    return np.array(mnc[1:])


def mnc2cum(mnc):
    '''convert non-central moments to cumulants
    recursive formula produces as many cumulants as moments

    http://en.wikipedia.org/wiki/Cumulant#Cumulants_and_moments
    '''
    mnc = [1] + list(mnc)
    kappa = [1]
    for nn,m in enumerate(mnc[1:]):
        n = nn+1
        kappa.append(m)
        for k in range(1,n):
            kappa[n] -= comb(n-1,k-1,exact=1) * kappa[k]*mnc[n-k]

    return np.array(kappa[1:])





def getOptions():
    parser = OptionParser()
    parser.add_option("-m", dest="mean", default=0.5, type="float",
                      help="mean")
    parser.add_option("-s", dest="sigma", default=-1, type="float",
                      help="sigma, if negative use uniform distribution") 
    parser.add_option("--pdf", dest="pdffile", default="", type="string",
                      help="PDF filename")
    parser.add_option("-n", dest="n", default=20, type="int",
                      help="# of leg polys")
    parser.add_option("--lplot", dest="lplot", default=9, type="int",
                      help="# of leg polys to plot")
    parser.add_option("-B", dest="B", default=50, type="int",
                      help="# bootstraps")
    parser.add_option("-N", dest="N", default=10000, type="int",
                      help="# samples per boostrap")
    parser.add_option("-H", dest="h", default=100, type="int",
                      help="# histogram points")
    parser.add_option("-E", dest="he", default=1, type="int",
                      help="# histogram expansion (1=-1..2)")

    parser.add_option("-x", dest="noise", default=0.0, type="float",
                      help="add this noise to points")
    parser.add_option("--show", dest="show", default=False,
                      action="store_true", help="plot shit")


    (o, args) = parser.parse_args()

    return o

def generateMeasurements(o):
    if (o.pdffile):
        tmp=np.random.uniform(0.,1.,o.N)
        Fs=pdfi(tmp)
    elif o.sigma<0:
        Fs=np.random.uniform(0.,1.,o.N)
    else:
        l=[]
        while len(l)<o.N:
            g=random.gauss(o.mean, o.sigma)
            if (g>=0) and (g<=1):
                l.append(g)
        Fs=np.array(l)


    if (o.noise>0):
        Fs+=np.random.normal(0.,o.noise,o.N)
    return Fs

def getMeanF(o):
    if (o.pdffile):
        print "Here?"
        stop()
    if (o.sigma<0):
        return 0.5
    else:
        return math.erf(0.1e1 / o.sigma * o.mean * math.sqrt(0.2e1) / 0.2e1) / 0.2e1 - math.erf(math.sqrt(0.2e1) * (-0.1e1 + o.mean) / o.sigma / 0.2e1) / 0.2e1

def getTrueF(o, nbins):
    xx=np.linspace(0.,1.0,nbins)
    print xx,'AA'
    if (o.pdffile):
        yy=pdft(xx)
    elif (o.sigma<0):
        yy=np.ones(len(xx))
    else:
        norm=math.erf(0.1e1 / o.sigma * o.mean * math.sqrt(0.2e1) / 0.2e1) / 0.2e1 - math.erf(math.sqrt(0.2e1) * (-0.1e1 + o.mean) / o.sigma / 0.2e1) / 0.2e1
        yy=1/np.sqrt(2*np.pi*o.sigma**2)*np.exp(-(xx-o.mean)**2/(2*o.sigma**2))/norm
    return xx,yy

def histX(o):
    return np.linspace(-o.he,1+o.he,o.h*(1+2*o.he))

def fixpdf(p):
    p=abs(p)
    return p/p.sum()

def deconvolveLogLike(bins,binr,thist,N):
    prop=np.zeros(len(thist))
    bins=fixpdf(bins)
    for b,r in zip(bins,binr):
        prop+=b*r
    #last bin
    #prop+=(1.0-bins.sum())*binr[-1]
    prop*=N
    err2=prop
    err2[err2<0]=1e-10
    err2[err2==0]=1e-3 ## to avoid division by zero
    
    #plt.plot(thist)
    #plt.plot(prop)
    #plt.show()
    chi2=((thist-prop)**2/err2).sum()
    print chi2

    return chi2


def deconvolvePDF(o,thist, work):
    if work is None:
        xx=histX(o)
        binr=[]
        for i in range(o.n):
           ## response of n-th bin
            lo,up=i*1.0/o.n,1.0*(i+1)/o.n
            print lo,up
            if (o.noise)>0:
                cbin=(-erf(math.sqrt(0.2e1) * (lo - xx) / o.noise / 0.2e1) / 0.2e1 +
                      erf(math.sqrt(0.2e1) * (up - xx) / o.noise / 0.2e1) / 0.2e1)
            else:
                cbin=np.ones(len(xx))
                cbin[xx<lo]=0.0
                cbin[xx>=up]=0.0
            cbin/=cbin.sum()
            binr.append(cbin)
        N=thist.sum()
        work=(binr,N)        
    else:
        binr,N=work
        
    ax,iguess=getTrueF(o,2*o.n+1)
    iguess=iguess[1::2] ## these are now centered!!
    iguess*=(1+np.random.normal(0,0.1,o.n))

    #iguess=np.random.uniform(0,1,o.n)
    #iguess/=iguess.sum()

    toret=minimize(deconvolveLogLike,iguess,args=(binr,thist,N),method='powell')
    toret=fixpdf(toret.x)
    #toret=np.array(list(toret)+[1.0-toret.sum()])
    return toret,work
    

def getBootStraps(o):
    p=LegendrePDF.Mom2Leg(o.n)
    resM=np.zeros((o.n, o.B))
    resC=np.zeros((o.n, o.B))
    resL=np.zeros((o.n, o.B))
    resD=np.zeros((o.n, o.B))
    hN=o.h*(1+o.he*2)
    hist=np.zeros(hN,int)
    dwork=None
    ## bootstrap
    for b in xrange(o.B):
        nums=generateMeasurements(o)
        chist=np.bincount((nums*o.h+o.he*o.h).astype(int),minlength=hN)
        hist+=chist
        moms=[1.0]
        for m in range(1,o.n):
            moms.append((nums**(m)).mean())
        moms=np.array(moms)
        kappa=np.ones(o.n)
        kappa[1:]=mnc2cum(moms[1:])
        if (o.noise>0):
            kappa[2]-=o.noise**2
            moms[1:]=cum2mnc(kappa[1:])
            
        resM[:,b]=moms
        resC[:,b]=kappa
        resL[:,b]=p.mom2leg(resM[:,b])
        dres,dwork=deconvolvePDF(o,chist,dwork)
        resD[:,b]=dres
        
    hist=hist/(1.*o.B*o.N/o.h)

    return hist, resM, resC, resL, resD


def analyzeBootStraps(o, li):
    oli=[]
    for res in li:
        mm=[]
        for i in range(0,o.n):
            mm_=res[i,:].mean()
            mm.append(mm_)

        mm=np.array(mm)
        cm=np.zeros((o.n,o.n))
        for i in range(o.B):
            vec=res[:,i]-mm
            cm+=np.outer(vec,vec)

        cm/=o.B**2 ## a factor of two to get to per sample
        oli.append((mm,cm))
    
    return oli

def cov2cor(m):
    d=np.sqrt(m.diagonal())
    return m/np.outer(d,d)
    
    cm/=np.outer(cmd,cmd)
    cld=np.sqrt(cl.diagonal())
    cl/=np.outer(cld,cld)


def plotResults(o,hist,oli):
    plt.figure(figsize=(16,16))
    plt.subplot(5,1,1)
    xv=histX(o)
    plt.plot (xv,hist,'r-', label='histo')
    ml=oli[2][0]
    ls=legser(ml[:o.lplot]) 
    xvp=np.linspace(0,1,o.h)
    plt.plot (xvp,ls(2*xvp-1),'b--',lw=2, label='leg')
    md=oli[3][0]
    xvp=np.linspace(0,1,o.n)
    plt.plot (xvp,md*o.n,'g--',lw=2, label='decon')
    xx,yy=getTrueF(o,100)
    plt.plot(xx,yy,'k-',lw=1,label=True)
    plt.xlabel('val')
    plt.ylabel('p(val)')
    plt.legend()

    s=0
    names=['moment','cumulant','leg coeff','decon']
    for me,cov in oli:
        plt.subplot(5,2,3+2*s)
        plt.plot ([0,o.n],[0,0],'k--')
        plt.errorbar(range(1,o.n),me[1:],np.sqrt(cov.diagonal()[1:]))
        if (o.sigma<0) and (s==0):
            plt.plot(range(1,o.n), 1./(1+np.arange(1,o.n)),'bo')
        plt.ylabel(names[s])

        plt.subplot(5,2,4+2*s)
        plt.imshow(cov2cor(cov), interpolation='nearest')
        plt.xlabel("cov matrix")
        plt.colorbar()
        s+=1
    plt.savefig("fig.pdf")
    plt.show()
        

if __name__=="__main__":
    main()
