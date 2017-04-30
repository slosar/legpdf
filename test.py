#!/usr/bin/env python

from optparse import OptionParser
from numpy.polynomial.legendre import Legendre as legser
import LegendrePDF
import random
import numpy as np
import pylab



def main():
    o=getOptions()
    hist,resM,resL=getBootStraps(o)
    mm,cm,ml,cl=analyzeBootStraps(o,resM,resL)
    if o.show:
        plotResults(o,hist,mm,cm,ml,cl)    


def getOptions():
    parser = OptionParser()
    parser.add_option("-m", dest="mean", default=0.5, type="float",
                      help="mean")
    parser.add_option("-s", dest="sigma", default=-1, type="float",
                      help="sigma, if negative use uniform distribution")
    parser.add_option("-n", dest="n", default=20, type="int",
                      help="# of leg polys")
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
    if o.sigma<0:
        Fs=np.random.uniform(0.,1.,o.N)
    else:
        l=[]
        while len(l)<o.N:
            g=random.gauss(o.mean, o.sigma)
            if (g>=0) and (g<=1):
                l.append(g)
        Fs=np.array(l)
    if (o.noise>0):
        Fs+=np.random.normal(0,o.noise,o.N)
    return Fs

def getBootStraps(o):
    p=LegendrePDF.Mom2Leg(o.n)
    resM=np.zeros((o.n, o.B))
    resL=np.zeros((o.n, o.B))
    hN=o.h*(1+o.he*2)
    hist=np.zeros(hN,int)
    ## bootstrap
    for b in xrange(o.B):
        nums=generateMeasurements(o)
        hist+=np.bincount((nums*o.h+o.he*o.h).astype(int),minlength=hN)
        moms=[]
        ## this is an idiotally inefficient way of producing one
        ## for m=0, but who cares
        for m in range(o.n):
            moms.append((nums**(m)).mean())
        resM[:,b]=np.array(moms)
        resL[:,b]=p.mom2leg(resM[:,b])
    hist=hist/(1.*o.B*o.N/o.h)
    print hist
    return hist, resM, resL


def analyzeBootStraps(o, resM,resL):

    mm=[1]
    ml=[1]
    ##means and errors
    for i in range(1,o.n):
        mm_=resM[i,:].mean()
        ml_=resL[i,:].mean()
        mm.append(mm_)
        ml.append(ml_)

    mm=np.array(mm)
    ml=np.array(ml)

    cm=np.zeros((o.n,o.n))
    cl=np.zeros((o.n,o.n))
    for i in range(o.B):
        vec=resM[:,i]-mm
        cm+=np.outer(vec,vec)
        vec=resL[:,i]-ml
        cl+=np.outer(vec,vec)
        
    cmd=np.sqrt(cm.diagonal())
    cm/=np.outer(cmd,cmd)
    cld=np.sqrt(cl.diagonal())
    cl/=np.outer(cld,cld)
    cm/=o.B
    cl/=o.B
    return mm,cm,ml,cl


def plotResults(o,hist,mm,cm,ml,cl):


    pylab.figure(figsize=(16,16))
    pylab.subplot(3,1,1)
    xv=np.linspace(-o.he,1+o.he,o.h*(1+2*o.he))
    pylab.plot (xv,hist,'r-', label='histo')
    ls=legser(ml) 
    xvp=np.linspace(0,1,o.h)
    pylab.plot (xvp,ls(2*xvp-1),'b-', label='leg')
    pylab.xlabel('val')
    pylab.ylabel('p(val)')
    pylab.legend()

    for s in range(2):
        pylab.subplot(3,2,3+2*s)
        pylab.plot ([0,o.n],[0,0],'k--')
        if (s==0):
            pylab.errorbar(range(1,o.n),mm[1:],np.sqrt(cm.diagonal()[1:]))
            pylab.ylabel('moment')
        else:
            pylab.errorbar(range(1,o.n),ml[1:],np.sqrt(cl.diagonal()[1:]))
            pylab.ylabel('leg. coeff.')
        pylab.xlabel('n')

        pylab.subplot(3,2,4+2*s)
        if (s==0):
            pylab.imshow(cm, interpolation='nearest')
        else:
            pylab.imshow(cl, interpolation='nearest')
        pylab.xlabel("cov matrix")
        pylab.colorbar()
    pylab.savefig("fig.pdf")
    pylab.show()
        

if __name__=="__main__":
    main()
