#!/usr/bin/env python

from optparse import OptionParser
from numpy.polynomial.legendre import Legendre as legser
import LegendrePDF
import random
import numpy as np
import pylab

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
parser.add_option("--show", dest="show", default=False,
                  action="store_true", help="plot shit")


(o, args) = parser.parse_args()

if o.sigma<0:
    ra=lambda : random.uniform(0,1)
else:
    def ra():
        while True:
            g=random.gauss(o.mean, o.sigma)
            if (g>=0) and (g<=1):
                return g

p=LegendrePDF.Mom2Leg(o.n)

resM=np.zeros((o.n, o.B))
resL=np.zeros((o.n, o.B))
hist=np.zeros(o.h)
for b in xrange(o.B):
    nums=np.array([ra() for i in xrange(o.N)])
    for v in nums:
        hist[int(v*o.h)]+=1
    # get moments
    moms=[]
    ## this is an idiotally inefficient way of producing one
    ## for m=0, but who cares
    for m in range(o.n):
        moms.append((nums**(m)).mean())
    resM[:,b]=np.array(moms)
    resL[:,b]=p.mom2leg(resM[:,b])
hist/=(1.*o.B*o.N/o.h)
mmv=[1]
lmv=[1]
mev=[0]
lev=[0]

for i in range(1,o.n):
    mm,ms=resM[i,:].mean(),np.sqrt(resM[i,:].var()/o.B)
    lm,ls=resL[i,:].mean(),np.sqrt(resL[i,:].var()/o.B)
    print i,"m= {0} +/- {1}".format(mm,ms)
    print "l= {0} +/- {1}".format(lm,ls)
    mmv.append(mm)
    mev.append(ms)
    lmv.append(lm)
    lev.append(ls)

mmv=np.array(mmv)
lmv=np.array(lmv)


if o.show:
    pylab.figure(figsize=(16,16))
    pylab.subplot(3,1,1)
    xv=np.linspace(0,1,o.h)
    pylab.plot (xv,hist,'r-', label='histo')
    ls=legser(lmv)
    pylab.plot (xv,ls(-1+2*xv),'b-', label='leg')
    pylab.xlabel('val')
    pylab.ylabel('p(val)')
    pylab.legend()

    for s in range(2):
        pylab.subplot(3,2,3+2*s)
        pylab.plot ([0,o.n],[0,0],'k--')
        if (s==0):
            pylab.errorbar(range(1,o.n),mmv[1:],mev[1:])
            pylab.ylabel('moment')
        else:
            pylab.errorbar(range(1,o.n),lmv[1:],lev[1:])
            pylab.ylabel('leg. coeff.')
        pylab.xlabel('n')

        pylab.subplot(3,2,4+2*s)
        c=np.zeros((o.n,o.n))
        for i in range(o.B):
            if (s==0):
                vec=resM[:,i]-mmv
            else:
                vec=resL[:,i]-lmv
            c+=np.outer(vec,vec)
        # now we have the full cov matrix
        cd=np.sqrt(c.diagonal())
        c/=np.outer(cd,cd)
        pylab.imshow(c, interpolation='nearest')
        pylab.xlabel("cov matrix")
        pylab.colorbar()
    pylab.savefig("fig.pdf")
    pylab.show()
        
