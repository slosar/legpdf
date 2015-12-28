#
# A quick class to covnert moments to legendre polynomials
# 
#

from math import *
import numpy as np
import numpy.linalg as la
import math
class Mom2Leg(object):
    def __init__ (self, N):
        self.N=N
        self.mat=self.leg2mom_mat(N)
        self.imat=la.inv(self.mat)

    def mom2leg_mat(self,N):
        """ Returns matrix that converts momets into leg polynomials """
        mat=self.leg2mon(N)
        return la.inv(mat)
    def leg2mom_mat(self,N):
        """ Returns matrix that converts legendrey poly into moments"""
        mat=np.zeros((N,N))
        for n in range(N):
            for l in range(n+1):
                ## note that our indices are one off
                mat[n,l]=float(math.factorial(n))**2/math.factorial(n-l)/math.factorial(n+l+1)
        return mat
    
    def mom2leg(self, vec):
        return np.dot(self.imat,vec)

