from numpy import array, dot, kron, zeros, ones, identity, outer, array
from scipy.fftpack import fft, ifft
from math import sqrt, log, ceil

def ket0():
	return array([1,0])

def ket1():
	return array([0,1])

def H2x2():
	a = array([[1,1],[1,-1]])
	h = (1/sqrt(2)) * a
	return h

def H2x2_kron_log_dim(dim):
	log_dim = int(ceil(log(dim,2)))
	# print 'log(' + str(dim) + ') = ' + str(log_dim)
	h = H2x2()
	for i in xrange(0,log_dim-1):
		h = kron(h,H2x2())
	return h

def QFT(a):
	dim = a.shape[0]
	factor = sqrt(1.0/dim)
	m = factor * fft(a)
	return m

def IQFT(a):
	dim = a.shape[0]
	factor = sqrt(1.0/dim)
	# print 'IQFT mult factor = ' + str(factor)
	# print 'a = ' + str(a)
	#NOTE: In IFFT, the elements of the matrix are multiplied by factor 1/dim
	#		so that DFT*IFFT = I. So, in the quantum setting, we need to compensate
	#		for the inherently included multiplicative factor 1/dim in each element of
	#		IFFT by factoring in dim. This will make our quantum version of IFFT,
	#		IQFT, be equal to (1/dim)*(IFFT) without the previous factor dim inside
	#		IFFT.
	m = factor * dim * ifft(a)
	return m

def Euclid_Norm(cplex):
	return sqrt((cplex.conjugate()*cplex).real)

def Euclid_Norm_Squared(cplex):
	return (cplex.conjugate()*cplex).real

