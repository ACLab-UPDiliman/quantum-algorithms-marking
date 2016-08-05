# This is a simple simulation of a quantum algorithm for approximate string matching.
# Please refer to conference paper 27 of Workshop on Computation: Theory and Practice, 2015.

# Outline of filtering subroutine
# Let Loc be a set of candidate solution indices in T. Let Loc={}
# For 1 to N/(q-d) do
# 	Prepare a superposition of indices i in T. //|i>|0> 
#	Given each symbol T[i] determine the index j_{T[i]} of first occurrence of each symbol T[i] in P. //|i>|j_{T[i]}> 
#	Subtract each j_{T[i]} from their corresponding index i. //|i>|i-j_{T[i]}>
#	Measure state of second register. //k = i-j_{T[i]}
#  	Add k to Loc.
# Return Loc.
from numpy import zeros, ones, identity, empty, array, dot, kron, tensordot, sort, matrix
from math import sqrt, ceil, log
from QOperators import H2x2, ket0, ket1, CNOT

def constructP_Sym(pattern):
# Given P, setify P to get all distinct symbols in P convert the set into a list by
# listifying it. Sortification is not necessary actually.
	return list(set(pattern))

def constructP_Loc(sym_set):
# there are many ways to construct P_Loc using Python's list comprehension feature but
# using map we explicitly say that there is a one-to-one and onto mapping between 
# P_Sym and P_Loc

	return map(getFirstLoc,sym_set)

def getFirstLoc(x):
	# print 'P.index('+x+'): ', P.index(x)
	return P.index(x)

# construct operators
def constructU_Loc(text,pattern,sym_set,loc_set):
	T_len = len(text)
	P_len = len(pattern)
	matrix_dim = T_len*P_len
	matrix_operator = zeros((matrix_dim,matrix_dim),int)
	# print "sym_set: ", sym_set
	for x in xrange(0,T_len):
		if text[x] in sym_set:
			# print "text["+str(x)+"]: ", text[x], " available in P_sym and is at index " + str(sym_set.index(text[x])) + "."
			index = sym_set.index(text[x])
			loc = loc_set[index]
			matrix_operator[x*P_len][(x*P_len)+loc] = 1
			matrix_operator[(x*P_len)+loc][x*P_len] = 1
			for y in xrange(0,P_len):
				if y != 0 and y != loc:
					matrix_operator[(x*P_len)+y][(x*P_len)+y] = 1
		else:
			# print "text["+str(x)+"]: ", text[x], " not available in P_sym"
			for y in xrange(0,P_len):
				matrix_operator[(x*P_len)+y][(x*P_len)+y] = 1

	return matrix_operator

def constructU_Sub():
	return

def constructU_Mis():
	return

def constructU_Ham():
	return

def constructU_Mark():
	return

def constructU_Amp():
	return

# initialize states
def constructPsi_init(dimension):
	num_qubits = int(ceil(log(dimension,2)))
	register1 = ket0()
	H = H2x2()
	for x in xrange(1,num_qubits):
		register1 = kron(register1,ket0())
		H = kron(H,H2x2())
	register1 = dot(H,register1)
	return kron(register1,ket0())

################ Filtering Phase ################
# Define T
T = "aabc"
# Define P
P = "ab"
print "T=",T, ", N=", str(len(T))
print "P=",P, ", M=", str(len(P))

# Construct P_Sym and P_Loc from P
P_Sym = constructP_Sym(P)
print 'P_Sym: ', P_Sym, ", q=", str(len(P_Sym))
P_Loc = constructP_Loc(P_Sym)
print 'P_Loc: ', P_Loc

######### Define unitary operator operators #########
# U_Loc
U_Loc = constructU_Loc(T,P,P_Sym,P_Loc)
print "U_Loc: "
print U_Loc
# U_Sub

####### Execute algorithm steps #######
# Prepare initial superposition state
psi_init = constructPsi_init(len(T))
print "|psi_init>: ", psi_init

# Apply U_Loc to |psi_init> to identify location of first occurence of each symbol T[i]
psi_loc = dot(U_Loc,psi_init)
print "|psi_loc>: ", psi_loc

# CNOT
control = 0
target = 0
result = CNOT(control,target)
print 'CNOT(',control,',',target,') = (',result[0],',',result[1],')'

# x_i = 1
# y_i = 1
# previous_carry = 1
# z_i = BitToBitFullAdder(x_i,y_i,previous_carry)
# print "x_i: ", x_i, "y_i: ", y_i, "previous carry: ", previous_carry
# print "Sum: ", z_i.sum, ", Carry: ", z_i.carry
################ Verification Phase ################








