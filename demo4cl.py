import numba

@numba.jit()
def _demo():
	x=0
	for i in xrange(100000):
		for j in xrange(100000):
			x+=1
	return x

#print _demo()
