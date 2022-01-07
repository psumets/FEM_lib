def solveBVP(bc,c,b,a,f,xl,x,d,w):
    m=numpy.size(x)-1
    n=numpy.size(xl) 
    N=(n-1)*m+1 # number of mesh points
    x_mesh=numpy.zeros((N,)) # all mesh points
    for l in range(n-1):
        ne=l*m+numpy.arange(m+1) # mesh point indices
        x_mesh[ne]=xl[l]+(xl[l+1]-xl[l])*x 
    A,F=assemblingAF(c,b,a,f,xl,x,d,w)
    A0,F0=assemblingBC(bc,A,F)
    u=solveFEM(bc,A0,F0)
    return u,x_mesh