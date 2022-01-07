def GaussQuad(n,s0,s1):
    nargin=len(locals())
    d=[]; w=[]
    if nargin==1:
        s0=-1; s1=1
    if n==1:
        d=(s0+s1)/2; w=s1-s0  
        return d,w 
    j=numpy.arange(1,n)
    bta=j/numpy.sqrt(4*(j**2)-1)   
    Q=numpy.diag(bta, k=-1)+numpy.diag(bta, k=1)
    L,U = numpy.linalg.eig(Q)
    k=numpy.argsort(L)
    d=L[k]; w=2*(U[0][k]**2)
    d=s0+0.5*(s1-s0)*(d+1)
    w=0.5*(s1-s0)*w
    return d,w