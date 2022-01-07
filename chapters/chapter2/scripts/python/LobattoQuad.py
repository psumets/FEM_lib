def LobattoQuad (n,s0,s1):
    nargin=len(locals())
    d=[]; w=[]
    if n<2:
        print('LobattoQuad: number of nodes is less then 2') 
        return d,w 
    if nargin==1:
        s0=-1; s1=1 
    if n==2:
        d=[s0,s1]; w=[(s1-s0)/2,(s1-s0)/2]  
        return d,w 
    k=numpy.arange(1,n-2)
    bta=numpy.sqrt(k*(k+2)/(2*k+1)/(2*k+3)) 
    Q=numpy.diag(bta,k=-1)+numpy.diag(bta,k=1) 
    L,U = numpy.linalg.eig(Q)
    k=numpy.argsort(L)
    d=L[k]; w=4/3*(U[0][k]**2); w=w/(1-d**2)
    w=numpy.append(w,2/(n**2-n))
    w=numpy.insert(w,0,2/(n**2-n))
    d=numpy.append(d,1)
    d=numpy.insert(d,0,-1)
    d=s0+0.5*(s1-s0)*(d+1);
    w=0.5*(s1-s0)*w;
    return d,w