def interpolatingMat(x,d):
    nx=numpy.size(x); nd=numpy.size(d)
    E=numpy.zeros((nd,nx)) 
    D=numpy.zeros((nd,nx))
    b=numpy.zeros((nx,))
    j_all=numpy.arange(nx)
    for i in range(nx):
        j=numpy.delete(j_all,i) 
        b[i]=1/numpy.prod(x[i]-x[j]) 
    for i in range(nx):
        j=numpy.delete(j_all,i)
        for k in range(nd):
            E[k,i]=b[i]*numpy.prod(d[k]-x[j])
            ds=0
            for l in j:
                i1=min(i,l); i2=max(i,l) 
                jj=numpy.delete(j_all,[i1,i2]) 
                ds=ds+numpy.prod(d[k]-x[jj]) 
            D[k,i]=b[i]*ds
    return E,D
