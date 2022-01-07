def getRefinedValues(u,xl,d,k):
    d_ref=numpy.linspace(0,1,k+1) # refined basis element
    E,D=interpolatingMat(d,d_ref)
    n=numpy.size(xl); m=numpy.size(d)-1
    x_ref=(n-1)*k+1; # number of all refined mesh points
    x=numpy.zeros((x_ref,))
    y=numpy.zeros(x.shape)
    dy=numpy.zeros(x.shape)
    for l in range(n-1):
        h=xl[l+1]-xl[l]  #size of element l
        i_c=l*m+numpy.arange(m+1) # coarse point indices
        i_ref=l*k+numpy.arange(k+1)  # refined points indices
        x[i_ref]=numpy.linspace(xl[l],xl[l+1],k+1)
        y[i_ref]=E @ u[i_c] 
        dy[i_ref]=D @ u[i_c]/h 
    return x,y,dy
