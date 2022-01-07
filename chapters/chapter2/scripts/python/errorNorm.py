def errorNorm(u,du,y,xl,d):
    m=numpy.size(d)-1
    d_new,w=GaussQuad(m+2,0,1)
    E,D=interpolatingMat(d,d_new) 
    e0=0; e1=0
    for l in range(numpy.size(xl)-1):
        h=xl[l+1]-xl[l]
        dl=xl[l]+h*d_new
        ic=l*m+numpy.arange(m+1)
        uh=E @ y[ic]; duh=D @ y[ic]/h 
        e=u(dl)-uh; de=du(dl)-duh
        e0=e0+h/2*sum(w*(e**2))
        e1=e1+h/2*sum(w*(de**2))
    e0=e0**(1/2); e1=e1**(1/2)
    return e0,e1
