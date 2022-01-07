def mainTest():
    # set BVP data:
    # integration domain
    s_0=-1; s_1=1 
    # exact solution
    def u(x): return numpy.sin(x) 
    def du(x): return numpy.cos(x) 
    # equation coefficients
    def c(x): return numpy.cos(x)
    def b(x): return -2*numpy.sin(x) 
    def a(x): return numpy.sin(x)
    def f(x): return numpy.sin(x)**2
    # boundary conditions
    sigma=0
    mu=c(s_1)*du(s_1)+sigma*u(s_1)
    bc=numpy.array([[1, 0, u(s_0)],
                     [3, sigma, mu] ])    
    # set FEM parameters
    m=3 # polynomial degree
    ne=6 # number of final elements
    M=m+1 # number of quadrature nodes
    xl=numpy.linspace(s_0,s_1,ne+1) # finite element nodes
    t,_ =LobattoQuad(m+1,0,1) # interpolation nodes
    d,w=LobattoQuad(M,0,1) # quadrature nodes and weighs
    # solve BVP
    start_time = timeit.default_timer()
    uh,x=solveBVP(bc,c,b,a,f,xl,t,d,w) 
    print("solving time =",timeit.default_timer() - start_time)
    # find error in L_2 and H^1 norms
    e0,e1=errorNorm(u,du,uh,xl,t)
    print("L2 norm error =",e0)
    print("H1 norm error =",e1)
    X,Uh,dUh=getRefinedValues(uh,xl,t,10*m) 
    # find error in various norms
    Z=u(X)-Uh; z=u(x)-uh; dZ=du(X)-dUh 
    ie=numpy.arange(0,numpy.size(x),m) 
    print("uh error inf norm:",numpy.linalg.norm(Z,numpy.inf))
    print("duh error inf norm:",numpy.linalg.norm(dZ,numpy.inf))
    print("error at f.e. nodes inf norm:", \
        numpy.linalg.norm(u(xl)-uh[ie],numpy.inf))
    plt.figure() 
    # plot solution
    plt.plot(X,Uh,'-b',x,uh,'rx',x[ie],uh[ie],'ro')
    plt.xlabel('x'); plt.ylabel('u_h')
    plt.show()
    #plot solution error
    plt.plot(X,Z,'-b',x,z,'rx',x[ie],z[ie],'ro')
    plt.xlabel('x'); plt.ylabel('u-u_h')
    plt.show()
    #plot derivatives error
    plt.plot(X,dZ,'-b', \ 
            x,numpy.zeros(numpy.size(x)),'rx', \ 
            x[ie],numpy.zeros(numpy.size(x[ie])),'ro')
    plt.xlabel('x'); plt.ylabel('u''-u''_h'); 
    plt.show()