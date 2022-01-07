import numpy
import scipy.sparse as sparse
def assemblingBC(bc,p,e,nargout):
    # assemblingBC - assembling boundary conditions
    # The following call is also allowed :
    # N,S = assemblingBC(bc,p,e,nargout)
    # N,S,UD = assemblingBC(bc,p,e,nargout)
    np=p.shape[1]
    S=[] 
    M=[]  
    N=sparse.eye(np)
    UD=[]
    if (bc["bsR"]==None) and \ 
       (bc["sg"]==None) and \ 
       (bc["mu"]==None):
        return # all boundaries are homogeneous Robin b.c.
    e=e[:,numpy.logical_or(e[5,:]==0,e[6,:]==0)] #boundary edges
    k=e[4,:] 
    bsg=numpy.transpose(numpy.nonzero( \
                        sparse.coo_matrix( \
                (numpy.ones(len(k)),(k,k)),shape=(np,np))))[:,0]
    nbsg=len(bsg) 
    # set local numeration of boundary segments
    local=numpy.zeros((nbsg)) 
    # set mixed boundary conditions
    if nargout>=2:
        bsR=bc["bsR"] 
        if len(bsR)>0:
            if bsR==None: # all BC are mixed 
                bsR=bsg
            local[bsR]=range(0,len(bsR)) 
            # find boudary edges with mixed b.c.
            eR=e[:,numpy.isin(k,bsR)]
            eR=eR[[0, 1, 4],:]
            k1=eR[0,:] 
            k2=eR[1,:] # indices of starting points and endpoints
            x1=0.5*(p[0,k1]+p[0,k2]) 
            x2=0.5*(p[1,k1]+p[1,k2]) 
            h=numpy.sqrt(numpy.square(p[0,k2]-p[0,k1])
                +numpy.square(p[1,k2]-p[1,k1])) # edges length
            sdl=local[eR[2,:]] 
            # evaluate sigma , mu on edges barycenter
            sf=bc["sg"](x1,x2,sdl) 
            mf=bc["mu"](x1,x2,sdl) 
            # 'exact' integration
            so=(sf/6)*h 
            sd=2*so
            # quadrature rule
            #so=(sf/4)*h 
            #sd = so  
            S=sparse.coo_matrix((so,(k1,k2)),shape=(np,np)) 
            S=S+sparse.coo_matrix((so,(k2,k1)),shape=(np,np)) 
            S=S+sparse.coo_matrix((sd,(k1,k1)),shape=(np,np)) 
            S=S+sparse.coo_matrix((sd,(k2,k2)),shape=(np,np)) 
            if nargout==4:
                mf=bc["mu"](x1,x2,sdl) 
                mf=(mf/2)*h 
                k_zero=numpy.zeros(k1.shape)
                M=sparse.coo_matrix((mf,(k1,k_zero)), \ 
                                    shape=(np,1)) 
                M=M+sparse.coo_matrix((mf,(k2,k_zero)), \ 
                                    shape=(np,1)) 
    # set Dirichlet boundary conditions
    bsD=bc["bsD"] 
    if len(bsD)>0:
        if bsD==None: # all b.c. are Dirichlet
            bsD=bsg 
        local[bsD]=range(0,len(bsD)) 
        if all(local==0):
            print('error. bsD+bsR~=number of boundary segments')
        eD=e[:,numpy.isin(k,bsD)] # boudary with Dirichlet BC
        eD=eD[[0, 1, 4],:]
        sdl=local[eD[2,:]] 
        # indices of start and end points
        k1=eD[0,:] 
        k2=eD[1,:] 
        iD=numpy.concatenate([k1, k2]) # Dirichlet points indices
        i_d=numpy.transpose(numpy.nonzero(
            sparse.coo_matrix(
                (numpy.ones(len(iD)),(iD,iD)), \ 
                shape=(np,np))))[:,0]
        # iN - indices of non−Dirichlet points
        iN=numpy.ones((np)) 
        iN[i_d]=numpy.zeros((len(i_d))) 
        iN=numpy.transpose(numpy.nonzero(iN)).flatten() 
        niN=len(iN) 
        N=sparse.coo_matrix((numpy.ones(niN), \ 
                    (iN,range(0,niN))),shape=(np,niN)) 
        if nargout>=3: # evaluate UD on Dirichlet points
            UD=sparse.csr_matrix((1,np))
            UD[0,k1]=bc["uD"](p[0,k1],p[1,k1],sdl)
            UD[0,k2]=bc["uD"](p[0,k2],p[1,k2],sdl) 
            UD=UD.reshape((np,1))
    if nargout==2:return N,S
    elif nargout==3:return N,S,UD
    elif nargout==4:return N,S,UD,M
    else:print('assemblingBC: Wrong number of output parameters.') 