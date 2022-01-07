import numpy
import scipy.sparse as sparse
from scipy.sparse.linalg import spsolve 
def assemblingPDE(bc,p,e,t,c,a,b1,b2,f,nargout):
    # A,F,N,UD,S,M=assemblingPDE(bc,p,e,t,c,a,b1,b2,f,nargout)  - 
    #returns FEM matrices of the PDE problem.
    # u=assemblingPDE(bc,p,e,t,c,a,b1,b2,f,nargout) - 
    #returns solution to the PDE problem
    A=[]; F=[] ; UD=[]; S=[]; M=[]
    if nargout==6:
        A,F=assemblingAF(p,t,c,a,b1,b2,f)
        N,S,UD,M=assemblingBC(bc,p,e,4)
        return A,F,N,UD,S,M
    elif (nargout == 1):
        A,F=assemblingAF(p,t,c,a,b1,b2,f) 
        N,S,UD,M=assemblingBC(bc,p,e,4)
        if N.shape[1]==p.shape[1]: # no Dirichlet BC
            A=A+S; F=F+M 
        else:
            Nt=sparse.coo_matrix.transpose(N)
            A=A+S
            F=Nt @ ((F+M)-A @ UD)  
            A=Nt @ A @ N 
        u=sparse.linalg.spsolve(A,F)
        u=(sparse.coo_matrix(u)).reshape((u.shape[0],1))
        return (N @ u+UD)
    else: print('assembpde: Wrong number of output parameters.') 