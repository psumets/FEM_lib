import numpy
import scipy.sparse as sparse
def assemblingAF(p,t,c,a,b1,b2,f=None):
    # assemblingAF -  vectorized algorithm of the assembling 
    #stiffness matrix A and forcing vector F for the G1 mesh.
    # Stiffness matrix A is associated with the PDE operator
    # -div(c*grad(u)) + b*grad(u)+a*u
    #The following call is allowed :
    # A=assemblingAF(p,t,c,a,b1,b2)
    #
    np=p.shape[1] 
    # mesh point indices
    k1=t[0,:] 
    k2=t[1,:]
    k3=t[2,:] 
    sdl=t[3,:] # subdomain labels
    # barycenter of the triangles
    x1=(p[0,k1]+p[0,k2]+p[0,k3])/3 
    x2=(p[1,k1]+p[1,k2]+p[1,k3])/3 
    # gradient of the basis functions, multiplied by J
    g1_x1=p[1,k2]-p[1,k3] 
    g1_x2=p[0,k3]-p[0,k2] 
    g2_x1=p[1,k3]-p[1,k1] 
    g2_x2=p[0,k1]-p[0,k3] 
    g3_x1=p[1,k1]-p[1,k2] 
    g3_x2=p[0,k2]-p[0,k1] 
    J=abs(g3_x2*g2_x1-g3_x1*g2_x2)  # J=2*area
    # evaluate c , b , a on triangles barycenter
    cf=c(x1,x2,sdl) 
    af=a(x1,x2,sdl) 
    b1f=b1(x1,x2,sdl) 
    b2f=b2(x1,x2,sdl) 
    # diagonal and off diagonal elements of mass matrix
    ao=(af/24)*J
    ad=4*ao  # 'exact' integration
    # ao=(af/18).*J ; ad=3*ao ; # quadrature rule
    # coefficients of the stiffness matrix
    cf =(0.5*cf)/J 
    a12=cf*(g1_x1*g2_x1+g1_x2*g2_x2)+ao 
    a23=cf*(g2_x1*g3_x1+g2_x2*g3_x2)+ao 
    a31=cf*(g3_x1*g1_x1+g3_x2*g1_x2)+ao
    A=[]
    F=[]
    if all(b1f==0) and all(b2f==0): # symmetric problem
        A=sparse.coo_matrix((a12,(k1,k2)),shape=(np,np)) 
        A=A+sparse.coo_matrix((a23,(k2,k3)),shape=(np,np)) 
        A=A+sparse.coo_matrix((a31,(k3,k1)),shape=(np,np)) 
        A=A+sparse.coo_matrix.transpose(A)
        A=A+sparse.coo_matrix((ad-a31-a12,(k1,k1)),shape=(np,np)) 
        A=A+sparse.coo_matrix((ad-a12-a23,(k2,k2)),shape=(np,np)) 
        A=A+sparse.coo_matrix((ad-a23-a31,(k3,k3)),shape=(np,np)) 
    else:
        # b contributions
        b1f=b1f/6 
        b2f=b2f/6 
        bg1=b1f*g1_x1+b2f*g1_x2 
        bg2=b1f*g2_x1+b2f*g2_x2 
        bg3=b1f*g3_x1+b2f*g3_x2 
        A=sparse.coo_matrix((a12+bg2,(k1,k2)),shape=(np,np)) 
        A=A+sparse.coo_matrix((a23+bg3,(k2,k3)),shape=(np,np)) 
        A=A+sparse.coo_matrix((a31+bg1,(k3,k1)),shape=(np,np)) 
        A=A+sparse.coo_matrix((a12+bg1,(k2,k1)),shape=(np,np)) 
        A=A+sparse.coo_matrix((a23+bg2,(k3,k2)),shape=(np,np)) 
        A=A+sparse.coo_matrix((a31+bg3,(k1,k3)),shape=(np,np)) 
        A=A+sparse.coo_matrix((ad-a31-a12+bg1,(k1,k1)), \ 
                                shape=(np,np)) 
        A=A+sparse.coo_matrix((ad-a12-a23+bg2,(k2,k2)), \ 
                                shape=(np,np)) 
        A=A+sparse.coo_matrix((ad-a23-a31+bg3,(k3,k3)), \ 
                                shape=(np,np)) 
    if f!=None:
        ff=f(x1,x2,sdl) 
        ff=(ff/6)*J 
        k_zero=numpy.zeros(k1.shape)
        F=sparse.coo_matrix((ff,(k1,k_zero)),shape=(np,1)) 
        F=F+sparse.coo_matrix((ff,(k2,k_zero)),shape=(np,1)) 
        F=F+sparse.coo_matrix((ff,(k3,k_zero)),shape=(np,1))
        return A,F
    else:return A