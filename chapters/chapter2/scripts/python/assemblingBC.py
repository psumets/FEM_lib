def assemblingBC(bc,A,F):
    # bc - boundary conditions matrix of size 2x3 :
    # bc(1,:)=[1,0,u_s0] , or bc(1,:)=[3,sig_s0,mu_s0]
    # bc(2,:)=[1,0,u_s1] , or bc(2,:)=[3,sig_s1,mu_s1]
    n=F.shape[0]-1 #number of all mesh points
    # set boundary conditions of 3-d type
    if bc[0,0]==3:
        A[0,0]=A[0,0]+bc[0,1]
        F[0]=F[0]+bc[0,2]
    if bc[1,0]==3:
        A[n,n]=A[n,n]+bc[1,1] 
        F[n,0]=F[n,0]+bc[1,2]
    #set boundary conditions of the first type
    ib=0; ie=n+3
    if bc[0,0]==1: 
        ib=1; u_s0=bc[0,2] 
        F=F-u_s0*A[:,0]
    if bc[1,0]==1: 
        ie=n+2; u_s1=bc[1,2] 
        F=F-u_s1*A[:,n]
    A0=A[ib:ie,ib:ie] 
    F0=F[ib:ie,0]
    return A0,F0