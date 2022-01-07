def solveFEM(bc,A0,F0):
    u_s0 =[]; u_s1 =[] 
    if bc[0,0]==1: u_s0=bc[0,2] 
    if bc[1,0]==1: u_s1=bc[1,2]
    u0=sparse.linalg.spsolve(A0,F0)
    u=u_s0
    u=numpy.append(u,u0)
    u=numpy.append(u,u_s1)
    return u