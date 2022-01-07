function A=stiffnessAssemb_4(p,t) 
  np=size(p,2) ; 
  nt=size(t,2) ; 
  %indices of local points  
  k1=t(1,:) ; 
  k2=t(2,:) ; 
  k3=t(3,:) ; 
  % aij(m)=Al(i,j) on element m
  a12=ones(nt,1) ; A=sparse(k1,k2,a12,np,np) ;
  a13=ones(nt,1) ; A=A+sparse(k1,k3,a13,np,np) ;
  a23=ones(nt,1) ; A=A+sparse(k2,k3,a23,np,np) ;
  % replace next 3 lines by A=A+A.' for symmetric A
  a21=ones(nt,1) ; A=A+sparse(k2,k1,a21,np,np) ;
  a31=ones(nt,1) ; A=A+sparse(k3,k1,a31,np,np) ;
  a32=ones(nt,1) ; A=A+sparse(k3,k2,a32,np,np) ;
  a11=ones(nt,1) ; A=A+sparse(k1,k1,a11,np,np) ;
  a22=ones(nt,1) ; A=A+sparse(k2,k2,a22,np,np) ;
  a33=ones(nt,1) ; A=A+sparse(k3,k3,a33,np,np) ;
end