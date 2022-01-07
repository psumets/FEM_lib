function A=stiffnessAssemb_5(p,t) 
  np=size(p,2) ; 
  nt=size(t,2) ; 
  ntep=3;
  %indices of local points
  k1=t(1,:) ; 
  k2=t(2,:) ; 
  k3=t(3,:) ; 
  Al=ones(nt,ntep,ntep) ; % all local stiffness matrices
  A=sparse(k1,k2,Al(:,1,2),np,np) ;
  A=A+sparse(k1,k3,Al(:,1,3),np,np) ;
  A=A+sparse(k2,k3,Al(:,2,3),np,np) ;
  % replace next 3 lines by A=A+A.'for symmetric A
  A=A+sparse(k2,k1,Al(:,2,1),np,np) ;
  A=A+sparse(k3,k1,Al(:,3,1),np,np) ;
  A=A+sparse(k3,k2,Al(:,3,2),np,np) ;
  A=A+sparse(k1,k1,Al(:,1,1),np,np) ;
  A=A+sparse(k2,k2,Al(:,2,2),np,np) ;
  A=A+sparse(k3,k3,Al(:,3,3),np,np) ;
end