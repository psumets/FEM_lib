function A=stiffnessAssemb_2(p,t) 
  np=size(p,2) ; 
  nt=size(t,2) ; 
  ntep=3;
  A=sparse(np,np) ; % global stiffness matrix
  Al=ones(ntep,ntep) ; % local stiffness matrix
  for k=1:nt
    I=t(1:ntep,k) ;   
    A(I,I)=A(I,I)+Al ; % assembling matrix
  end
end