function A=stiffnessAssemb_1(p,t)
  np=size(p,2) ; 
  nt=size(t,2) ; 
  ntep=3; % number of triangular element points
  A=sparse(np,np) ; % global stiffness matrix
  Al=ones(ntep,ntep) ; % local stiffness matrix
  for k=1:nt
    I=t(1:ntep,k) ; % global point indices on element k
    for i=1:ntep  
      for j=1:ntep 
        A(I(i),I(j))=A(I(i),I(j))+Al(i,j) ;
      end 
    end
  end
end