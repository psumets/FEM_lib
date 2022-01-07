function A=stiffnessAssembGeneric(ntep,p,t)
  np=size(p,2) ; 
  nt=size(t,2) ; 
  Al=zeros(ntep,ntep,nt) ; % all local stiffness matrices
  for k=1:nt
    Al(k,:,:)=calcLocalStiffness(k,p,t) ; % local stiffness matrix 
  end
  A=sparse(np,np) ; % global matrix
  for i=1:ntep  
    for j=1:ntep
      A=A+sparse(t(i,:),t(j,:),Al(:,i,j),np,np) ;
    end 
  end
end