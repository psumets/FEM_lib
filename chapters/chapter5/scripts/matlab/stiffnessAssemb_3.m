function A=stiffnessAssemb_3(p,t)
  np=size(p,2) ; 
  nt=size(t,2) ; 
  ntep=3;
  ntep2=ntep^2; % number of local stiffness matrix elements
  m=ntep2*nt ; % number of elemnts of all local matrices
  % coordinates representation of matrix A
  i=zeros(m,1) ; 
  j=zeros(m,1) ; 
  v=ones(m,1) ;
  Al=rand(ntep,ntep) ; % local stiffness matrix
  % assembling matrix A in coordinate form (i,j,v)
  for k=1:nt
    I=t(1:ntep,k) ; 
    il=repmat(I(:),1,ntep) ; 
    jl=il';
    m=ntep2*(k-1)+(1:ntep2) ;
    i(m)=il(:) ;
    j(m)=jl(:) ;
    v(m)=Al(:) ;
  end
  A=sparse(i,j,v,np,np) ; 
end