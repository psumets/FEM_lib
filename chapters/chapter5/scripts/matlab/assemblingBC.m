function [N,S,UD,M]=assemblingBC(bc,p,e)
  % The following call is also allowed :
  % [N,S]=assemblingBC(bc,p,e)
  % [N,S,UD]=assemblingBC(bc,p,e)
  np=size(p,2) ;
  S=sparse(np,np) ; M=sparse(np,1) ; N=speye(np,np) ;
  UD=sparse(np,1) ;
  if all(bc.bsR==inf) && (numel(bc.sg)==0) && (numel(bc.mu)==0)
      return % all boundary are homogeneous Robin BC
  end
  e=e(:,e(6,:)==0 | e (7,:)==0) ; % all boundary edges
  k=e(5,:) ;
  [bsg,~]=find(sparse(k,k,1,np,np)) ; % boundary segment indices
  nbsg=numel(bsg) ;
  % set localal enumeration of boundary segments
  local=zeros(1,nbsg) ;
  % set Robin boundary conditions
  if nargout>=2
  bsR=bc.bsR ;
  if numel(bsR)>0
    if bsR==inf , bsR=bsg ; end % all BC are mixed
      local(bsR)=1:numel(bsR) ;
      % find boudary edges with Robin b.c.
      eR=e([1 2 5],ismember(k,bsR)) ;
      k1=eR(1,:) ; k2=eR(2,:) ; % corner point indices
      x1=0.5*(p(1,k1)+p(1,k2)) ; x2=0.5*(p(2,k1)+p(2,k2)) ;
      h=sqrt((p(1,k2)-p(1,k1)).^2+(p(2,k2)-p(2,k1)).^2) ; 
      sdl=local(eR(3,:)) ;
      % evaluate sigma and mu on edges barycenter
      sf=feval(bc.sg,x1,x2,sdl) ;
      mf=feval(bc.mu,x1,x2,sdl) ;
      % diagonal and off diagonal elements 
      so=(sf/6).*h ; sd=2*so ; % 'exact' integration
      %so=(sf/4).*h ; sd = so ; % quadrature rule
      S=sparse(k1,k2,so,np,np) ;
      S=S+sparse(k2,k1,so,np,np) ;
      S=S+sparse(k1,k1,sd,np,np) ;
      S=S+sparse(k2,k2,sd,np,np) ;
      if nargout==4
        mf=feval(bc.mu,x1,x2,sdl) ;
        mf=(mf/2).*h ;
        M=sparse(k1,1,mf,np,1) ;
        M=M+sparse(k2,1,mf,np,1) ;
      end
    end
  end
  % set Dirichlet boundary conditions
  bsD=bc.bsD ;
  if numel (bsD)>0
    if bsD==inf  
      bsD=bsg ; 
    end % all b.c. are the Dirichlet one
    local(bsD)=1:numel(bsD) ;
    if all(local==0)
      disp ('error. bsD+bsR~=number of boundary segments')
    end
    eD=e([1 2 5],ismember(k,bsD)) ; % Dirichlet BC edges
    sdl=local(eD(3,:)) ;
    k1=eD(1,:) ; k2=eD(2,:) ; % indices of corner points
    iD=[k1 k2] ; 
    [id,~]=find(sparse(iD,iD,1,np,np)) ; % Dirichlet point indices
    iN=ones(1,np) ; iN(id)=zeros(1,numel(id)) ;
    iN=find(iN) ; % indices of non-Dirichlet points 
    niN=numel(iN) ;
    N=sparse(iN,1:niN,1,np,niN) ;
    if nargout>=3 % evaluate UD on Dirichlet points
      UD(k1)=feval(bc.uD,p(1,k1),p(2,k1),sdl) ;
      UD(k2)=feval(bc.uD,p(1,k2),p(2,k2),sdl) ;
    end
  end
end
