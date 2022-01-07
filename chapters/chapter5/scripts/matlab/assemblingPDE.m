function [A,F,N,UD,S,M]=assemblingPDE(bc,p,e,t,c,a,b1,b2,f)
  %[A,F,N,UD,S,M]=assemblingPDE(bc,p,e,t,c,a,b1,b2,f)  - returns 
  % FEM matrices.
  % u=assemblingPDE(bc,p,e,t,c,a,b1,b2,f) - returns 
  % FEM solution represented as a column vector.
  if nargout==6
    [A,F]=assemblingAF(p,t,c,a,b1,b2,f) ;
    [N,S,UD,M]=assemblingBC(bc,p,e) ;
  elseif nargout==1
    [A,F]=assemblingAF(p,t,c,a,b1,b2,f) ;
    [N,S,UD,M]=assemblingBC(bc,p,e) ;
    if size(N,2)==size(p,2) % no Dirichlet boundary conditions
      A=A+S; 
      F=F+M; 
    else
      Nt=N.';
      A=A+S;
      F=Nt*((F+M)-A*UD) ; 
      A=Nt*A*N; 
    end
    u=A\F;
    A=N*u+UD; % returns A=u
  else
    error('Wrong number of output parameters.' ) ;
  end
end
