function [u,x_mesh]=solveBVP(bc,c,b,a,f,xl,x,d,w)
  m=numel(x)-1; n=numel(xl);
  N=(n-1)*m+1; % number of mesh points
  x_mesh=zeros(N,1); % mesh points
  for l=1:n-1
    ne=(l-1)*m+(1:m+1); % mesh point indices
    x_mesh(ne)=xl(l)+(xl(l+1)-xl(l))*x;
  end
  [A,F]=assemblingAF(c,b,a,f,xl,x,d,w);
  [A0,F0]=assemblingBC(bc,A,F);
  u=solveFEM(bc,A0,F0);
end