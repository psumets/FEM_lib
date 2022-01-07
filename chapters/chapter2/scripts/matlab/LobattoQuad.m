function [d,w] = LobattoQuad (n,s0,s1)
  if n<2
    disp ('LobattoQuad: number of nodes is less then 2' ); 
    return;
  end
  if nargin==1
    s0=-1; s1=1; 
  end
  if n==2
    d=[s0;s1]; w=[(s1-s0)/2; (s1-s0)/2]; 
    return; 
  end
  j=1:n-3;
  bta=sqrt(j.*(j+2)./(2*j+1)./(2*j+3));
  Q=diag(bta,-1)+diag(bta,1);
  [U,L]=eig(Q);
  [d,k]=sort(diag(L));
  w=4/3*(U(1,k).^2)';
  w=[2/(n^2-n); w./(1-d.^2); 2/(n^2-n)];
  d=[-1; d ; 1 ];
  d=s0+0.5*(s1-s0)*(d+1);
  w=0.5*(s1-s0)*w;
end