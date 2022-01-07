function [A,F]=assemblingAF(c,b,a,f,xl,x,d,w)
  [E,D]=interpolatingMat(x,d);
  ET=E'; DT=D';
  m=numel(x)-1; n=numel(xl); n_total=(n-1)*m+1;
  A=sparse(n_total,n_total); F=zeros(n_total,1); 
  for l=1:n-1
    h=xl(l+1)-xl(l); % size of element l
    ne=(l-1)*m+(1:m+1); % point indices
    xd=xl(l)+h*d; % quadrature nodes on element l
    c_loc=diag((w/h).*c(xd));
    b_loc=diag(w.*b(xd));
    a_loc=diag((h*w).*a(xd));
    Al=(DT*c_loc+ET*b_loc)*D+ET*a_loc*E;
    Fl=(h*w).*f(xd);
    A(ne,ne)=A(ne,ne)+Al;
    F(ne)=F(ne)+ET*Fl(:);
  end
end