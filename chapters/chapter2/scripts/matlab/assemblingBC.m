function [A0,F0]=assemblingBC(bc,A,F)
  %   bc - boundary conditions matrix of size 2x3 :
  %   bc(1,:)=[1,0,u_s0] , or bc(1,:)=[3,sig_s0,mu_s0]
  %   bc(2,:)=[1,0,u_s1] , or bc(2,:)=[3,sig_s1,mu_s1]
   n=numel(F);% number of all mesh points
  if bc(1,1)==3
    A(1,1)=A(1,1)+bc(1,2);
    F(1)=F(1)+bc(1,3);
  end
  if bc(2,1)==3
    A(n,n)=A(n,n)+bc(2,2);
    F(n)=F(n)+bc(2,3);
  end
  ib=1; ie=n;
  if bc(1,1)==1 
    ib=2; u_s0=bc(1,3); 
    F=F-u_s0*A(:,1); 
  end
  if bc(2,1)==1 
    ie=n-1; u_s1=bc(2,3); 
    F=F-u_s1*A(:,n); 
  end
  I=(ib:ie)';
  A0=A(I,I);
  F0=F(I);
end