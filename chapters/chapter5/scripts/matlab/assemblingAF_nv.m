function [A,F]=assemblingAF_nv(p,t,c,a,b1,b2,f)
  % Non-vectorized assembling algorithm.
  % The following call is allowed :
  % A=assemblingAF_nv(p,t,c,a,b1,b2)
  np=size(p,2) ; nt=size(t,2) ;
  Al=zeros(nt,3,3) ; % all local stiffness matrix
  F=zeros(np,1) ; % global force vector
  for k=1:nt
    % mesh point indices
    k1=t(1,k) ; k2=t(2,k) ; k3=t(3,k) ;
    sdl=t(4,k) ; % subdomain labels
    % barycenter of the triangle
    x1=(p(1,k1)+p(1,k2)+p(1,k3))/3 ;
    x2=(p(2,k1)+p(2,k2)+p(2,k3))/3 ;
    % gradient of the basis functions, multiplied by J
    g1_x1=p(2,k2)-p(2,k3) ; g1_x2=p(1,k3)-p(1,k2) ;
    g2_x1=p(2,k3)-p(2,k1) ; g2_x2=p(1,k1)-p(1,k3) ;
    g3_x1=p(2,k1)-p(2,k2) ; g3_x2=p(1,k2)-p(1,k1) ;
    J=abs(g3_x2.*g2_x1-g3_x1.*g2_x2) ; % J=2*area
    % evaluate c , b , a on triangles barycenter
    cf=feval(c,x1,x2,sdl) ;
    af=feval(a,x1,x2,sdl) ;
    b1f=feval(b1,x1,x2,sdl) ;
    b2f=feval(b2,x1,x2,sdl) ;
    % diagonal and off-diagonal elements of mass matrix
    ao=(af/24).*J ; ad=4*ao ; % 'exact' integration
    % ao=(af/18).*J ; ad=3*ao ; % quadrature rule
    % b contributions
    b1f=b1f/6 ; b2f=b2f/6 ;
    bg1=b1f.*g1_x1+b2f.*g1_x2 ;
    bg2=b1f.*g2_x1+b2f.*g2_x2 ;
    bg3=b1f.*g3_x1+b2f.*g3_x2 ;
    % coefficients of the stiffness matrix
    cf=(0.5*cf)./J ;
    a12=cf.*(g1_x1.*g2_x1+g1_x2.*g2_x2)+ao ;
    a23=cf.*(g2_x1.*g3_x1+g2_x2.*g3_x2)+ao ;
    a31=cf.*(g3_x1.*g1_x1+g3_x2.*g1_x2)+ao ;
    Al(k,:,:)=[ad-a31-a12+bg1 a12+bg2 a31+bg3 
               a12+bg1 ad-a12-a23+bg2 a23+bg3
               a31+bg1 a23+bg2 ad-a23-a31+bg3] ;
    if nargout==2
      ff=feval(f,x1,x2,sdl) ;
      I=[k1;k2;k3] ;
      F(I)=F(I)+ff*J/6 ;
    end
  end
  A=sparse(np,np) ; % global stiffness matrix
  for i=1:3  
    for j=1:3
      A=A+sparse(t(i,:),t(j,:),Al(:,i,j),np,np) ;
    end  
  end
  if nargout==1 
    F=[] ; 
  end
end
