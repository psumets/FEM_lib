function [A,F]=assemblingAF(p,t,c,a,b1,b2,f)
  % Vectorized assembling algorithm.
  % The following call is allowed :
  % A=assemblingAF(p,t,c,a,b1,b2)
  np=size(p,2) ;
  k1=t(1,:) ; k2=t(2,:) ; k3=t(3,:) ; % mesh point indices
  sdl=t(4,:) ; % subdomain labels
  % barycenter of the triangles
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
  % diagonal and off diagonal elements of mass matrix
  ao=(af/24).*J ; ad=4*ao ; % 'exact' integration
  % ao=(af/18).*J ; ad=3*ao ; % quadrature rule
  % coefficients of the stiffness matrix
  cf =(0.5*cf)./J ;
  a12=cf.*(g1_x1.*g2_x1+g1_x2.*g2_x2)+ao ;
  a23=cf.*(g2_x1.*g3_x1+g2_x2.*g3_x2)+ao ;
  a31=cf.*(g3_x1.*g1_x1+g3_x2.*g1_x2)+ao ;
  if all(b1f==0) && all( b2f==0) % symmetric problem
    A=sparse(k1,k2,a12,np,np) ;
    A=A+sparse(k2,k3,a23,np,np) ;
    A=A+sparse(k3,k1,a31,np,np) ;
    A=A+A.';
    A=A+sparse(k1,k1,ad-a31-a12,np,np) ;
    A=A+sparse(k2,k2,ad-a12-a23,np,np) ;
    A=A+sparse(k3,k3,ad-a23-a31,np,np) ;
  else
    % b contributions
    b1f=b1f/6 ; b2f=b2f/6 ;
    bg1=b1f.*g1_x1+b2f.*g1_x2 ;
    bg2=b1f.*g2_x1+b2f.*g2_x2 ;
    bg3=b1f.*g3_x1+b2f.*g3_x2 ;
    A=sparse(k1,k2,a12+bg2,np,np) ;
    A=A+sparse(k2,k3,a23+bg3,np,np) ;
    A=A+sparse(k3,k1,a31+bg1,np,np) ;
    A=A+sparse(k2,k1,a12+bg1,np,np) ;
    A=A+sparse(k3,k2,a23+bg2,np,np) ;
    A=A+sparse(k1,k3,a31+bg3,np,np) ;
    A=A+sparse(k1,k1,ad-a31-a12+bg1,np,np) ;
    A=A+sparse(k2,k2,ad-a12-a23+bg2,np,np) ;
    A=A+sparse(k3,k3,ad-a23-a31+bg3,np,np) ;
  end
  if nargout==2
    ff=feval(f,x1,x2,sdl) ;
    ff=(ff/6).*J ;
    F=sparse(k1,1,ff,np,1) ;
    F=F+sparse(k2,1,ff,np,1) ;
    F=F+sparse(k3,1,ff,np,1) ;
  else
    F=[] ;
  end
end