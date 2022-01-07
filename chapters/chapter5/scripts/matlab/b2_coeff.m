function f=b2_coeff(x1,x2,sdl)
  f=zeros(size(x1)) ;
  % 1-st subdomain
  I=find(sdl==1); f(I)=1;
  % 2-nd subdomain
  I=find(sdl==2); f(I)=1 ;
  % 3-d subdomain
  I=find(sdl==3); f(I)=1 ;
end