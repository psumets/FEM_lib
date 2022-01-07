function f=b1_coeff(x1,x2,sdl)
  f=zeros(size(x1)) ;
  % 1-st subdomain
  I=find(sdl==1); f(I)=x1(I)+x2(I);
  % 2-nd subdomain
  I=find(sdl==2); f(I)=x1(I)-x2(I) ;
  % 3-d subdomain
  I=find(sdl==3); f(I)=x1(I).*x2(I) ;
end