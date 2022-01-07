function f=c_coeff_def(x1,x2,sdl)
  f=zeros(size(x1)) ;
  % 1-st subdomain
  I=find(sdl==1); f(I)=1;
  % 2-nd subdomain
  I=find(sdl==2); f(I)=x1(I).^2+x2(I) ;
  % 3-d subdomain
  I=find(sdl==3); f(I)=sin(x1(I)+x2(I)) ; 
end