function f=u_D(x1,x2,sdl) 
  f=zeros(size(x1)) ;
  % segments in local numeration:
  % 1-st segment 
  I=find(sdl==1); f(I)=x1(I).^2+x2(I).^2;
  % 2-nd segment 
  I=find(sdl==2); f(I)=x1(I).^2+x2(I).^2 ;
end