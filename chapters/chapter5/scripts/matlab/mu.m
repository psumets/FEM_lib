function f=mu(x1,x2,sdl) 
  f=zeros(size(x1)) ;
  % segments in local numeration:
  % 1-st segment 
  I=find(sdl==1);
  f(I)=2*(2+x1(I)+x2(I)).*(x1(I).^2+x2(I).^2) ;
  % 2-nd segment 
  I=find(sdl==2); 
  f(I)=2*(2+x1(I)+x2(I)+1).*(x1(I).^2+x2(I).^2) ;
end