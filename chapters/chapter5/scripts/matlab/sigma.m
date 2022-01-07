function f=sigma(x1,x2,sdl) 
  f=zeros(size(x1)) ;
  % segments in local numeration:
  % 1-st segment 
  I=find(sdl==1); f(I)=0;
  % 2-nd segment 
  I=find(sdl==2); f(I)=2 ;
end