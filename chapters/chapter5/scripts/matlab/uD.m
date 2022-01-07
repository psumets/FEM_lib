function f=uD(x1,x2,sdl) 
  f=zeros(size(x1)) ;
  % 1-st segment (local enumeration)
  I=find(sdl==1); f(I)=1;
  % 2-nd segment
  I=find(sdl==2); f(I)=x1(I)+x2(I) ;
  % 3-d segment
  I=find(sdl==3); f(I)=1;
end