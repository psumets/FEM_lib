function [x1,x2]= geometryFunction (b_index, s)
  % Geometry data 
  s_param=[ 0 0 0 0 0 0 0 0 0 0 0     pi 3/2*pi    0 pi/2
         1 1 1 1 1 1 1 1 1 1 1 3/2*pi   2*pi pi/2   pi
         0 3 3 1 3 3 0 0 0 3 3      0      0    0    0
         1 2 2 0 0 2 1 3 3 2 1      1      1    1    1 ];
  x1x2=[2 1.25 1.25  2 2 0.75  0 0 0 0.75 0 
        0 1.25 0.75  2 2 1.25  0 0 2 0.75 2 
       -1 0.75 0.25 -1 0 0.75 -1 0 1 0.25 0 
       -1 0.25 0.25  0 1 0.75  0 1 1 0.75 0 ]; 
  if nargin==0, x1=size(s_param,2); return; end
  if nargin==1, x1=s_param(:,b_index); return; end
  x1=zeros(size(s)); 
  x2=zeros(size(s)); 
  [m,n]=size(b_index);
  if m==1 && n==1, b_index=b_index*ones(size(s)); end
  if ~isempty(s),
    % line segments
    for k=1:11 
      i=find(b_index==k);
      if ~isempty(i)
        x1(i)=x1x2(1,k)+(x1x2(2,k)-x1x2(1,k))*s(i);
        x2(i)=x1x2(3,k)+(x1x2(4,k)-x1x2(3,k))*s(i);
      end
    end
    % circle segments
    for k=12:15 
      i=find(b_index==k);
      if ~isempty(i)
        c_x1=1; c_x2=-0.5; rc=0.25;
        x1(i)=rc*cos(s(i))+c_x1;
        x2(i)=rc*sin(s(i))+c_x2;
      end
    end
  end
end