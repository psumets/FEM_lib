function [x,y,dy]=getRefinedValues(u,xl,d,k)
  d_ref=linspace(0,1,k+1); % refined basis element
  [E,D]=interpolatingMat(d,d_ref);
  n=numel(xl); m=numel(d)-1;
  x_ref=(n-1)*k+1; % number of all refined mesh points
  x=zeros(x_ref,1); y=zeros(size(x)); dy=zeros(size(x));
  % Iterate over finite elements and refine them 
  for l=1:n-1
    h=xl(l+1)-xl(l); % size of element l
    i_c=(l-1)*m+(1:m+1); % coarse point indices
    i_ref=(l-1)*k+(1:k+1); % refined points indices
    x(i_ref)=linspace(xl(l),xl(l+1),k+1)';
    y(i_ref)=E*u(i_c);
    dy(i_ref)=D*u(i_c)/h;
  end
end