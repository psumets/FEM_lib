% main script
clc ; clear all ; close all
% set BVP data:
% integration domain
s_0=-1; s_1=1; 
% exact solution
u=@(x) sin(x); 
du=@(x) cos(x);
% equation coefficients
c=@(x) cos(x);
b=@(x) -2*sin(x);
a=@(x) sin(x);
f=@(x) sin(x).^2;
% boundary conditions
sigma=0;
mu=c(s_1)*du(s_1)+sigma*u(s_1);
bc=[1 0 u(s_0)
    3 sigma mu ];
% set FEM parameters
m=3; % polynomial degree
ne=6; % number of final elements
M=m+1; % number of quadrature nodes
xl=linspace(s_0,s_1,ne+1); % finite element nodes
t=LobattoQuad(m+1,0,1); % interpolation nodes on basis element
[d,w]=LobattoQuad(M,0,1); % quadrature nodes and weighs
% solve BVP
tic
[uh,x]=solveBVP(bc,c,b,a,f,xl,t,d,w);
timeFemSol=toc
% find error in various norms
[e0,e1]=errorNorm(u,du,uh,xl,t)
[X,Uh,dUh]=getRefinedValues(uh,xl,t,10*m);
Z=u(X)-Uh; z=u(x)-uh; dZ=du(X)-dUh; 
ie=1:m:numel(x);
einfh=norm(Z,inf)
einfdh=norm(dZ,inf)
einfbh=norm(u(xl(:))-uh(ie),inf)
figure; % plot solution
plot(X,Uh,'-b',x,uh,'rx',x(ie),uh(ie),'ro')
xlabel('x'); ylabel('u_h')
figure;% plot solution error
plot(X,Z,'-b',x,z,'rx',x(ie),z(ie),'ro')
xlabel('x'); ylabel('u-u_h')
figure;% plot derivatives error
plot(X,dZ,'-b',x,zeros(size(x)),'rx',x(ie),zeros(size(x(ie))),'ro')
xlabel('x'); ylabel('u''-u''_h')