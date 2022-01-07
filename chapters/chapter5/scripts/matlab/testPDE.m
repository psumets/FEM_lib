function testPDE
  clear all ;close all ; clc
  g='circleg' ; np=[] ; nt=[] ; errh2=[] ; err=[] ;
  t_assemblingAF=[] ; t_assemblingBC=[] ; t_solving=[] ;
  for h=[0.5 0.1 0.05 0.02]
    [p,e,t]=initmesh(g,'hmax',h) ;
    exact=@(x1,x2,sdl) x1.^2+x2.^2; % exact solution
    x1=p(1,:) ; x2=p(2,:) ;
    u_exact=exact(x1,x2,1)';
    % PDE coefficients
    c=@(x1,x2,sdl) 2+x1+x2 ;
    a=@(x1,x2,sdl) x1+x2 ; 
    f=@(x1,x2,sdl) -8-6*(x1+x2)+(2+x1+x2).*(x1.^2+x2.^2) ;
    b1=@(x1,x2,sdl) x1 ; 
    b2=@(x1,x2,sdl) x2 ;
    % boundary conditions
    bc.bsD=[1 3] ; bc.uD=@u_D ;
    bc.bsR=[2 4] ; bc.sg=@sigma ; bc.mu=@mu ;
    tic
    [A,F]=assemblingAF(p,t,c,a,b1,b2,f) ;
    t_assemblingAF=[t_assemblingAF toc] ;
    tic 
    [N,S,UD,M]=assemblingBC(bc,p,e) ; 
    t_assemblingBC=[t_assemblingBC toc] ;
    tic
    u=assemblingPDE(bc,p,e,t,c,a,b1,b2,f) ;
    t_solving =[t_solving toc] ;
    np=[np size(p,2)] ; nt=[nt size(t,2)] ; 
    err=[err norm(u_exact-u,inf)] ;
    errh2=[errh2 norm(u_exact-u,inf)/h^2] ;
  end
  disp ('np   assemblingAF    assemblingBC    solving')
  disp ( [np' t_assemblingAF' t_assemblingBC' t_solving'] )
  err
  errh2
end