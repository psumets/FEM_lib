function testAssemblingAF
  time=[] ; time_nv=[] ; % assembling time
  nt=[] ; np=[] ;
  g=geometryMatrix;
  % PDE coefficients;
  c=@c_coeff ; 
  a=@a_coeff; 
  b1=@b1_coeff ; 
  b2=@b2_coeff; 
  f=@f_coeff ; 
  for h=[0.1 0.05 0.02]
    [p,e,t]=initmesh(g,'hmax',h) ;
    np=[np size(p,2)] ; nt=[nt size(t,2)] ;
    tic ; 
    [A,F]=assemblingAF(p,t,c,a,b1,b2,f) ; 
    time=[time toc];
    tic ; 
    [A,F]=assemblingAF_nv(p,t,c,a,b1,b2,f) ; 
    time_nv=[time_nv toc];
  end 
  disp (' np nt assemblingAF assemblingAF_nv')
  disp([np' nt' time' time_nv'])  
end