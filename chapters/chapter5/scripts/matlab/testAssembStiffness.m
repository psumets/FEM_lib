function testAssembStiffness
  clc
  % assembling time
  time_1=[] ; time_2=[] ; time_3=[] ; time_4=[] ; time_5=[] ; 
  nt=[] ; np=[] ;
  for nx=[10 50 100]
    % set a mesh of rectangular domain
    [p,e,t] =poimesh('squareg',nx,nx) ;
    np=[np size(p,2)] ; 
    nt=[nt size(t,2)] ; 
    tic , stiffnessAssemb_1(p,t) ; time_1=[time_1 toc] ;
    tic , stiffnessAssemb_2(p,t) ; time_2=[time_2 toc] ;
    tic , stiffnessAssemb_3(p,t) ; time_3=[time_3 toc] ;
    tic , stiffnessAssemb_4(p,t) ; time_4=[time_4 toc] ;
    tic , stiffnessAssemb_5(p,t) ; time_5=[time_5 toc] ;
  end
  disp('np  nt  ver_1   ver_2   ver_3   ver_4   ver_5')
  disp([np' nt' time_1' time_2' time_3' time_4' time_5'])
end