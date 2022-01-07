function F=forcingAssemb_1(p,t)
  np=size(p,2) ; nt=size(t,2) ;
  k1=t(1,:) ; k2=t(2,:) ; k3=t(3,:) ;
  % Fl - local forcing vector
  %fj(m)= Fl(j) on element m, 
  f1=ones(nt,1) ; F=sparse(k1,1,f1,np,1) ;
  f2=ones(nt,1) ; F=F+sparse(k2,1,f2,np,1) ;
  f3=ones(nt,1) ; F=F+sparse(k3,1,f3,np,1) ;
end