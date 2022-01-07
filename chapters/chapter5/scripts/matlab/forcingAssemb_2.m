function F=forcingAssemb_2(p,t)
  np=size(p,2) ; nt=size(t,2) ;
  k1=t(1,:) ; k2=t(2,:) ; k3=t(3,:) ;
  Fl=ones(nt,3) ; % all local forcing vectors
  F=sparse(k1,1,Fl(:,1),np,1);
  F=F+sparse(k2,1,Fl(:,2),np,1) ;
  F=F+sparse(k3,1,Fl(:,3),np,1) ;
end