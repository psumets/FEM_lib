function F=forcingAssemb_3(p,t)
  np=size(p,2); nt=size(t,2) ; ntep=3;
  F=zeros(np,1) ; 
  for k=1:nt
    I=t(1:ntep,k) ;
    Fl=ones(ntep,1) ; 
    F(I)=F(I)+Fl ;
  end
end