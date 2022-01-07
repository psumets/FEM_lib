function [it,pt]=extendT(p,t)
  np=size(p,2); 
  nt=size(t,2);
  j =[1:nt;1:nt;1:nt];
  T_ext=sparse(t(1:3,:),j,1,np,nt); % extended connectivity matrix
  % coordinate form of T_ext
  [it,j]=find(T_ext'); 
  pt=find(diff([0;j;np+1]))';
end