function [e_c,t_c]=complMeshG1(p,e,t)
  ne=size(e,2); nt=size(t,2);
  tt=(t(1:3,:))';
  ee = [tt(:,[1,2]);tt(:,[2,3]);tt(:,[3,1])]; % not unique edges
  [ee,~,j]=unique(sort(ee,2),'rows'); % unique mesh edges
  ee=ee'; j=j'; nee=size(ee,2);
  mt=1:nt;
  t_c=[j(mt);j(mt+nt);j(mt+2*nt)]; % elements edges
  t_c=[t_c;t(4,:)];
  et=zeros(2,nee); % edge to element connectivity matrix
  for k=1:nt
    for j=1:3
      it=t_c(j,k);
      if et(1,it)==0, et(1,it)=k;
      else et(2,it)=k; end
    end
  end
  % sort first 2 rows of e as in ee
  for i=1:ne
    if e(1,i)>e(2,i);
      c=e(1,i); e(1,i)=e(2,i); e(2,i)=c;
      c=e(3,i); e(3,i)=e(4,i); e(4,i)=c;
      c=e(6,i); e(6,i)=e(7,i); e(7,i)=c;
    end
  end
  [~,ib,it]=intersect(ee',e(1:2,:)','rows');
  % set 3:5 rows of e_c
  e_c=[ee;zeros(5,nee)];
  e_c(3:5,ib)=e(3:5,it);
  % set 6 ,7 rows of e_c
  for ie=1:nee % ie=edge
    t1=et(1,ie); t2=et(2,ie); % neighbour triangles
    i=ee(1,ie); j=ee(2,ie); k=setdiff(tt(t1,:),[i j]); 
    % oriented area of triangle (i,j,k)
    d12=p(:,j)-p(:,i);
    d13=p(:,k)-p(:,i);
    S = d12(1,:).*d13(2,:)-d12(2,:).*d13(1,:);
    if S>0, e_c(6,ie)=t1; e_c(7,ie)=t2;
    else e_c(6,ie)=t2; e_c(7,ie)=t1; end
  end
end