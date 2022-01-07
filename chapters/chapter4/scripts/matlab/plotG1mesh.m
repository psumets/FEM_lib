function h=plotG1mesh (p,e,t,opt,fs)
  h=figure;
  pdemesh(p,e,t); axis equal; hold on
  plot(p(1,:),p(2,:),'.k','MarkerSize',8)
  xlabel('x_1'), ylabel('x_2')
  if opt(1)==1 % output point labels
    for k=1:size(p,2)
      text(p(1,k),p(2,k),[' ' int2str(k)], ...
          'FontSize',fs,'Color','r');
    end
  end
  if opt(2)==1 % output boundary labels
    for k=1:size(e,2)
      i=e(1,k); j=e(2,k);
      x1=(p(1,i)+p(1,j))/2; x2=(p(2,i)+p(2,j))/2;
      text(x1,x2+0.05,[' ' int2str(e(5,k))], ... 
          'FontSize',fs,'Color','k');
    end
  end
  if opt(3)==1 % output element labels
    for k=1:size(t,2)
      I=t(1:3,k); x1x2=sum(p(:,I)')/3;
      text(x1x2(1),x1x2(2),int2str(k),'FontSize',fs,'Color','b');
    end
  end
end