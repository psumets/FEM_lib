function h=plotComplMeshG1(p,e,t,e_c,fs)
  h=figure; pdemesh(p,e,t); axis equal; hold on
  plot(p(1,:),p(2,:),'.k','MarkerSize',8);
  xlabel('x_1'), ylabel('x_2')
  for i = 1:size(p,2)
    text(p(1,i),p(2,i),[' ' int2str(i)], ...
        'FontSize',fs-1,'Color','k');
  end
  for it=1:size(e_c,2) % edge labels
    i=e_c(1,it); j=e_c(2,it);
    text((p(1,i)+p(1,j))/2-0.06,(p(2,i)+p(2,j))/2, ...
    [' ' int2str(it)],'FontSize',fs-2,'Color','r');
    if e_c(5,it);
      text((p(1,i)+p(1,j))/2+0.05,(p(2,i)+p(2,j))/2+0.05, ...
      [' ' int2str(e_c(5,it))],'FontSize',fs-2,'Color','k');
    end
  end
  for i=1:size(t,2) % element labels
    I=t(1:3,i);
    text(sum(p(1,I))/3,sum(p(2,I))/3,int2str(i), ...
        'FontSize',fs-1,'Color','b');
  end
end