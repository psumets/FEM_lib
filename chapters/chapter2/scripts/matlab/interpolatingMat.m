function [E,D]=interpolatingMat(x,d)
  nx=numel(x); nd=numel(d);
  E=zeros(nd,nx); D=zeros(nd,nx); b=zeros(1,nx);
  for i=1:nx
    j=[1:(i-1),(i+1):nx];
    b(i)=1/prod(x(i)-x(j));
  end
  for i=1:nx
    j=[1:(i-1),(i+1):nx];
    for k=1:nd
      E(k,i)=b(i)*prod(d(k)-x(j));
      ds=0;
      for l=j
        i1=min(i,l); i2=max(i,l);
        jj=[1:(i1-1),(i1+1):(i2-1),(i2+1):nx];
        ds=ds+prod(d(k)-x(jj));
      end
      D(k,i)=b(i)*ds;
    end
  end
end