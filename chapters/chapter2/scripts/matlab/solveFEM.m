function u = solveFEM(bc,A0,F0)
  u_s0 =[]; u_s1=[];
  if bc(1,1)==1 , u_s0=bc(1,3) ; end
  if bc(2,1)==1 , u_s1=bc(2,3) ; end
  u0=A0\F0;
  u=[u_s0;u0;u_s1];
end