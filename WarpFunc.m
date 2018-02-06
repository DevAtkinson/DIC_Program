function Wout = WarpFunc(dx,dy,in3)
%WARPFUNC
%    WOUT = WARPFUNC(DX,DY,IN3)

%    This function was generated by the Symbolic Math Toolbox version 7.2.
%    05-Feb-2018 11:26:07

P1 = in3(:,1);
P2 = in3(:,2);
P3 = in3(:,3);
P4 = in3(:,4);
P5 = in3(:,5);
P6 = in3(:,6);
Wout = [P1+P3.*dy+dx.*(P2+1.0),P4+P5.*dx+dy.*(P6+1.0)];
