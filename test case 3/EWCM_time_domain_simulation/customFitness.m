function [fit,kv] = customFitness(x,X0,fv,L,Pipe,Fluid,g,Node,Pump,PRV,TCV,...
    R,M_p,M_f,M_h,A1,A2,h2,Q_dem,Leak,denum,P_P,C_Nor,N_h)
TCV.Kv=abs(x);
X0 = fsolve(@(x) ODEequations(0,x,fv,L,Pipe,Fluid,g,Node,Pump,PRV,TCV,...
    R,M_p,M_f,M_h,A1,A2,h2,Q_dem,Leak,denum,P_P,C_Nor,N_h), X0);
q0=X0(1:M_p,1);
h0=X0(M_p+1:M_p+M_h-Node.N_Res,1);
S0=X0(M_p+M_h-Node.N_Res+1:M_p+M_h-Node.N_Res+Pump.num,1);
Xm0=X0(M_p+M_h-Node.N_Res+Pump.num+1:M_p+M_h-Node.N_Res+Pump.num+PRV.num,1);
fit=abs(X0(1)-0.01073)*3000+abs(X0(10)-64.275)*10+abs(X0(14)-0.0066)*3000;
var=x;
end