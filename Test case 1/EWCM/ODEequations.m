function dxdt = ODEequations(t,x,fv,L,Pipe,Fluid,g,Node,Pump,PRV,TCV,...
    R,M_p,M_f,M_h,A1,A2,h2,Q_dem,Leak,denum,P_P,C_Nor,N_h,d_dem)
%% prepairing value for friction term and 
q(:,1)=x(1:M_p);
h1(:,1)=x(M_p+1:M_p+M_h-Node.N_Res);
[f_q]=friction(Pipe,g,q,fv,Fluid,R);

%% main eq
% equations for flow of independent pipes
dxdt(1:M_p)=(-C_Nor.*(f_q).*q.*abs(q)./L+(A1'*h1./L+A2'*h2./L));
dxdt(M_p+1:M_p+M_h-Node.N_Res) = (-A1*q-Leak.*h1.^0.5-Q_dem-d_dem.*(h1.^0.5))./(denum);
dxdt = dxdt(:);
end