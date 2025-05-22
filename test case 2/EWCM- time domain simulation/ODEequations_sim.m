function dxdt = ODEequations_sim(t,x,fv,L,Pipe,Fluid,g,Node,...
    R,M_p,M_h,A1,A2,h2,Q_dem,Leak,denum,C_Nor,freq,d_dem,input_dem)
%% prepairing value for friction term and 
q(:,1)=x(1:M_p);
h1(:,1)=x(M_p+1:M_p+M_h-Node.N_Res);
[f_q]=friction(Pipe,g,q,fv,Fluid,R);
Q_dem(input_dem)=0.01+0.05*0.01*sin(freq*t);
%% main eq
if t>10
    dd=0;
end
% equations for flow of independent pipes
dxdt(1:M_p)=(-C_Nor.*(f_q).*q.*abs(q)./L+(A1'*h1./L+A2'*h2./L));
dxdt(M_p+1:M_p+M_h-Node.N_Res) = (-A1*q-Leak.*h1.^0.5-Q_dem-d_dem.*(h1.^0.5))./(denum);
dxdt = dxdt(:);
end