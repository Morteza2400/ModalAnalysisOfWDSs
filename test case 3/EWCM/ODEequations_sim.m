function dxdt= ODEequations_sim(t,x,fv,L,Pipe,Fluid,g,Node,Pump,PRV,TCV,...
    R,M_p,M_f,M_h,A1,A2,h2,Q_dem,Leak,denum,P_P,C_Nor,N_h,Q_dem_down,node_active,deltaQ)
%% prepairing value for friction term and
q(:,1)=x(1:M_p);
h1(:,1)=x(M_p+1:M_p+M_h-Node.N_Res);

[f_q]=friction(Pipe,g,q,fv,Fluid,R);
ts = 10;  % Start time of the rise
td = 6;   % Duration of the rise
Q_dem_up=Q_dem_down+deltaQ;
if t<= ts
    Q_dem(node_active) = Q_dem_down;  % Before the rise starts, Kv is constant at 50
elseif t > ts && t <= ts + td
    % During the rise, smoothly transition from 50 to 400
    Q_dem(node_active) = (Q_dem_up- Q_dem_down)* (sin(pi * (t - ts) / td - pi / 2) + 1) / 2 + Q_dem_down;
elseif t > ts + td
    Q_dem(node_active) = Q_dem_up;  % After the rise ends, Kv is constant at 400
end
% Q_dem_out=Q_dem(node_active);
%% main eq
% equations for flow of independent pipes
dxdt(1:M_p)=(-C_Nor.*(f_q).*q.*abs(q)./L+(A1'*h1./L+A2'*h2./L));
dxdt(M_p+1:M_p+M_h-Node.N_Res) = (-A1*q-Leak.*h1.^0.5-Q_dem)./(denum);
dxdt = dxdt(:);
end