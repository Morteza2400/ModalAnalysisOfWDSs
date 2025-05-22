clc
clear
close all;
load 'input.mat'
[t_total,deltaT,Fluid,Node,Pipe,g,PRV,TCV,Pump]=read_data(data);

% Pre-process of data
[Pipe,Node,Fluid,A,A1,A2,Res_row,B,B1,q0,h0,h1,...
    h2,q,L,R,C,M_p,M_h,P_P,N_h,C_PRV,C_pump,C_TCV,C_Nor,RR,A_save]=PreProcess(Pipe,Fluid,Node,g);
%% initial condition
M_f=M_p-Pipe.ctr;
X0=zeros(M_p+M_h-Node.N_Res+Pump.num+PRV.num,1);
X0(1:M_p,1)=q0(1:M_p,1);
X0(M_p+1:M_p+M_h-Node.N_Res)=h1;
Leak=zeros(M_h-Node.N_Res,1);
C_T=zeros(M_h-Node.N_Res,1);
Q_dem=zeros(M_h-Node.N_Res,1);
Q_dem(1:Node.N-Node.N_Res)=Node.Q0_demand(1:Node.N-Node.N_Res,1);
% [A1,A2,B1,P_P,N_h]=make_spars(A1,A2,B1,P_P,N_h)
%% display netwrok
edgeValues=zeros(M_p,1);
nodeValues=zeros(M_h,1);

%%
dH(1:M_f,1)=0;
vec=1:Node.N-Node.N_Res;
A1=sparse(A1);
A2=sparse(A2);
B1=sparse(B1);
denum=(0.5*B1*C')+C_T;
fv = zeros (sum(Pipe.reach),1);
num_samples = 1;
num_nodes=M_h-Node.N_Res;

%%  steady state
options = optimoptions('fsolve', 'MaxFunctionEvaluations', 50000);
X0 = fsolve(@(x) ODEequations(0,x,fv,L,Pipe,Fluid,g,Node,Pump,PRV,TCV,...
    R,M_p,M_f,M_h,A1,A2,h2,Q_dem,Leak,denum,P_P,C_Nor,N_h), X0, options);
q0=X0(1:M_p,1);
h0=X0(M_p+1:M_p+M_h-Node.N_Res,1);
%% Transient simulation
%%%%%%%%%%%%%%%%%  changing set points and solving trasient condition
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);
tspan = 0:0.5:t_total; % time span with fixed time step
deltaQ=+0.5;
margin_percent = 0.01; % 10% margin
node_active1=17;
node_active1=node_active1-1;
Q_dem_down=Q_dem(node_active1);
[t, x1] = ode15s(@(t, x) ODEequations_sim(t, x, fv, L, Pipe, Fluid, g, Node, Pump, PRV, TCV, ...
    R, M_p, M_f, M_h, A1, A2, h2, Q_dem, Leak, denum, P_P, C_Nor, N_h,Q_dem_down,node_active1,deltaQ), tspan, X0, options);
Risk1 = calculateRisk(X0, x1, M_p, margin_percent);
%%
node_active2=18;
node_active2=node_active2-1;
Q_dem_down=Q_dem(node_active2);
[t, x2] = ode15s(@(t, x) ODEequations_sim(t, x, fv, L, Pipe, Fluid, g, Node, Pump, PRV, TCV, ...
    R, M_p, M_f, M_h, A1, A2, h2, Q_dem, Leak, denum, P_P, C_Nor, N_h,Q_dem_down,node_active2,deltaQ), tspan, X0, options);
Risk2 = calculateRisk(X0, x2, M_p, margin_percent);
%%
node_active3=19;
node_active3=node_active3-1;
Q_dem_down=Q_dem(node_active3);
[t, x3] = ode15s(@(t, x) ODEequations_sim(t, x, fv, L, Pipe, Fluid, g, Node, Pump, PRV, TCV, ...
    R, M_p, M_f, M_h, A1, A2, h2, Q_dem, Leak, denum, P_P, C_Nor, N_h,Q_dem_down,node_active3,deltaQ), tspan, X0, options);
Risk3 = calculateRisk(X0, x3, M_p, margin_percent);
%%
node_active4=10;
node_active4=node_active4-1;
Q_dem_down=Q_dem(node_active4);
[t, x4] = ode15s(@(t, x) ODEequations_sim(t, x, fv, L, Pipe, Fluid, g, Node, Pump, PRV, TCV, ...
    R, M_p, M_f, M_h, A1, A2, h2, Q_dem, Leak, denum, P_P, C_Nor, N_h,Q_dem_down,node_active4,deltaQ), tspan, X0, options);
Risk4 = calculateRisk(X0, x4, M_p, margin_percent);
%%
node_active5=489;
node_active5=node_active5-1;
Q_dem_down=Q_dem(node_active5);
[t, x5] = ode15s(@(t, x) ODEequations_sim(t, x, fv, L, Pipe, Fluid, g, Node, Pump, PRV, TCV, ...
    R, M_p, M_f, M_h, A1, A2, h2, Q_dem, Leak, denum, P_P, C_Nor, N_h,Q_dem_down,node_active5,deltaQ), tspan, X0, options);
Risk5 = calculateRisk(X0, x5, M_p, margin_percent);
%%
node_active6=4;
node_active6=node_active6-1;
Q_dem_down=Q_dem(node_active6);
[t, x6] = ode15s(@(t, x) ODEequations_sim(t, x, fv, L, Pipe, Fluid, g, Node, Pump, PRV, TCV, ...
    R, M_p, M_f, M_h, A1, A2, h2, Q_dem, Leak, denum, P_P, C_Nor, N_h,Q_dem_down,node_active6,deltaQ), tspan, X0, options);
Risk6 = calculateRisk(X0, x6, M_p, margin_percent);
%%
x1=x1-X0';
x2=x2-X0';
x3=x3-X0';
x4=x4-X0';
x5=x5-X0';
x6=x6-X0';
% Create figure with adjustable size
fig = figure('Units', 'centimeters', 'Position', [5, 5, 16, 12]); % Adjust as needed

hold on

% Marker incidence
markerInc = floor(length(t)/40);

% Plot each line with markers


% plot(t, x5(:,947), '-v', 'MarkerIndices', 1:markerInc:length(t), 'LineWidth', 1.5, 'DisplayName', 'Node 2');

% plot(t, x4(:,M_p+17-1), '-<', 'MarkerIndices', 1:markerInc:length(t), 'LineWidth', 1.5, 'DisplayName', 'Node 17 at 10');
plot(t, x6(:,M_p+4-1), '->', 'MarkerIndices', 1:markerInc:length(t), 'LineWidth', 1.5, 'DisplayName', 'Node 4');
plot(t, x5(:,M_p+489-1), '-^', 'MarkerIndices', 1:markerInc:length(t), 'LineWidth', 1.5, 'DisplayName', 'Node N_A');
plot(t, x4(:,M_p+10-1), '-h', 'MarkerIndices', 1:markerInc:length(t), 'LineWidth', 1.5, 'DisplayName', 'Node 10');
plot(t, x1(:,M_p+17-1), '--o', 'MarkerIndices', 1:markerInc:length(t), 'LineWidth', 1.5, 'DisplayName', 'Node 17');
plot(t, x2(:,M_p+18-1), '-.s', 'MarkerIndices', 1:markerInc:length(t), 'LineWidth', 1.5, 'DisplayName', 'Node 18');
plot(t, x3(:,M_p+19-1), 'k:d', 'MarkerIndices', 1:markerInc:length(t), 'LineWidth', 1.5, 'DisplayName', 'Node 19');

% Labels and limits
xlabel('Time (s)', 'FontSize', 12, 'FontName', 'Times New Roman')
ylabel('Pressure Head Change (m)', 'FontSize', 12, 'FontName', 'Times New Roman')
xlim([0 150])
ylim([-32 18])
% Grid and legend
grid on
legend('Location', 'best', 'FontSize', 10, 'FontName', 'Times New Roman')

% Set axis font properties
set(gca, 'FontSize', 12, 'FontName', 'Times New Roman')

% Save the figure in high-quality formats
saveas(fig, 'Pressure_Head_Change.pdf')  % High-quality PDF
saveas(fig, 'Pressure_Head_Change.bmp')  % High-quality BMP

hold off

