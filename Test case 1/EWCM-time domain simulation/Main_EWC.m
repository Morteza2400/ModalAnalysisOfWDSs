clc
clear
close all;
%$$$$$$$$$$$$$$$$$$$$$
load 'input.mat'
[t_total,deltaT,Fluid,Node,Pipe,g,PRV,TCV,Pump]=read_data(data);
% Pre-process of data
[Pipe,Node,Fluid,A,A1,A2,Res_row,B,B1,q0,h0,h1,...
    h2,q,L,R,C,M_p,M_h,P_P,N_h,C_PRV,C_pump,C_TCV,C_Nor,RR,A_save,N_h_total]=PreProcess(Pipe,Fluid,Node,g);
%% initial condition
M_f=M_p-Pipe.ctr;
X0=zeros(M_p+M_h-Node.N_Res+Pump.num+PRV.num,1);
X0(1:M_p,1)=q0(1:M_p,1);
X0(M_p+1:M_p+M_h-Node.N_Res)=h1;

Leak=zeros(M_h-Node.N_Res,1);
C_T=zeros(M_h-Node.N_Res,1);
Q_dem=zeros(M_h-Node.N_Res,1);
d_dem=zeros(M_h-Node.N_Res,1);
Q_dem(1:Node.N-Node.N_Res)=Node.Q0_demand(1:Node.N-Node.N_Res,1);
d_dem(1:Node.N-Node.N_Res)=Node.d(1:Node.N-Node.N_Res,1);

%% display netwrok
edgeValues=zeros(M_p,1);
nodeValues=zeros(M_h,1);
% plotGraphWithValues(Node, A_save, nodeValues, edgeValues,0,0,0);
%%
h_out=0;
dH(1:M_f,1)=0;
vec=1:Node.N-Node.N_Res;
denum=(0.5*B1*C')+C_T;
fv = zeros (sum(Pipe.reach),1);
input_dem=[4,8]; % excitation points node 4 and 8. if we add reservior nodes they become 6 and 10
for kk=1:2
demand(1:M_h-Node.N_Res) =0;
demand(input_dem(kk)) =0.01;
% demand=Q_dem;
Q_dem=demand';
%% steady state
options = optimoptions('fsolve', 'MaxFunctionEvaluations', 10000);
X0 = fsolve(@(x) ODEequations(0,x,fv,L,Pipe,Fluid,g,Node,Pump,PRV,TCV,...
    R,M_p,M_f,M_h,A1,A2,h2,Q_dem,Leak,denum,P_P,C_Nor,N_h,d_dem), X0, options);
q0=X0(1:M_p,1);
h0=X0(M_p+1:M_p+M_h-Node.N_Res,1);
%% Transient simulation
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
tspan = 0:0.01:200; % time span with fixed time step
[t, x] = ode15s(@(t, x) ODEequations_sim(t, x, fv, L, Pipe, Fluid, g, Node, Pump, PRV, TCV, ...
    R, M_p, M_f, M_h, A1, A2, h2, Q_dem, Leak, denum, P_P, C_Nor, N_h,d_dem,input_dem(kk)), tspan, X0, options);
%% plot
fontSize = 14;  % Modify this value for text size
% Create figure with specific size
figure('Units', 'inches', 'Position', [1, 1, 10, 6]); % Adjust figure width and height
%%node 10
t_s=17085;
t_e=t_s+200;
% Define markers and line styles
markers = {'o', 's', 'd', '^', 'v', '<', '>', 'p', 'h', 'x', '+', '*', 'd', 'P', 'H'};
lines = {'-', '--', '-.', ':', '-', '--', '-.', ':', '-', '--', '-.', ':', '-', '--', '-.'};
% Labels for the legend
labels = {'Node 3', 'Node 4', 'Node 5', 'Node 6', 'Node 7', ...
          'Node 8', 'Node 9', 'Node 10', 'Node 11', 'Node 12', ...
          'Node 13', 'Node 14', 'Node 15', 'Node 16', 'Node 17'};

% Marker indices for controlling the number of markers (e.g., every 5th point)
num_points = length(t(t_s:t_e));
marker_spacing = 15;  % Adjust this value to regulate marker frequency
marker_indices = 1:marker_spacing:num_points;

% Plot the data
hold on;
for i = 1:15
    plot(t(t_s:t_e)-t(t_s), x(t_s:t_e, i+16)-X0(i+16), ...
        'LineStyle', lines{i},'LineWidth',1, 'Marker', markers{i}, 'MarkerIndices', marker_indices, ...
        'DisplayName', labels{i});
end
hold off;

% Set axis labels and title with Times New Roman font
xlabel('Time (s)', 'FontSize', fontSize, 'FontName', 'Times New Roman');
ylabel('Head Difference', 'FontSize', fontSize, 'FontName', 'Times New Roman');
title('Head at Different Nodes Over Time', 'FontSize', fontSize + 2, 'FontName', 'Times New Roman');

% Add legend with width equal to figure width
legendHandle = legend('show', 'Location', 'eastoutside', 'NumColumns', 1);
legendHandle.FontSize = fontSize;
legendHandle.FontName = 'Times New Roman';
% legendHandle.Position = [0.1, 0.02, 0.8, 0.05]; % Adjust legend width and placement
% Improve layout
grid on;
set(gca, 'FontSize', fontSize, 'FontName', 'Times New Roman');  % Adjust font size for readability
end
