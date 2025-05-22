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
freq=[3.5271  3.71];
for kk=1:2  % run for two frequncy

num_samples = 1;
demand_min = 0.001;
demand_max = 0.001;
num_nodes=M_h-Node.N_Res;
demand = demand_min + (demand_max - demand_min) * rand(num_samples, num_nodes);
demand =demand*0;
input_dem=21;
demand(input_dem) =0.01;
Q_dem=demand';
%% steady state
options = optimoptions('fsolve', 'MaxFunctionEvaluations', 10000);
X0 = fsolve(@(x) ODEequations(0,x,fv,L,Pipe,Fluid,g,Node,Pump,PRV,TCV,...
    R,M_p,M_f,M_h,A1,A2,h2,Q_dem,Leak,denum,P_P,C_Nor,N_h,d_dem), X0, options);
q0=X0(1:M_p,1);
h0=X0(M_p+1:M_p+M_h-Node.N_Res,1);
%% Transient simulation
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
tspan = 0:0.01:300; % time span with fixed time step
[t, x] = ode15s(@(t, x) ODEequations_sim(t, x, fv, L, Pipe, Fluid, g, Node, ...
    R, M_p, M_h, A1, A2, h2, Q_dem, Leak, denum, C_Nor,freq(kk),d_dem,input_dem), tspan, X0, options);

% Example plot for your data
figure;
fig = figure('Color', 'w', 'Units', 'inches', 'Position', [1, 1, 10, 6]);
if kk==1
t_s = 27000+24+54;
t_e = 27000 + 338/2+24+54;
else
t_s = 27769+0;
t_e = 27769 + 338/2+0;
end
% Define the list of specific nodes to plot (e.g., Node 3, Node 5, Node 10)
nodes_to_plot = [3, 4, 5, 6, 8, 12, 16, 21, 23, 25, 27];  % You can modify this list with any nodes you want
% nodes_to_plot =nodes_to_plot -2;
% Define markers and line styles (repeated if needed)
markers = {'o', '*', 'd', '^', 'v', '<', '>', 'p', 'h', 'x', '+', 's', 'd', 'P', 'H'};
lines = {'-', '--', '-.', ':', '-', '--', '-.', ':', '-', '--', '-.', ':', '-', '--', '-.'};

% Repeat markers and lines to cover 25 nodes
markers = repmat(markers, 1, ceil(25 / length(markers))); % Repeat to get enough markers
lines = repmat(lines, 1, ceil(25 / length(lines))); % Repeat to get enough lines

% Labels for the legend (based on the selected nodes)
labels = arrayfun(@(i) sprintf('Node %d', i), nodes_to_plot, 'UniformOutput', false);
set(gca, 'FontName', 'Times New Roman');
% Marker indices for controlling the number of markers (e.g., every 5th point)
num_points = length(t(t_s:t_e));
marker_spacing = 15;  % Adjust this value to regulate marker frequency
marker_indices = 1:marker_spacing:num_points;

% Plot the data
hold on;
for i = 1:length(nodes_to_plot)  % Loop over the selected nodes
    node_index = nodes_to_plot(i)-2;  % Get the current node from the list
    plot(t(t_s:t_e)-t(t_s), x(t_s:t_e, node_index+27)-X0(node_index+27), ...
        'LineStyle', lines{i},'LineWidth',1, 'Marker', markers{i}, 'MarkerIndices', marker_indices, ...
        'DisplayName', labels{i});
end
hold off;

% Set axis labels and title
xlabel('Time (s)');
ylabel('Head Difference (m)');
% title('Head change at Selected Nodes Over Time');
xlim([0 1.7])
% Add legend and grid
legend('show', 'Location', 'eastoutside', 'NumColumns', 1);
grid on;

% Improve layout
set(gca, 'FontSize', 20);  % Adjust font size for readability
filename = sprintf('head_change%d.png', kk);
exportgraphics(fig, filename, 'Resolution', 600, 'BackgroundColor', 'white');
end