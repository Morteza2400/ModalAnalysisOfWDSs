clc
clear
close all;
%$$$$$$$$$$$$$$$$$$$$$
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
%% display netwrok
edgeValues=zeros(M_p,1);
nodeValues=zeros(M_h,1);
%%
h_out=0;
dH(1:M_f,1)=0;
vec=1:Node.N-Node.N_Res;
denum=(0.5*B1*C')+C_T;
fv = zeros (sum(Pipe.reach),1);
num_samples = 1;
demand_min = 0.001;
demand_max = 0.001;
num_nodes=M_h-Node.N_Res;

demand = demand_min + (demand_max - demand_min) * rand(num_samples, num_nodes);

demand =demand*0;
% demand(1) =0.01;
Q_dem=demand(1,:)';
%% main body
%%  steady state
options = optimoptions('fsolve', 'MaxFunctionEvaluations', 10000);
X0 = fsolve(@(x) ODEequations(0,x,fv,L,Pipe,Fluid,g,Node,Pump,PRV,TCV,...
    R,M_p,M_f,M_h,A1,A2,h2,Q_dem,Leak,denum,P_P,C_Nor,N_h), X0, options);
q0=X0(1:M_p,1);
h0=X0(M_p+1:M_p+M_h-Node.N_Res,1);
%% Transient simulation
%%%%%%%%%%%%%%%%%  changing set points and solving trasient condition
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);
tspan = 0:t_total; % time span with fixed time step
[t, x] = ode15s(@(t, x) ODEequations(t, x, fv, L, Pipe, Fluid, g, Node, Pump, PRV, TCV, ...
    R, M_p, M_f, M_h, A1, A2, h2, Q_dem, Leak, denum, P_P, C_Nor, N_h), tspan, X0, options);
%% modal analysis
M = zeros (M_p+M_h-Node.N_Res,M_p+M_h-Node.N_Res);
temp = -(RR.*R)'./L.*2.*abs(X0(1:M_p));
for i = 1:M_p
    M(i,i) = temp(i);
end
M (1:M_p,M_p+1:M_p+M_h-Node.N_Res) = diag(L)^(-1)*A1';
M(M_p+1:M_p+M_h-Node.N_Res,1:M_p) = - diag(denum)^(-1)*A1;
eigenvalues = eig(M);

[V,D] = eig(M);
e = diag (D);
eigenvalues_magnitude = abs(e);

W = inv(V)';
P = W.*V;
modes=real(e);

PF = abs(real(P));
for i = 1:size(D,1)
    PF(:,i) = PF(:,i)/max(PF(:,i));
end
LL=Pipe.l./Pipe.reach;
L_critical=max(Pipe.l./Pipe.reach);
F_c=2*pi*Pipe.a(1)/10/(L_critical);

% Filter
valid_modes_indices = find(abs(imag(e)) < F_c);
valid_modes = e(valid_modes_indices);
valid_PF = PF(:, valid_modes_indices);

stateNames = 1:M_p+M_h-Node.N_Res;

numModes = size(PF, 2); 
modeLabels = arrayfun(@num2str, 1:numModes, 'UniformOutput', false); 

figure;

plot(real(e), imag(e), 'x', 'MarkerSize', 5, 'LineWidth', 2 , 'DisplayName', 'EWCM');
hold on;
xLimits=xlim;
line(xLimits, [F_c F_c], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 2, 'DisplayName', 'Critical Frequency');
legend('Location','northwest')
ylim([-1 inf])
xlabel('Decay Rate');
ylabel('Frequency rad/s');
title('polzes Plot H3');
% Create the heatmap 
figure;
heatmap(modeLabels, stateNames, PF, 'Colormap', jet, 'Title', 'Participation Factor Heatmap');
xlabel('Modes');
ylabel('States');

%  the heatmap for valid modes
valid_modeLabels = arrayfun(@num2str, valid_modes_indices, 'UniformOutput', false);

figure;
heatmap(valid_modeLabels, stateNames, valid_PF, 'Colormap', jet, 'Title', 'Participation Factor Heatmap for Valid Modes');
xlabel('Valid Modes');
ylabel('States');


% Initialize damping ratios array
damping_ratios = zeros(size(e));

% Calculate damping ratio for each eigenvalue
for i = 1:length(e)
    real_part(i) = real(e(i));
    imag_part(i) = imag(e(i));
    damping_ratios(i) = -real_part(i) / sqrt(real_part(i)^2 + imag_part(i)^2)*100;
end

k=0;
for i=1:length(valid_modes_indices)
ModeNumber=valid_modes_indices(i);
if imag_part(ModeNumber)>0 
    k=k+1;
plot_m_number(k)=ModeNumber;
end
if imag_part(ModeNumber)==0
    k=k+1;
plot_m_number(k)=ModeNumber;
end
end

for i=1:length(plot_m_number)
ModeNumber=plot_m_number(i);
plot_mode=PF(:,ModeNumber);
edgeValues=plot_mode(1:M_p);
node_Values=plot_mode(M_p+1:end);
nodeValues=zeros(M_h,1);
nodeValues(Node.N_Res+1:end)=node_Values;
plotGraphWithValues(Node, A_save, nodeValues, edgeValues,ModeNumber,real_part(ModeNumber),imag_part(ModeNumber));
end
