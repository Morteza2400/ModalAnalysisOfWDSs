clc
clear
close all;
%$$$$$$$$$$$$$$$$$$$$$
%% INPUTS
load 'inputs.mat'
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
%%
h_out=0;
dH(1:M_f,1)=0;
vec=1:Node.N-Node.N_Res;
denum=(0.5*B1*C')+C_T;
fv = zeros (sum(Pipe.reach),1);

%% steady state
options = optimoptions('fsolve', 'MaxFunctionEvaluations', 10000);
X0 = fsolve(@(x) ODEequations(0,x,fv,L,Pipe,Fluid,g,Node,Pump,PRV,TCV,...
    R,M_p,M_f,M_h,A1,A2,h2,Q_dem,Leak,denum,P_P,C_Nor,N_h,d_dem), X0, options);
q0=X0(1:M_p,1);
h0=X0(M_p+1:M_p+M_h-Node.N_Res,1);
%% modal analysis

M = zeros (M_p+M_h-Node.N_Res,M_p+M_h-Node.N_Res);
temp = -RR'./L.*2.*abs(X0(1:M_p))';
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
P_re=real(P);
PF = abs(real(P));
for i = 1:size(D,1)
    PF(:,i) = PF(:,i)/max(PF(:,i));
end
LL=Pipe.l./Pipe.reach;
L_critical=max(Pipe.l./Pipe.reach);
F_c=2*pi*Pipe.a(1)/10/(L_critical);

% Filter
valid_modes_indices = find(imag(e) < F_c & (imag(e)) >= 0);
valid_modes = e(valid_modes_indices);
valid_PF = PF(:, valid_modes_indices);

figureHandle = figure;
set(figureHandle, 'Position', [100, 50, 800, 600]);

numModes = size(PF, 2); 
modeLabels = arrayfun(@num2str, 1:numModes, 'UniformOutput', false); 
plot(real(e), imag(e), 'x', 'MarkerSize', 8,'MarkerFaceColor', 'b', 'LineWidth', 2 , 'DisplayName', 'EWCM');
hold on;
xLimits=xlim;
line(xLimits, [F_c F_c], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 2, 'DisplayName', 'Critical Frequency');
legend('Location','northwest')
xlabel('Decay Rate');
ylabel('Frequency rad/s');
title('Location of Poles');
hold off


% Create the heatmap 
figure;
heatmap(modeLabels, 1:M_p+M_h-Node.N_Res, PF, 'Colormap', jet, 'Title', 'Participation Factor Heatmap');
xlabel('Modes');
ylabel('States');
%%
%  the heatmap for valid modes
valid_modeLabels = arrayfun(@num2str, 1:length(valid_modes_indices), 'UniformOutput', false);

valid_PF = round(valid_PF, 3);
valid_PF(valid_PF < 1e-10) = 0;

% Create the figure
figure('Color', 'w', 'Position', [100, 100, 800, 800]); % Set white background and size

% Create the heatmap
h = heatmap(valid_modeLabels, 1:length(valid_PF), valid_PF, ...
    'Colormap', jet, ...
    'Title', 'PF Heatmap for Valid Modes');

xlabel('Mode number');
ylabel('States');

% Adjust the figure for high-quality saving
set(gca, 'FontSize', 14); % Adjust font size for better readability
set(gca, 'FontName', 'Times New Roman'); % Use a readable font

% Save the figure in high quality
saveas(gcf, 'Participation_Factor_Heatmap.pdf'); % Save as PDF
% saveas(gcf, 'Participation_Factor_Heatmap.jpg'); % Save as JPG

% If you want even higher quality for the JPG, use:
print(gcf, 'Participation_Factor_Heatmap', '-djpeg', '-r600'); % Save JPG with 600 dpi resolution


%%
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


% Define adjustable font size
fontSize = 14;  % Modify this value for text size

% Define subplot height based on number of plots
numPlots = length(plot_m_number) - 1;  
subplotHeight = 0.20;  
totalHeight = numPlots * subplotHeight + 0.1; 

% Create figure with dynamic height
figure('Units', 'inches', 'Position', [1, 1, 10, totalHeight * 8]);

for i = 2:length(plot_m_number)
    ModeNumber = plot_m_number(i);
    plot_mode = PF(:, ModeNumber);
    edgeValues = plot_mode(1:M_p);
    node_Values = plot_mode(M_p+1:end);
    nodeValues = zeros(M_h, 1);
    nodeValues(1:M_h-Node.N_Res) = node_Values;
    nodeValues = N_h_total * nodeValues;
    
    subplot(numPlots, 1, i-1);  

    plotGraphWithValues(Node, A_save, nodeValues, edgeValues, ModeNumber, real_part(ModeNumber), imag_part(ModeNumber),i);
    colorbar('off');  % Remove individual colorbars
    
    % Ensure proper box and spacing
    set(gca, 'XColor', 'k', 'YColor', 'k', 'LineWidth', 1);
    ylim([min(nodeValues) - 50, max(nodeValues) + 50]); % Adjust padding
end

h = colorbar;
h.Position = [0.92, 0.15, 0.02,  1.1*totalHeight]; 
h.Label.String = 'PF Magnitude';
h.Label.FontSize = 16; 
h.Label.FontName = 'Times New Roman';
h.TickLabelInterpreter = 'latex'; 
h.FontSize = 14; 
h.FontName = 'Times New Roman'; 
saveas(gcf, 'PF_plot.jpg'); % Save as JPG
saveas(gcf, 'PF_plot.pdf'); % Save as PDF
print('PF_plot', '-dpng', '-r300'); % Save high-resolution PNG (300 DPI)

