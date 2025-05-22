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
t = 0:0.1:20000; % time span with fixed time step
node_active=223;
node_active=node_active-1;
Q_dem_down=Q_dem(node_active);
deltaQ=+0.5;
%% input signal

ts = 10;  % Start time of the rise
td = 2;   % Duration of the rise
Q_dem_up=Q_dem_down+deltaQ;
Q_dem_input=zeros(length(t),1);
hold_time = 2900;  % Time to hold the peak

for i = 1:length(t)
    if t(i) <= ts
        Q_dem_input(i) = Q_dem_down;  % Before the rise starts
    elseif t(i) > ts && t(i) <= ts + td
        % Smooth transition up
        Q_dem_input(i) = (Q_dem_up - Q_dem_down) * (sin(pi * (t(i) - ts) / td - pi / 2) + 1) / 2 + Q_dem_down;
    elseif t(i) > ts + td && t(i) <= ts + td + hold_time
        Q_dem_input(i) = Q_dem_up;  % Hold the peak
    elseif t(i) > ts + td + hold_time && t(i) <= ts + 2*td + hold_time
        % Smooth transition back down
        Q_dem_input(i) = (Q_dem_up - Q_dem_down) * (sin(pi * (t(i) - (ts + td + hold_time)) / td + pi / 2) + 1) / 2 + Q_dem_down;
    else
        Q_dem_input(i) = Q_dem_down;  % Return to original value
    end
end

% Q_dem_input = Q_dem_input - mean(Q_dem_input);


% Q_dem_input=x(:,M_p+19);
figure
plot(t,Q_dem_input,"-r")
xlim([0 30]);
% Sampling frequency (assumes uniform sampling)
Fs = 1 / (t(2) - t(1));  % Compute the sampling frequency from time step

% Compute FFT
N = length(Q_dem_input);  % Number of points
Y = fft(Q_dem_input);  % Compute FFT
P2 = abs(Y/N);  % Normalize
P1 = P2(1:N/2+1);  % Take single-sided spectrum
P1(2:end-1) = 2*P1(2:end-1);  % Double the amplitude for non-DC components
fprintf('Sampling Frequency: %.3f Hz\n', Fs);
fprintf('FFT Frequency Resolution: %.5f Hz\n', Fs/N);
fprintf('Nyquist Frequency: %.3f Hz\n', Fs/2);
% Frequency vector
f = Fs * (0:(N/2)) / N;

threshold= max(P1) / sqrt(2);  % -3 dB threshold (70.7% of max amplitude)
valid_frequencies = f(P1 >= threshold);  % Find frequencies above threshold

% Compute bandwidth
if ~isempty(valid_frequencies)
    bandwidth = valid_frequencies(end) - valid_frequencies(1);
else
    bandwidth = 0;  % No significant frequencies found
end
[max_amp, max_idx] = max(P1);
dominant_freq = f(max_idx);
fprintf('Dominant Frequency: %.3f Hz\n', dominant_freq);
% Display bandwidth
fprintf('The 3 dB bandwidth of the signal is %.3f Hz\n', bandwidth);
Fs = 1 / (t(2) - t(1));  % Compute the sampling frequency from time step
omega_s = 2 * pi * Fs;   % Convert to rad/s

% Compute FFT
N = length(Q_dem_input);  % Number of points
Y = fft(Q_dem_input);  % Compute FFT
P2 = abs(Y/N);  % Normalize
P1 = P2(1:N/2+1);  % Take single-sided spectrum
P1(2:end-1) = 2*P1(2:end-1);  % Double the amplitude for non-DC components

fprintf('Sampling Frequency: %.3f rad/s\n', omega_s);
fprintf('FFT Frequency Resolution: %.5f rad/s\n', omega_s/N);
fprintf('Nyquist Frequency: %.3f rad/s\n', omega_s/2);

% Frequency vector in rad/s
omega = omega_s * (0:(N/2)) / N;

threshold = max(P1) / sqrt(2);  % -3 dB threshold (70.7% of max amplitude)
valid_frequencies = omega(P1 >= threshold);  % Find frequencies above threshold

% Compute bandwidth
if ~isempty(valid_frequencies)
    bandwidth = valid_frequencies(end) - valid_frequencies(1);
else
    bandwidth = 0;  % No significant frequencies found
end

[max_amp, max_idx] = max(P1);
dominant_freq = omega(max_idx);
fprintf('Dominant Frequency: %.3f rad/s\n', dominant_freq);
fprintf('The 3 dB bandwidth of the signal is %.3f rad/s\n', bandwidth);

% Plot amplitude spectrum up to 10 rad/s
figure;
plot(omega, P1, 'b', 'LineWidth', 1.5);
xlim([0 0.1]); % Limit frequency range to 10 rad/s
xlabel('Frequency (rad/s)');
ylabel('Amplitude');
title('Amplitude Spectrum of Q_{dem\_input}');
xlim([0.0005 1]) % use a small positive number instead of 0
set(gca, 'XScale', 'log')
grid on
%% modal analysis


M = zeros (M_p+M_h-Node.N_Res,M_p+M_h-Node.N_Res);
temp = -RR'./L.*2.*abs(X0(1:M_p))';
for i = 1:M_p
    M(i,i) = temp(i);
end
M (1:M_p,M_p+1:M_p+M_h-Node.N_Res) = diag(L)^(-1)*A1';
M(M_p+1:M_p+M_h-Node.N_Res,1:M_p) = - diag(denum)^(-1)*A1;
eigenvalues = eig(M);
j=1;
[V,D] = eig(M);
e = diag (D);
eigenvalues_magnitude = abs(e);
[strongest_mode_magnitude(j), strongest_mode_index(j)] = max(eigenvalues_magnitude);
W = inv(V)';
P = W.*V;
modes=real(e);

PF = abs(real(P));
LL=Pipe.l./Pipe.reach;
L_critical=max(Pipe.l./Pipe.reach);
F_c=2*pi*Pipe.a(1)/10/(L_critical);

% Filter
valid_modes_indices = find((imag(e)) < F_c*1 & (imag(e)) > 0);
% valid_modes_indices = find(abs(imag(e)) < F_c*1.0);
valid_modes = e(valid_modes_indices);
valid_PF = PF(:, valid_modes_indices);

numModes = size(PF, 2);
modeLabels = arrayfun(@num2str, 1:numModes, 'UniformOutput', false);

figure;

plot(real(e), imag(e), 'x', 'MarkerSize', 5, 'LineWidth', 2 , 'DisplayName', 'EWCM');
hold on;
xLimits=xlim;
line(xLimits, [F_c F_c], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 2, 'DisplayName', 'Critical Frequency');
legend('Location','northwest')
ylim([0 F_c*2])
xlabel('Decay Rate');
ylabel('Frequency rad/s');
title('polzes Plot H3');


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
end

for i=1:length(plot_m_number)
    ModeNumber=plot_m_number(i);
    plot_mode=PF(:,ModeNumber);
    edgeValues=plot_mode(1:M_p);
    node_Values=plot_mode(M_p+1:end);
    nodeValues=zeros(M_h,1);
    nodeValues(Node.N_Res+1:end)=node_Values;
    % subplot(length(plot_m_number), 1, i);
    % plotGraphWithValues_un_leg(Node, A_save, nodeValues, edgeValues,ModeNumber,real_part(ModeNumber),imag_part(ModeNumber));
end

weights = interp1(f, P1, abs(imag(valid_modes)), 'linear', 0); % Interpolate amplitudes as weights
weights = weights / sum(weights);
max_pf = valid_PF * weights;
max_pf_norm = (max_pf - min(max_pf)) / (max(max_pf) - min(max_pf));
max_pf = max_pf_norm;

% Separate edge and node values
edgeValues = max_pf(1:M_p);
node_Values = max_pf(M_p+1:end);
nodeValues = zeros(M_h,1);
nodeValues(Node.N_Res+1:end) = node_Values;
plotGraphWithValues_un_leg2(Node, A_save, nodeValues, edgeValues, 0, 0, 0);
h = colorbar;
h.Label.String = 'PF Magnitude';
h.Label.FontSize = 16; % Adjust font size
h.Label.FontName = 'Times New Roman';
h.TickLabelInterpreter = 'latex'; % Optional: Use LaTeX formatting
h.FontSize = 14; % Adjust tick label size
h.FontName = 'Times New Roman'; % Ensure consistent font
saveas(gcf, 'Av_PF_plot.jpg'); % Save as JPG
saveas(gcf, 'Av_PF_plot.pdf'); % Save as PDF
print('Av_PF_plot', '-dpng', '-r300'); % Save high-resolution PNG (300 DPI)

