clear; clc; close all;
rng('shuffle');

%% ================= IEEE figure config =================
figW = 6.5; figH = 4.8; margin = 0.05;
paperW = figW + 2*margin; paperH = figH + 2*margin;

set(0,'DefaultAxesFontSize',11);
set(0,'DefaultTextFontSize',11);
set(0,'DefaultLineLineWidth',1.5);
set(0,'DefaultAxesLineWidth',0.9);
set(0,'DefaultFigureUnits','inches');
set(0,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontName','Times New Roman');
set(0,'DefaultAxesTickLabelInterpreter','latex');
set(0,'DefaultLegendInterpreter','latex');
set(0,'DefaultFigureColor','white');

%% ================= Common params =================
K_reuse    = 7;
num_trials = 200000;

% 6 first-tier interferers at angles 30, 90, ..., 330 deg (as in many classes)
theta_inter = (30:60:330) * pi/180;   % length=6

%% ============================================================
% FIGURE 1 — C/I CDF (R = 3000 m)
%% ============================================================
R_cell = 3000;
D = sqrt(3*K_reuse) * R_cell;

BS_inter_x = D*cos(theta_inter);
BS_inter_y = D*sin(theta_inter);

n_vals     = [2 4];
sigma_vals = [0 6];      % dB

fig1 = figure; hold on; grid on; box on;

color_map = {[0 0.447 0.741], [0.85 0.325 0.098]};  % n=2 blue, n=4 orange
line_map  = {'-', '--'};                            % sigma=0 solid, sigma=6 dashed

for ni = 1:length(n_vals)
    n = n_vals(ni);

    for si = 1:length(sigma_vals)
        sigma = sigma_vals(si);

        CI_lin = zeros(1, num_trials);

        for t = 1:num_trials
            % MS uniformly in disk radius R_cell
            r_ms = R_cell*sqrt(rand);
            phi  = 2*pi*rand;
            ms_x = r_ms*cos(phi);
            ms_y = r_ms*sin(phi);

            % Desired link shadowing (dB) -> linear gain
            Xc_dB = sigma*randn;
            Gc    = 10^(Xc_dB/10);

            % Carrier power (relative, Pt cancels in C/I)
            C = (r_ms^(-n)) * Gc;

            % Interference
            I_sum = 0;
            for j = 1:6
                d_ij = sqrt((ms_x-BS_inter_x(j))^2 + (ms_y-BS_inter_y(j))^2);
                Xi_dB = sigma*randn;
                Gi    = 10^(Xi_dB/10);
                I_sum = I_sum + (d_ij^(-n)) * Gi;
            end

            CI_lin(t) = C / I_sum;
        end

        CI_dB = 10*log10(CI_lin);
        CI_dB = sort(CI_dB);
        prob  = (1:num_trials)/num_trials;

        plot(CI_dB, prob, 'Color', color_map{ni}, 'LineStyle', line_map{si});
    end
end

set(gca,'YScale','log'); ylim([1e-4 1]);
xlabel('$C/I$ (dB)','Interpreter','latex');
ylabel('CDF','Interpreter','latex');
legend({'$n=2,\sigma=0$','$n=2,\sigma=6$',...
        '$n=4,\sigma=0$','$n=4,\sigma=6$'},...
        'Location','best');

export_ieee_pdf(fig1,'IEEE_CI_CDF',figW,figH,paperW,paperH);

%% ============================================================
% FIGURE 2 — C/(I+N) CDF (n=4, sigma=6 dB, R = 300/3000/30000 m)
%% ============================================================
n_b     = 4;
sigma_b = 6;         % dB
radii   = [300 3000 30000];

% ---- Noise power: -120 dBm -> W ----
N_dBm = -120;
N_W   = 10^((N_dBm - 30)/10);

% ---- Calibration: Pr at 100 m is -15 dBm (nominal) ----
d_ref       = 100;        % m
Pr_ref_dBm  = -15;        % dBm
Pr_ref_W    = 10^((Pr_ref_dBm - 30)/10);   % W

fig2 = figure; hold on; grid on; box on;
radius_colors = {[0 0.447 0.741], [0.85 0.325 0.098], [0.466 0.674 0.188]};

for rr = 1:length(radii)
    R_b = radii(rr);

    % Reuse distance for this cell radius
    D_b = sqrt(3*K_reuse) * R_b;
    BS_ix = D_b*cos(theta_inter);
    BS_iy = D_b*sin(theta_inter);

    CINR_lin = zeros(1, num_trials);

    for t = 1:num_trials
        % MS uniformly in disk radius R_b
        r_ms = R_b*sqrt(rand);
        phi  = 2*pi*rand;
        ms_x = r_ms*cos(phi);
        ms_y = r_ms*sin(phi);

        % Desired received power using reference at 100 m:
        % Pr(d) = Pr(100) * (d/100)^(-n) * 10^(X/10)
        Xc_dB = sigma_b*randn;
        Gc    = 10^(Xc_dB/10);
        C_W   = Pr_ref_W * (r_ms/d_ref)^(-n_b) * Gc;

        % Interference sum (each interferer same nominal power at 100m)
        I_W = 0;
        for j = 1:6
            d_ij = sqrt((ms_x-BS_ix(j))^2 + (ms_y-BS_iy(j))^2);
            Xi_dB = sigma_b*randn;
            Gi    = 10^(Xi_dB/10);
            I_W   = I_W + Pr_ref_W * (d_ij/d_ref)^(-n_b) * Gi;
        end

        CINR_lin(t) = C_W / (I_W + N_W);
    end

    CINR_dB = 10*log10(CINR_lin);
    CINR_dB = sort(CINR_dB);
    prob    = (1:num_trials)/num_trials;

    plot(CINR_dB, prob, 'Color', radius_colors{rr});
end

set(gca,'YScale','log'); ylim([1e-4 1]);
xlabel('$C/(I+N)$ (dB)','Interpreter','latex');
ylabel('CDF','Interpreter','latex');
legend({'$R=300$m','$R=3000$m','$R=30000$m'}, 'Location','best');

export_ieee_pdf(fig2,'IEEE_CINR_CDF',figW,figH,paperW,paperH);

disp('IEEE PDFs saved successfully: IEEE_CI_CDF.pdf and IEEE_CINR_CDF.pdf');

%% ========================= Export function =========================
function export_ieee_pdf(figHandle,fileName,figW,figH,paperW,paperH)
    set(figHandle,'PaperUnits','inches');
    set(figHandle,'PaperSize',[paperW paperH]);
    set(figHandle,'PaperPosition',[(paperW-figW)/2, (paperH-figH)/2, figW, figH]);
    print(figHandle, fileName, '-dpdf', '-painters', '-r300');
end
