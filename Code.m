% ============================================================
%   IEEE Double-Column Style CDF Plot (Pure Line Version)
% ============================================================
clear; clc; close all;
rng('shuffle');

%% ================= IEEE 图配置 =================
figW = 6.5;
figH = 4.8;
margin = 0.05;
paperW = figW + 2*margin;
paperH = figH + 2*margin;

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

%% ================= 通用参数 =================
R_cell = 3000;
K = 7;
D = sqrt(3*K) * R_cell;
num_trials = 200000;
theta_inter = (30:60:330) * pi/180;

BS_inter_x = D*cos(theta_inter);
BS_inter_y = D*sin(theta_inter);

%% ============================================================
% FIGURE 1 — C/I CDF
%% ============================================================
n_vals = [2 4];
sigma_vals = [0 6];

fig1 = figure; hold on; grid on; box on;

% n 不同颜色
color_map = {[0 0.447 0.741], [0.85 0.325 0.098]}; % 蓝 / 橙

% sigma 不同线型
line_map = {'-', '--'}; % 实线 / 虚线

for ni=1:length(n_vals)
    n = n_vals(ni);

    for si=1:length(sigma_vals)
        sigma = sigma_vals(si);

        CI_dB = zeros(1,num_trials);

        for i=1:num_trials
            r_ms = R_cell*sqrt(rand);
            phi = 2*pi*rand;
            ms_x = r_ms*cos(phi);
            ms_y = r_ms*sin(phi);

            gain_C = (r_ms^-n)*10^(sigma*randn/10);

            I_total=0;
            for j=1:6
                d_ij = sqrt((ms_x-BS_inter_x(j))^2 + ...
                            (ms_y-BS_inter_y(j))^2);
                I_total = I_total + ...
                    (d_ij^-n)*10^(sigma*randn/10);
            end

            CI_dB(i)=10*log10(gain_C/I_total);
        end

        sorted_CI = sort(CI_dB);
        prob = (1:num_trials)/num_trials;

        plot(sorted_CI,prob,...
            'Color',color_map{ni},...
            'LineStyle',line_map{si});
    end
end

set(gca,'YScale','log');
ylim([1e-4 1]);

xlabel('$C/I$ (dB)','Interpreter','latex');
ylabel('CDF','Interpreter','latex');

legend({'$n=2,\sigma=0$','$n=2,\sigma=6$',...
        '$n=4,\sigma=0$','$n=4,\sigma=6$'},...
        'Location','best');

export_ieee_pdf(fig1,'IEEE_CI_CDF',figW,figH,paperW,paperH);

%% ============================================================
% FIGURE 2 — C/(I+N)
%% ============================================================
n_b=4; sigma_b=6;
radii=[300 3000 30000];

N_dBm=-120;
N_linear=10^(N_dBm/10);

d_ref=100;
Pr_ref=10^(-15/10);
Pt=Pr_ref/(d_ref^-n_b);

fig2=figure; hold on; grid on; box on;

% 三种半径三种颜色
radius_colors = {[0 0.447 0.741],...
                 [0.85 0.325 0.098],...
                 [0.466 0.674 0.188]};

for rr=1:length(radii)
    R_b=radii(rr);
    CINR_dB=zeros(1,num_trials);

    D_b=sqrt(3*K)*R_b;
    BS_ix=D_b*cos(theta_inter);
    BS_iy=D_b*sin(theta_inter);

    for i=1:num_trials
        r_ms=R_b*sqrt(rand);
        phi=2*pi*rand;
        ms_x=r_ms*cos(phi);
        ms_y=r_ms*sin(phi);

        C=Pt*(r_ms^-n_b)*10^(sigma_b*randn/10);

        I_sum=0;
        for j=1:6
            d_ij=sqrt((ms_x-BS_ix(j))^2 + ...
                      (ms_y-BS_iy(j))^2);
            I_sum=I_sum+...
                Pt*(d_ij^-n_b)*10^(sigma_b*randn/10);
        end

        CINR_dB(i)=10*log10(C/(I_sum+N_linear));
    end

    sorted_CINR=sort(CINR_dB);
    prob=(1:num_trials)/num_trials;

    plot(sorted_CINR,prob,...
        'Color',radius_colors{rr});
end

set(gca,'YScale','log');
ylim([1e-4 1]);

xlabel('$C/(I+N)$ (dB)','Interpreter','latex');
ylabel('CDF','Interpreter','latex');

legend({'$R=300$m','$R=3000$m','$R=30000$m'},...
       'Location','best');

export_ieee_pdf(fig2,'IEEE_CINR_CDF',figW,figH,paperW,paperH);

disp('IEEE pure-line PDFs saved successfully.');

%% ========================= 导出函数 ==============================
function export_ieee_pdf(figHandle,fileName,figW,figH,paperW,paperH)
    set(figHandle,'PaperUnits','inches');
    set(figHandle,'PaperSize',[paperW paperH]);
    set(figHandle,'PaperPosition',[(paperW-figW)/2,...
        (paperH-figH)/2,figW,figH]);
    print(figHandle,fileName,'-dpdf','-painters','-r300');
end
