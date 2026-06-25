%% Figure 3 Code
% Generates Figure 3: model mean isotope trajectories for selected root
% orders and mean root longevity across orders.
%
% Required File(s): 'face_dataset.csv'
%
% Required toolboxes: None
%
clear; clc;

%% Load Data Into Table and Arrays
dtab = readtable('face_dataset.csv'); % read data table
tshft = 18; % time offset in days
t0 = min(dtab.DateSampled) - days(tshft);

dat.t = days(dtab.DateSampled - t0)/365;
dat.Order = dtab.Order; % list of orders 
dat.d13C = dtab.d13C; % list of d13C values

dat.Cind0 = strcmp(dtab.CO2Treatment,'Ambient'); % indices CO2 treatment, ambient
dat.Cind1 = strcmp(dtab.CO2Treatment,'Elevated'); % indices CO2 treatment, elevated

dat.Nind0 = strcmp(dtab.NTreatment,'Unfertilized'); % indices for unfertilized
dat.Nind1 = strcmp(dtab.NTreatment,'Fertilized'); % indices for fertilized

%% Model parameters from Table 1
% q = 1 unfertilized, q = 2 fertilized
% model fields: unc, con, hyp
xeq = -28.85; % equilibrium value for mean d13C

% Unfertilized parameters 
% Parameters from Table 1 for unconstrained model
par(1).unc.x0 = -38.18; par(1).unc.k = 1.42;
par(1).unc.lambda = [2.54 1.26 1.01 1.23 1.49];

% Parameters from Table 1 for constrained model
par(1).con.x0 = -38.18;  par(1).con.k = 1.42;
par(1).con.lambda = [2.52 1.26 1.21 1.21 1.21];

% Parameters from Table 1 for hypoexponential model
par(1).hyp.x0 = -38.13; par(1).hyp.k = 1.40; par(1).hyp.sd = 1.40;
par(1).hyp.lambda = [2.27 1.84 3.46 100 100.001];

% Fertilized parameters 
% Parameters from Table 1 for unconstrained model
par(2).unc.x0 = -37.80;  par(2).unc.k = 1.21;
par(2).unc.lambda = [1.06 1.91 0.89 0.93 0.73];

% Parameters from Table 1 for constrained model
par(2).con.x0 = -37.79;  par(2).con.k = 1.24;
par(2).con.lambda = [1.33 1.33 0.90 0.90 0.72];

% Parameters from Table 1 for hypoexponential model
par(2).hyp.x0 = -37.73;  par(2).hyp.k = 1.29; par(2).hyp.sd = 1.88;
par(2).hyp.lambda = [1.27 100 1.41 100.001 1.72];

%% Plotting

orders_show = [1 5]; % Choose orders to display
order_titles = {'Order 1','Order 5'};

treat_names = {'Unfertilized','Fertilized'};
Ninds = {dat.Nind0, dat.Nind1};

cols = lines(4); 
ms = 10;
lw = 4;

figure('Color','w','Position',[100 100 1450 750]);

tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

gobjects(2,3);

% Left and Middle Panels
for q = 1:2  % treatment row
    for j = 1:2  % order column
        ord = orders_show(j);
        nexttile((q-1)*3 + j);
        hold on;

        ind0 = dat.Cind1 & dat.Order == ord & Ninds{q};

        tvals = unique(dat.t(ind0));
        ymean = NaN(size(tvals));
        yerr  = NaN(size(tvals));

        for k = 1:length(tvals)
            ii = ind0 & dat.t == tvals(k);
            x = dat.d13C(ii);

            ymean(k) = mean(x,'omitnan');
            yerr(k) = std(x,'omitnan'); % SD across roots
        end

        % Data
        errorbar(tvals,ymean,yerr,'ko',...
            'MarkerFaceColor','w',...
            'MarkerSize',ms,...
            'LineWidth',lw,...
            'CapSize',6);

        tgrid = linspace(0,4,200);
        % Means of models
        % Unconstrained exponential mean
        m_unc = mean_exp_model(tgrid, xeq, par(q).unc.x0, ...
            par(q).unc.k, par(q).unc.lambda(ord));

        % Constrained exponential mean
        m_con = mean_exp_model(tgrid, xeq, par(q).con.x0, ...
            par(q).con.k, par(q).con.lambda(ord));

        % Hypoexponential mean
        [m_hyp, sd_hyp] = hypo_moments(tgrid, xeq, par(q).hyp.x0, ...
            par(q).hyp.k, par(q).hyp.lambda(1:ord), par(q).hyp.sd);
        
        fill([tgrid fliplr(tgrid)], ...
            [(m_hyp-sd_hyp)' fliplr((m_hyp+sd_hyp)')], ...
            cols(3,:), ...
            'FaceAlpha',0.15, ...
            'EdgeColor','none', ...
            'HandleVisibility','off'); 

        plot(tgrid,m_unc,'-','LineWidth',lw,'Color',cols(1,:));
        plot(tgrid,m_con,'--','LineWidth',lw,'Color',cols(2,:));
        plot(tgrid,m_hyp,':','LineWidth',lw+0.5,'Color',cols(3,:));

        title(sprintf('%s, %s',order_titles{j},treat_names{q}));
        xlabel('Time since FACE shutdown (years)');
        ylabel('Mean \delta^{13}C (‰)');
        box on;
        set(gca,'FontSize',11);

        if q == 1 && j == 1
            legend({'Data mean \pm SD','Unconstrained','Constrained','Hypoexponential'},...
                'Location','best');
        end
    end
end

% Right panel - longevity plot 
nexttile(3,[2 1]);
hold on;

root_orders = 1:5;

% Values from Table 1
lon_N0_unc = 1./par(1).unc.lambda;
lon_N0_con = 1./par(1).con.lambda;
lon_N0_hyp = cumsum(1./par(1).hyp.lambda);

lon_N1_unc = 1./par(2).unc.lambda;
lon_N1_con = 1./par(2).con.lambda;
lon_N1_hyp = cumsum(1./par(2).hyp.lambda);

cN0 = [0.0000 0.4470 0.7410];
cN1 = [0.8500 0.3250 0.0980];

plot(root_orders,lon_N0_unc,'o-','Color',cN0,'LineWidth',lw,'MarkerFaceColor','w');
plot(root_orders,lon_N0_con,'s--','Color',cN0,'LineWidth',lw,'MarkerFaceColor','w');
plot(root_orders,lon_N0_hyp,'^:','Color',cN0,'LineWidth',lw+0.5,'MarkerFaceColor','w');

plot(root_orders,lon_N1_unc,'o-','Color',cN1,'LineWidth',lw,'MarkerFaceColor','w');
plot(root_orders,lon_N1_con,'s--','Color',cN1,'LineWidth',lw,'MarkerFaceColor','w');
plot(root_orders,lon_N1_hyp,'^:','Color',cN1,'LineWidth',lw+0.5,'MarkerFaceColor','w');

xlabel('Root order');
ylabel('Mean longevity (years)');
title('Mean root longevity by order');
set(gca,'FontSize',11,'XTick',root_orders);
box on;

legend({...
    'Unfertilized unconstrained',...
    'Unfertilized constrained',...
    'Unfertilized hypoexponential',...
    'Fertilized unconstrained',...
    'Fertilized constrained',...
    'Fertilized hypoexponential'},...
    'Location','northwest');


%% Functions for Moment Computation
function m = mean_exp_model(t, xeq, x0, k, lambda)
% Mean isotope trajectory for exponential lifespan model.
%
% k       : carbon-pool relaxation rate
% lambda  : root replacement/turnover rate
%
% Formula:
% m(t) = xeq + (x0 - xeq) * (lambda*exp(-k*t) - k*exp(-lambda*t))/(lambda-k)
% Reference: Poll (2026, Tree Physiology, Supplemental Material).

    t = t(:);
    A = x0 - xeq;

    tol = 1e-10;

    if abs(lambda - k) < tol
        m = xeq + A.* exp(-k*t).* (1 + k*t);
    else
        m = xeq + A.*(lambda*exp(-k*t) - k*exp(-lambda*t)) ./ (lambda - k);
    end
end


function [m1, sd] = hypo_moments(t, xeq, x0, k, lambdas, sig)
% Numerically computes mean and SD of isotope distribution for hypo model.
% sig is fitted isotope SD.

    t = t(:);
    lambdas = lambdas(:);
    mu = sum(1 ./ lambdas);

    mB = @(s) xeq + (x0 - xeq).*exp(-k*s);
    mB2 = @(s) sig^2 + mB(s).^2;

    m1 = zeros(size(t));
    m2 = zeros(size(t));

    for ii = 1:length(t)
        ti = t(ii);

        tau = linspace(0,ti,300);
        S = hypo_survival(ti - tau, lambdas);
        R = hypo_residual_survival(ti, lambdas);

        m1(ii) = R*x0 + (1/mu)*trapz(tau, S(:).*mB(tau(:)));
        m2(ii) = R*(sig^2 + x0^2) + (1/mu)*trapz(tau, S(:).*mB2(tau(:)));
    end

    varx = max(m2 - m1.^2, 0);
    sd = sqrt(varx);
end

function S = hypo_survival(t, lambdas)
% Survival for sequential hypoexponential phases.
% Handles repeated/near-repeated rates.

    t = t(:);
    lambdas = lambdas(:);
    n = length(lambdas);

    T = zeros(n,n);
    for i = 1:n
        T(i,i) = -lambdas(i);
        if i < n
            T(i,i+1) = lambdas(i);
        end
    end

    alpha = zeros(1,n);
    alpha(1) = 1;

    one = ones(n,1);
    S = zeros(size(t));

    for ii = 1:length(t)
        S(ii) = alpha*expm(T*t(ii))*one;
    end
end

function R = hypo_residual_survival(t, lambdas)
% Residual-life survival contribution:
% R(t) = (1/mu) * integral_t^infty S(u) du

    lambdas = lambdas(:);
    n = length(lambdas);
    mu = sum(1 ./ lambdas);

    T = zeros(n,n);
    for i = 1:n
        T(i,i) = -lambdas(i);
        if i < n
            T(i,i+1) = lambdas(i);
        end
    end

    alpha = zeros(1,n);
    alpha(1) = 1;

    one = ones(n,1);

    R = (alpha*expm(T*t)*(-T\one)) / mu;
end