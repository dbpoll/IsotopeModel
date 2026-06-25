%% Figure 2: Stochastic Convergence to Continuum Distribution
% This script generates Figure 2 from Poll (2026).
%
% The script compares a finite stochastic root replacement model against
% the corresponding continuum isotope-distribution model. Stochastic
% simulations are shown for two values of N, and convergence is quantified
% using RMSE over a range of N and isotope-bin widths.
% 
% Required toolboxes:
%     - Statistics and Machine Learning Toolbox
%     - Curve Fitting Toolbox
% 
clear; clc; close all;

rng(1); % seed for reproducibility 

%% Base parameters
s = make_base_settings(); % parameters for simulation

N_small = 1e3; % N = 1000 roots
N_large = 1e4; % N = 10000 roots

out_small = run_rootdist_sim(s,N_small,true,true);
out_large = run_rootdist_sim(s,N_large,true,false);

% RMSE sweep
N_list = [20 50 100 200 500 1000]; % number of roots for stochastic sim
dC_width_list = [0.25 0.10 0.05]; % isotope-bin widths delta X

RMSE = NaN(length(dC_width_list),length(N_list));

for a = 1:length(dC_width_list)

    dx = dC_width_list(a);
    s2 = s;
    s2.grid.dC = -45:dx:-23;

    % compute model once per grid
    out_model = run_rootdist_sim(s2,N_list(1),false,true);
    P_model = out_model.P_model;

    for b = 1:length(N_list)
        fprintf('RMSE sweep: dx = %.3f, N = %d\n',dx,N_list(b));
        out = run_rootdist_sim(s2,N_list(b),true,false);

        RMSE(a,b) = sqrt(mean((out.P_stoch_raw(:)-P_model(:)).^2,'omitnan'));
    end
end

%% Plots for Figure 2

fs = 22; % font-size
lw = 4; % line-width
ms = 12; % markersize

figure('Color','w','Position',[100 100 1200 950]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

% Shared axis settings
xlims = [0 4];
ylims = [-45 -25];
xtks = 0:1:4;
ytks = -45:5:-25;
clims = [0 0.2];

% Top Left Panel - Small N
[Xsmall,Tsmall] = meshgrid(out_small.s.grid.dC,out_small.s.grid.t);

nexttile;
pcolor(Tsmall,Xsmall,out_small.P_stoch_raw');
shading flat; colormap hot;
axis([xlims ylims]);
clim(clims);
xlabel('t (years)');
ylabel('\delta^{13}C');
set(gca,'FontSize',fs,'XTick',xtks,'YTick',ytks);
box on;

% Top Right Panel - Large N
[Xlarge,Tlarge] = meshgrid(out_large.s.grid.dC,out_large.s.grid.t);

nexttile;
pcolor(Tlarge,Xlarge,out_large.P_stoch_raw');
shading flat; colormap hot;
axis([xlims ylims]);
clim(clims);
xlabel('t (years)');
ylabel('\delta^{13}C');
set(gca,'FontSize',fs,'XTick',xtks,'YTick',ytks);
box on;

% Bottom Left Panel - Continuum Model (small grid)
nexttile;
pcolor(Tsmall,Xsmall,out_small.P_model');
shading flat; colormap hot;
axis([xlims ylims]);
clim(clims);
xlabel('t (years)');
ylabel('\delta^{13}C');
set(gca,'FontSize',fs,'XTick',xtks,'YTick',ytks);
box on;

% Optional shared colorbar for pcolor plots
cb = colorbar;
cb.Layout.Tile = 'east';
cb.Label.String = 'Probability density';
cb.Label.FontSize = fs;
cb.FontSize = fs;

% Bottom Right Panel - RMSE convergence
nexttile;
hold on;

line_styles = {'-','--',':'};
markers = {'o','s','^'};

for a = 1:length(dC_width_list)

    dx = dC_width_list(a);

    good = isfinite(RMSE(a,:)) & RMSE(a,:) > 0;
    pfit = polyfit(log(N_list(good)),log(RMSE(a,good)),1);

    alpha = pfit(1);
    Cfit = exp(pfit(2));

    Nfit = logspace(log10(min(N_list)),log10(max(N_list)),200);
    RMSE_fit = Cfit*Nfit.^alpha;

    ls = line_styles{mod(a-1,length(line_styles))+1};
    mk = markers{mod(a-1,length(markers))+1};

    plot(N_list,RMSE(a,:),mk,...
        'LineStyle','none',...
        'MarkerSize',ms,...
        'MarkerFaceColor','auto',...
        'LineWidth',lw,...
        'DisplayName',sprintf('\\Delta X = %.2f data',dx));

    plot(Nfit,RMSE_fit,ls,...
        'LineWidth',lw,...
        'DisplayName',sprintf('\\Delta X = %.2f fit, slope %.2f',dx,alpha));
end

set(gca,'XScale','log','YScale','log',...
    'FontSize',fs,...
    'LineWidth',1.5);
xlabel('N');
ylabel('RMSE');
legend('Location','southwest','FontSize',14);
box on;


%% Functions
% Runs stochastic sim and model, stores approximate distributions
function out = run_rootdist_sim(s,N,run_stoch,run_model)

s.n_sample.stoch = N;

out.s = s;
out.P_stoch_raw = [];
out.P_model = [];

if run_stoch
    out.P_stoch_raw = sim_stochastic_eq(s);
end

if run_model
    out.P_model = solve_distribution_eq(s);
end
end

% Parameters for simulation 
function s = make_base_settings()
s.par.r_iso = 1; % turnover of carbon pool
s.par.mdC_0 = -38; % initial mean d13C
s.par.mdC_eq = -28.85; % equilibrium mean d13C
s.par.sd = 1.5; % standard deviation of Gaussian for pool

s.par.lam_life = 1.46; % lifespan rate for root replacement

s.grid.t = linspace(0,5,501); % time grid
s.grid.dC = linspace(-45,-23,201); % default d13C width

s.n_sample.stoch = 100; % default samples
s.n_sample.model = 1e6; % default samples

s.Tol.Pmodel = 1e-4; % tolerances
s.Tol.rel = 1e-6;
s.Tol.abs = 1e-12;
end

% Define and generate random samples from chosen distribution L(X)
function X = gen_distribution_samples(par,n_sample)
X = exprnd(1/par.lam_life,n_sample,1); % generate lifespan samples
end

% Define distribution for d13C values within the tree pool
function phi_vals = gen_distribution_phi(dC_grid,t_grid,par)

mphi = par.mdC_eq + (par.mdC_0 - par.mdC_eq).*exp(-par.r_iso*t_grid); % mean of carbon pool

phi_vals = zeros(length(dC_grid),length(t_grid));
for i = 1:length(t_grid)
    phi_vals(:,i) = normpdf(dC_grid,mphi(i),par.sd); % B(X,t) - carbon pool distribution
end

end

% Define Stochastic Simulation Function
function P_stoch_raw = sim_stochastic_eq(s)

dC_counts = zeros(length(s.grid.dC),length(s.grid.t));
% Given distribution L(X), get samples
% Generate lifespans and initialize residual lifetimes from the stationary
% renewal distribution using size-biased resampling.
X = gen_distribution_samples(s.par,s.n_sample.stoch); % samples from L(X)
w = X/sum(X); % weights
inds = randsample(s.n_sample.stoch,s.n_sample.stoch,'true',w); %size-weighted resampling

% Initialize residual ages, d13C values at t0
Y = rand(s.n_sample.stoch,1).*X(inds); % initial residual ages
dC = normrnd(s.par.mdC_0,s.par.sd,s.n_sample.stoch,1); % initial d13C values
store_dC = zeros(s.n_sample.stoch,length(s.grid.t));

% Run stochastic simulation
t = s.grid.t(1); T = s.grid.t(end); % beginning and end time
i = 1;
while t < T
    [~,idx] = min(Y); % find the next earliest death
    t = t + Y(idx); % update time to earliest event

    while i <= length(s.grid.t) && s.grid.t(i) < t  % get d13 values at time grid
        store_dC(:,i) = dC;
        i = i+1;
    end
    Y = Y - Y(idx); % decrement remaining time alive

    Y_renew = gen_distribution_samples(s.par,1);
    Y(idx) = Y_renew;

    md13C = s.par.mdC_eq + (s.par.mdC_0 - s.par.mdC_eq).*exp(-s.par.r_iso*t);
    dC(idx) = normrnd(md13C,s.par.sd,1,1);

end

% Count d13C values into bins on grid s.grid.dC
edge_first = s.grid.dC(1) - (s.grid.dC(2)-s.grid.dC(1))/2;
edge_last = s.grid.dC(end) + (s.grid.dC(end)-s.grid.dC(end-1))/2;
edge_int = (s.grid.dC(2:end) + s.grid.dC(1:end-1))/2;
edges = [edge_first, edge_int, edge_last];


for i = 1:length(s.grid.t)
    dC_counts(:,i) = histcounts(store_dC(:,i),edges)';
end
P_stoch_raw = dC_counts ./ (sum(dC_counts, 1) .* diff(edges)');

end

% Define Model Integration Function
% Note: The lifespan density L is estimated from samples rather than specified
% analytically so that the same solver can be used for more general lifespan
% distributions.
function P_model = solve_distribution_eq(s)

% Generate approximate distribution of lifespans L(X) for model
X = gen_distribution_samples(s.par,s.n_sample.model); % samples from L(X)
mu = mean(X); % mean of lifespan L(X)

X_edges = linspace(s.grid.t(1),s.grid.t(end),round(5*sqrt(s.n_sample.model))+1);
X_centers = (X_edges(1:end-1) + X_edges(2:end))/2;
dX = X_centers(2) - X_centers(1);
X_counts = histcounts(X,X_edges);
L_Xmass = X_counts./(sum(X_counts.*diff(X_edges))); % estimated probability mass
L_Xmass_sm = smooth(X_centers,L_Xmass,0.05,'loess'); % smooth the mass function
L_Xmass_sm = max(L_Xmass_sm,0); % catch any small errors
L_X = L_Xmass_sm/trapz(X_centers,L_Xmass_sm); % estimated probability density

% Define Distribution of Tree Carbon Pools
% -> we require phi(dC,0) = phi_0(dC) for continuity/stationary
phi = gen_distribution_phi(s.grid.dC,X_centers,s.par);
phi0 = gen_distribution_phi(s.grid.dC,0,s.par);

% Compute Survival Probability Function: S_L(t)
L_interp = @(t) interp1(X_centers, L_X, t, 'linear', 0);
S_L = arrayfun(@(t) integral(L_interp, t, X_centers(end), 'RelTol',s.Tol.rel,'AbsTol',s.Tol.abs), X_centers);

% Compute Initial Residual Survival Function: S_R(t)
% -> Can be computed as S_R(t) = int_t^\infty R_0(X) dX,
% -> R_0(X) is the initial residual distribution (time alive left)
% -> If stationary (well-mixed), then S_R(X) = (1/mu) int_t^inf S_L(X) dX
S_L_interp =  @(t) interp1(X_centers, S_L, t, 'linear', 0);
S_R = (1/mu)*arrayfun(@(t) integral(S_L_interp, t, X_centers(end), 'RelTol',s.Tol.rel,'AbsTol',s.Tol.abs), s.grid.t);

% Compute distribution P(dC,t) = S_R*phi(dC,0) + (1/mu)*(conv(S_L,phi))
P_model = zeros(length(s.grid.dC),length(s.grid.t));
for i = 1:length(s.grid.dC)
    conv_temp_full = conv(phi(i,:),S_L,'full').*dX;
    conv_temp = conv_temp_full(1:length(X_centers));
    P_model(i,:) = phi0(i)*S_R + (1/mu)*interp1(X_centers,conv_temp,s.grid.t,'pchip','extrap');
end

er =  sum((1 - trapz(s.grid.dC,P_model)).^2 )/length(s.grid.t);
P_model = P_model./trapz(s.grid.dC,P_model); % normalize for consistency

if er > s.Tol.Pmodel
    disp('Distribution Poor - Consider More Samples/Finer Grid')
end

end
