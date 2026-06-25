%% Main Script - Run Isotope Model Under Various Assumptions
%
% Generates fitted parameter values and fit statistics used in Tables 1 and 2.
% Set model.type and model.treatment below, then rerun the script for each
% model/treatment combination.
%
% Required File(s): 'face_dataset.csv'
%
% Required toolboxes: 
%      - Optimization Toolbox
%      - Statistics and Machine Learning Toolbox
% 
% Note:
% The hypoexponential model is slower to fit than the exponential models.
% For numerical stability, its survival and residual-survival functions are
% evaluated using matrix exponentials.
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

%% Define Parameters
% Fixed Parameters 
pars.xeq = -28.85; % equilibrium d13C mean
pars.ngrid = 1000; % # of grid points for integral

model.type = 'unc'; % 'unc', 'con', or 'hypo'
model.treatment = 'unfert'; %'unfert', 'fert', or 'comb'

% Initial guesses for parameters to fit
theta.x0 = -38.0;
theta.sd = 1.6; 
theta.k = 1.16;
theta.lambda = [1 2 4 10 10]; % initial lambda choices

%% Data for Fitting
switch lower(model.treatment)
    case 'unfert'
        fitInd = dat.Cind1 & dat.Nind0; % elevated CO2, unfertilized
        nLabel = "Unfertilized";
    case 'fert'
        fitInd = dat.Cind1 & dat.Nind1; % elevated CO2, fertilized
        nLabel = "Fertilized";
    case 'comb'
        fitInd = dat.Cind1 & (dat.Nind0 | dat.Nind1); % elevated CO2, both N treatments
        nLabel = "Combined";
    otherwise
        error('Unknown treatment option. Use unfert, fert, or comb.');
end

% Put selected data into structure
datfit.t = dat.t(fitInd);
datfit.Order = dat.Order(fitInd);
datfit.d13C = dat.d13C(fitInd);


%% Run Optimization Algorithm from Data & Model
theta0_vec = [theta.x0, theta.sd, theta.k, theta.lambda];

% Bounds
lb = [-45, 0.05, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001];
ub = [-30, 5.00, 10.000, 100.000, 100.000, 100.000, 100.000, 100.000];

% Inequality Constraints - No constraints on lifespan
A = [];
b = [];

if strcmpi(model.type,'con')
    % Inequality Constraints - constraints by order
    A = zeros(4,8);
    b = zeros(4,1);
    A(1,5) =  1;
    A(1,4) = -1;
    A(2,6) =  1;
    A(2,5) = -1;
    A(3,7) =  1;
    A(3,6) = -1;
    A(4,8) =  1;
    A(4,7) = -1;
end

% No linear equality constraints
Aeq = [];
beq = [];

% Define optimization options
options = optimoptions('fmincon', ...
    'Display','iter', ...
    'Algorithm','interior-point', ...
    'MaxFunctionEvaluations',1e5, ...
    'MaxIterations',1e4);

[theta_hat_vec, nll_hat, exitflag, output] = fmincon(@(x) negloglik(x, datfit, pars, model), ...
    theta0_vec, A, b, Aeq, beq, lb, ub, [], options);

%% Display Optimal Parameter Values
theta_hat = unpackTheta(theta_hat_vec);

switch model.type
    case 'unc'
        typeLabel = "Unconstrained";
    case 'con'
        typeLabel = "Constrained";
    case 'hypo'
        typeLabel = "Hypoexponential";
    otherwise
        typeLabel = "Unknown model type";
end


fprintf('\n%s Model - %s\n', typeLabel, nLabel);
fprintf('----------------------------------------\n');
fprintf('\nParameter values:\n');
fprintf('x0:     %.4f\n', theta_hat.x0);
fprintf('sd:     %.4f\n', theta_hat.sd);
fprintf('k:      %.4f\n', theta_hat.k);

% Turnover rates/average lifespans
fprintf('\nTurnover rates by order:\n');
fprintf('----------------------------------------\n');

if strcmpi(model.type,'hypo')
    fprintf('Hypoexponential model: reporting effective order-level turnover rates.\n');
    fprintf('Parentheses show fitted intrinsic stage parameter lambda_i.\n\n');
    for order = 1:5
        mu_order = sum(1 ./ theta_hat.lambda(1:order));
        lambda_eff = 1 / mu_order;
        lambda_param = theta_hat.lambda(order);
        fprintf('Order %d: %.4f  (lambda_%d = %.4f), mean lifespan = %.4f years\n', ...
            order, lambda_eff, order, lambda_param, mu_order);
    end
else
    fprintf('Exponential model: fitted lambda_i are order-level turnover rates.\n\n');
    for order = 1:5
        lambda_eff = theta_hat.lambda(order);
        mu_order = 1 / lambda_eff;

        fprintf('Order %d: %.4f  (lambda_%d = %.4f), mean lifespan = %.4f years\n', ...
            order, lambda_eff, order, theta_hat.lambda(order), mu_order);
    end
end
fprintf('\n');

% likelihood
fprintf('\nNegative log-likelihood: %.4f\n', nll_hat);
fprintf('\n');


%% Functions
% Initial Distribution phi0 - Assumed to be N(m0,sig0)
function y = phi0_pdf(X, theta)
mu0 = theta.x0;
sig0 = theta.sd;
y = normpdf(X,mu0,sig0);
end

% Carbon-pool Distribution B(X,t)
function y = B_pdf(X,t,pars,theta)
    xeq = pars.xeq;
    x0 = theta.x0;
    sig = theta.sd;
    k = theta.k;

    muB = xeq + (x0-xeq).*exp(-k*t);
    y = normpdf(X,muB,sig);
end

function [R, S, renewrate] = buildFunctions(t,tau,theta,order,model)
    % Hypoexponential Model
    if strcmpi(model.type,'hypo')
        rates = theta.lambda(1:order);
        rates = rates(:);
        n = length(rates);
        % Build phase-type generator matrix
        Q = zeros(n,n);
        for i = 1:n
            Q(i,i) = -rates(i);
            if i < n
                Q(i,i+1) = rates(i);
            end
        end

        alpha = zeros(1,n);
        alpha(1) = 1;

        onevec = ones(n,1);

        % Mean lifespan
        v = (-Q)\onevec;
        mu = alpha*v;

        % Residual survival R(t) at time t
        R = (alpha * expm(Q*t) * v) / mu;

        % Survival S(t - tau)
        age = t - tau;
        S = zeros(size(age));

        for j = 1:length(age)
            S(j) = alpha * expm(Q*age(j)) * onevec;
        end
        renewrate = 1/mu;
    
    % Exponential Model
    else
        lamb = theta.lambda(order);
        R = exp(-lamb*t);
        S = exp(-lamb*(t - tau));
        renewrate = lamb;
    end
end

% Full continuum distribution P(X,t)
function y = P_pdf(X,t,pars,theta,order,model)
    tau = linspace(0,t,pars.ngrid); % time grid for integral
    B = B_pdf(X,tau,pars,theta); % Evaluate points for B(X,tau);

    [R, S, renewrate] = buildFunctions(t,tau,theta,order,model); % compute survival/residual

    P0 = R .* phi0_pdf(X,theta); % build initial distribution
    Pbirth = renewrate*trapz(tau, S.*B ); % build birthed distribution
    y = P0 + Pbirth; % full distribution
end 

% Function for computing negative loglikelihood
function nll = negloglik(theta_vec, datfit, pars, model)
    theta = unpackTheta(theta_vec); % put back into structure
    
    LL = 0;
    for i = 1:length(datfit.d13C)

        X = datfit.d13C(i);
        t = datfit.t(i);
        order = datfit.Order(i);

        P = P_pdf(X, t, pars, theta, order, model);

        if ~isfinite(P) || P <= 0
            nll = 1e100;
            return
        end
        LL = LL + log(P);
    end

    nll = -LL;
end

% transform parameter vector back into structure
function theta = unpackTheta(theta_vec)

    theta.x0 = theta_vec(1);
    theta.sd = theta_vec(2);
    theta.k  = theta_vec(3);

    theta.lambda = theta_vec(4:8);

end