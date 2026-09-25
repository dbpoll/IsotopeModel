clear; clc
%%  Supplemental Script - Run Fits for Isotope Dataset
%
% Generates fitted parameter values and fit statistics used in the 
% Supplemental Material table.
%
% Required File(s): 'face_dataset.csv'
%
% Required toolboxes: 
%      - Optimization Toolbox
%      - Statistics and Machine Learning Toolbox
%
% 
% Note: the optimal values and negative loglikelihood reported in the 
% Supplemental Material are contained in variable structures 
% 'theta_hat' and 'nll_hat'. 
%

% Fitting Options
model.cat = 'all'; % 'all' or 'diam'
model.type = 'exp'; % 'exp', 'hypo', 'lognormal'
model.treatment = 'unfert'; %'unfert', 'fert', or 'comb'


%% Parameters
pars.xeq = -28.85; % equilibrium d13C mean
pars.ngrid = 1000; % # of grid points for integral

% Initial guesses for parameters to fit
theta.x0 = -38.0;
theta.sd = 1.8; 
theta.k = 1.2;

%% Load Data Into Table and Arrays
dtab = readtable('face_dataset_diamm.csv'); % read data table
tshft = 18; % time offset in days
t0 = min(dtab.DateSampled) - days(tshft);

% Select nitrogen treatment
ind = strcmp(dtab.CO2Treatment,'Elevated');
switch model.treatment
    case 'unfert'
        ind = ind & strcmp(dtab.NTreatment,'Unfertilized');
    case 'fert'
        ind = ind & strcmp(dtab.NTreatment,'Fertilized');
    case 'comb'
        % no changes needed
end

% Select diameter category/categories
switch model.cat
    case 'all'
         model.ngroup = 1;
        ind = ind & dtab.AvgDiam_mm < 2;
        dat.Group = ones(sum(ind),1);
    case 'diam'
        model.ngroup = 3;
        ind = ind & dtab.AvgDiam_mm < 2;
        diam = dtab.AvgDiam_mm(ind);
        dat.Group = zeros(sum(ind),1);
        dat.Group(diam >= 0   & diam < 0.5) = 1;
        dat.Group(diam >= 0.5 & diam < 1.0) = 2;
        dat.Group(diam >= 1.0 & diam < 2.0) = 3;
        theta.lambda = [1 2 3]; % initial lambda choices
end

%% Initial lifespan parameters

switch model.type

    case 'exp'
        theta.lambda = ones(1,model.ngroup);

    case 'hypo'
        % two sequential exponential rates
        theta.lambda1 = 1.8*ones(1,model.ngroup);
        theta.lambda2 = 2.2*ones(1,model.ngroup);

    case 'lognormal'
        % log(T) ~ N(logmu, logsd^2)
        theta.logmu = -0.3*ones(1,model.ngroup);
        theta.logsd = 0.8*ones(1,model.ngroup);

end



% Store only data actually used in fit
dat.t    = days(dtab.DateSampled(ind) - t0)/365;
dat.d13C = dtab.d13C(ind);

%% Run NLL Optimizer
switch model.type
    case 'exp'
        theta0_vec = [theta.x0, theta.sd, theta.k, ...
                      theta.lambda];
    case 'hypo'
        theta0_vec = [theta.x0, theta.sd, theta.k, ...
                      theta.lambda1, theta.lambda2];
    case 'lognormal'
        theta0_vec = [theta.x0, theta.sd, theta.k, ...
                      theta.logmu, theta.logsd];
end

% Bounds
switch model.type
    case 'exp'
        lb = [-45, 0.05, 0.001, ...
              0.001*ones(1,model.ngroup)];
        ub = [-30, 5, 10, ...
              100*ones(1,model.ngroup)];
    case 'hypo'
        lb = [-45, 0.05, 0.001, ...
              0.001*ones(1,2*model.ngroup)];
        ub = [-30, 5, 10, ...
              100*ones(1,2*model.ngroup)];
    case 'lognormal'
        lb = [-45, 0.05, 0.001, ...
              -5*ones(1,model.ngroup), ...      % logmu
              0.05*ones(1,model.ngroup)];      % logsd
        ub = [-30, 5, 10, ...
              5*ones(1,model.ngroup), ...
              5*ones(1,model.ngroup)];
end

% Optimization options
options = optimoptions('fmincon', ...
    'Display','iter', ...
    'Algorithm','interior-point', ...
    'MaxFunctionEvaluations',1e5, ...
    'MaxIterations',1e4);

% Fit model
[theta_hat_vec,nll_hat,exitflag,output] = fmincon( ...
    @(x) negloglik(x,dat,pars,model), ...
    theta0_vec,[],[],[],[],lb,ub,[],options);

theta_hat = unpackTheta(theta_hat_vec,model);


%% Print Out Stuff
switch model.type
    case 'exp'
        meanlife = 1 ./ theta_hat.lambda;
    case 'hypo'
        meanlife = 1./theta_hat.lambda1 + 1./theta_hat.lambda2;
    case 'lognormal'
        meanlife = exp(theta_hat.logmu + theta_hat.logsd.^2/2);
end

fprintf('NLL = %.4f\n',nll_hat);
disp('Mean lifespan (years):')
disp(meanlife)



%% Helper Functions
% Initial Distribution phi0 - Assumed to be N(m0,sig0)
function y = phi0_pdf(X, theta)
mu0 = theta.x0;
sig0 = theta.sd;
y = normpdf(X,mu0,sig0);
end

% Carbon-pool Distribution B(X,t), also normal with mean muB
function y = B_pdf(X,t,pars,theta)
    xeq = pars.xeq;
    x0 = theta.x0;
    sig = theta.sd;
    k = theta.k;

    muB = xeq + (x0-xeq).*exp(-k*t);
    y = normpdf(X,muB,sig);
end


function [R,S,renewrate] = buildFunctions(t,tau,theta,group,model)

    age = t - tau;

    switch model.type

        %% Exponential
        case 'exp'

            lambda = theta.lambda(group);

            S = exp(-lambda*age);
            R = exp(-lambda*t);

            renewrate = lambda;


        %% Two-stage hypoexponential
        case 'hypo'

            l1 = theta.lambda1(group);
            l2 = theta.lambda2(group);

            mu = 1/l1 + 1/l2;
            renewrate = 1/mu;

            % Erlang-2 limit if rates are nearly equal
            if abs(l1-l2) < 1e-6

                lambda = 0.5*(l1+l2);

                S = exp(-lambda*age).*(1 + lambda*age);

                R = exp(-lambda*t).*(1 + lambda*t/2);

            else

                S = (l2*exp(-l1*age) - ...
                     l1*exp(-l2*age)) / (l2-l1);

                intS = (l2/l1*exp(-l1*t) - ...
                        l1/l2*exp(-l2*t)) / (l2-l1);

                R = intS / mu;

            end


        %% Lognormal
        case 'lognormal'

            muLog = theta.logmu(group);
            sigLog = theta.logsd(group);

            % Mean lifespan
            mu = exp(muLog + sigLog^2/2);
            renewrate = 1/mu;

            % Ordinary survival
            S = 1 - logncdf(age,muLog,sigLog);

            % Equilibrium residual survival
            if t == 0
                R = 1;
            else
                St = 1 - logncdf(t,muLog,sigLog);

                Etail = mu * normcdf( ...
                    (muLog + sigLog^2 - log(t))/sigLog);

                R = (Etail - t*St)/mu;
            end

    end
end

% Construct Isotope Distribution Model
function P = P_pdf(X,t,pars,theta,group,model)

    tau = linspace(0,t,pars.ngrid); % integration grid
    B = B_pdf(X,tau,pars,theta); % values for B(X,s)

    [R,S,renewrate] = buildFunctions( ...
        t,tau,theta,group,model); % survival, residual, mu

    P0 = R .* phi0_pdf(X,theta); % initial distribution
    Pbirth = renewrate * trapz(tau,S.*B); % replacement distribution

    P = P0 + Pbirth;
end

% Function for computing negative log-likelihood
function nll = negloglik(theta_vec,dat,pars,model)
    theta = unpackTheta(theta_vec,model);
    LL = 0;

    for i = 1:length(dat.d13C)
        X = dat.d13C(i);
        t = dat.t(i);
        group = dat.Group(i);

        P = P_pdf(X,t,pars,theta,group,model);

        if ~isfinite(P) || P <= 0
            nll = 1e100;
            return
        end
        LL = LL + log(P);
    end
    nll = -LL;
end


function theta = unpackTheta(theta_vec,model)

theta.x0 = theta_vec(1);
theta.sd = theta_vec(2);
theta.k  = theta_vec(3);

n = model.ngroup;

switch model.type

    case 'exp'
        theta.lambda = theta_vec(4:3+n);

    case 'hypo'
        theta.lambda1 = theta_vec(4:3+n);
        theta.lambda2 = theta_vec(4+n:3+2*n);

    case 'lognormal'
        theta.logmu = theta_vec(4:3+n);
        theta.logsd = theta_vec(4+n:3+2*n);

end

end