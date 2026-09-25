clear; clc
%%  Supplemental Script - Run Fits for Minirhizotron Dataset
%
% Generates fitted parameter values and fit statistics used in the 
% Supplemental Material table.
%
% Required File(s): 'root_ages.csv'
%
% Required toolboxes: 
%      - Optimization Toolbox
%      - Statistics and Machine Learning Toolbox
%
% 
% Note: the values reported in the Supplemental Material are contained
% in variable structures 'fit' and 'nll'. 
%


%% Modeling Parameters
model.cat = 'all'; % 'all or 'diam'
model.treatment = 'comb'; % 'fert', 'unfert', or 'comb'

mtab = readtable('root_ages.csv'); % read data table

% Select Elevated CO2 Data
ind = strcmpi(mtab.co2treat,'ELEVATED');
ind = ind & mtab.session > 63; % check on this?

% Nitrogen treatment
switch model.treatment
    case 'unfert'
        ind = ind & mtab.ntreat == 0;
    case 'fert'
        ind = ind & mtab.ntreat == 1;
    case 'all'
        % keep both nitrogen treatments
end

% Convert diameter to mm

% Diameter categories
diam_mm = 1000 * mtab.maxdiam;
ind = ind & diam_mm < 2;
switch model.cat
    case 'all'
        dat = mtab(ind,:);
        dat.Group = ones(height(dat),1);
    case 'diam'
        dat = mtab(ind,:);
        d = 1000 * dat.maxdiam;
        dat.Group = zeros(height(dat),1);
        dat.Group(d >= 0   & d < 0.5) = 1;
        dat.Group(d >= 0.5 & d < 1.0) = 2;
        dat.Group(d >= 1.0 & d < 2.0) = 3;
end
dat.maxdiam = mtab.maxdiam(ind)*1000; 
dat.rightcens = strcmpi(mtab.rightcens(ind),'TRUE');

%% Fit all survival models and plot

models = {'exp','hypo','lognormal'};
ngroup = max(dat.Group);

figure;

for g = 1:ngroup

    indg = dat.Group == g;
    age  = dat.ageyears(indg);
    cens = dat.rightcens(indg);

    % Kaplan-Meier
    [Skm,tkm] = ecdf(age, ...
        'Censoring',cens, ...
        'Function','survivor');

    tfit = linspace(0,max(age),300);

    subplot(1,ngroup,g)
    stairs(tkm,Skm,'LineWidth',1.5); hold on


    %% Exponential
    model.type = 'exp';

    [fit.exp{g},nll.exp(g)] = fminbnd( ...
        @(x) nll_surv(x,age,cens,model), ...
        0.001,100);

    lambda = fit.exp{g};

    Sfit = exp(-lambda*tfit);
    plot(tfit,Sfit,'LineWidth',1.5)

    meanlife.exp(g) = 1/lambda;


    %% Hypoexponential
    model.type = 'hypo';

    theta0 = [1 2];
    lb = [0.001 0.001];
    ub = [100 100];

    fit.hypo{g} = fmincon( ...
        @(x) nll_surv(x,age,cens,model), ...
        theta0,[],[],[],[],lb,ub);

    nll.hypo(g) = nll_surv( ...
        fit.hypo{g},age,cens,model);

    l1 = fit.hypo{g}(1);
    l2 = fit.hypo{g}(2);

    if abs(l1-l2) < 1e-6
        l = 0.5*(l1+l2);
        Sfit = exp(-l*tfit).*(1+l*tfit);
    else
        Sfit = (l2*exp(-l1*tfit) - ...
                l1*exp(-l2*tfit))/(l2-l1);
    end

    plot(tfit,Sfit,'LineWidth',1.5)

    meanlife.hypo(g) = 1/l1 + 1/l2;


    %% Lognormal
    model.type = 'lognormal';

    theta0 = [mean(log(age)),std(log(age))];
    lb = [-10 0.001];
    ub = [10 10];

    fit.lognormal{g} = fmincon( ...
        @(x) nll_surv(x,age,cens,model), ...
        theta0,[],[],[],[],lb,ub);

    nll.lognormal(g) = nll_surv( ...
        fit.lognormal{g},age,cens,model);

    mu  = fit.lognormal{g}(1);
    sig = fit.lognormal{g}(2);

    Sfit = 1-logncdf(tfit,mu,sig);
    plot(tfit,Sfit,'LineWidth',1.5)

    meanlife.lognormal(g) = exp(mu + sig^2/2);


    %% Plot formatting
    xlabel('Age (years)')
    ylabel('Survival probability')
    ylim([0 1])
    box on

    legend('KM','Exponential','Hypoexponential','Lognormal')

    if ngroup == 1
        title('All roots < 2 mm')
    else
        title(sprintf('Diameter group %d',g))
    end

    hold off
end

%% Helper Functions
function nll = nll_surv(theta_vec,age,cens,model)

    switch model.type
        case 'exp'
            lambda = theta_vec(1);
            f = lambda .* exp(-lambda.*age);
            S = exp(-lambda.*age);

        case 'hypo'
            lambda1 = theta_vec(1);
            lambda2 = theta_vec(2);

            % Handle nearly equal rates using Erlang-2 limit
            if abs(lambda1-lambda2) < 1e-6
                lambda = 0.5*(lambda1+lambda2);
                f = lambda.^2 .* age .* exp(-lambda.*age);
                S = exp(-lambda.*age) .* (1 + lambda.*age);
            else
                f = (lambda1.*lambda2./(lambda2-lambda1)) .* ...
                    (exp(-lambda1.*age) - exp(-lambda2.*age));
                S = (lambda2.*exp(-lambda1.*age) - ...
                     lambda1.*exp(-lambda2.*age)) ./ ...
                    (lambda2-lambda1);
            end

        case 'lognormal'
            mu = theta_vec(1);
            sigma = theta_vec(2);
            f = lognpdf(age,mu,sigma);
            S = 1 - logncdf(age,mu,sigma);
    end

    % Avoid log(0)
    f = max(f,realmin);
    S = max(S,realmin);

    % cens = true means right-censored
    LL = sum((~cens).*log(f) + cens.*log(S));

    nll = -LL;

end