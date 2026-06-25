%% Figure 4 - C:N ratio visualizations
% Generates Figure 4: mean C:N ratios plotted against root order and
% against time for order 5
%
% Required File(s): 'face_dataset.csv'
%
% Required toolboxes: None
%
clear; clc;

%% Load data
dtab = readtable('face_dataset.csv'); % read data table
tshft = 18; % time offset in days
t0 = min(dtab.DateSampled) - days(tshft);

dat.t = days(dtab.DateSampled - t0)/365;
dat.Cind0 = strcmp(dtab.CO2Treatment,'Ambient');
dat.Cind1 = strcmp(dtab.CO2Treatment,'Elevated');

dat.Nind0 = strcmp(dtab.NTreatment,'Unfertilized');
dat.Nind1 = strcmp(dtab.NTreatment,'Fertilized');

%% Settings
base_ind = dat.Cind1; % elevated CO2 only

orders = 1:5;
Ninds = {dat.Nind0, dat.Nind1};
treat_names = {'Unfertilized','Fertilized'};

cols = lines(2);
ms = 9;
lw = 2;

%% Compute C:N summaries by order and treatment
CN_mean = NaN(2,length(orders));
CN_sd   = NaN(2,length(orders));

for q = 1:2
    for j = 1:length(orders)
        ind = base_ind & Ninds{q} & dtab.Order == orders(j);
        x = dtab.C_N(ind);

        CN_mean(q,j) = mean(x,'omitnan');
        CN_sd(q,j)   = std(x,'omitnan');
    end
end

%% Plotting
figure('Color','w','Position',[100 100 1400 500]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

% Left panel - Mean C:N against root order
axA = nexttile;
hold(axA,'on');

for q = 1:2
    errorbar(orders,CN_mean(q,:),CN_sd(q,:),'o-',...
        'Color',cols(q,:),...
        'MarkerFaceColor','w',...
        'MarkerEdgeColor',cols(q,:),...
        'MarkerSize',ms,...
        'LineWidth',lw,...
        'CapSize',6);
end

xlabel(axA,'Root Order');
ylabel(axA,'C:N Ratio');
xlim([0.5 5.5])
set(axA,'FontSize',22,'XTick',orders);
legend(axA,treat_names,'Location','best');
box(axA,'on');

% Right panel - C:N of order 5 roots through time
axB = nexttile;
hold(axB,'on');

ord_show = 5;

for q = 1:2
    ind = base_ind & Ninds{q} & dtab.Order == ord_show;

    % Raw points
    scatter(dat.t(ind),dtab.C_N(ind),55,...
        'MarkerFaceColor',cols(q,:),...
        'MarkerEdgeColor','k',...
        'MarkerFaceAlpha',0.45);

    % Time-specific means +/- SD
    tvals = unique(dat.t(ind));
    ym = NaN(size(tvals));
    ys = NaN(size(tvals));

    for k = 1:length(tvals)
        ii = ind & dat.t == tvals(k);
        x = dtab.C_N(ii);

        ym(k) = mean(x,'omitnan');
        ys(k) = std(x,'omitnan');
    end

    errorbar(tvals,ym,ys,'o-',...
        'Color',cols(q,:),...
        'MarkerFaceColor','w',...
        'MarkerEdgeColor',cols(q,:),...
        'MarkerSize',ms,...
        'LineWidth',lw,...
        'CapSize',6);
end

xlabel(axB,'t (years)');
ylabel(axB,'C:N Ratio');
set(axB,'FontSize',22);
box(axB,'on');
