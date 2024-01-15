% HP97 Shock recovery
% SSF: ---------------------------------------------------------------------------------------------
%   Z(t) = D1*X(t) + D2*X(t-1) + R*ε(t),      X(t) = latent States
%   X(t) = A*X(t-1)            + C*ε(t),      ε(t) ~ MN(0,I)
% --------------------------------------------------------------------------------------------------
clear; clc; tic;
% set plotting defaults
set(groot,'defaultLineLineWidth',2); set(groot,'defaultAxesFontSize',14)
set(groot,'defaultAxesXTickLabelRotationMode','manual')
set(groot,'defaultAxesFontName','Times New Roman')
addpath(genpath('../../functions'))
addpath(genpath('../../utility.Functions'))               % set path to db functions
% addpath(genpath('D:/matlab.tools/db.toolbox/db'))   % set path to db functions
% CALL: get_all_db_toolbox_function_calls.m from Directory of code to be shared

% Sample size and seed for random number generator in simulation
Ts = 1e5; rng(123);   % takes about 1 sec for 1e5, 10 secs.for 1e6, 90 secs. for 1e7. --> does not change correlations from sims much
PLOT_STATES     = 1;  % set to 1 to plot ε(t) states
ESTIMATE_HP_US  = 0;  % set to 1 to plot the comparison with HP-filter function

% --------------------------------------------------------------------------------------------------    
% PARAMETERS: standard deviation sqrt(lambda = 1600)
phi = 40;





% --------------------------------------------------------------------------------------------------    
% offset for the first k latent state variables before the shocks.
k = 0;
% DEFINE SSM INPUT MATRICES ------------------------------------------------------------------------
dim_Z = 1;       % rows Z(t)
dim_X = 3;       % rows X(t)
dim_R = 2;       % rows ε(t)
% --------------------------------------------------------------------------------------------------    
% Define D1
D1 = zeros(dim_Z,dim_X); 
D1(1,1) = 1;  D1(1,2) = phi; D1(1,3) = -2*phi;
% Define D2
D2 = zeros(dim_Z,dim_X); 
D2(1,3) = phi;
% Define R
R  = zeros(dim_Z,dim_R);
% --------------------------------------------------------------------------------------------------
% Define A
A = zeros(dim_X); 
A(3,2) = 1; 
% Define C
C = [eye(dim_R); zeros(1,2); ];

% --------------------------------------------------------------------------------------------------

% CALL TO THE KURZ_SSM FUNCTION --------------------------------------------------------------------
P = Kurz_steadystate_P(D1, D2, R, A, C);
ss = k+1:dim_X;
% make display names % row_names = make_table_names('ε',1:dim_R,'(t)');          
row_names = {'ε1(t)','ε2(t)','ε2(t-1)'} ; 

Pstar = array2table([ diag(P.tT(ss,ss)) diag(P.tt(ss,ss)) ], ...
        'VariableNames',{'P(t|T)','P(t|t)'}, 'RowNames', row_names);
% select what to print to screen
sep; print_table(Pstar(1:dim_R,:),4,1,0)

% SIMULATE DATA FROM THE MODEL --> compute 'theoretical' properites of states
[Zs, Xs, Us] = Kurz_simulate_SSF(D1, D2, R, A, C, Ts);

% --------------------------------------------------------------------------------------------------
% CALL TO FUNCTIONS FROM KURZ's GITHUB PAGE, MILDLY MODIFIED TO SIMPLIFY INPUT AND COMPARABILTY WITH 
% MY CODE ABOVE AND USE OF PINV IN AM SMOOTHER OTHERWISE NON-SINGULARITY ISSUES.
% --------------------------------------------------------------------------------------------------
% INITIALIZE FILTER 
% Note: errors will always be N(0,1), but latent states may need more careful initialization.
a00 = zeros(dim_X, 1); P00 = eye(dim_X);
% Filter
[~, Kurz_KF] = Kurz_Filter(Zs, D1, D2, R, A, C, a00, P00);
% Smoothers 
% --------------------------------------------------------------------------------------------------
% Modified de Jong (1988, 1989) and Kohn and Ansley (1989) smoother (Eq. (4.11) in Kurz (2018))
KFS_deJ = Kurz_DeJongKohnAnsley_Smoother(D1, D2, A, Kurz_KF); % Contains KF and KS output. NO INV, NO INITVALS FOR STATES
% --------------------------------------------------------------------------------------------------

%% PLOT THE KF/KS ESTIMATES OF THE STATES 
% --------------------------------------------------------------------------------------------------
PLOT_KF = 0; % set to 1 to use KF output, otherwise use KS
if PLOT_STATES
  clf; tiledlayout(4,2,TileSpacing = "compact", Padding = "compact");
  % plot_names = make_table_names('$\epsilon_{', 1:k, 't}$');
  % loop through plots
  if PLOT_KF
    state_t = KFS_deJ.att;
    state_name = 'Estimate:$\,\hat{X}_{t|t}$'; 
  else 
    state_t = KFS_deJ.atT;
    state_name = 'Estimate:$\,\hat{X}_{t|T}$'; 
  end

  for ii = k+1:dim_R
    nexttile([1 2])
    hold on;
      plot(Xs(:,ii), 'LineWidth',3);                              % 'true' simulated state X
      plot(state_t(:,ii),'--','Color',clr(3),'LineWidth',2.5);    % Filtered or smoothed estimate of state X
    hold off;
    hline(0)
    box on; grid on;
    set(gca,'GridLineStyle',':' ,'GridAlpha',1/3, 'LineWidth',5/5);
    add2yaxislabel;
    addlegend({'True State',state_name},1)
    % addsubtitle(row_names(ii),-1.115)
    addsubtitle(['$\varepsilon_{' num2str(ii) 't}$'],-1.115,16)
  end
  nexttile; % scatter(true, estimate)
  dims = [-5 5]; i = 1; 
  scatter(Xs(:,i),state_t(:,i),'x');
    ylim(dims); setxticklabels([dims(1):1:dims(2)]);
    line(dims, dims, 'Color', 'k', 'LineWidth', 1); hline(0)
    addsubtitle(['$\varepsilon_{' num2str(i) 't}$ (true)'],-1.115,16)
    ylabel(['$E_T\varepsilon_{' num2str(i) 't}$'],'Interpreter','latex')
    addgrid(4/5)
  nexttile; i = 2; 
  scatter(Xs(:,i),state_t(:,i),'x');
    ylim(dims); setxticklabels([dims(1):1:dims(2)]); 
    line(dims, dims, 'Color', 'k', 'LineWidth', 1); hline(0)
    addsubtitle(['$\varepsilon_{' num2str(i) 't}$ (true)'],-1.115,16)
    ylabel(['$E_T\varepsilon_{' num2str(i) 't}$'],'Interpreter','latex')
    addgrid(4/5)
  % print2pdf('HP97_plots_KS',2);
end
%% --------------------------------------------------------------------------------------------------

% CORRELATIONS:
% NOTE: 𝜙²yᶜ(t) = ε2(t), so the correlation between the true and estimated ε2(t) is equivalent to the correlation between the true and estimated HP output gap. 
corr_table = array2table( diag(corr(Xs(:,ss), KFS_deJ.atT(:,ss))), ...
               'RowNames', row_names, 'VariableNames', {'Corr(.)'});
% print correlations simulated and KS shocks
print_table(corr_table(1:dim_R,:),4,1,'Correlation between True X(t) and (estimated) Kalman Smoothed States X(t|T)');sep

% Correlation matrix from KS estimates, Truth is uncorrelated
corr_XtT = array2table( corr(KFS_deJ.atT), 'RowNames', row_names, 'VariableNames', row_names);
print_table(corr_XtT(1:dim_R,1:dim_R),4,1,'Correlation Matrix of (estimated) Kalman Smoothed States X̂(t|T)');sep
% Correlation matrix from KF estimates, Truth is uncorrelated
corr_Xtt = array2table( corr(KFS_deJ.att), 'RowNames', row_names, 'VariableNames', row_names);
print_table(corr_Xtt(1:dim_R,1:dim_R),4,1,'Correlation Matrix of (estimated) Kalman Filtered States X̂(t|t)',[],0);sep

% Define/Make eta(i) = eps(i) 
for jj = 1:dim_R
  eval(['ETe' num2str(jj) 't = KFS_deJ.atT(:,k+' num2str(jj) ');']);
  eval(['Ete' num2str(jj) 't = KFS_deJ.att(:,k+' num2str(jj) ');']);
end

% IDENTITIES: Run (dynamic) OLS regressions: ie., ∆ETη5t = 0.107∆ETη3t − 0.028ETη4t and the like
sep(133,'=',1); fprintf('Filter Identity on page 17 in EER(2022). Dependent variable: Etε2(t) \n')
Xnames_ID1 = {'Etε1(t)'};  % Xnames_ID1 = [];
ID1 = ols(Ete2t, [ Ete1t ], 1, Xnames_ID1);

sep(133,'=',1); fprintf('Smoother Identity on page 12, eq. 10 in Stars(2023). Dependent variable: Δ²ETε1(t) \n')
sep;fprintf('NOTE:    Identity should be: Δ²ETε1(t) = 1/40 ETε2(t-2) (and not ETε2(t) as stated)\n');sep
Xnames_ID2 = {'ETε2(t-2)'};  % Xnames_ID2 = [];
ID2 = ols( delta(ETe1t, 2) , [ lag(ETe2t, 2)], 1, Xnames_ID2);

%%
if ESTIMATE_HP_US
% -------------------------------------------------------------------------------------------------
% USE HP Filter routine to back out trend and cycle shocks and compare from US REAL GDP DATA
% --------------------------------------------------------------------------------------------------
load('US-GDP.mat')
y   = log(usGDP);
D2y = delta(y, 2);    % don't use diff due to y being a TimeTable, --> use delta instead
% make Z(t) variable for 'shock recovery' SSM form: ie. Z(t)=Δ²y(t)
ZZ  = D2y.("Δ²gdp");
HP  = hp_filter(y, phi^2); % call to 'standard' HP filter code (see 
dHP_trend = delta(HP.trend,1,1); 
[~, Kurz_KF_US] = Kurz_Filter(removenans(ZZ), D1, D2, R, A, C, a00, P00);
% construct the cycle from the shock recovery SSM
KS_deJ_US = Kurz_DeJongKohnAnsley_Smoother(D1, D2, A, Kurz_KF_US);
SSM.cycle = phi*addnans(KS_deJ_US.atT(:,2),2,0);

% Reconstruct HP.trend from KS_deJ_US.atT(:,1) using initVals from HP.trend (∆2y*(t) = ε1(t))
% cumsum([HP.trend(1); cumsum([ΔHP.trend(1); ETε1(t)])])
SSM.trend = cumsum([HP.trend(1); cumsum([dHP_trend(1); KS_deJ_US.atT(:,1)])]);
% Δ⁴HP-trend(t) =  1/φ² HP-cycle(t-2) 
sep(133,'=',1); fprintf('Filter Identity on page 17 in EER(2022) for US data from HP-Filter . Dependent variable: Δ⁴HP-trend(t) \n')
Xnames_ID3 = {'HP-Cycle(t-2)'};
ID3 = ols( delta(SSM.trend, 4) , [ lag(SSM.cycle, 2)], 1, Xnames_ID3);

% PLOT HP ON US REAL GDP DATA
  figure(2); clf;
  % tiledlayout(3,1)
  tiledlayout(4,1, TileSpacing = 'loose', Padding = 'compact');
  nexttile
  hold on;
    plot(HP.trend)
    plot(SSM.trend,'--')
    hline(0)
  hold off; 
  box on; addgrid;
  % ylim([7 11])
  setyticklabels(7:.5:10.5,1)
  setdateticks(HP.Date,25)
  addlegend({'HP-Filter','Shock-Recovery SSM'},1)
  addsubtitle('Trend in US-GDP data', -1.16)
  % PLOT cycles 
  nexttile
  hold on;
    plot(HP.cycle)
    plot(SSM.cycle,'--')
    hline(0)
  hold off; 
  box on; addgrid;
  setyticklabels(-.10:.02:.06)
  setdateticks(HP.Date,25)
  addlegend({'HP-Filter','Shock-Recovery SSM'})
  addsubtitle('Cycle in US-GDP data', -1.16)
end










toc































%EOF