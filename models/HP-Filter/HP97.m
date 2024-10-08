% HP97 Shock recovery SSF
% --------------------------------------------------------------------------------------------------
%   Observed: Z(t) = D1*X(t) + D2*X(t-1) + R*ε(t),    X(t) = latent States
%   State:    X(t) =  A*X(t-1)           + Q*ε(t),    ε(t) ~ MN(0,I)
% --------------------------------------------------------------------------------------------------
clear; clc; tic;
% set plotting defaults
set(groot,'defaultLineLineWidth',2); set(groot,'defaultAxesFontSize',14)
set(groot,'defaultAxesXTickLabelRotationMode','manual')
set(groot,'defaultAxesFontName','Times New Roman')
addpath('../../functions', '../../utility.Functions')         % addpath to functions used
% addpath(genpath('D:/matlab.tools/db.toolbox/db')) % set path to db functions folder (including all subfolders)
% CALL: get_all_db_toolbox_function_calls.m from Directory of code to be shared

% Sample size and seed for random number generator in simulation
Ts = 5e2; rng(10);    % takes about 1 sec for 1e5, 10 secs. for 1e6, 90 secs. for 1e7. --> does not change correlations from sims much
PLOT_STATES     = 1;  % set to 1 to plot ε(t) states
ESTIMATE_HP_US  = 0;  % set to 1 to plot the comparison with HP-filter function
PLOTS2PDF       = 0;  % set to 1 to print plots to PDF.

% --------------------------------------------------------------------------------------------------    
% PARAMETERS: standard deviation sqrt(lambda = 1600)
% psi = sqrt(100000)        % very large   
psi = 40;                   % sqrt(1600)
% psi = 0.6375;             % gives approximate 50% R2 for each shock
% psi = 0.1;                % very low Noise 2 signal
N2S = psi^2;                % Noise-to-Signal ratio (NOTE: Harvey uses q = Signal-to-Noise ratio)
Lambda = psi^2;             % HP(Lambda)

% DEFINE SSF INPUT MATRICES ------------------------------------------------------------------------
dim_Z = 1;                  % rows Z(t)
dim_X = 3;                  % rows X(t)
dim_R = 2;                  % rows ε(t)
% offset for the first k latent state variables before the shocks.
k = 0;

% --------------------------------------------------------------------------------------------------    
% Define D1
D1 = zeros(dim_Z,dim_X); 
D1(1,1) = 1;  D1(1,2) = psi; D1(1,3) = -2*psi;

% Define D2
D2 = zeros(dim_Z,dim_X); 
D2(1,3) = psi;
% Define R
R  = zeros(dim_Z,dim_R);
% --------------------------------------------------------------------------------------------------
% Define A
A = zeros(dim_X); 
A(3,2) = 1; 
% Define Q
Q = [eye(dim_R); zeros(1,2); ];
% --------------------------------------------------------------------------------------------------

% CALL TO FUNCTIONS FROM KURZ's GITHUB PAGE, MILDLY MODIFIED TO SIMPLIFY INPUT AND COMPARABILTY WITH 
% MY CODE ABOVE AND USE OF PINV IN AM SMOOTHER OTHERWISE NON-SINGULARITY ISSUES.
% --------------------------------------------------------------------------------------------------
% SIMULATE DATA FROM THE MODEL --> compute 'theoretical' properites of states
[Zs, Xs, Us] = Kurz_simulate_SSF(D1, D2, R, A, Q, Ts);
% NOW JUST DO ONE CALL TO Kurz_FilterSmoother RATHER THAN MULTIPLE CALLS
% Note: errors will always be N(0,1), but latent states may need more careful initialization.
a00 = zeros(dim_X, 1); P00 = eye(dim_X);
[~,KFS] = Kurz_FilterSmoother(Zs, D1, D2, R, A, Q, a00, P00);
% % USE THIS IF YOU DON'T HAVE DATA TO PUT INTO THE FULL Kurz_FilterSmoother() FUNCTION
% % Pstar0 = Kurz_Pstar(D1, D2, R, A, Q);
% % [~, Kurz_KF] = Kurz_Filter(Zs, D1, D2, R, A, Q, a00, P00);
% % KFS = Kurz_Smoother(D1, D2, R, A, Q, Kurz_KF); % Contains KF and KS output. 

Neps = k+1:dim_X-1;    % shock index in States X(t)
row_names = make_table_names('ε',1:dim_R,'(t)');              % make display names 

Pstar = array2table([ diag(KFS.Pstar.tT(Neps,Neps)) diag(KFS.Pstar.tt(Neps,Neps)) ], 'VariableNames',{'P*(t|T)','P*(t|t)'}, 'RowNames', row_names);
% select what to print to screen
sep; print_table(Pstar,4,1,0)

% CORRELATIONS:
% --------------------------------------------------------------------------------------------------
% Correlation between the true and estimated states. 
% R2 of Plagborg-Møller and Wolf (2022) 
R2 = [];
for jj = Neps
  pwR2.( ['e' num2str(jj)]) = ols(KFS.atT(:,jj),Xs(:,jj),1,[],[],[],0); % set last 0 to 1 to print to screen
  R2(jj-k,:) = eval((['pwR2.e' num2str(jj) '.R2']));
end

% compute the theoretical correlations implied by formula (10)
STDs  = [ones(dim_R,1)]; % theoretical/model stdevs.
rho_theory = corr_theory(STDs, std(KFS.atT(:,Neps))', Pstar.('P*(t|T)'));

corr_table = array2table( [ diag(corr(Xs(:,Neps),KFS.atT(:,Neps)))  rho_theory  R2], 'RowNames', row_names, 'VariableNames', {'ρ(Sim)','ρ(Theory)','R²(Sim)'});
% print correlations simulated and KS shocks
print_table(corr_table(1:dim_R,:),4,1,'Correlation between True X(t) and (estimated) Kalman Smoothed States ETX(t)',[],0);

% Correlation matrix from KS estimates, Truth is uncorrelated
corr_XtT = array2table( corr(KFS.atT(:,Neps)), 'RowNames', row_names, 'VariableNames', row_names);
print_table(corr_XtT(1:dim_R,1:dim_R),4,1,'Correlation Matrix of (estimated) Kalman Smoothed States ETX(t)',[],0);

% Correlation matrix from KF estimates, Truth is uncorrelated
corr_Xtt = array2table( corr(KFS.att(:,Neps)), 'RowNames', row_names, 'VariableNames', row_names);
print_table(corr_Xtt(1:dim_R,1:dim_R),4,1,'Correlation Matrix of (estimated) Kalman Filtered States EtX(t)',[],0); sep

% DISPLAY RECOVEY DIAGNOSTICS ALL IN ONE MATRIX TO PRINT TO LATEX
matRowNames = { 'P*(t|T)  ';'R²(Sim)  ';'ρ(Theory)'};
fprintf('Recovery Measures (Order is)\n'); sep
mat2latex([Pstar.("P*(t|T)")'; corr_table.("R²(Sim)")'; corr_table.("ρ(Theory)")'], 4, matRowNames); sep
toc

%% PLOT THE KF/KS ESTIMATES OF THE STATES 
% --------------------------------------------------------------------------------------------------
% make some plotting variables
WHICH_STATE_2_PLOT  = 0;                % set to 1 to use KF output, otherwise use KS
state_t = KFS.atT;
shock_names = {'Trend';'Cycle'};
if WHICH_STATE_2_PLOT; state_t = KFS.att; end
xgrd = linspace(-5,5,100)';             % make xgrd for plotting
STL = -1.21;                  % subtitle location
OST = -.00;
dims = [-6:2:6]; FNS = 12; XOS = 11;
PtT1 = sqrt(squeeze(KFS.PtT(1,1,:))); % needed for CI computation

if PLOT_STATES 
  clf; TL = tiledlayout(7,2); if ~verLessThan('matlab','9.9');TL.TileSpacing='compact';TL.Padding='loose'; end
  % loop through plots
  for ii = k+1:dim_R
    nexttile
    hold on;
      plot(Xs(:,ii), 'LineWidth',3);                              % 'true' simulated state X
      plot(state_t(:,ii),'--','Color',clr(3),'LineWidth',2.5);    % Filtered or smoothed estimate of state X
    hold off; 
    xlim([-XOS length(Xs(:,ii))+XOS]); 
    setyticklabels(dims,0);
    addgrid(5/5); hline(0); 
    addlegend({'True','Estimate'},1,FNS)
    addsubtitle(['$\varepsilon_{' num2str(ii-k) 't}$'],STL,FNS)
    add2yaxislabel(1)

    nexttile
    hold on; 
      scatter(Xs(:,ii),state_t(:,ii),'x');
      plot(xgrd, xgrd*pwR2.( ['e' num2str(ii)]).bhat)
    hold off; 
    setxticklabels([dims(1):1:dims(end)]);
    setyticklabels(dims,0);
    addgrid(5/5); hline(0); 
    line(dims, dims, 'Color', 'k', 'LineWidth', 1); 
    ylabel(['$E_T\varepsilon_{' num2str(ii-k) 't}$ (Estimate)'],'Interpreter','latex','FontSize',FNS)
    addsubtitle(['$\varepsilon_{' num2str(ii-k) 't}$ (True ' char(shock_names(ii)) ' Shock)'],STL,FNS)
    addlegend({['$R^2=' num2str(pwR2.( ['e' num2str(ii)]).R2,'%2.4f') '$']},1,FNS)
    add2yaxislabel(1)
  end

  % plot Observed and latent state
  nexttile(5,[1 2])
  hold on;
    plot(Zs(:,1) , 'Color',clr(3), 'LineWidth',3);                              % Observed series
    plot(Xs(:,1) , 'Color',clr(1), 'LineWidth',2);                              % Latent state X
    % plot(KFS.atT(:,1),'--','Color',clr(3),'LineWidth',2.0);    % Filtered or smoothed estimate of state X
    % plot(KFS.att(:,1),'--','Color',clr(4),'LineWidth',1.0);    % Filtered or smoothed estimate of state X
  hold off; 
  xlim([-XOS length(Xs(:,ii))+XOS]); 
  % setyticklabels(dims,0,FNS);
  addgrid(5/5); hline(0); 
  addlegend({ 'Observed series: $y_t$', ...
              ['True state: $\mu_t$ (Noise-to-signal = ' num2str(N2S^2,4) ')'], ...
            } ...
            ,1,FNS-1)
  addsubtitle(['Observed vs Latent' ],STL+OST,FNS)
  add2yaxislabel(1)
  % plot states with CIs and estimates
  nexttile(7,[1 2]); LH = [];
  hold on;
    addCIs( [KFS.atT(:,1) - 2*PtT1, KFS.atT(:,1) + 2*PtT1 ], 'r' );
    plot(Xs(:,1), 'LineWidth', 2);                              % 'true' simulated state X
    plot(KFS.atT(:,1),'--','Color',clr(2),'LineWidth', 2);    % Filtered or smoothed estimate of state X
  hold off; 
  xlim([-XOS length(Xs(:,ii))+XOS]); 
  % setyticklabels(dims,0,FNS);
  addgrid(5/5); hline(0); 
  legend({'$95\%$ CI' , ...
          'True state $\mu_t$', ...
          'KS Estimate'}, ...
          'Interpreter','latex','FontSize', FNS-1, 'Location','northwest')
  addsubtitle(['True and KS estimate of $\mu$'],STL+OST,FNS)
  add2yaxislabel(1)

  % TO PRINT TO PDF
  if PLOTS2PDF   
  % fig_name = ['LL_plots_KS_Noise2Signal_', num2str(N2S^2,4) '_T_' num2str(Ts,'%d') '.pdf'];
  fig_name = 'HP97_plots_KS.pdf';
  % print2pdf(fig_name); % super slow here
    exportgraphics(gcf, fig_name, 'ContentType','vector')
  end
end
% --------------------------------------------------------------------------------------------------

% IDENTITIES REGRESSIONS:
% --------------------------------------------------------------------------------------------------
% define/make: ETεi(t) or ETηi(t) as needed
for jj = 1:dim_R
  eval(['ETe' num2str(jj) 't = KFS.atT(:,k+' num2str(jj) ');']);
  eval(['Ete' num2str(jj) 't = KFS.att(:,k+' num2str(jj) ');']);
end

sep(133,'=',1); fprintf('Filter Identity on page 17 in EER(2022). Dependent variable: Etε2(t) \n')
Xnames_ID1 = {'Etε1(t)'};  % Xnames_ID1 = [];
ID1 = ols(Ete2t, [ Ete1t ], 1, Xnames_ID1);

sep(133,'=',1); fprintf('Smoother Identity on page 12, eq. 10 in Stars(2023). Dependent variable: Δ²ETε1(t) \n')
sep;fprintf('NOTE:    Identity should be: Δ²ETε1(t) = 1/40 ETε2(t-2) (and not ETε2(t) as stated)\n');sep
Xnames_ID2 = {'ETε2(t-2)'};  % Xnames_ID2 = [];
ID2 = ols( delta(ETe1t, 2) , [ lag(ETe2t, 2)], 1, Xnames_ID2);


%
if ESTIMATE_HP_US
% -------------------------------------------------------------------------------------------------
% USE HP Filter routine to back out trend and cycle shocks and compare from US REAL GDP DATA
% --------------------------------------------------------------------------------------------------
load('US-GDP.mat')
y = usGDP;
y.gdp = log(usGDP.gdp);
D2y   = delta(y, 2);    % don't use diff due to y being a TimeTable, --> use delta instead
% make Z(t) variable for 'shock recovery' SSF form: ie. Z(t)=Δ²y(t)
ZZ  = D2y.("Δ²gdp");
HP  = hp_filter(y, psi^2); % call to 'standard' HP filter code
dHP_trend = delta(HP.trend,1,1); 
[~, Kurz_KF_US] = Kurz_Filter(removenans(ZZ), D1, D2, R, A, Q, a00, P00);
% construct the cycle from the shock recovery SSM
KFS_deJ_US = Kurz_DeJongKohnAnsley_Smoother(D1, D2, A, Kurz_KF_US);
SSM.cycle  = psi*addnans(KFS_deJ_US.atT(:,2),2,0);

% Reconstruct HP.trend from KS_deJ_US.atT(:,1) using initVals from HP.trend (∆2y*(t) = ε1(t))
% cumsum([HP.trend(1); cumsum([ΔHP.trend(1); ETε1(t)])])
SSM.trend = cumsum([HP.trend(1); cumsum([dHP_trend(1); KFS_deJ_US.atT(:,1)])]);
% Δ⁴HP-trend(t) =  1/φ² HP-cycle(t-2) 
sep(133,'=',1); fprintf('Filter Identity on page 17 in EER(2022) for US data from HP-Filter . Dependent variable: Δ⁴HP-trend(t) \n')
Xnames_ID3 = {'HP-Cycle(t-2)'};
ID3 = ols( delta(SSM.trend, 4) , [ lag(SSM.cycle, 2)], 1, Xnames_ID3);

% PLOT HP ON US REAL GDP DATA
  figure(2); clf;
  TL = tiledlayout(4,1); if ~verLessThan('matlab','9.9');TL.TileSpacing='loose';TL.Padding='loose'; end
  nexttile
  hold on;
    plot(HP.trend)
    plot(SSM.trend,'--')
    hline(0)
  hold off; addgrid;
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
  hold off; addgrid;
  setyticklabels(-.10:.02:.06)
  setdateticks(HP.Date,25)
  addlegend({'HP-Filter','Shock-Recovery SSM'})
  addsubtitle('Cycle in US-GDP data', -1.16)
end









% toc































%EOF