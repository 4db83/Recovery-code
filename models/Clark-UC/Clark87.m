% Clark87 Shock recovery
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
ESTIMATE_CLARK  = 0;  % set to 1 to estimate the parameters given below from US GDP Data

% --------------------------------------------------------------------------------------------------    
% PARAMETERS: taken from the model fitted to US-GDP data from 1947:Q1 to 2019:Q4 (which can be estimated with the code below)
a1  =  1.51023433332729;   % a1   
a2  = -0.56787952465929;   % a2   
% standard deviations                                                                             
s1  =  0.54396737899273;   % sigma(ystar)                                                                     
s2  =  0.02093523340402;   % sigma(g)                                                                        
s3  =  0.59796738263493;   % sigma(ytilde)                                                                         
% --------------------------------------------------------------------------------------------------    
% offset for the first k latent state variables before the shocks.
k = 0;
% DEFINE SSM INPUT MATRICES ------------------------------------------------------------------------
dim_Z = 1;       % rows Z(t)
dim_X = 8;       % rows X(t)
dim_R = 3;       % rows ε(t)
% --------------------------------------------------------------------------------------------------    
% Define D1
D1 = zeros(dim_Z,dim_X); 
D1(1,[4:5 8]) = [s1 -a1*s1 s3];
% Define D2
D2 = zeros(dim_Z,dim_X); 
D2(1,[2 5:8]) = [s2 -a2*s1 -a1*s2 -a2*s2 -s3];
% Define R
R  = zeros(dim_Z,dim_R);
% --------------------------------------------------------------------------------------------------
% Define A
A = zeros(dim_X);   
A(4,1) = -1; A(5,4) = 1; A(6,2) = 1; A(7,6) = 1; A(8,3) = -1;
% Define C
C = [eye(dim_R); zeros(dim_X-dim_R,dim_R)];
C(4,1) = 1; C(8,3) = 1;
% --------------------------------------------------------------------------------------------------

%% CALL TO THE KURZ_SSM FUNCTION --------------------------------------------------------------------
P = Kurz_steadystate_P(D1, D2, R, A, C);
ss = k+1:dim_R;
% make display names % row_names = make_table_names('ε',1:dim_R,'(t)');          
row_names = {'ε1(t)','ε2(t)','ε3(t)'} ; 

Pstar = array2table([ diag(P.tT(ss,ss)) diag(P.tt(ss,ss)) ], ...
        'VariableNames',{'P(t|T)','P(t|t)'}, 'RowNames', row_names);
% select what to print to screen
sep; print_table(Pstar,4,1,0)

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

% PLOT THE KF/KS ESTIMATES OF THE STATES: % STATE vector = [trend(t); g(t); cycle(t); cycle(t-1)]; 
% --------------------------------------------------------------------------------------------------
PLOT_KF = 0; % set to 1 to use KF output, otherwise use KS
if PLOT_STATES
  clf; tiledlayout(5,1,TileSpacing = "compact", Padding = "compact");
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
    nexttile
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
    addsubtitle(['$\varepsilon_{' num2str(ii) 't}$'],-1.17,16)
  end
  % print2pdf('Clark87_plots_KS',2);
end
% --------------------------------------------------------------------------------------------------

% CORRELATIONS: 
corr_table = array2table( diag(corr(Xs(:,ss), KFS_deJ.atT(:,ss))), ...
               'RowNames', row_names, 'VariableNames', {'Corr(.)'});
% print correlations simulated and KS shocks
print_table(corr_table(1:dim_R,:),4,1,'Correlation between True X(t) and (estimated) Kalman Smoothed States X(t|T)');sep

% Correlation matrix from KS estimates, Truth is uncorrelated
corr_XtT = array2table( corr(KFS_deJ.atT(:,ss)), 'RowNames', row_names, 'VariableNames', row_names);
print_table(corr_XtT,4,1,'Correlation Matrix of (estimated) Kalman Smoothed States X̂(t|T)');sep
% Correlation matrix from KF estimates, Truth is uncorrelated
corr_Xtt = array2table( corr(KFS_deJ.att(:,ss)), 'RowNames', row_names, 'VariableNames', row_names);
print_table(corr_Xtt,4,1,'Correlation Matrix of (estimated) Kalman Filtered States X̂(t|t)',[],0);sep
fprintf('NOTE: Filtered trend growth shocks Etε2(t) == 0 for all t, hence we get a NaN in correlations\n')
% Define/Make eta(i) = eps(i)
for jj = 1:dim_R
  eval(['ETe' num2str(jj) ' = KFS_deJ.atT(:,k+' num2str(jj) ');']);
  eval(['Ete' num2str(jj) ' = KFS_deJ.att(:,k+' num2str(jj) ');']);
end

%% IDENTITIES: Run (dynamic) OLS regressions: ie., ∆ETη5t = 0.107∆ETη3t − 0.028ETη4t and the like
clc
fprintf('\n');sep(133,'=');
fprintf('Filter Identity similar to HP Filter. Dependent variable: Etε1(t)\n')
% sep;fprintf('NOTE:    Identity should be: Δ²ETε1(t) = 1/40 ETε2(t-2) (and not ETε2(t) as stated)\n');sep
Xnames_ID1 = {'Etε3(t)'};  
ID1 = ols( Ete1,  [ Ete3  ], 1, Xnames_ID1);

sep(133,'=',1); fprintf('Filter Identity similar to HP Filter. Dependent variable: ∆ETε3(t)\n')
Xnames_ID2 = {'∆ETε1(t-1)','ETε3(t-1)','ETε3(t-2)'};
ID2 = ols( delta(ETe3,1),  [ mlag(delta(ETe1),1)  mlag(ETe3,2) ], 1, Xnames_ID2);

sep(133,'=',1); fprintf('Filter Identity similar to HP Filter. Dependent variable: ∆ETε3(t)\n')
Xnames_ID3 = {'∆ETε2(t-1)','∆ETε2(t-1)','ETε3(t-1)','ETε3(t-2)'};
ID3 = ols( delta(ETe3,1),  [ mlag(delta(ETe2),2) mlag(ETe3,2)], 1, Xnames_ID3);

%
% sep(133,'=',1); fprintf('Filter Identity similar to HP Filter. Dependent variable: ∆Etε3(t)\n')
% Xnames_ID4 = {'∆Etε1(t-1)','Etε3(t-1)','Etε3(t-2)'};
% Xnames_ID4 = []
% ID4 = ols( delta(ETe1,1),  [ mlag(ETe3,3) mlag(ETe2,3) ], 1, Xnames_ID4);


%% -------------------------------------------------------------------------------------------------
if ESTIMATE_CLARK
% NOW ESTIMATE THE DOUBLE DRIFT UNOBSERVED COMPONENT (UC) MODEL OF CLARK US REAL GDP DATA
% --------------------------------------------------------------------------------------------------
% STATE vector alpha = [trend(t); g(t); cycle(t); cycle(t-1)];
% --------------------------------------------------------------------------------------------------
load('US-GDP.mat')
Y  = log(usGDP);
% set sample period
TT = timerange('Q1-1947','Q4-2019','closed');
y  = 100.*Y(TT,:);
% HP Filter on Y to get initial estiamtes
HP = hp_filter(y,1600);		% Lambda = 36000 is what is used by them in the paper
% ols on HP_cycle to get intial estimates AR(2) parameter estimates
AR2 = estimate_armax(HP.cycle,0,[1 2],[]);% print_arma_results(AR2);

% INITIVAL VALUES
sig_Dg		= nanstd(delta(HP.trend));
sig_trend = 1;
sig_cycle = AR2.SE_reg;
% INITIVAL VALUE AND PRIORS FOR UC0 
Clark_initVals	= [-AR2.aL(2:3)'; sig_trend; sig_Dg; sig_cycle];
Clark_prior.a00 = [HP.trend(1); diff(HP.trend(1:2)); 0; 0];
Clark_prior.P00 = 1e6;
% RHO CORRELATION PARAMETER
Clark_prior.rho = 0;
% DATA TO BE PARSED IN;
Y_in = y.gdp';
LL_Clark_initVals = LogLike_Clark_UC0(Clark_initVals, Y_in, Clark_prior,1);

% ESTIMATE USING UNCORRELATED ERRORS MODEL WITH FMINUNC
% opts_fminunc = optimoptions(@fminunc,'FunctionTolerance',1e-8,'FiniteDifferenceType', 'central');
opts_fminunc = optimoptions(@fminunc,'Display','none');
[theta_Clark, LL_Clark,~,~,~,H_Clark] = fminunc(@LogLike_Clark_UC0, Clark_initVals, opts_fminunc, Y_in, Clark_prior);
% EXTRACT NOW THE SMOOTHED SERIES OF TREND GROWTH
[~, KFS_Clark] = LogLike_Clark_UC0(theta_Clark, Y_in, Clark_prior, 1);
sep; fprintf(' Clark (1987) Model parameter estimates \n')
model_names = { 'Fminunc'; 'Stderr'; 'InitVals' };
par_names	  = {'AR(1)'; 'AR(2)'; 'sigma_y*'; 'sigma_g'; 'sigma_y~'; 'Log-Like'};
par_est	    = [ [theta_Clark sqrt(diag(inv(H_Clark)))  Clark_initVals ];
                [-LL_Clark    nan                      -LL_Clark_initVals ] ];
% print to screen
print2screen( par_est, par_names, model_names, '%17.4f'); %, '/_output/Table5.Clarks.UC.xls');
sep

% make Z(t) variable for 'shock recovery' SSM form: ie. Z(t)=a(L)Δ²y(t)
DY = delta(y.gdp,2);
ZZ = DY - a1*lag(DY) -a2*lag(DY,2);
[~, Kurz_KF_US] = Kurz_Filter(removenans(ZZ), D1, D2, R, A, C, a00, P00);
% Smoothers 
% --------------------------------------------------------------------------------------------------
% Modified de Jong (1988, 1989) and Kohn and Ansley (1989) smoother (Eq. (4.11) in Kurz (2018))
KS_deJ_US  = Kurz_DeJongKohnAnsley_Smoother(D1, D2, A, Kurz_KF_US);

% PLOT TREND GROWTH ETC FOR US GDP DATA
% Shock-recovery SSM
FILTERED = 0;
  eps_1 =  s1*KS_deJ_US.atT(:,1); eps_2 =  s2*KS_deJ_US.atT(:,2); eps_3 =  s3*KS_deJ_US.atT(:,3);
  ystr  = KFS_Clark.KFS.atT(:,1); g     = KFS_Clark.KFS.atT(:,2); ytld  = KFS_Clark.KFS.atT(:,3);
  Ktype = 'Smoothed~ ';
% filtered comparison
if FILTERED
  eps_1 =  s1*KS_deJ_US.att(:,1); eps_2 =  s2*KS_deJ_US.att(:,2); eps_3 =  s3*KS_deJ_US.att(:,3);
  ystr  = KFS_Clark.KFS.att(:,1); g     = KFS_Clark.KFS.att(:,2); ytld  = KFS_Clark.KFS.att(:,3);
  Ktype = 'Filtered~ ';
end
% now compute the shocks from Clark's SSM to be compatible with the shock-recovery SSM's shocks
eps_ystr  = delta(ystr) - lag(g);
eps_g     = delta(g);
eps_ytld  = ytld - a1*lag(ytld) - a2*lag(ytld,2);

% plot the different shocks
figure(2);clf; ST = -1.18;
tiledlayout(4,1, Padding="compact");
nexttile
plot([eps_ystr  [nan(4,1); eps_1]]); hline(0); setyticklabels(-2:0.5:2,1)
setdateticks(HP.Date,20)
addsubtitle(strcat(Ktype, '$\eta_1 = \sigma_1\varepsilon_{1}$'),ST)
addgrid(3/4)
nexttile
plot([eps_g     [nan(4,1); eps_2]]); hline(0); setyticklabels([-0.5:.1:0.5]/1,1)
setdateticks(HP.Date,20)
addlegend({'Clark SSM','Shock Recovery SSM'},3)
addsubtitle(strcat(Ktype, '$\eta_2 = \sigma_2\varepsilon_{2}$'),ST)
addgrid(3/4)
nexttile
plot([eps_ytld  [nan(4,1); eps_3]]); hline(0); setyticklabels(-2:.5:2,1)
setdateticks(HP.Date,20)
addsubtitle(strcat(Ktype, '$\eta_3 = \sigma_3\varepsilon_{3}$'),ST)
addgrid(3/4)
% print2pdf('Clark_SSM_Filtered')
% if SMOOTHED; print2pdf('Clark_SSM_Smoothed',0); else print2pdf('Clark_SSM_Filtered',0); end
% Smoothed

% plot the states from Clarks model
figure(3);clf; ST = -1.18;
tiledlayout(4,1, Padding="compact");
nexttile
plot([KFS_Clark.KFS.att(:,1)  KFS_Clark.KFS.atT(:,1)]); % hline(0); 
setdateticks(HP.Date,20)
ylim([700 1100])
addsubtitle('Filtered and Smoothed estimates of $y_t^\ast$',ST)
addgrid(3/4)
nexttile
plot([4*KFS_Clark.KFS.att(:,2)  4*KFS_Clark.KFS.atT(:,2)]); hline(0); 
setdateticks(HP.Date,20)
ylim([-2 8])
addsubtitle('Filtered and Smoohted estimates of $g_t$',ST)
addlegend({'Filtered','Smoothed'},3)
addgrid(3/4)
nexttile
plot([KFS_Clark.KFS.att(:,3)  KFS_Clark.KFS.atT(:,3)]); hline(0); 
setdateticks(HP.Date,20)
ylim([-6 6])
addsubtitle('Filtered and Smoothed estimates of $\tilde{y}_t$',ST)
addgrid(3/4)
% print2pdf('Clark_SSM',0)


end














































%EOF