% LW03 Shock recovery (https://www.newyorkfed.org/research/policy/rstar)
% SSF: ---------------------------------------------------------------------------------------------
%   Z(t) = D1*X(t) + D2*X(t-1) + R*ε(t),      X(t) = latent States
%   X(t) = A*X(t-1)            + C*ε(t),      ε(t) ~ MN(0,I)
% --------------------------------------------------------------------------------------------------
clear; clc;
% set plotting defaults
set(groot,'defaultLineLineWidth',2); set(groot,'defaultAxesFontSize',14)
set(groot,'defaultAxesXTickLabelRotationMode','manual')
set(groot,'defaultAxesFontName','Times New Roman')
addpath(genpath('./functions'))
addpath(genpath('./utility.Functions'))               % set path to db functions
% addpath(genpath('D:/matlab.tools/db.toolbox/db'))   % set path to db functions
% CALL: get_all_db_toolbox_function_calls.m from Directory of code to be shared

% Sample size and seed for random number generator in simulation
Ts = 1e4; rng(123);
PLOT_STATES = 0;

% DEFINE SSM INPUT MATRICES ------------------------------------------------------------------------
dim_Z = 1;       % rows Z(t)
dim_X = 3;       % rows X(t)
dim_R = 2;       % rows ε(t)
% offset for the first k latent state variables before the shocks.
k = dim_X - dim_R; 
k = 0;
% --------------------------------------------------------------------------------------------------    
% standard deviation sqrt(lambda = 1600)
phi = 40;
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
[PtT, Ptt] = Kurz_steadystate_P(D1, D2, R, A, C);
ss = k+1:dim_X;
% row_names = make_table_names('ε',1:dim_R,'(t)');          % make display names 
row_names = {'ε1(t)','ε2(t)','ε2(t-1)'} ; 

Pstar = array2table([ diag(PtT(ss,ss)) diag(Ptt(ss,ss)) ], ...
        'VariableNames',{'P(t|T)','P(t|t)'}, 'RowNames', row_names);
sep; print_table(Pstar,4,1,0)

% SIMULATE DATA FROM THE MODEL --> compute 'theoretical' properites of states
[Zs, Xs, Us] = Kurz_simulate_SSF(D1, D2, R, A, C, dim_Z, dim_X, dim_R, Ts);
Z = Zs; 

% --------------------------------------------------------------------------------------------------
% CALL TO FUNCTIONS FROM KURZ's GITHUB PAGE, MILDLY MODIFIED TO SIMPLIFY INPUT AND COMPARABILTY WITH 
% MY CODE ABOVE AND USE OF PINV IN AM SMOOTHER OTHERWISE NON-SINGULARITY ISSUES.
% --------------------------------------------------------------------------------------------------
% INITIALIZE FILTER 
% Note: errors will always be N(0,1), but latent states may need more careful initialization.
a00 = zeros(dim_X, 1); P00 = eye(dim_X);
% Filter
[~, Kurz_KF] = Kurz_Filter(removenans(Z), D1, D2, R, A, C, a00, P00);
% Smoothers 
% --------------------------------------------------------------------------------------------------
% Modified de Jong (1988, 1989) and Kohn and Ansley (1989) smoother (Eq. (4.11) in Kurz (2018))
KS_deJ  = Kurz_DeJongKohnAnsley_Smoother(D1, D2, A, Kurz_KF); % NO INV, NO INITVALS FOR STATES
% KS_deJ  = Kurz_AndersonMoore_Smoother(   D1, D2, A, Kurz_KF); % uses inv() --> changed to pinv(), needs initial values 
% --------------------------------------------------------------------------------------------------

% PLOT THE KF/KS ESTIMATES OF THE STATES 
% --------------------------------------------------------------------------------------------------
if PLOT_STATES
  clf; tiledlayout(4,2,TileSpacing="compact",Padding="compact");
  % make plot names
  % plot_names = make_table_names('$\epsilon_{', 1:k, 't}$');
  % loop through plots
  for ii = k+1:dim_X
    nexttile
    hold on;
      plot(Xs(:,ii), 'LineWidth',3); 
      plot(KS_deJ.att(:,ii),'--','Color',clr(3),'LineWidth',2.5);   % Filtered States
      % plot(KS_deJ.atT(:,ii),'--','Color',clr(3),'LineWidth',2.5);   % Smoothed States
    hold off;
    hline(0)
    box on; grid on;
    set(gca,'GridLineStyle',':' ,'GridAlpha',1/3, 'LineWidth',5/5);
    add2yaxislabel;
    addlegend({'True','Estimate:$\,a_{t|T}$'},1)
    % addsubtitle(plot_names(ii-k),-1.10)
  end
end
% % UNCOMMENT TO SHOW SCATTER PLOTS 
% % for ii = 6:size(KS_deJ.atT,2)
% %   nexttile
% %   hold on;
% %     scatter(Xs(:,ii),KS_deJ.att(:,ii)); 
% %   hold off;
% %   box on; grid on;
% %   set(gca,'GridLineStyle',':' ,'GridAlpha',1/3, 'LineWidth',5/5);
% %   add2yaxislabel;
% %   addlegend({'True vs. Estimate'},1)
% %   addsubtitle(row_names{ii-5})
% % end
% --------------------------------------------------------------------------------------------------

% CORRELATIONS (can also read off directly from the corr_table below -------------------------------
% print simple correlations
corr_table = array2table( corr(Xs(:,ss), KS_deJ.atT(:,ss)), ...
             'RowNames', row_names, 'VariableNames',row_names);
print_table(corr_table,4,1,'Correlation matrix of True and estimated smoothed States');sep

% % CALCULATE CORRELATION between ∆r*(t) and ET∆r*(t) [ie., actual vs. estimate)
% if ADD_Drstar 
%   PHI = PtT(end,end);
%   % Theoretical stdev 
%   sig_KS_Drstar = std( KS_deJ.atT(:,end) );
%   % Theoretical from model
%   sig_Drstar = sqrt( (4*c*s5)^2 + s3^2 );
%   % correlation between smoothed dr* and actual 
%   rho = 0.5*(sig_Drstar^2 + sig_KS_Drstar^2 - PHI) / (sig_Drstar*sig_KS_Drstar);
%   fprintf('Corr(ET(∆r*(t)),∆r*(t) from simulation) Analytical formulas: %4.4f\n', rho);
% end

% Define/Make eta(i) = eps(i)
for jj = 1:dim_R
  eval(['etT_' num2str(jj) ' = KS_deJ.atT(:,k+' num2str(jj) ');']);
  eval(['ett_' num2str(jj) ' = KS_deJ.att(:,k+' num2str(jj) ');']);
end

% Check some identities by running ols regressions: ie., ∆ETη5t = 0.107∆ETη3t − 0.028ETη4t 
fprintf('\n');sep('=');fprintf('Filter Identity on page 17 in EER(2022). Dependent variable: Etε2(t) \n')
Xnames_ID1 = {'Etε1(t)'};  % Xnames_ID1 = [];
ID1 = ols(ett_2, [ ett_1 ], 1, Xnames_ID1);

fprintf('\n');sep('=');fprintf('Smoother Identity on page 10 in Stars(2023). Dependent variable: Δ²ETε1(t) \n')
sep;fprintf('NOTE:    Identity should be: Δ²ETε1(t) = 1/40 ETε2(t-2) (and not ETε2(t) as stated)\n');sep
Xnames_ID2 = {'ETε2(t-2)'};  % Xnames_ID2 = [];
ID2 = ols( delta(etT_1, 2) , [ lag(etT_2, 2)], 1, Xnames_ID2);

% % ------------------------------------------------------------------------------------------------
% UNCOMMENT TO USE US REAL GDP DATA
% USE HP Filter routine to back out trend and cycle shocks and compare
% --------------------------------------------------------------------------------------------------
load('US-GDP.mat')
y   = 100.*log(usGDP);
D2y = delta(y, 2);    % don't use diff due to y being a TimeTable, --> use delta instead
% make Z(t) variable for 'shock recovery' SSM form: ie. Z(t)=Δ²y(t)
ZZ  = D2y.("Δ²gdp");
HP  = hp_filter(y, phi^2);
[~, Kurz_KF_US] = Kurz_Filter(removenans(ZZ), D1, D2, R, A, C, a00, P00);
% construct the cycle from the shock recovery SSM
KS_deJ_US = Kurz_DeJongKohnAnsley_Smoother(D1, D2, A, Kurz_KF_US);
SSM.cycle = phi*addnans(KS_deJ_US.atT(:,2),2,0);
SSM.trend = HP.trend(1)+cumsum(cumsum(KS_deJ_US.atT(:,1)));

% PLOT trend
clf;
% tiledlayout(3,1)
tiledlayout(3,1, TileSpacing = 'loose', Padding = 'compact');
nexttile
hold on;
  plot(HP.trend)
  % plot(y.gdp-SSM.cycle,'--')
  plot(SSM.trend,'--')
  hline(0)
hold off; 
box on; addgrid;
ylim([7 11]*100)
setdateticks(HP.Date,25)
addlegend({'HP-Filter','Shock-Recovery SSM'},1)
addsubtitle('Trend in US-GDP data')

% PLOT cycles 
% subplot(2,1,2)
nexttile
hold on;
  plot(HP.cycle)
  plot(SSM.cycle,'--')
  hline(0)
hold off; 
box on; addgrid;
% ylim([-.10 .06])
setdateticks(HP.Date,25)
addlegend({'HP-Filter','Shock-Recovery SSM'})
addsubtitle('Cycle in US-GDP data')































%EOF