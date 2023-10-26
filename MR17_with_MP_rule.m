% MR17_with_MP_rule Shock recovery (https://www.newyorkfed.org/research/policy/rstar)
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
% set to 1 if wanting to add ∆r*(t) to State vector X(t)
ADD_Drstar  = 1;
PLOT_STATES = 0;

% DEFINE SSM INPUT MATRICES ------------------------------------------------------------------------
dim_Z = 4;                % rows Z(t)
dim_X = 15 + ADD_Drstar;  % rows X(t)
dim_R = 8;                % rows ε(t)
% offset for the first k latent state variables before the shocks.
k = dim_X - dim_R - ADD_Drstar; 
% --------------------------------------------------------------------------------------------------
% Posterior Means (from Table A2: Parameter Estimates on page 17 in Bulletin).    
% MP rule equation is from MARTIN Model equation 36 converted to real rate via i(t) = r(t) + pi(t).
ay1   =  1.48;      % ay1   paper only gives the sum as 0.945                                         
ay2   =  -.53;      % ay2                                                                            
ar    =   .06;      % ar -> NOTE: MR17 use -ar/2 in SSM instead of + ar/2 as in LW03 who have a negative estimate. It is the same in the end, but be careful with the construction and comparison
beta2 = -0.33;      % beta_2   --> unemployment gap
beta  =  0.64;      % beta     --> output gap in Okun's law
b1    =  0.3;
% standard deviations                                                                               
s1 = 0.37;          % s(ytild) --> IS curve                                                                        
s2 = 0.80;          % s(pi)    --> Phillips curve                                                                       
s3 = 0.34;          % s(z)     --> Unexplained r*                                                                           
s4 = 0.55;          % s(ystar) --> Trend output
s5 = 0.05;          % s(g)     --> Trend growth, quarterly rate, not annualized ???
s6 = 0.15;          % s(u*)    --> NAIRU
s7 = 0.07;          % s(u)     --> Unemployment 
s8 = 1.19;          % s(r)     --> MP rule shock
% Structural parameters -------- Mode   Mean  ------------------------------------------------------
% IS curve – y!t−1               1.53   1.48      Normal 1.10 1.50
% IS curve – y!t−2              –0.54  –0.53      Normal –0.20 1.50
% IS curve – rt(L)−rt∗(L)        0.05   0.06      Inverse gamma 0.15 1.00
% Phillips curve – t (L)         0.39   0.41      Beta 0.50 0.25
% Phillips curve – ut−1−ut−1∗   –0.32  –0.33      Normal –0.50 0.30
% Okun’s law – y!t (L)           0.62   0.64      Normal 0.50 0.3
% Shock processes ----------------------------------------------------------------------------------
% IS curve                       0.38   0.37      Inverse gamma 1.00 1.00
% Phillips curve                 0.79   0.80      Inverse gamma 1.00 1.00
% Unemployment                   0.07   0.07      Inverse gamma 0.25 0.25
% Trend output                   0.54   0.55      Inverse gamma 1.00 1.00
% NAIRU                          0.15   0.15      Inverse gamma 0.40 0.25
% Trend growth                   0.05   0.05      Inverse gamma 0.25 0.50
% Unexplained r*                 0.22   0.34      Inverse gamma 0.40 0.25

% --------------------------------------------------------------------------------------------------
% Define D1
D1 = zeros(dim_Z,dim_X); 
D1(1,1:3) = [1 -ay1 -ay2];    D1(1,8)   = s1;
                              D1(2,9)   = s2;
D1(3,1:3) = -beta*[.4 .3 .2]; D1(3,14)  = s7;
D1(4,[5 7]) = [b1 2*b1];      D1(4,15)  = s8;
% Define D2
D2 = zeros(dim_Z,dim_X); 
D2(1,5:6) = -ar/2;
D2(2,7)   = -beta2;
D2(3,3)   = -beta*.1;
% Define R
R  = zeros(dim_Z,dim_R);
% --------------------------------------------------------------------------------------------------
% Define A
A = zeros(dim_X); 
A([1 2],1) = 1; A(3,2) = 1; A([1 4],4) = 1; 
A(5:6,5) = 1; A([1 4],4) = 1; A(7,7) = 1;
% Define C
C = zeros(dim_X,dim_R); C(k+1:(dim_X-ADD_Drstar),:) = eye(dim_R);
C(1,4) = s4; C(4,5) = s5; C(5,[3 5]) = [s3 4*s5]; C(7,6) = s6;
if ADD_Drstar; C(end,[3 5]) = [s3 4*s5]; end
% --------------------------------------------------------------------------------------------------

% CALL TO THE KURZ_SSM FUNCTION --------------------------------------------------------------------
[PtT, Ptt] = Kurz_steadystate_P(D1, D2, R, A, C);
ss = k+1:dim_X;
row_names = make_table_names('ε',1:dim_R,'(t)');              % make display names 
if ADD_Drstar; row_names = [row_names; {'∆r*(t)'}]; end   % add ∆r* to display names 

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
[~, Kurz_KF] = Kurz_Filter(Z, D1, D2, R, A, C, a00, P00);
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
  plot_names = make_table_names('$\epsilon_{', 1:k, 't}$');
  if ADD_Drstar; plot_names = [plot_names; '$\Delta r^{\ast}_{t}$']; end
  % loop through plots
  for ii = k+1:dim_X
    nexttile
    hold on;
      plot(Xs(:,ii), 'LineWidth',3); 
      % plot(KS_deJ.att(:,ii),'--','Color',clr(3),'LineWidth',2.5);   % Filtered States
      plot(KS_deJ.atT(:,ii),'--','Color',clr(3),'LineWidth',2.5);   % Smoothed States
    hold off;
    hline(0)
    box on; grid on;
    set(gca,'GridLineStyle',':' ,'GridAlpha',1/3, 'LineWidth',5/5);
    add2yaxislabel;
    addlegend({'True','Estimate:$\,a_{t|T}$'},1)
    addsubtitle(plot_names(ii-k),-1.10)
  end
end
% % % UNCOMMENT TO SHOW SCATTER PLOTS 
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

%% CORRELATIONS (can also read off directly from the corr_table below -------------------------------
% print simple correlations
corr_table = array2table( corr(Xs(:,k+1:end), KS_deJ.atT(:,k+1:end)), ...
             'RowNames', row_names, 'VariableNames',row_names);
print_table(corr_table,4,1,'Correlation matrix of True and estimated smoothed States');sep

% CALCULATE CORRELATION between ∆r*(t) and ET∆r*(t) [ie., actual vs. estimate)
if ADD_Drstar 
  PHI = PtT(end,end);
  % Theoretical stdev 
  sig_KS_Drstar = std( KS_deJ.atT(:,end) );
  % Theoretical from model
  sig_Drstar = sqrt( (4*s5)^2 + s3^2 );
  % correlation between smoothed dr* and actual 
  rho = 0.5*(sig_Drstar^2 + sig_KS_Drstar^2 - PHI) / (sig_Drstar*sig_KS_Drstar);
  fprintf('Corr(ET(∆r*(t)),∆r*(t) from simulation) Analytical formulas: %4.4f\n', rho);
end

% Define/Make eta(i) = sig(i)*eps(i)
for jj = 1:dim_R
  eval(['n' num2str(jj) ' = s' num2str(jj) '*KS_deJ.atT(:,k+' num2str(jj) ');']);
end

% % Check some identities by running ols regressions: ie., ∆ETη5t = 0.107∆ETη3t − 0.028ETη4t 
% fprintf('\n');sep('=');fprintf('Identity (20). Dependent variable: ∆ETη5t \n')
% Xnames_ID1 = {'∆ETη3(t)','ETη4(t)'};
% % Xnames_ID1 = [];
% ID1 = ols( delta(n1), [lag(delta(n1)) lag(delta(n3)) lag(n4,2)], 1, Xnames_ID1);

% % --------------------------------------------------------------------------------------------------
% % OTHER IDENTITIES 
% % --------------------------------------------------------------------------------------------------
% % ET∆r*(t) = ETη5t + ETη3t (21) --> this should be: ET∆r*(t) = 4*ETη5t + ETη3t
% % make ET∆r*(t) from X(:,4)-X(:,5) or alternatively from KS_deJ.atT(:,11);
% % --------------------------------------------------------------------------------------------------
% Drstar = KS_deJ.atT(:,4)-KS_deJ.atT(:,5);
% 
% fprintf('\n');sep('=');fprintf('Identity (21). Dependent variable: ET∆r*(t) \n')
% ID2 = ols(Drstar, [4*n5 n3], 1, {'4*ETη5(t)','ETη3(t)'});
% 
% fprintf('\n');sep('=');fprintf('Identity (22). Dependent variable: ET∆r*(t) \n')
% ID3 = ols(Drstar, [lag(Drstar) n1 n2 n4 lag(n1) lag(n4)], 1, ...
%       {'ET∆r*(t-1)','ETη1(t)','ETη2(t)','ETη4(t)','ETη1(t-1)','ETη2(t-4)'});
% 
% 



































%EOF