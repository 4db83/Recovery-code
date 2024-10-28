% HLW17 Shock recovery SSF (https://www.newyorkfed.org/research/policy/rstar)
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

% IN PAPER USE 1e5: Sample size and seed for random number generator in simulation
Ts = 1e5; rng(10);  % takes about 1 sec for 1e5, 10 secs. for 1e6, 90 secs. for 1e7. --> does not change correlations from sims much
PLOT_STATES = 0;    % set to 1 to plot ε(t) states
ADD_Drstr   = 1;    % set to 1 if wanting to add ∆r*(t) to State vector X(t)
PLOTS2PDF   = 0;    % set to 1 to print plots to PDF.

% ----------------------------------------------------------------------------- % THIS IS WHAT YOU GET WHEN RUNNING THEIR LW CODE                                          
% PARAMETERS: (from the published paper, Table 1 on page S60.                                           
ay1 =  1.530;       % paper only gives the sum as 0.942
ay2 = -0.588;
ar  = -0.071;
by  =  0.079;
bpi =  0.668;       % not needed
c   =  1.0;         % set to 1 in HLW17
% standard deviations      (NOTE the order in Table 1 is different, here I use the same order as in the LW03 code)                                                                         
s1  =  0.354;       % sigma(ytild)                                  % 0.354;                                         
s2  =  0.791;       % sigma(pi)                                     % 0.791;                                         
s3  =  0.150;       % sigma(z)                                      % 0.150;                                         
s4  =  0.575;       % sigma(ystar)                                  % 0.575;                                         
s5  =  0.122/4;     % sigma(g) --> reported as annualized rate      % 0.122/4;                                  
% divide by 4 to express in quarterly rate (annualized later in the C Matrix)                         
sDr = sqrt( (c*4*s5)^2+s3^2 ); % stdev(∆r*(t))
















% DEFINE SSF INPUT MATRICES ------------------------------------------------------------------------
dim_Z = 2;                  % rows Z(t)
dim_X = 10 + ADD_Drstr;     % rows X(t)
dim_R = 5;                  % rows ε(t)
% offset for the first k latent state variables before the shocks.
k = dim_X - dim_R - ADD_Drstr; 

% --------------------------------------------------------------------------------------------------    % over the sample.start <- c(1961,1) sample.end   <- c(2002,2) with data vintage 2018      
% Define D1
D1 = zeros(dim_Z,dim_X); 
D1(1,1) = 1;    D1(1,6) = s1;
D1(2,2) = -by;  D1(2,7) = s2;
% Define D2
D2 = zeros(dim_Z,dim_X); 
D2(1,1:2) = [-ay1 -ay2]; D2(1,4:5) = -ar/2;
% Define R
R  = zeros(dim_Z,dim_R);
% --------------------------------------------------------------------------------------------------
% Define A
A = zeros(dim_X); 
A(1:2,1) = 1; A(4:5,4) = 1; A([1 3],3) = 1;
% Define Q
Q = zeros(dim_X,dim_R); Q(k+1:(dim_X-ADD_Drstr),:) = eye(dim_R);
Q(1,4) = s4; Q(3,5) = s5; Q(4,[3 5]) = [s3 4*c*s5];
% USE g(t) instead of g(t-1), ie., ∆y*(t) = g(t) + sigma_4*ε4(t), add sigma_5*ε5(t) to the baseline ∆y*(t) = g(t-1) + sigma_4*ε4(t)
% C(1,5) = s5;
if ADD_Drstr;  Q(end,[3 5]) = [s3 4*c*s5]; end
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

Neps = k+1:dim_X;    % shock index in States X(t)
row_names = make_table_names('ε',1:dim_R,'(t)');              % make display names 
% add ∆r* to display names 
if ADD_Drstr; row_names = [row_names; {'∆r*(t)'}]; end        
% make Pstar table
Pstar = array2table([ diag(KFS.Pstar.tT(Neps,Neps)) diag(KFS.Pstar.tt(Neps,Neps)) ], 'VariableNames',{'P*(t|T)','P*(t|t)'}, 'RowNames', row_names);

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
STDs  = [ones(dim_R,1)]; if ADD_Drstr; STDs  = [ones(dim_R,1); sDr]; end    % theoretical/model stdevs.
rho_theory = corr_theory(STDs, std(KFS.atT(:,Neps))', Pstar.('P*(t|T)'));

if ADD_Drstr 
  % normalize ∆r*(t) by Var(∆r*(t)) from theoretical model to make comparable to other measures
  KFS.Pstar.tT(end,end) = KFS.Pstar.tT(end,end)/sDr^2; 
  KFS.Pstar.tt(end,end) = KFS.Pstar.tt(end,end)/sDr^2; 
  % make/update Pstar to be normalized by  var(∆r*) for the last entry
  Pstar = array2table([ diag(KFS.Pstar.tT(Neps,Neps)) diag(KFS.Pstar.tt(Neps,Neps)) ], 'VariableNames',{'P*(t|T)','P*(t|t)'}, 'RowNames', row_names);
end 
% print Pstar to screen
sep; print_table(Pstar,4,1,0)

corr_table = array2table( [ diag(corr(Xs(:,Neps),KFS.atT(:,Neps)))  rho_theory  R2], 'RowNames', row_names, 'VariableNames', {'ρ(Sim)','ρ(Theory)','R²(Sim)'});
% print correlations simulated and KS shocks
print_table(corr_table(1:dim_R+ADD_Drstr,:),4,1,'Correlation between True X(t) and (estimated) Kalman Smoothed States ETX(t)',[],0);

% Correlation matrix from KS estimates, Truth is uncorrelated
corr_XtT = array2table( corr(KFS.atT(:,Neps)), 'RowNames', row_names, 'VariableNames', row_names);
print_table(corr_XtT(1:dim_R+ADD_Drstr,1:dim_R+ADD_Drstr),4,1,'Correlation Matrix of (estimated) Kalman Smoothed States ETX(t)',[],0);

% Correlation matrix from KF estimates, Truth is uncorrelated
corr_Xtt = array2table( corr(KFS.att(:,Neps)), 'RowNames', row_names, 'VariableNames', row_names);
print_table(corr_Xtt(1:dim_R+ADD_Drstr,1:dim_R+ADD_Drstr),4,1,'Correlation Matrix of (estimated) Kalman Filtered States EtX(t)',[],0);sep

% DISPLAY RECOVEY DIAGNOSTICS ALL IN ONE MATRIX TO PRINT TO LATEX
matRowNames = { 'P*(t|T)  ';'R²(Sim)  ';'ρ(Theory)'};  matColNames = row_names; %matColNames = []
fprintf('Recovery Measures (Order is)\n'); sep
% sprintf('% s \n', row_names')
mat2latex([Pstar.("P*(t|T)")'; corr_table.("R²(Sim)")'; corr_table.("ρ(Theory)")'], 4, matRowNames, matColNames); sep
toc

%% PLOT THE KF/KS ESTIMATES OF THE STATES 
% --------------------------------------------------------------------------------------------------
% make some plotting variables
WHICH_STATE_2_PLOT  = 0;      % set to 1 to use KF output, otherwise use KS
state_t = KFS.atT;
if WHICH_STATE_2_PLOT; state_t = KFS.att; end
xgrd = linspace(-5,5,100)';   % make xgrd for plotting
STL = -1.26;                  % subtitle location
dims = [-6:2:6]; FNS = 11; XOS = 11;

if PLOT_STATES 
  clf; TL = tiledlayout(9,2); if ~verLessThan('matlab','9.9');TL.TileSpacing='compact';TL.Padding='loose'; end
  % loop through plots
  for ii = Neps(1:end-ADD_Drstr)
    nexttile
    hold on;
      plot(Xs(:,ii), 'LineWidth',3);                              % 'true' simulated state X
      plot(state_t(:,ii),'--','Color',clr(3),'LineWidth',2.5);    % Filtered or smoothed estimate of state X
    hold off; 
    xlim([-XOS length(Xs(:,ii))+XOS]); 
    setyticklabels(dims,0,FNS);
    addgrid(5/5); hline(0); 
    addlegend({'True','Estimate'},1,FNS-1)
    addsubtitle(['$\varepsilon_{' num2str(ii-k) 't}$'],STL,FNS)
    add2yaxislabel(1)

    nexttile
    hold on; 
      scatter(Xs(:,ii),state_t(:,ii),'x');
      plot(xgrd, xgrd*pwR2.( ['e' num2str(ii)]).bhat)
    hold off; 
    setxticklabels([dims(1):1:dims(end)]);
    setyticklabels(dims,0,FNS);
    addgrid(5/5); hline(0); 
    line(dims, dims, 'Color', 'k', 'LineWidth', 1); 
    ylabel(['$E_T\varepsilon_{' num2str(ii-k) 't}$ (Estimate)'],'Interpreter','latex','FontSize',FNS)
    addsubtitle(['$\varepsilon_{' num2str(ii-k) 't}$ (True)'],STL,FNS)
    addlegend({['$R^2=' num2str(pwR2.( ['e' num2str(ii)]).R2,'%2.4f') '$']},1,FNS-1)
    add2yaxislabel(1)
  end
  if ADD_Drstr
    dims = [-3:1:3]; xgrd = linspace(-2,2,100)';                   % adjust xgrd for plotting 
    nexttile
    hold on;
      plot(Xs(:,end), 'LineWidth',3);                              % 'true' simulated state X
      plot(state_t(:,end),'--','Color',clr(3),'LineWidth',2.5);    % Filtered or smoothed estimate of state X
    hold off; 
    xlim([-XOS length(Xs(:,ii))+XOS]); 
    setyticklabels(dims,0,FNS);
    addgrid(5/5); hline(0); 
    addlegend({'True','Estimate'},1,FNS-1)
    addsubtitle(['$\Delta r^{\ast}_{t}$'],STL,FNS)
    add2yaxislabel(1)

    nexttile
    hold on; 
      scatter(Xs(:,end),state_t(:,end),'x');
      plot(xgrd, xgrd*pwR2.( ['e' num2str(ii+1)]).bhat)
    hold off; 
    setxticklabels([dims(1):1:dims(end)]);
    setyticklabels(dims,0,FNS);
    addgrid(5/5); hline(0); 
    line(dims, dims, 'Color', 'k', 'LineWidth', 1); 
    ylabel(['$E_T\Delta r^{\ast}_{t}$ (Estimate)'],'Interpreter','latex','FontSize',FNS)
    addsubtitle(['$\Delta r^{\ast}_{t}$ (True)'],STL,FNS)
    addlegend({['$R^2=' num2str(pwR2.( ['e' num2str(ii+1)]).R2,'%2.4f') '$']},1,FNS-1)
    add2yaxislabel(1)
    % print the normalized steady-state P(t|T) for ∆r*, as it has a much smaller variance than the
    % unit variance of the shocks ε(t).
    sep;fprintf('Normalized P*(t|T)(∆r*) = %2.4f \n', Pstar.('P*(t|T)')(end)/var(Xs(:,end)))
  end    
  % TO PRINT TO PDF
  if PLOTS2PDF   
  % fig_name = ['LL_plots_KS_Noise2Signal_', num2str(N2S^2,4) '_T_' num2str(Ts,'%d') '.pdf'];
  fig_name = 'HLW17_plots_KS_0.pdf';
  % print2pdf(fig_name); % super slow here
    exportgraphics(gcf, fig_name, 'ContentType','vector')
  end
end
% --------------------------------------------------------------------------------------------------

%% IDENTITIES REGRESSIONS:
clc
% --------------------------------------------------------------------------------------------------
% define/make: ETεi(t) or ETηi(t) as needed
for jj = 1:dim_R
  eval(['ETn' num2str(jj) 't = KFS.atT(:,k+' num2str(jj) ')*s' num2str(jj) ';']);
  eval(['Etn' num2str(jj) 't = KFS.att(:,k+' num2str(jj) ')*s' num2str(jj) ';']);
end
Drstar = KFS.atT(:,4) - KFS.atT(:,5);

sep(133,'=',1); fprintf('Identity (20). Dependent variable: ∆ETη5(t) \n')
Xnames_ID1 = {'∆ETη3(t)','ETη4(t)'};
ID1 = ols( delta(ETn5t), [ delta(ETn3t) ETn4t], 1, Xnames_ID1);

sep(133,'=',1); fprintf('Identity (21). Dependent variable: ET∆r*(t) \n')
Xnames_ID2 = {'ETη5(t)','ETη3(t)'};
ID2 = ols(Drstar, [ETn5t ETn3t], 1, Xnames_ID2);

sep(133,'=',1); fprintf('Identity (22). Dependent variable: ET∆r*(t) \n')
Xnames_ID3 = {'ET∆r*(t-1)','ETη1(t)','ETη2(t)','ETη4(t)','ETη1(t-1)','ETη4(t-1)'};
ID3 = ols(Drstar, [lag(Drstar) ETn1t ETn2t ETn4t lag(ETn1t) lag(ETn4t)], 1, Xnames_ID3 );

sep(133,'=',1); fprintf('Dependent variable: ETη3(t) \n')
Xnames_ID4 = {'ETη5(t)','ET∆r*(t-1)','ETη1(t)','ETη2(t)','ETη4(t)','ETη1(t-1)','ETη4(t-1)'};
ID4 = ols(ETn3t, [ ETn5t lag(Drstar) ETn1t ETn2t ETn4t lag(ETn1t) lag(ETn4t)], 1, Xnames_ID4 );































%EOF









