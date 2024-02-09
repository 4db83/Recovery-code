% MR17_with_MP_rule Shock recovery McCririck R and D Rees (2017),  The Neutral Interest Rate’, RBA Bulletin, September, pp 9–18.
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
% addpath(genpath('D:/matlab.tools/db.toolbox/db'))       % set path to db functions
% CALL: get_all_db_toolbox_function_calls.m from Directory of code to be shared

% Sample size and seed for random number generator in simulation
Ts = 1e5; rng(10);    % takes about 1 sec for 1e5, 10 secs. for 1e6, 90 secs. for 1e7. --> does not change correlations from sims much
PLOT_STATES     = 1;  % set to 1 to plot ε(t) states
ADD_Drstr       = 1;  % set to 1 if wanting to add ∆r*(t) to State vector X(t)

% ----------------------------------------------------------------------------- % THIS IS WHAT YOU GET WHEN RUNNING THEIR LW CODE                                          
% PARAMETERS: Posterior MEANS (from Table A2: Parameter Estimates on page 17 in Bulletin).                         
ay1   =  1.48;      % ay1   
ay2   =  -.53;      % ay2                                                                            
ar    =   .06;      % ar -> NOTE: MR17 use -ar/2 in SSM instead of + ar/2 as in LW03 who have a negative estimate. It is the same in the end, but be careful with the construction and comparison
beta2 = -0.33;      % beta_2   --> unemployment gap
beta  =  0.64;      % beta     --> output gap in Okun's law
b1    =  0.3;       % b1       --> interest rate smoothing parameter
% standard deviations                                                                               
s1 = 0.37;          % s(ytild) --> IS curve                                                                        
s2 = 0.80;          % s(pi)    --> Phillips curve                                                                       
s3 = 0.34;          % s(z)     --> Unexplained r*                                                                           
s4 = 0.55;          % s(ystar) --> Trend output
s5 = 0.05;          % s(g)     --> Trend growth, QUARTERLY RATE, NOT ANNUALIZED ???
s6 = 0.15;          % s(u*)    --> NAIRU
s7 = 0.07;          % s(u)     --> Unemployment 
s8 = 1.19;          % s(r)     --> MP rule shock
sDr = sqrt( (4*s5)^2+s3^2 ); % stdev(∆r*(t))
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
% Unexplained r*                 0.22   0.34      Inverse gamma 0.40 0.25
% Trend output                   0.54   0.55      Inverse gamma 1.00 1.00
% Trend growth                   0.05   0.05      Inverse gamma 0.25 0.50
% NAIRU                          0.15   0.15      Inverse gamma 0.40 0.25
% Unemployment                   0.07   0.07      Inverse gamma 0.25 0.25

% DEFINE SSF INPUT MATRICES ------------------------------------------------------------------------
dim_Z = 4;                  % rows Z(t)
dim_X = 15 + ADD_Drstr;     % rows X(t)
dim_R = 8;                  % rows ε(t)
% offset for the first k latent state variables before the shocks.
k = dim_X - dim_R - ADD_Drstr; 

% --------------------------------------------------------------------------------------------------    % over the sample.start <- c(1961,1) sample.end   <- c(2002,2) with data vintage 2018      
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
C = zeros(dim_X,dim_R); C(k+1:(dim_X-ADD_Drstr),:) = eye(dim_R);
C(1,4) = s4; C(4,5) = s5; C(5,[3 5]) = [s3 4*s5]; C(7,6) = s6;
% USE g(t) instead of g(t-1), ie., ∆y*(t) = g(t) + sigma_4*ε4(t), add sigma_5*ε5(t) to the baseline ∆y*(t) = g(t-1) + sigma_4*ε4(t)
C(1,5) = s5;
if ADD_Drstr; C(end,[3 5]) = [s3 4*s5]; end
% --------------------------------------------------------------------------------------------------

% CALL TO THE KURZ_SSF FUNCTION --------------------------------------------------------------------
Pstar = Kurz_steadystate_P(D1, D2, R, A, C);
Neps  = k+1:dim_X;    % shock index in States X(t)
row_names = make_table_names('ε',1:dim_R,'(t)');              % make display names 
if ADD_Drstr; row_names = [row_names; {'∆r*(t)'}]; end        % add ∆r* to display names 

Pstar = array2table([ diag(Pstar.tT(Neps,Neps)) diag(Pstar.tt(Neps,Neps)) ], 'VariableNames',{'P*(t|T)','P*(t|t)'}, 'RowNames', row_names);
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
% Smoother 
% --------------------------------------------------------------------------------------------------
% Modified de Jong (1988, 1989) and Kohn and Ansley (1989) smoother (Eq. (4.11) in Kurz (2018))
KFS_deJ = Kurz_DeJongKohnAnsley_Smoother(D1, D2, A, Kurz_KF); % Contains KF and KS output. NO INV, NO INITVALS FOR STATES
% --------------------------------------------------------------------------------------------------

% CORRELATIONS:
% --------------------------------------------------------------------------------------------------
% Correlation between the true and estimated states. 
% R2 of Plagborg-Møller and Wolf (2022) 
R2 = [];
for jj = Neps
  pwR2.( ['e' num2str(jj)]) = ols(KFS_deJ.atT(:,jj),Xs(:,jj),1,[],[],[],0); % set last 0 to 1 to print to screen
  R2(jj-k,:) = eval((['pwR2.e' num2str(jj) '.R2']));
end

% compute the theoretical correlations implied by formula (10)
STDs  = [ones(dim_R,1)]; if ADD_Drstr; STDs  = [ones(dim_R,1); sDr]; end    % theoretical/model stdevs.
rho_theory = corr_theory(STDs, std(KFS_deJ.atT(:,Neps))', Pstar.('P*(t|T)'));

corr_table = array2table( [ diag(corr(Xs(:,Neps),KFS_deJ.atT(:,Neps))) rho_theory R2], ...
  'RowNames', row_names, 'VariableNames', {'ρ(Sim)','ρ(Theory)','R²(Sim)'});
% print correlations simulated and KS shocks
print_table(corr_table(1:dim_R+ADD_Drstr,:),4,1,'Correlation between True X(t) and (estimated) Kalman Smoothed States ETX(t)');sep

% Correlation matrix from KS estimates, Truth is uncorrelated
corr_XtT = array2table( corr(KFS_deJ.atT(Neps,Neps)), 'RowNames', row_names, 'VariableNames', row_names);
print_table(corr_XtT(1:dim_R+ADD_Drstr,1:dim_R+ADD_Drstr),4,1,'Correlation Matrix of (estimated) Kalman Smoothed States ETX(t)');sep
% Correlation matrix from KF estimates, Truth is uncorrelated
corr_Xtt = array2table( corr(KFS_deJ.att(Neps,Neps)), 'RowNames', row_names, 'VariableNames', row_names);
print_table(corr_Xtt(1:dim_R+ADD_Drstr,1:dim_R+ADD_Drstr),4,1,'Correlation Matrix of (estimated) Kalman Filtered States EtX(t)',[],0);sep

% DISPLAY RECOVEY DIAGNOSTICS ALL IN ONE MATRIX TO PRINT TO LATEX
% mat2latex([Pstar.("P*(t|T)")'; corr_table.("ρ(Sim)")'; corr_table.("R²(Sim)")']);

%% PLOT THE KF/KS ESTIMATES OF THE STATES 
% --------------------------------------------------------------------------------------------------
% make some plotting variables
WHICH_STATE_2_PLOT  = 0;      % set to 1 to use KF output, otherwise use KS
state_t = KFS_deJ.atT;
if WHICH_STATE_2_PLOT; state_t = KFS_deJ.att; end
xgrd = linspace(-5,5,100)';   % make xgrd for plotting
STL = -1.26;                  % subtitle location
dims = [-6:2:6]; FNS = 11; XOS = 11;

if PLOT_STATES
  TL = tiledlayout(9,2); if ~verLessThan('matlab','9.9');TL.TileSpacing='compact';TL.Padding='loose'; end
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
    fprintf('     Normalized P*(t|T)(∆r*) = %2.4f \n', Pstar.('P*(t|T)')(end)/var(Xs(:,end)))
  end    

  % UNCOMMENT TO PRINT TO PDF
  % print2pdf('MR17_MP_plots_KS',2); % super slow here
  % exportgraphics(gcf,'MR17_MP_plots_KS.pdf','ContentType','vector')
end
% --------------------------------------------------------------------------------------------------

% IDENTITIES:  
% --------------------------------------------------------------------------------------------------
% ET∆r*(t) = ETη5t + ETη3t (21) --> this should be: ET∆r*(t) = 4*c*ETη5t + ETη3t
% make ET∆r*(t) from X(:,4)-X(:,5) or alternatively from KS_deJ.atT(:,11) if ADD_rstr == 1;
% --------------------------------------------------------------------------------------------------
% define/make: ETεi(t) or ETηi(t) as needed
for jj = 1:dim_R
  eval(['ETn' num2str(jj) 't = KFS_deJ.atT(:,k+' num2str(jj) ')*s' num2str(jj) ';']);
  eval(['Etn' num2str(jj) 't = KFS_deJ.att(:,k+' num2str(jj) ')*s' num2str(jj) ';']);
end
Drstar = KFS_deJ.atT(:,4) - KFS_deJ.atT(:,5);

sep(133,'=',1); fprintf('Identity (20). Dependent variable: ∆ETη5(t) \n')
Xnames_ID1 = {'∆ETη3(t)','ETη4(t)'};
ID1 = ols( delta(ETn5t), [ delta(ETn3t) ETn4t], 1, Xnames_ID1);

sep(133,'=',1); fprintf('Identity (21). Dependent variable: ET∆r*(t) \n')
Xnames_ID2 = {'4*ETη5(t)','ETη3(t)'};
ID2 = ols(Drstar, [4*ETn5t ETn3t], 1, Xnames_ID2);

sep(133,'=',1); fprintf('Identity (22). Dependent variable: ET∆r*(t) \n')
Xnames_ID3 = {'ET∆r*(t-1)','ETη1(t)','ETη2(t)','ETη4(t)','ETη1(t-1)','ETη4(t-1)'};
ID3 = ols(Drstar, [lag(Drstar) ETn1t ETn2t ETn4t lag(ETn1t) lag(ETn4t)], 1, Xnames_ID3 );





toc































%EOF









