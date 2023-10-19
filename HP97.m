% HP97 Shock recovery
% SSF: ------------------------------------------------------------------------------------
%   Z(t) = D1*X(t) + D2*X(t-1) + R*ε(t),      X(t) = latent States
%   X(t) = A*X(t-1)            + C*ε(t),      ε(t) ~ MN(0,I)
% -----------------------------------------------------------------------------------------
clear; clc;
% set plotting defaults
set(groot,'defaultLineLineWidth',2); set(groot,'defaultAxesFontSize',14)
set(groot,'defaultAxesXTickLabelRotationMode','manual')
set(groot,'defaultAxesFontName','Times New Roman')

addpath(genpath('./functions'))
addpath(genpath('D:/matlab.tools/db.toolbox/db'))           % set path to db functions
% get US GDP data for empirical application
% usGDP = as_timetable(getFredData('GDPC1', '1947-01-01', '2023-12-31','lin','q'),'gdp');  % 'https://fred.stlouisfed.org/data/GDPC1.txt'
% save('usGDP','usGDP')

% standard deviations
phi = 40;
% -----------------------------------------------------------------------------------------
dim_Z = 1;       % rows Z(t)
dim_X = 3;       % rows X(t)
dim_R = 2;       % rows ε(t)
% -----------------------------------------------------------------------------------------
% Define D1
D1 = zeros(dim_Z,dim_X); 
D1(1,1) = 1;  D1(1,2) = phi; D1(1,3) = -2*phi;
% Define D2
D2 = zeros(dim_Z,dim_X); 
D2(1,3) = phi;
% Define R
R  = zeros(dim_Z,dim_R);
% -----------------------------------------------------------------------------------------
% Define A
A = zeros(dim_X); 
A(3,2) = 1; 
% Define C
C = [eye(dim_R); zeros(1,2); ];
% C(1,1) = 1; C(2,2) = 1;
% --------------------------------------------------------------------------------------------------
% outerjoin(d2y,dy)

% -----------------------------------------------------------------------------------------
% call to the KURZ_SSM function
[PtT, Ptt] = Kurz_steadystate_P(D1, D2, R, A, C);
ss = 1:3;
row_names = {'ε1(t)',' ε2(t)','ε2(t-1)'} ; 
% row_names = {'ε*(t)','εᶜ(t)'} ; 
Pstar = array2table([ diag(Ptt(ss,ss)) diag(PtT(ss,ss)) ], ...
        'VariableNames',{'P(t|t)','P(t|T)'}, 'RowNames', row_names(ss));
print_table(Pstar,4,1)

%% USE SIMULATE DATA (FIX seed for random number generator)
TT = 5e4; rng(1);
[Zs, Xs, Us] = Kurz_simulate_SSF(D1, D2, R, A, C, dim_Z, dim_X, dim_R, TT);
Z = Zs; 
% [HPc, HPt] = hp_filter( cumsum(cumsum(Z)), phi^2);

% USE US REAL GDP DATA
% load('US-GDP.mat')
% y   = log(usGDP);
% d2y = delta( delta(y,1) ,1,1); % don't use diff due to y being a TimeTable
% Z   = d2y.("ΔΔgdp");
% [HPc, HPt] = hp_filter(y.gdp, phi^2);

% --------------------------------------------------------------------------------------------------
% FUNCTIONS FROM KURZ's GITHUB PAGE, MILDLY MODIFIED TO SIMPLIFY INPUT AND COMPARABILTY WITH MY CODE 
% ABOVE AND USE OF PINV IN AM SMOOTHER OTHERWISE NON-SINGULARITY ISSUES.
% --------------------------------------------------------------------------------------------------
% INITIALIZE FILTER 
% Note: errors will always be N(0,1), but latent states may need more careful initialization.
a00 = zeros(dim_X, 1); P00 = eye(dim_X);
[~, Kurz_KF] = Kurz_Filter(Z, D1, D2, R, A, C, a00, P00);
% SMOOTHERS 
% --------------------------------------------------------------------------------------------------
% Modified de Jong (1988, 1989) and Kohn and Ansley (1989) smoother (Eq. (4.11) in Kurz (2018))
KS_deJ  = Kurz_DeJongKohnAnsley_Smoother(D1, D2, A, Kurz_KF); % NO INV, NO INITVALS FOR STATES
% KS_deJ  = Kurz_AndersonMoore_Smoother(   D1, D2, A, Kurz_KF); % uses inv() --> changed to pinv(), needs initial values 
% --------------------------------------------------------------------------------------------------

% PLOT THE KF/KS ESTIMATES OF THE STATES
clf;
k = 1;
hold on;
  plot(Xs(:,k))
  plot(KS_deJ.atT(:,k))
  % plot(Us(:,k))
  % plot(HPc,'-','Color',clr(2),'LineWidth',3); 
  % plot(addnans(phi*KS_deJ.atT(:,2),2,0),'--','Color',clr(3),'LineWidth',2);
  % plot(addnans(phi*KS_AM.atT(:,2),2,0) ,'- ','Color',[0 0 0],'LineWidth',1);
hold off;
% setrotatedateticks(y.Date,7,'m')
% ylim([-.1 .06]); hline(0);
box on; grid on;
set(gca,'GridLineStyle',':' ,'GridAlpha',1/3, 'LineWidth',5/5);
setoutsideTicks; add2yaxislabel
disp( 'Correlation matrix of True and estimated smoothed States')
disp(corr(Xs,KS_deJ.atT))

%% CALCULATE CORRELATION 
% phi = PtT(end,end);
% % Theoretical stdev 
% KS_K_drr = KS_K.a_t_T(:,1);
% disp('std dr* estimate 2');
% std_KS_AM_drr = std(KS_AM_drr)
% 
% %stdev_estdrstar_sm =0.04222116394588856;
% stdev_estdrstar_sm = std_KS_AM_drr;
% % Theoretical stdev_rstar
% stdev_rstar = sqrt((cc*s5)^2 + s3^2);
% %correlation between smoothed dr* and actual 
% disp('correlation smoothed dr* estimate and dr*');
% rho = 0.5*(stdev_estdrstar_sm^2 - phi + stdev_rstar^2)/(stdev_rstar*stdev_estdrstar_sm);
% disp(rho)

%% Estimate the identity
%∆ETη5t = 0.107∆ETη3t − 0.028ETη4t.
% static correlation from KF: ε1(t) = 40*ε2(t)
% row_names = {'ε1(t)',' ε2(t)','ε2(t-1)'} ; 
eps1 = Kurz_KF.att(:,1);
eps2 = Kurz_KF.att(:,2);
% clc
ols_KF = ols(eps1, eps2, 1, 'ε2(t)');
clf; scatter(eps1, eps2)

% dynamic correlation 
eps1 = KS_deJ.atT(:,1);
eps2 = KS_deJ.atT(:,2);
ols_KF = ols(eps1, eps2, 1, 'ε2(t)');

ols_KS = ols(delta(delta(eps1)), eps2, 1, 'ε2(t)');
% clf; scatter(delta(delta(ep_tT)), ec_tT)

%Create epsilons (as in paper)
% Use KS_K.a_t_T
% epsilon1 = s1*KS_K.a_t_T(:,6);
% epsilon2 = s2*KS_K.a_t_T(:,7);
% epsilon3 = s3*KS_K.a_t_T(:,8);
% epsilon4 = s4*KS_K.a_t_T(:,9);
% epsilon5 = s5*KS_K.a_t_T(:,10);	
% 
% %Create differences
% depsilon3 = diff(epsilon3);
% depsilon5 = diff(epsilon5); 
% 
% %Estimate equation
% %ddnrr depsilon5 depsilon3
% Xmat1 = [depsilon3,epsilon4(2:end)];
% disp('coefficients first identity');
% b1 = regress(depsilon5,Xmat1)




































%EOF