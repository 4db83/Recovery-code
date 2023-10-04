% LW(2003) recovery check
% SSF: ------------------------------------------------------------------------------------
%   Z(t) = D1*S(t) + D2*S(t-1) + R*ε(t),      S(t) = latent States
%   S(t) = A*S(t-1)            + C*ε(t),      ε(t) ~ MN(0,I)
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

%set seed for random number generator
rng(1);

% standard deviations
phi = 40;
% -----------------------------------------------------------------------------------------
dim_Z = 1;       % rows Z(t)
dim_S = 3;       % rows S(t)
dim_R = 2;       % rows ε(t)
% -----------------------------------------------------------------------------------------
% Define D1
D1 = zeros(dim_Z,dim_S); 
D1(1,1) = 1;  D1(1,2) = phi; D1(1,3) = -2*phi;
% Define D2
D2 = zeros(dim_Z,dim_S); 
D2(1,3) = phi;
% Define R
R  = zeros(dim_Z,dim_R);
% -----------------------------------------------------------------------------------------
% Define A
A = zeros(dim_S); 
A(3,2) = 1; 
% Define C
C = [eye(dim_R); zeros(1,2); ];
% C(1,1) = 1; C(2,2) = 1;
% --------------------------------------------------------------------------------------------------

% outerjoin(d2y,dy)

% -----------------------------------------------------------------------------------------
% call to the KURZ_SSM function
[PtT, Ptt] = Kurz_steadystate_P(D1,D2,R,A,C);
DD = 1:2;
% row_names = {'ε*(t)','εᶜ(t)','εᶜ(t-1)'} ; 
row_names = {'ε*(t)','εᶜ(t)'} ; 
Pstar = array2table([ diag(PtT(DD,DD)) diag(Ptt(DD,DD))],'VariableNames',{'P(t|T)','P(t|t)'},...
        'RowNames',row_names);
print_table(Pstar,4,1)

%% SIMULATE DATA
% 
TT = 1e4; 
UU = randn(TT,1);
[Zs, Xs] = Kurz_simulate_SSF(D1, D2, R, A, C, dim_Z, dim_S, dim_R, TT);
Z = Zs;

% USE REAL GDP DATA
load('usGDP.mat')
y   = log(usGDP);
dy  = delta(y,1,1);
d2y = delta(dy,1,1);
Z   = d2y.("ΔΔgdp");
[HPc, HPt] = hp_filter(y.gdp, phi^2);
% --------------------------------------------------------------------------------------------------
% FUNCTIONS FROM KURZ'S GITHUB PAGE, MILDLY MODIFIED TO NOT REQUIRE INITVAL COMP. 
% AND USE OF PINV IN AM SMOOTHER OTHERWISE NON-SINGULARITY ISSUES.
% --------------------------------------------------------------------------------------------------
% FILTER
% initialize filter. Note that errors will always be N(0,1), but latent states may need to be more carefully initialized.
a00 = zeros(dim_S, 1); P00 = eye(dim_S);
[~, Kurz_KF] = Kurz_Filter(Z, D1, D2, A, C, R, a00, P00);
% SMOOTHERS 
% Modified Anderson and Moore (1979) smoother (Eq. (4.3) in Kurz (2018))
% KS_AM   = Kurz_AndersonMoore_Smoother(   D1, D2, A, Kurz_KF.Z_tilde, Kurz_KF.Finv, Kurz_KF.K, Kurz_KF.a_t_t, Kurz_KF.P_t_t);
% Modified de Jong (1988, 1989) and Kohn and Ansley (1989) smoother (Eq. (4.11) in Kurz (2018))
% KS_deJ  = Kurz_DeJongKohnAnsley_Smoother(D1, D2, A, Kurz_KF.Z_tilde, Kurz_KF.Finv, Kurz_KF.K, Kurz_KF.a_t_t, Kurz_KF.P_t_t);
% --------------------------------------------------------------------------------------------------
% Modified Koopman (1993) smoother (Eq. (4.14)-(4.15) in Kurz (2018))
KS_K    = Kurz_Koopman_Smoother(D1, D2, A, C, R, Kurz_KF.Z_tilde, Kurz_KF.Finv, Kurz_KF.K);

%% plot the KF/KS estimates of the states
clf;
hold on;
% plot(y.gdp);
  % plot(HPt,'--','Color',clr(2)); 
  plot(HPc,'-','Color',clr(2)); 
  plot(addnans(KS_K.a_t_T(:,2),2,0)*phi,'--','Color',clr(1)); 
hold off;
setrotatedateticks(y.Date,7,'m')
ylim([-.1 .1]); hline(0);
box on; grid on;
set(gca,'GridLineStyle',':' ,'GridAlpha',1/3, 'LineWidth',5/5);
setoutsideTicks; add2yaxislabel

lst([HPc addnans(KS_K.a_t_T(:,2),2,0)*phi],10)


% open setplot



% for ii = 1:dim_S
%   subplot(1,3,ii);
%   hold on;
%   plot(Xs(:,ii))
%   plot(KS_AM.a_t_T(:,ii), '--')
%   plot(KS_deJ.a_t_T(:,ii), '--')
%   plot(KS_K.a_t_T(:,ii), '--')
%   plot(Kurz_KF.a_t_t(:,ii), '--')
%   hold off;
%   box on; grid on;
% end



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
%∆ET η5t = 0.107∆ET η3t − 0.028ET η4t.

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