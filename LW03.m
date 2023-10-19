% LW03 Shock recovery (https://www.newyorkfed.org/research/policy/rstar)
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

ADD_dr = 1; % set to 1 if 
% -----------------------------------------------------------------------------------------
dim_Z = 2;            % rows Z(t)
dim_X = 10 + ADD_dr;  % rows X(t)
dim_R = 5;            % rows ε(t)
% -----------------------------------------------------------------------------------------
% parameters (from the published paper, 'Baseline' in Table 1 on page 1065.
a1 =  1.517;      % ay1   paper only gives the sum as 0.945
a2 =  -.572;      % ay2
a3 =  -.098;      % ar
b3 =   .043;      % by 
c  =  1.068;      % c
% standard deviations
s1 =  0.387;      % sig_ytild 
s2 =  0.731;      % sig_pi    
s3 =  0.323;      % sig_z     
s4 =  0.605;      % sig_ystar 
s5 =  0.102/4;    % sig_g, annualized rate devided by 4 to express in quarterly rate --> annualized later in the relevant matrices
% this is what you get when running their LW code 
% D:\_research\_current\stars\recovery-code\_notes\LW_replication\kalman.standard.errors.R
%  [1,] "1.56495477703249"    "a_1"      
%  [2,] "-0.614537663213145"  "a_2"      
%  [3,] "-0.0614762410516162" "a_3"      
%  [4,] "0.560395442735451"   "b_1"      
%  [5,] "0.377932050558007"   "b_2"      
%  [6,] "0.0440540347304658"  "b_3"      
%  [7,] "0.00228697742079232" "b_4"      
%  [8,] "0.0379650747803849"  "b_5"      
%  [9,] "1.24070752805134"    "c"        
% [10,] "0.335593534431049"   "sigma_1"  
% [11,] "0.760431178504621"   "sigma_2"  
% [12,] "0.592270881950106"   "sigma_4"  
% over the sample.start <- c(1961,1) sample.end   <- c(2002,2) with data vintage 2018

% -----------------------------------------------------------------------------------------
% Define D1
D1 = zeros(dim_Z,dim_X); 
D1(1,1) = 1;    D1(1,6) = s1;
D1(2,2) = -b3;  D1(2,7) = s2;
% Define D2
D2 = zeros(dim_Z,dim_X); 
D2(1,1:2) = [-a1 -a2]; D2(1,4:5) = -a3/2;
% Define R
R  = zeros(dim_Z,dim_R);
% -----------------------------------------------------------------------------------------
% Define A
A = zeros(dim_X); 
A(1:2,1) = 1; A(4:5,4) = 1; A([1 3],3) = 1;
% Define C
C = [zeros(5,5); eye(dim_R)];
C(1,3) = s3; C(3,4) = s4; C(4,[4 5]) = [4*c*s4 s5];
% --------------------------------------------------------------------------------------------------

% -----------------------------------------------------------------------------------------
% call to the KURZ_SSM function
[PtT, Ptt] = Kurz_steadystate_P(D1, D2, R, A, C);
ss = 6:10;
row_names = {'ε1(t)';'ε2(t)';'ε3(t)';'ε4(t)';'ε5(t)'}; 
% row_names = {'ε*(t)','εᶜ(t)'} ; 
Pstar = array2table([ diag(Ptt(ss,ss)) diag(PtT(ss,ss)) ], ...
        'VariableNames',{'P(t|t)','P(t|T)'}, 'RowNames', row_names);
sep
print_table(Pstar,4,1,0)

%% USE SIMULATE DATA (FIX seed for random number generator)
TT = 1e4; rng(1);
[Zs, Xs, Us] = Kurz_simulate_SSF(D1, D2, R, A, C, dim_Z, dim_X, dim_R, TT);
Z = Zs; 

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
KS_deJ  = Kurz_AndersonMoore_Smoother(   D1, D2, A, Kurz_KF); % uses inv() --> changed to pinv(), needs initial values 
% --------------------------------------------------------------------------------------------------

%% PLOT THE KF/KS ESTIMATES OF THE STATES
clf; tiledlayout(3,2,TileSpacing="compact",Padding="compact");

% UNCOMMENT TO PLOT TRUE VS ESTIMATED
for k = 6:10
  nexttile
  hold on;
    plot(Xs(:,k), 'LineWidth',3); 
    % plot(KS_deJ.att(:,k),'--','Color',clr(3),'LineWidth',3); 
    plot(KS_deJ.atT(:,k),'--','Color',clr(3),'LineWidth',3); 
  hold off;
  box on; grid on;
  set(gca,'GridLineStyle',':' ,'GridAlpha',1/3, 'LineWidth',5/5);
  add2yaxislabel;
  addlegend({'True','Estimate:$\,a_{t|T}$'},1)
  addsubtitle(['$\varepsilon_{' num2str(k-5) 't}$'])
end

% UNCOMMENT TO SHOW SCATTER PLOTS 
% for k = 6:10
%   nexttile
%   hold on;
%     scatter(Xs(:,k),KS_deJ.att(:,k)); 
%   hold off;
%   box on; grid on;
%   set(gca,'GridLineStyle',':' ,'GridAlpha',1/3, 'LineWidth',5/5);
%   add2yaxislabel;
%   addlegend({'True vs. Estimate'},1)
%   addsubtitle(row_names{k-5})
% end

%% --------------------------------------------------------------------------------------------------
% print the correlations
sep
disp( 'Correlation matrix of True and estimated smoothed States')
sep
lst(corr(Xs(:,6:end), KS_deJ.atT(:,6:end)), 2)

% Estimate the identities 
% ∆ETη5t = 0.107∆ETη3t − 0.028ETη4t, where eta(i) = sig(i)*eps(i)
n1 = s1*KS_deJ.atT(:,5+1);
n2 = s2*KS_deJ.atT(:,5+2);
n3 = s3*KS_deJ.atT(:,5+3);
n4 = s4*KS_deJ.atT(:,5+4);
n5 = s5*KS_deJ.atT(:,5+5);

% run ols regressions
ols_KS = ols(delta(n3), [delta(n5) n4], 1, {'∆ε5(t)','ε4(t)'});




% CALCULATE CORRELATION 
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

































%EOF