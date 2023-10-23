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
% set to 1 if wanting to add ∆r*(t) to State vector X(t)
ADD_Drstar = 1;

% DEFINE SSM INPUT MATRICES ------------------------------------------------------------------------
dim_Z = 2;                % rows Z(t)
dim_X = 10 + ADD_Drstar;  % rows X(t)
dim_R = 5;                % rows ε(t)
% --------------------------------------------------------------------------------------------------    % this is what you get when running their LW code                                          
% parameters (from the published paper, 'Baseline' in Table 1 on page 1065.                             % D:\_research\_current\stars\recovery-code\_notes\LW_replication\kalman.standard.errors.R 
a1 =  1.517;      % ay1   paper only gives the sum as 0.945                                             %  [1,] "1.56495477703249"    "a_1"                                                        
a2 =  -.572;      % ay2                                                                                 %  [2,] "-0.614537663213145"  "a_2"                                                        
a3 =  -.098;      % ar                                                                                  %  [3,] "-0.0614762410516162" "a_3"                                                        
b3 =   .043;      % by                                                                                  %  [4,] "0.560395442735451"   "b_1"                                                        
c  =  1.068;      % c                                                                                   %  [5,] "0.377932050558007"   "b_2"                                                        
% standard deviations                                                                                   %  [6,] "0.0440540347304658"  "b_3"                                                        
s1 =   .387;      % s(ytild)                                                                            %  [7,] "0.00228697742079232" "b_4"                                                        
s2 =   .731;      % s(pi)                                                                               %  [8,] "0.0379650747803849"  "b_5"                                                        
s3 =   .323;      % s(z)                                                                                %  [9,] "1.24070752805134"    "c"                                                          
s4 =   .605;      % s(ystar)                                                                            % [10,] "0.335593534431049"   "sigma_1"                                                    
s5 =   .102/4;    % s(g) reportd as annualized rate -->                                                 % [11,] "0.760431178504621"   "sigma_2"                                                    
% devided by 4 to express in quarterly rate (annualized later in the C Matrix)                          % [12,] "0.592270881950106"   "sigma_4"                                                    
% --------------------------------------------------------------------------------------------------    % over the sample.start <- c(1961,1) sample.end   <- c(2002,2) with data vintage 2018      
% Define D1
D1 = zeros(dim_Z,dim_X); 
D1(1,1) = 1;    D1(1,6) = s1;
D1(2,2) = -b3;  D1(2,7) = s2;
% Define D2
D2 = zeros(dim_Z,dim_X); 
D2(1,1:2) = [-a1 -a2]; D2(1,4:5) = -a3/2;
% Define R
R  = zeros(dim_Z,dim_R);
% --------------------------------------------------------------------------------------------------
% Define A
A = zeros(dim_X); 
A(1:2,1) = 1; A(4:5,4) = 1; A([1 3],3) = 1;
% Define C
C = zeros(dim_X,dim_R); C(dim_R+1:10,:) = eye(dim_R);
C(1,4) = s4; C(3,5) = s5; C(4,[3 5]) = [s3 4*c*s5];
if ADD_Drstar; C(end,[3 5]) = [s3 4*c*s5]; end
% --------------------------------------------------------------------------------------------------

% CALL TO THE KURZ_SSM FUNCTION --------------------------------------------------------------------
[PtT, Ptt] = Kurz_steadystate_P(D1, D2, R, A, C);
ss = 6:dim_X;
row_names = {'ε1(t)';'ε2(t)';'ε3(t)';'ε4(t)';'ε5(t)'}; 
if ADD_Drstar 
  row_names = {'ε1(t)';'ε2(t)';'ε3(t)';'ε4(t)';'ε5(t)';'∆r*(t)'}; 
end
% row_names = {'ε*(t)','εᶜ(t)'} ; 
Pstar = array2table([ diag(Ptt(ss,ss)) diag(PtT(ss,ss)) ], ...
        'VariableNames',{'P(t|t)','P(t|T)'}, 'RowNames', row_names);
sep
print_table(Pstar,4,1,0)

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

% % UNCOMMENT TO PLOT TRUE VS ESTIMATED STATES -----------------------------------------------------
% % PLOT THE KF/KS ESTIMATES OF THE STATES 
% clf; tiledlayout(3,2,TileSpacing="compact",Padding="compact");
% for k = 6:size(KS_deJ.atT,2)
%   nexttile
%   hold on;
%     plot(Xs(:,k), 'LineWidth',3); 
%     plot(KS_deJ.att(:,k),'--','Color',clr(3),'LineWidth',2);  % Filtered States
%     plot(KS_deJ.atT(:,k),'--','Color',clr(3),'LineWidth',2);    % Smoothed States
%   hold off;
%   hline(0)
%   box on; grid on;
%   set(gca,'GridLineStyle',':' ,'GridAlpha',1/3, 'LineWidth',5/5);
%   add2yaxislabel;
%   addlegend({'True','Estimate:$\,a_{t|T}$'},1)
%   if k == 11; addsubtitle(['$\Delta r^{\ast}_{t}$']); 
%   else;       addsubtitle(['$\varepsilon_{' num2str(k-5) 't}$']); end
% end
% % UNCOMMENT TO SHOW SCATTER PLOTS 
% for k = 6:size(KS_deJ.atT,2)
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
% --------------------------------------------------------------------------------------------------

% CORRELATIONS (can also read off directly from the corr_table below -------------------------------
% print simple correlations
corr_table = array2table( corr(Xs(:,6:end), KS_deJ.atT(:,6:end)), ...
             'RowNames', row_names, 'VariableNames',row_names);
print_table(corr_table,4,1,'Correlation matrix of True and estimated smoothed States');sep

% Define eta(i) = sig(i)*eps(i)
n1 = s1*KS_deJ.atT(:,5+1);
n2 = s2*KS_deJ.atT(:,5+2);
n3 = s3*KS_deJ.atT(:,5+3);
n4 = s4*KS_deJ.atT(:,5+4);
n5 = s5*KS_deJ.atT(:,5+5);

% Estimate the identities: ∆ETη5t = 0.107∆ETη3t − 0.028ETη4t by running ols regressions:
fprintf('\n');sep('=');fprintf('Identity (20). Dependent variable: ∆ETη5t \n')
ID1 = ols(delta(n5), [delta(n3) n4], 1, {'∆ETη3(t)','ETη4(t)'});

% CALCULATE CORRELATION between ∆r*(t) and ET∆r*(t) [ie., actual vs. estimate)
if ADD_Drstar 
  PHI = PtT(end,end);
  % Theoretical stdev 
  sig_KS_Drstar = std( KS_deJ.atT(:,end) );
  % Theoretical from model
  sig_Drstar = sqrt( (4*c*s5)^2 + s3^2 );
  % correlation between smoothed dr* and actual 
  rho = 0.5*(sig_Drstar^2 + sig_KS_Drstar^2 - PHI) / (sig_Drstar*sig_KS_Drstar);
  fprintf('Corr(ET(∆r*(t)),∆r*(t) from simulation) Analytical formulas: %4.4f\n', rho);
end

% OTHER IDENTITIES ---------------------------------------------------------------------------------
% ET∆r*(t) = ETη5t + ETη3t (21) --> this should be: ET∆r*(t) = 4*c*ETη5t + ETη3t
% make ET∆r*(t) from X(:,4)-X(:,5) or alternatively from KS_deJ.atT(:,11);
Drstar = KS_deJ.atT(:,4)-KS_deJ.atT(:,5);

fprintf('\n');sep('=');fprintf('Identity (21). Dependent variable: ET∆r*(t) \n')
ID2 = ols(Drstar, [4*c*n5 n3], 1, {'4*c*ETη5(t)','ETη3(t)'});
 
fprintf('\n');sep('=');fprintf('Identity (22). Dependent variable: ET∆r*(t) \n')
ID3 = ols(Drstar, [lag(Drstar) n1 n2 n4 lag(n1) lag(n4)], 1, ...
      {'ET∆r*(t-1)','ETη1(t)','ETη2(t)','ETη4(t)','ETη1(t-1)','ETη2(t-4)'});





































%EOF