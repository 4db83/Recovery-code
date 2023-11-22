% Clark87 Shock recovery (https://www.newyorkfed.org/research/policy/rstar)
% SSF: ---------------------------------------------------------------------------------------------
%   Z(t) = D1*X(t) + D2*X(t-1) + R*ε(t),      X(t) = latent States
%   X(t) = A*X(t-1)            + C*ε(t),      ε(t) ~ MN(0,I)
% --------------------------------------------------------------------------------------------------
clear; clc;
% set plotting defaults
set(groot,'defaultLineLineWidth',2); set(groot,'defaultAxesFontSize',15)
set(groot,'defaultAxesXTickLabelRotationMode','manual')
set(groot,'defaultAxesFontName','Times New Roman')
addpath(genpath('./functions'))
addpath(genpath('./utility.Functions'))               % set path to db functions
% addpath(genpath('D:/matlab.tools/db.toolbox/db'))   % set path to db functions
% CALL: get_all_db_toolbox_function_calls.m from Directory of code to be shared

% Sample size and seed for random number generator in simulation
Ts = 1e4; rng(123);
PLOT_STATES = 0;
% --------------------------------------------------------------------------------------------------    
% parameters are from the model fitted to US-GDP data from 1947:Q1 to 2019:Q4 (which can be estimated with the code below)
% a1  =  1.51023433332729;   % a1   
% a2  = -0.56787952465929;   % a2   
% % standard deviations                                                                             
% s1  =  0.54396737899273;   % sigma(ystar)                                                                     
% s2  =  0.02093523340402;   % sigma(g)                                                                        
% s3  =  0.59796738263493;   % sigma(ytilde)                                                                         
% --------------------------------------------------------------------------------------------------    
% offset for the first k latent state variables before the shocks.
k = 2;
% DEFINE SSM INPUT MATRICES ------------------------------------------------------------------------
dim_Z = 1;       % rows Z(t)
dim_X = 4;       % rows X(t)
dim_R = 2;       % rows ε(t)
% --------------------------------------------------------------------------------------------------    
% Define D1
D1 = zeros(dim_Z,dim_X); 
D1(1,2) = 1;
% Define D2
D2 = zeros(dim_Z,dim_X); 
% D2(1,[2 5:8]) = [s2 -a2*s1 -a1*s2 -a2*s2 -s3];
% Define R
R  = zeros(dim_Z,dim_R);
R(1,1) = 1;
% --------------------------------------------------------------------------------------------------
% Define A
A = zeros(dim_X); 
A([1,2],1) = 1;
% Define C
C = zeros(dim_X,dim_R);
C(1,2) = 1; C(3,1) = 1; C(4,2) = 1;
% --------------------------------------------------------------------------------------------------
% SIMULATE DATA FROM THE MODEL --> compute 'theoretical' properites of states
[Zs, Xs, Us] = Kurz_simulate_SSF(D1, D2, R, A, C, dim_Z, dim_X, dim_R, Ts); 
% --------------------------------------------------------------------------------------------------
% INITIALIZE FILTER 
% Note: errors will always be N(0,1), but latent states may need more careful initialization.
a00 = zeros(dim_X, 1); P00 = eye(dim_X);
% Filter NOTE: this uses the simulated Zs from above
[~, Kurz_KF] = Kurz_Filter(Zs, D1, D2, R, A, C, a00, P00);
% Smoother Modified de Jong (1988, 1989) and Kohn and Ansley (1989) smoother (Eq. (4.11) in Kurz (2018))
KS_deJ0  = Kurz_DeJongKohnAnsley_Smoother(D1, D2, A, Kurz_KF); % NO INV, NO INITVALS FOR STATES
% --------------------------------------------------------------------------------------------------
% CALL TO THE KURZ_SSM FUNCTION 
P0 = Kurz_steadystate_P(D1, D2, R, A, C);
ss = k+1:4;
% make display names 
row_names = make_table_names('ε',1:dim_R,'(t)');          
Pstar = array2table([ diag(P0.tT(ss,ss)) diag(P0.tt(ss,ss)) ], ...
         'VariableNames',{'P(t|T)','P(t|t)'}, 'RowNames', row_names);
sep; print_table(Pstar,4,1,0)
sep; fprintf('           Head Filtered                         Tail Filtered \n'); sep
disp([head(KS_deJ0.att) tail(KS_deJ0.att)])
sep; fprintf('           Head Smoothed                         Tail Smoothed \n'); sep
disp([head(KS_deJ0.atT) tail(KS_deJ0.atT)])

%% **************************************************************************************************
% NOW REFORMULATE THE MODEL TO THE SHOCK RECOVERY SSM
% **************************************************************************************************
sep('*')
disp('              SHOCK RECOVERY SSM')
sep('*')
% offset for the first k latent state variables before the shocks.
k = 0;
% DEFINE SSM INPUT MATRICES ------------------------------------------------------------------------
dim_Z = 1;       % rows Z(t)
dim_X = 3;       % rows X(t)
dim_R = 2;       % rows ε(t)
% --------------------------------------------------------------------------------------------------
% Define D1
D1 = zeros(dim_Z,dim_X); 
D1(1,1) = 1; D1(1,2) =  0;
% Define D2
D2 = zeros(dim_Z,dim_X); 
D2(1,1) = -1; D2(1,2) = 1;
% Define R
R  = zeros(dim_Z,dim_R);
% --------------------------------------------------------------------------------------------------
% Define A
A = zeros(dim_X); 
A(3,1) = 1;
% Define C
C = zeros(dim_X,dim_R);
C(1,1) = 1; C(2,2) = 1;
% --------------------------------------------------------------------------------------------------
% SIMULATE DATA FROM THE MODEL --> compute 'theoretical' properites of states
[Zs, Xs, Us] = Kurz_simulate_SSF(D1, D2, R, A, C, dim_Z, dim_X, dim_R, Ts); 
% --------------------------------------------------------------------------------------------------
% INITIALIZE FILTER 
% Note: errors will always be N(0,1), but latent states may need more careful initialization.
a00 = zeros(dim_X, 1); P00 = eye(dim_X);
% Filter NOTE: this uses the simulated Zs from above
[~, Kurz_KF] = Kurz_Filter(Zs, D1, D2, R, A, C, a00, P00);
% Smoother Modified de Jong (1988, 1989) and Kohn and Ansley (1989) smoother (Eq. (4.11) in Kurz (2018))
KS_deJ  = Kurz_DeJongKohnAnsley_Smoother(D1, D2, A, Kurz_KF); % NO INV, NO INITVALS FOR STATES
% --------------------------------------------------------------------------------------------------
% CALL TO THE KURZ_SSM FUNCTION 
P = Kurz_steadystate_P(D1, D2, R, A, C);
ss = k+1:2;
% make display names 
row_names = make_table_names('ε',1:dim_R,'(t)');          % make display names 
Pstar = array2table([ diag(P.tT(ss,ss)) diag(P.tt(ss,ss)) ], ...
         'VariableNames',{'P(t|T)','P(t|t)'}, 'RowNames', row_names);
sep; print_table(Pstar,4,1,0)
sep; fprintf('           Head Filtered                         Tail Filtered \n'); sep
disp([head(KS_deJ.att) tail(KS_deJ.att)])
sep; fprintf('           Head Smoothed                         Tail Smoothed \n'); sep
disp([head(KS_deJ.atT) tail(KS_deJ.atT)])



























 
%EOF