% HLW_after_Covid19(2023) recovery check
%This version adds the change in r* to the state (at the end), so as the correlation can
%be computed.

clear; clc;

%set seed for random number generator
rng(1);

% SSF: ------------------------------------------------------------------------------------
%   Z(t) = D1*X(t) + D2*X(t-1) + R*u(t)
%   X(t) = A*X(t-1)            + C*u(t),      u(t) ~ MN(0,I)
% -----------------------------------------------------------------------------------------
% pars
ay1 =  1.3872;
ay2 = -0.4507;    
by =   0.0733;
ar =  -0.0790;
cc =   1.1283;
% -----------------------------------------------------------------------------------------
% standard deviations
sig_ytild = 0.4516;     % sig_ytild = 0.33378693;
sig_pi    = 0.7873;     % sig_pi    = 0.78620285;
sig_ystar = 0.5000;     % sig_ystar = 0.57390968;
% sig_galready annualized, ie. multiplied by 4
sig_g     = 0.1453;     % sig_g     = 0.03073865; 
sig_z     = 0.1181;     % sig_z     = 0.17417263;

% crisis dummies that amplifiy the variance
kappa_2020 = 9.0326;
kappa_2021 = 1.7908;
kappa_2022 = 1.6760;

% choose which kappa value, ie, crisis year to use. if kappa = 1, pre-covid
%kappa = kappa_2021;
kappa = 1;

% -----------------------------------------------------------------------------------------
dimZ = 2;       % rows Z(t)
%dimX = 10;      % rows X(t)
%size of the state increased by 1 so as to include dr*
dimX = 11;
dimR = 5;       % rows eps(t)
% -----------------------------------------------------------------------------------------
% Define D1
D1 = zeros(dimZ,dimX); 
D1(1,1) = 1;  D1(1,2) = -ay1; D1(1,5) = -ar/2; D1(1,6) = kappa*sig_ytild;
D1(2,2) = by; D1(2,7) = kappa*sig_pi; 
% Define D2
D2 = zeros(dimZ,dimX); 
D2(1,2) = -ay2;  D2(1,5) = -ar/2; 
% Define R
R  = zeros(dimZ,dimR);
% -----------------------------------------------------------------------------------------
% Define A
A = zeros(dimX); 
A([1 2],1) = 1; A([1 3],3) = 1; A([4 5],4) = 1;

%% Define C
%C = [zeros(dimR); eye(dimR)];
%new - to include dr*
C = zeros(dimX,dimR);
C(dimR+1:2*dimR,1:dimR)=eye(dimR);
C(1,3) = sig_ystar; C(3,4) = sig_g; C([4],[4 5]) = [cc*sig_g; sig_z];
%dr*
C([11],[4 5]) = [cc*sig_g; sig_z];

% -----------------------------------------------------------------------------------------
% call to the KURZ_SSM function
[PtT,Ptt] = kurz_SSM(D1,D2,R,A,C);
disp('    Pt|T');
disp(diag(PtT(6:end,6:end)))
disp('    Pt|t');
disp(diag(Ptt(6:end,6:end)))

phi = PtT(11,11);

%% Calculate correlation
%Theoretical stdev = from EViews
stdev_estdrstar_sm =0.007599727673354755;

%Theoretical stdev_rstar
stdev_rstar = ((cc*sig_g)^2 + sig_z^2)^0.5;

%correlation between smoothed dr* and actual 
disp('correlation smoothed dr* estimate and dr*');
rho = 0.5*(stdev_estdrstar_sm^2 - phi + stdev_rstar^2)/(stdev_rstar*stdev_estdrstar_sm)


% % % COMPUTE THE CORRELATIONS
% % % Var of change in Smoothed r* from file Laubach_Williams_current_estimates.xlsx up to 2019:Q4
% % var_Drstar_LW03 = 0.000647776314276209;
% % % theoretical as per quations (c*s5)^2+s3^2
% % kappa   = ((cc*sig_g)^2 + sig_z^2)^0.5;
% % PtT_11  = PtT(end,end);
% % rho     = (var_Drstar_LW03 + kappa^2 - PtT_11)/(2*sqrt(var_Drstar_LW03)*kappa)

%% Simulate data for EViews
 nobs = 10000;
 dataout =modelsim(D1,D2,R,A,C,dimZ,dimX,dimR,nobs);
 writematrix(dataout,'HLWafterdata.csv') 


%% KURZ SSM FUNCTION
function [PtT_out,Ptt_out] = kurz_SSM(D1,D2,R,A,C,TT)
% TT = default number of iterations, optional
dimX = size(D1,2);  % rows X(t)
dimR = size(R,2);   % rows eps(t)

% set default number of iterations if not supplied
if nargin < 6 || isempty(TT); TT  = 1e4; end  

% how far back to go to check convergence
T0  = 1e2;    
% pre-compute some quantitities
Lam = (D1*C+R);
G   = (D1*A+D2);
CCT = C*C';
% -----------------------------------------------------------------------------------------
% FILTERING
% storage matrix if needed
Ptt = zeros(dimX,dimX,TT); 
% initialize state vector MSE 
P00 = eye(dimX);
Pt  = P00;  
% FORWARD RECURSIONS
for t = 1:TT
  Ft =  G*Pt*G' + Lam*Lam';
  Kt = (A*Pt*G' + CCT*D1' + C*R') / Ft;
  Pt =  A*Pt*A' + CCT - Kt*Ft*Kt';
  % STORE SOME PTT FOR CONVERGENCE CHECKING
  Ptt(:,:,t) = Pt;
end
% CONVERGENCE CHECK: compute the difference between Last entry and Last-T0 entry
Ptt_conv = Ptt(:,:,end) - Ptt(:,:,end-T0);
fprintf('Covergence of Steady-State P(t|t): %d\n', norm(Ptt_conv))
% -----------------------------------------------------------------------------------------
% SMOOTHING 
NT = zeros(dimX);
Nt = NT;
LL = (A-Kt*G);
FF = Ft;
GinvFG = G'/FF*G;
% storage matrix if needed
Ntt = zeros(dimX,dimX,TT); 
% BACKWARD RECURSIONS
for t = (TT-1):-1:1
  % THIS USES THE STEADY-STATE Kt AND Ft VALUES FROM ABOVE
  % Nt = G'/Ft*G + (A-Kt*G)'*Nt*(A-Kt*G);
  Nt = GinvFG + LL'*Nt*LL;
  % store if needed
  Ntt(:,:,t) = Nt;
end
% nt(:,:,j)=nt1+lt'*nt(:,:,j-1)*lt;
% CONVERGENCE CHECK: compute the difference between first entry and 10th entry
Ntt_conv = Ntt(:,:,1) - Ntt(:,:,T0);
fprintf('Covergence of Steady-State N(t):   %d\n', norm(Ntt_conv))

% P(t|T) = smoothed steady-state MSE
PtT = Pt - Pt*Nt*Pt; 

% RETURN THE FOLLOWING
Ptt_out = Pt;
PtT_out = PtT; 
end

%% Function: modelsim
%simulates observed data from the model
function [dataout] = modelsim(D1,D2,R,A,C,dimZ,dimX,dimR,nobs)

%declare the matrix to store the simulated data
dataout = zeros(dimX,dimZ);
%each row is one time period

%initialise the state, and its lag
X=zeros(dimX,1);
Xlag=zeros(dimX,1);

%simulate the shocks. Save them in UMAT
UMAT = mvnrnd(zeros(dimR,1),eye(dimR),nobs);

for i=1:nobs
    %   Z(t) = D1*X(t) + D2*X(t-1) + R*u(t)
    %   X(t) = A*X(t-1)            + C*u(t), u(t) ~ MN(0,I)
    
    %1. select the shocks
    u=UMAT(i,:)';
    %2. Simulate state
    X = A*Xlag + C*u;
    %3Simulate observed data
    Z = D1*X + D2*Xlag + R*u;
    dataout(i,:) = Z';
    %4. Update Xlag
    Xlag = X;
end

end






























%EOF