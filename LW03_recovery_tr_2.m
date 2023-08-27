% LW(2003) recovery check
clear; clc;

%set seed for random number generator
rng(1);

% SSF: ------------------------------------------------------------------------------------
%   Z(t) = D1*X(t) + D2*X(t-1) + R*u(t)
%   X(t) = A*X(t-1)            + C*u(t),      u(t) ~ MN(0,I)
% -----------------------------------------------------------------------------------------
% pars
a1 = 1.51;
a2 = -.57;    
by =  .043;
ar = -.098;
cc = 1.068;
% -----------------------------------------------------------------------------------------
% standard deviations
s1 = 0.387;     % sig_ytild = 0.33378693;
s2 = 0.731;     % sig_pi    = 0.78620285;
s3 = 0.323;     % sig_z     = 0.17417263;
s4 = 0.605;     % sig_ystar = 0.57390968;
s5 = 0.102;     % sig_g     = 0.03073865; 
% -----------------------------------------------------------------------------------------
dimZ = 2;       % rows Z(t)
dimX = 10;      % rows X(t)
dimR = 5;       % rows eps(t)
% -----------------------------------------------------------------------------------------
% Define D1
D1 = zeros(dimZ,dimX); 
D1(1,1) = 1;  D1(1,2) = -a1; D1(1,5) = -ar/2; D1(1,6) = s1;
D1(2,2) = by; D1(2,7) = s2; 
% Define D2
D2 = zeros(dimZ,dimX); 
D2(1,2) = -a2;  D2(1,5) = -ar/2; 
% Define R
R  = zeros(dimZ,dimR);
% -----------------------------------------------------------------------------------------
% Define A
A = zeros(dimX); 
A([1 2],1) = 1; A([1 3],3) = 1; A([4 5],4) = 1;
% Define C
C = [zeros(dimR); eye(dimR)];
C(1,4) = s4; C(4,3) = s3; C([3 4],5) = [s5; cc*s5];

% *****************************************************************************************
% ADD EXTRA ROWS AND COLS FOR DELTA_R* TO D1 AND D2
% *****************************************************************************************
D1 = [D1 zeros(2,1)];
D2 = [D2 zeros(2,1)];
A0 = zeros(11,11); 
A0(1:10,1:10) = A;
A = A0;
C = [C; C(4,:)];

%update dimX to take account of this
dimX = dimX + 1;

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
stdev_estdrstar_sm =0.03678712582719936;

%Theoretical stdev_rstar
stdev_rstar = ((cc*s5)^2+s3^2)^0.5;

%correlation between smoothed dr* and actual 
disp('correlation smoothed dr* estimate and dr*');
rho = 0.5*(stdev_estdrstar_sm^2 - phi + stdev_rstar^2)/(stdev_rstar*stdev_estdrstar_sm)

% COMPUTE THE CORRELATIONS
% Var of change in Smoothed r* from file Laubach_Williams_current_estimates.xlsx up to 2019:Q4
%var_Drstar_LW03 = 0.000647776314276209;
% theoretical as per quations (c*s5)^2+s3^2
%kappa   = ((cc*s5)^2 + s3^2)^0.5;
%PtT_11  = PtT(end,end);
%rho     = (var_Drstar_LW03 + kappa^2 - PtT_11)/(2*sqrt(var_Drstar_LW03)*kappa)

%% Simulate data for EViews
 nobs = 10000;
 dataout =modelsim(D1,D2,R,A,C,dimZ,dimX,dimR,nobs);
 writematrix(dataout,'LWdata.csv') 


%% KURZ FUNCTION
function [PtT_out,Ptt_out] = kurz_SSM(D1,D2,R,A,C,TT)
% TT = default number of iterations, optional
dimX = size(D1,2);  % rows X(t)
dimR = size(R,2);   % rows eps(t)

% set default number of iterations if not supplied
if nargin < 6 || isempty(TT); TT  = 5e3; end  

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