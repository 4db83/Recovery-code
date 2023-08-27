% HLW(2023) recovery check. Parameter values are from
% Holston_Laubach_Williams_current_estimates.xlsx see also Table 1 in their new paper.
clear; clc;
addpath(genpath('./Kurz-SSMwLS-github'))
rng(123)
% SSF: ------------------------------------------------------------------------------------
%   Z(t) = D1*X(t) + D2*X(t-1) + R*u(t)
%   X(t) = A*X(t-1)            + C*u(t),      u(t) ~ MN(0,I)
% -----------------------------------------------------------------------------------------
% parameters rounded to 4 decimal points in the Excel file
ay1 =  1.5399;
ay2 = -0.5986;    
by =   0.0756;
ar =  -0.0679;
cc =   1.0000;
bpi =  0.6708; % not used
% phi = -0.0854; % not needed for Kurz SSF
% -----------------------------------------------------------------------------------------
% standard deviations %CODE %PAPER  
sig_ytild = 0.3338;   %1    1  
sig_pi    = 0.7862;   %2    2 
sig_ystar = 0.5739;   %3    4   
% HLW reort sig_4 annualized, ie. multiplied by 4, devide it by 4 to get quarterly rate
sig_g     = 0.1230/4; %4    5  
sig_z     = 0.1742;   %5    3  
% crisis dummies that amplifiy the variance
kappa_2020 = 9.0326;
kappa_2021 = 1.7908;
kappa_2022 = 1.6760;
%--> CHOOSE WHICH kappa VALUE, IE, CRISIS YEAR TO USE. IF kappa = 1, PRE-COVID
kappa = kappa_2020;
kappa = 1;

% -----------------------------------------------------------------------------------------
dimR = 5;       % rows eps(t)
dimZ = 2;       % rows Z(t)
%--> CHANGE THIS TO 11 TO ADD \Delta r*(t) ROW AT END
dimX = 11;      % rows X(t)

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
% Define C
C = [zeros(dimR); eye(dimR)];
C(1,3) = sig_ystar; C(3,4) = sig_g; C([4],[4 5]) = [cc*4*sig_g; sig_z];
% -----------------------------------------------------------------------------------------
%--> ADD \Delta r* ROW TO STATE VECTOR
%C([11],[4 5]) = [cc*4*sig_g; sig_z];

%Normalised d(nrr)
%don't normalise if doing the correlation of the change in the neutral rate
%stddnrr = (sig_z^2+(cc*4*sig_g)^2)^0.5;
stddnrr = 1;
C(11,4)=cc*4*sig_g/stddnrr;
C(11,5)=sig_z/stddnrr;


% -----------------------------------------------------------------------------------------
% CALL to the KURZ_SSM function
[PtT,Ptt] = Kurz_SSM(D1,D2,R,A,C);
disp('    Pt|T');
disp(diag(PtT(6:end,6:end)))
disp('    Pt|t');
disp(diag(Ptt(6:end,6:end)))

% SIMULATE DATA
% 
TT = 1e4;
[Zs, Xs] = sim_Kurz_SSF(D1, D2, R, A, C, dimZ, dimX, dimR, TT);

% --------------------------------------------------------------------------------------------------
% FUNCTIONS TAKEN FROM KURZ'S GITHUB PAGE, MILDLY MODIFIED TO NOT REQUIRE INITVAL COMP. 
% AND USE OF PINV IN AM SMOOTHER OTHERWISE NON-SINGULARITY ISSUES. PLEASE HAVE A LOOK.
% --------------------------------------------------------------------------------------------------
% FILTER
[~, Kurz_KF] = modifiedFilter(Zs, D1, D2, A, C, R);
% SMOOTHERS 
% Modified Anderson and Moore (1979) smoother (Eq. (4.3) in Kurz (2018))
KS_AM   = modifiedAndersonMooreSmoother(   D1, D2, A, Kurz_KF.Z_tilde, Kurz_KF.Finv, Kurz_KF.K, Kurz_KF.a_t_t, Kurz_KF.P_t_t);
% Modified de Jong (1988, 1989) and Kohn and Ansley (1989) smoother (Eq. (4.11) in Kurz (2018))
KS_deJ  = modifiedDeJongKohnAnsleySmoother(D1, D2, A, Kurz_KF.Z_tilde, Kurz_KF.Finv, Kurz_KF.K, Kurz_KF.a_t_t, Kurz_KF.P_t_t);
% Modified Koopman (1993) smoother (Eq. (4.14)-(4.15) in Kurz (2018))
KS_K    = modifiedKoopmanSmoother(D1, D2, A, C, R, Kurz_KF.Z_tilde, Kurz_KF.Finv, Kurz_KF.K);

%% plot the KF/KS estimates of the states
clf;
for ii = 1:dimX
  subplot(4,3,ii);
  hold on;
  plot(Xs(:,ii))
  plot(KS_AM.a_t_T(:,ii), '--')
  plot(KS_deJ.a_t_T(:,ii), '--')
  plot(KS_K.a_t_T(:,ii), '--')
  plot(Kurz_KF.a_t_t(:,ii), '--')
  hold off;
  box on; grid on;
end

% check some output from the smoothers
%ii = 11
%head([Kurz_KF.a_t_t(:,ii) KS_AM.a_t_T(:,ii) KS_deJ.a_t_T(:,ii) KS_K.a_t_T(:,ii)])
%tail([Kurz_KF.a_t_t(:,ii) KS_AM.a_t_T(:,ii) KS_deJ.a_t_T(:,ii) KS_K.a_t_T(:,ii)])

% -----------------------------------------------------------------------------------------
% CALCULATE CORRELATION 
phi = PtT(end,end);
% Theoretical stdev 
%estimate 1
KS_AM_drr = KS_AM.a_t_T(:,11);
disp('std dr* estimate 1');
std_KS_AM_drr = std(KS_AM_drr)

%estimate 2
KS_K_drr = KS_K.a_t_T(:,11);
disp('std dr* estimate 2');
std_KS_AM_drr = std(KS_AM_drr)

%estimate 3
KS_deJ_drr = KS_deJ.a_t_T(:,11);
disp('std dr* estimate 3');
std_KS_deJ_drr = std(KS_deJ_drr)


%stdev_estdrstar_sm =0.04222116394588856;
stdev_estdrstar_sm = std_KS_AM_drr;
% Theoretical stdev_rstar
stdev_rstar = sqrt((cc*4*sig_g)^2 + sig_z^2);
%correlation between smoothed dr* and actual 
disp('correlation smoothed dr* estimate and dr*');
rho = 0.5*(stdev_estdrstar_sm^2 - phi + stdev_rstar^2)/(stdev_rstar*stdev_estdrstar_sm);
disp(rho)
% -----------------------------------------------------------------------------------------

%Create epsilons (as in paper)
% Use KS_K.a_t_T
epsilon1 =  sig_ytild*KS_K.a_t_T(:,6);
epsilon2 = sig_pi*KS_K.a_t_T(:,7);
epsilon3 = sig_z*KS_K.a_t_T(:,10);
epsilon4 = sig_ystar*KS_K.a_t_T(:,8);
epsilon5 = 4*sig_g*KS_K.a_t_T(:,9);	%note *4

%Create differences
dnrr = KS_K.a_t_T(:,11);
ddnrr = diff(KS_K.a_t_T(:,11)); %'11 is dnrr
depsilon5 = diff(epsilon5); 
depsilon3 = diff(epsilon3);

%Estimate equation - second difference
%ddnrr depsilon5 depsilon3
Xmat1 = [depsilon5,depsilon3];
disp('coefficients first identity - second difference');
[param1,bint1,r1,rint1,stats1] = regress(ddnrr,Xmat1);
bint1
stats1

%Estimate equation - first difference
%ddnrr depsilon5 depsilon3
Xmat1a = [epsilon5,epsilon3];
disp('coefficients first identity - first difference');
[param1a,bint1a,r1a,rint1a,stats1a] = regress(dnrr,Xmat1a);
bint1a
stats1a

%Identity 2 Second diff in nrr against eps 1,2,4 and 1(-1) and 4(-1)
%ddnrr epsilon1 epsilon2 epsilon4 epsilon1(-1) epsilon4(-1)
epsilon1lag = lagmatrix(epsilon1,1);
epsilon4lag = lagmatrix(epsilon4,1);
Xmat2 = [epsilon1,epsilon2,epsilon4,epsilon1lag,epsilon4lag];
Xmat2 = Xmat2(2:end,:);
disp('coefficients second identity - second difference');
[param2,bint2,r2,rint2,stats2] = regress(ddnrr,Xmat2);
bint2
stats2

%Identity 2 Second diff in nrr against eps 1,2,4 and 1(-1) and 4(-1)
%ddnrr epsilon1 epsilon2 epsilon4 epsilon1(-1) epsilon4(-1)
epsilon1lag = lagmatrix(epsilon1,1);
epsilon4lag = lagmatrix(epsilon4,1);
dnrrlag = lagmatrix(dnrr,1);
Xmat2a = [dnrrlag,epsilon1,epsilon2,epsilon4,epsilon1lag,epsilon4lag];
Xmat2a = Xmat2a(2:end,:);
disp('coefficients second identity - first difference');
[param2a,bint2a,r2a,rint2a,stats2a] = regress(dnrr(2:end),Xmat2a);
bint2a
stats2a


%% FUNCTIONS 
% (1) P*(t|t) and P*(t|T) from Kurz SSM
function [PtT_out,Ptt_out] = Kurz_SSM(D1,D2,R,A,C,TT)
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

% (2) Simulate from Kurz SSF
function [Zt,Xt] = sim_Kurz_SSF(D1,D2,R,A,C,dimZ,dimX,dimR,TT)
  Xt = zeros(dimX,TT);  % states initialized at 0;
  Zt = zeros(dimZ,TT); 
  Ut = randn(dimR,TT);  % random draw from U(t)
  
  for t = 2:TT
      Xt(:,t) =  A * Xt(:,t-1) + C  * Ut(:,t);
      Zt(:,t) = D1 * Xt(:,t)   + D2 * Xt(:,t-1) + R * Ut(:,t);
  end
  % return (TxK) matrices of Z(t) and X(t), dropping initial values X(1) and 
  Zt = Zt'; Xt = Xt';
end

























%EOF