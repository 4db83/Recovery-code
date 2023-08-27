% LW(2003) recovery check
clear; clc;
addpath(genpath('./Kurz-SSMwLS-github'))

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
s3 = 0.323;     % s3     = 0.17417263;
s4 = 0.605;     % sig_ystar = 0.57390968;
s5 = 0.102;     % s5     = 0.03073865; %annualized - could /4 and *4 everwhere s5 appears
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
%C = [C; C(4,:)];

%Normalised d(nrr)
%don't normalize if doing the correlation of dr*
%stddnrr = (s3^2+(cc*s5)^2)^0.5;
stddnrr=1;
C = [C; zeros(1,dimR)];
C(end,3)=s3/stddnrr;
C(end,5)=s5/stddnrr;

%update dimX to take account of this
dimX = dimX + 1;

% -----------------------------------------------------------------------------------------
% call to the KURZ_SSM function
[PtT,Ptt] = kurz_SSM(D1,D2,R,A,C);
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


%% CALCULATE CORRELATION 
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
stdev_rstar = sqrt((cc*s5)^2 + s3^2);
%correlation between smoothed dr* and actual 
disp('correlation smoothed dr* estimate and dr*');
rho = 0.5*(stdev_estdrstar_sm^2 - phi + stdev_rstar^2)/(stdev_rstar*stdev_estdrstar_sm);
disp(rho)

%% Estimate the identity
%∆ET η5t = 0.107∆ET η3t − 0.028ET η4t.

%Create epsilons (as in paper)
% Use KS_K.a_t_T
epsilon1 = s1*KS_K.a_t_T(:,6);
epsilon2 = s2*KS_K.a_t_T(:,7);
epsilon3 = s3*KS_K.a_t_T(:,8);
epsilon4 = s4*KS_K.a_t_T(:,9);
epsilon5 = s5*KS_K.a_t_T(:,10);	

%Create differences
depsilon3 = diff(epsilon3);
depsilon5 = diff(epsilon5); 

%Estimate equation
%ddnrr depsilon5 depsilon3
Xmat1 = [depsilon3,epsilon4(2:end)];
disp('coefficients first identity');
b1 = regress(depsilon5,Xmat1)

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