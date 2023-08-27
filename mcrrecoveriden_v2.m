%MR recover
%This version estimates the identities in the paper.
%This version assumes output growth is the observed
%This verison also calculates the correlation
% 
% 
clear; clc;
addpath(genpath('./Kurz-SSMwLS-github'))

%set seed for random number generator
rng(1);

% SSF: ------------------------------------------------------------------------------------
%   Z(t) = D1*X(t) + D2*X(t-1) + R*u(t)
%   X(t) = A*X(t-1)            + C*u(t),      u(t) ~ MN(0,I)
% -----------------------------------------------------------------------------------------
% pars
%now uses MEAN*****
alp1=1.48; %1.53; %IS lag alpha1
alp2=-0.54;%-.54; %IS lag2 alpha2
by=-.33;%-.32; %PC ugap beta2
b1=.64;%.62; %Okuns beta
alpr=-0.06;%-.05; %IS r ar

s1=0.37;%IS %.38;
s2=0.8;%PC %.79;
s3=0.34;%Unexplained %.22;
s4=0.55; %Trend output %.54;
s5=.05; %Trend growth
s6=.15; %NAIRU
s7=.07; %Unemployment

%states will be 
%y* y*(-1) y*(-2) r* r*(-1) u* u*(-1) g and 7 shocks
%will be using the Nimark version of SSF to reduce states. Done by using
%the D1 and D2 matrices


%state equations

m=zeros(15,15);

m(1,1)=1;
m(1,8)=1;
m(2,1)=1;
m(3,2)=1;

m(4,4)=1;

m(5,4)=1;
m(6,6)=1;
m(7,6)=1;
m(8,8)=1;

%Relabel
A=m;

%there are 15 states and 7 shocks. shocks are states 9-15 so need mapping 

c=zeros(15,7);
c(1,4)=s4;
c(1,5)=s5;
c(8,5)=s5;
c(4,5)=4*s5;
c(4,3)=s3;
c(6,6)=s6;
c(9,1)=1;
c(10,2)=1;
c(11,3)=1;
c(12,4)=1;
c(13,5)=1;
c(14,6)=1;
c(15,7)=1;


%ons equations


% disp('c');
% disp(c);

%observables are psi1, psi2

%need to map them into states
%formulas are below

%states will be 
% 1 y* 
% 2 y*(-1) 
% 3 y*(-2)  
% 4 r* 
% 5 r*(-1) 
% 6 u*
% 7 u*(-1)
% 8 g 
% 9 eps1 
% 10 eps2 
% 11 eps3 
% 12 eps4 
% 13 eps5 
% 14 eps6 
% 15 eps7 

%set up d1 and d2

d1=zeros(3,15);
d2=zeros(3,15);

%curent 

%d1(1,1)=1; %level of contemp potential output no longer appears
d1(1,8) = 1; %g now appears

d1(1,9)=s1;
d1(1,12)=s4; %s4*eps4 now appears

%lag 1 
%d1(1,2)=-alp1; %-1.53;
d1(1,2) = -(alp1-1); %To take into account lag when observable growth
d1(1,5)=-alpr/2;%.05/2;


%lag 2 
d2(1,2)=-alp2;%.54;
d2(1,5)=-alpr/2;%.05/2;



%2nd obs equation

% current 
d1(2,10)=s2;

%lag 1 shocks 

d2(2,6)=-by;%.32;

%3rd obs equation
%current

d1(3,6)=1;
%d1(3,1)=-(.62*.4);
d1(3,1)=-(b1*.4);
d1(3,15)=s7;
%lag 1
%d1(3,2)=-(.62*.3);
d1(3,2)=-(b1*.3);
%lag 2
%d1(3,3)=-(.62*.2);
d1(3,3)=-(b1*.2);
%lag 3
%d2(3,3)=-(.62*.1);
d2(3,3)=-(b1*.1);




r=zeros(3,7);

%Relabel
D1=d1;
D2=d2;
R=r;
C=c;
 

%define number states, number observables. number shocks, how far to get
%the steady state KF

nos=15;
nobs=3;
nshk=7;
nt=1000;

%relabel
dimZ = nobs;       % rows Z(t)
dimX = size(D1,2);      % rows X(t)
dimR = nshk;       % rows eps(t)



% *****************************************************************************************
% ADD EXTRA ROWS AND COLS FOR DELTA_R* TO D1 AND D2
% *****************************************************************************************
D1 = [D1 zeros(3,1)];
D2 = [D2 zeros(3,1)];
A0 = zeros(dimX+1,dimX+1); 
A0(1:dimX,1:dimX) = A;
A = A0;
%C = [C; C(4,:)];

%Normalised d(nrr)
%don't normalize if doing the correlation of dr*
%stddnrr = (s3^2+(cc*s5)^2)^0.5;
stddnrr=1;
C = [C; zeros(1,dimR)];
C(end,3)=s3/stddnrr;
C(end,5)=4*s5/stddnrr;

%update dimX to take account of this
dimX = dimX + 1;



% -----------------------------------------------------------------------------------------
% call to the KURZ_SSM function
[PtT,Ptt] = kurz_SSM(D1,D2,R,A,C);
disp('    Pt|T');
disp(diag(PtT(9:end,9:end)))
disp('    Pt|t');
disp(diag(Ptt(9:end,9:end)))

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


%% CALCULATE CORRELATION 
phi = PtT(end,end);
% Theoretical stdev 
%estimate 1
KS_AM_drr = KS_AM.a_t_T(:,end);
disp('std dr* estimate 1');
std_KS_AM_drr = std(KS_AM_drr)

%estimate 2
KS_K_drr = KS_K.a_t_T(:,end);
disp('std dr* estimate 2');
std_KS_AM_drr = std(KS_AM_drr)

%estimate 3
KS_deJ_drr = KS_deJ.a_t_T(:,end);
disp('std dr* estimate 3');
std_KS_deJ_drr = std(KS_deJ_drr)


%stdev_estdrstar_sm =0.04222116394588856;
stdev_estdrstar_sm = std_KS_AM_drr;
% Theoretical stdev_rstar
stdev_rstar = sqrt((4*s5)^2 + s3^2);
%correlation between smoothed dr* and actual 
disp('correlation smoothed dr* estimate and dr*');
rho = 0.5*(stdev_estdrstar_sm^2 - phi + stdev_rstar^2)/(stdev_rstar*stdev_estdrstar_sm);
disp(rho)


%% Estimate identities

%states will be
% 1 y* 
% 2 y*(-1) 
% 3 y*(-2)  
% 4 r* 
% 5 r*(-1) 
% 6 u* 
% 7 u*(-1) 
% 8 g 
% 9 eps1 
% 10 eps2 
% 11 eps3 
% 12 eps4 
% 13 eps5
% 14 eps6 
% 15 eps7 

%The epsilons are defined using the names in the paper

%Create epsilons (as in paper)
% Use KS_K.a_t_T
epsilon1 = s1*KS_K.a_t_T(:,9);
epsilon2 = s2*KS_K.a_t_T(:,10);
epsilon3 = s3*KS_K.a_t_T(:,11);
epsilon4 = s4*KS_K.a_t_T(:,12);
epsilon5 = s5*KS_K.a_t_T(:,13);	
epsilon6 = s6*KS_K.a_t_T(:,14);
epsilon7 = s7*KS_K.a_t_T(:,15);

%export epsilons to csv file
dataout=[epsilon1,epsilon2,epsilon3,epsilon4,epsilon5,epsilon6,epsilon7];
writematrix(dataout,'mcreesdata.csv')

%Create differences
depsilon1 = diff(epsilon1);
depsilon2 = diff(epsilon2);
depsilon3 = diff(epsilon3);
depsilon4 = diff(epsilon4); 
depsilon5 = diff(epsilon5);
depsilon6 = diff(epsilon6); 
depsilon7 = diff(epsilon7);

%Create lags
epsilon1lag = lagmatrix(epsilon1,1);
epsilon2lag = lagmatrix(epsilon2,1);
epsilon3lag = lagmatrix(epsilon3,1);
epsilon4lag = lagmatrix(epsilon4,1);
epsilon5lag = lagmatrix(epsilon5,1);
epsilon6lag = lagmatrix(epsilon6,1);
epsilon7lag = lagmatrix(epsilon7,1);

depsilon1lag = lagmatrix(depsilon1,1);
depsilon2lag = lagmatrix(depsilon2,1);
depsilon3lag = lagmatrix(depsilon3,1);
depsilon4lag = lagmatrix(depsilon4,1);
depsilon5lag = lagmatrix(depsilon5,1);
depsilon6lag = lagmatrix(depsilon6,1);
depsilon7lag = lagmatrix(depsilon7,1);

%Create lag 2
epsilon1lag2 = lagmatrix(epsilon1,2);
epsilon2lag2 = lagmatrix(epsilon2,2);
epsilon3lag2 = lagmatrix(epsilon3,2);
epsilon4lag2 = lagmatrix(epsilon4,2);
epsilon5lag2 = lagmatrix(epsilon5,2);
epsilon6lag2 = lagmatrix(epsilon6,2);
epsilon7lag2 = lagmatrix(epsilon7,2);

%Estimate equation


Xmat1 = [epsilon1lag(2:end),depsilon5lag,epsilon4lag2(2:end),epsilon7(2:end),epsilon7lag2(2:end)];
Xmat1 = Xmat1(1000:end,:);
[param1,bint1,r1,rint1,stats1] = regress(epsilon1(1001:end),Xmat1);
disp('coefficients first identity');
param1
stats1

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








