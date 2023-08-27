%MR recover
%This version estimates the identities in the paper.
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

alp1=1.53;
alp2=-.54;
by=-.32;
b1=.62;
alpr=-.05;


s1=.38;
s2=.79;
s3=.22;
s4=.54;
s5=.05;
s6=.15;
s7=.07;

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

%states will be y* y*(-1) y*(-2)  r* r*(-1) u* g eps1 eps2 eps3 eps4 eps5 eps6 eps7 eps

%set up d1 and d2

d1=zeros(3,15);
d2=zeros(3,15);

%curent 

d1(1,1)=1;
d1(1,9)=s1;

%lag 1 
d1(1,2)=-1.53;
d1(1,5)=.05/2;


%lag 2 
d2(1,2)=.54;
d2(1,5)=.05/2;


%2nd obs equation

% current 
d1(2,10)=s2;

%lag 1 shocks 

d2(2,6)=.32;

%3rd obs equation
%current

d1(3,6)=1;
d1(3,1)=-(.62*.4);
d1(3,15)=s7;
%lag 1
d1(3,2)=-(.62*.3);
%lag 2
d1(3,3)=-(.62*.2);
%lag 3
d2(3,3)=-(.62*.1);




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


% prec=zeros(nos,nos,nt);
% krec=zeros(nos,nobs,nt);
% frec=zeros(nobs,nobs,nt);
% capp=d1*m+d2;
% %lam=d1*c+r;
% lam=d1*c+r;
% %compute filter using Nimark
% 
% prec(:,:,1)=eye(nos);
% 
% 
% i=2;
% while i<=nt;
%     frec(:,:,i)=capp*prec(:,:,i-1)*capp'+lam*lam';
%     
%    
%     krec(:,:,i)=(m*prec(:,:,i-1)*capp'+c*c'*d1'+c*r')*pinv(frec(:,:,i));
%     
%     prec(:,:,i)=m*prec(:,:,i-1)*m'+c*c'-krec(:,:,i)*(capp*prec(:,:,i-1)*capp'+lam*lam')*krec(:,:,i)';
%     
%     i=i+1;
% end;
% %define gain etc at end. This should be SS
% k=zeros(nos,nobs);
% k=krec(:,:,nt);
% p=zeros(nos,nos);
% p=prec(:,:,nt);
% 
%  pbar=m*p*m'+c*c'-k*(capp*p*capp'+lam*lam')*k';
%  perr=pbar-p;
%  disp('ssf p error')
%  disp(p);
%  disp(trace(perr));
% 
% f=zeros(nobs,nobs);
% f=frec(:,:,nt);
% 
% a=m*(eye(nos)-p*capp'*inv(f)*capp);
% b=m*p*capp'*inv(f);
% %now smooth using Kurz
% nt=zeros(nos,nos,1000);
% finv=inv(f);
% nt1=capp'*finv*capp;
% nt(:,:,1)=nt1;
% j=2;
% lt=m-k*(d1*m+d2);
% 
% disp('d1, d2 and M matrices')
% disp(d1);
% disp(d2);
% disp(m);
% 
% %number of baclward steps
% nrep=3000;
% while j<=nrep;
% 
% nt(:,:,j)=nt1+lt'*nt(:,:,j-1)*lt;
% 
% j=j+1;
% end;
% ntend=nt(:,:,nrep);
% psm=p-p*ntend*p';
% disp('k filter gain');
% disp(k);
% disp('pt and pT for shocks');
% disp('filtered');
% disp(p(9:15,9:15));
% disp('smoothed');
%  disp(psm(9:15,9:15));
 

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

%% plot the KF/KS estimates of the states
%clf;
%for ii = 1:dimX
%  subplot(4,3,ii);
%  hold on;
%  plot(Xs(:,ii))
%  plot(KS_AM.a_t_T(:,ii), '--')
%  plot(KS_deJ.a_t_T(:,ii), '--')
%  plot(KS_K.a_t_T(:,ii), '--')
%  plot(Kurz_KF.a_t_t(:,ii), '--')
%  hold off;
%  box on; grid on;
%end

% check some output from the smoothers
%ii = 11
%head([Kurz_KF.a_t_t(:,ii) KS_AM.a_t_T(:,ii) KS_deJ.a_t_T(:,ii) KS_K.a_t_T(:,ii)])
%tail([Kurz_KF.a_t_t(:,ii) KS_AM.a_t_T(:,ii) KS_deJ.a_t_T(:,ii) KS_K.a_t_T(:,ii)])

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

%CODE TO PAPER MAPPING OF SHOCKS
%STATE  CODE    PAPER
% 9     1       eps1     
% 10    2       eps7 
% 11    3       eps2 
% 12    4       eps4 
% 13    5       eps5
% 14    6       eps6 
% 15    7       eps3 

%The epsilons are defined using the names in the paper

%Create epsilons (as in paper)
% Use KS_K.a_t_T
epsilon1 = s1*KS_K.a_t_T(:,9);
epsilon2 = s2*KS_K.a_t_T(:,11);
epsilon3 = s3*KS_K.a_t_T(:,15);
epsilon4 = s4*KS_K.a_t_T(:,12);
epsilon5 = s5*KS_K.a_t_T(:,13);	
epsilon6 = s6*KS_K.a_t_T(:,14);
epsilon7 = s7*KS_K.a_t_T(:,10);

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

%ET ∆η3t = −.136ET η5t + .3ET η5t−1 (39)
%ET η2t = 86.7ET ∆η6t + 398ET ∆η7t. 

%Estimate equation

%Equation 1
%ET ∆η3t = −.136ET η5t + .3ET η5t−1 (39)

Xmat1 = [epsilon5,epsilon5lag];
Xmat1 = Xmat1(2:end,:);
[b1,bint1,r1,rint1,stats1] = regress(depsilon3,Xmat1);
disp('coefficients first identity');
b1
stats1

%Equation 2
%ET η2t = 86.7ET ∆η6t + 398ET ∆η7t. 

Xmat2 = [depsilon6,depsilon7];
[b2,bint2,r2,rint2,stats2] = regress(epsilon2(2:end),Xmat2);
disp('coefficients second identity');
b2
stats2

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








