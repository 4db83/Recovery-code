%schmidt-grohe uribe model

%V3

%will be using the Nimark version of SSF to reduce states. Done by using
%the D1 and D2 matrices

%State-space form

%z = D1*x + D2*x(-1) + R*u
%x = M*x(-1) + C*u

%states, x, will be: dxm zm dx z dxr ytil pitil itil plus 8 shocks

%8 shocks: 1...5, y, pi, i

%SGU Model

%Gap variables are Phi. 
%Phi = y, pi, i (all tilde)

%Phi follow a VAR(1)
%Phi = B*Phi(-1) + CHAT*Xi

%Xi are AR(1). Summarise as:

%Xi = E*Xi(-1) + S*eps

%eps  = shocks 1...5 (i.e. excludes measurement error)

%% A. Parameters

%Measurement equation parameters
alpha=0.0;

delta=2;
%delta=8.3292;
%delta=4;

sigmay=sqrt(1.2304);
sigmapi=sqrt(.4862);
sigmai=sqrt(.3208);

%State Equation Parameters
B=[.2627, .0187, -.5031; 
.3129 , .3292 , -.1170; 
.2268 , -.0977 , .5048];

CHAT=[
-.0956 , 0 , -.2603 , 1 ,-.0051 ; 
-.4892 , 0 , .5632 , .8727 , .3651 ; 
1.3964 , 1.0 , -.0309 , .2579 , -.2184];

rho1=.2426;
rho2=.3298;
rho3=.2619;
rho4=.4254;
rho5=.3110;

sigma1=.4824;
sigma2=.6250;
sigma3=1.3624;
sigma4=1.0913;
sigma5=.4723;


%% B. State equations

%Form implied matrices

E = zeros(5,5);
E(1,1)=rho1;
E(2,2)=rho2;
E(3,3)=rho3;
E(4,4)=rho4;
E(5,5)=rho5;

S = zeros(5,5);
S(1,1) = sigma1;
S(2,2) = sigma2;
S(3,3) = sigma3;
S(4,4) = sigma4;
S(5,5) = sigma5;


%Substitute equation for Xi into Phi.

%Phi = B*Phi(-1) + CHAT*(E*Xi(-1) + S*eps) 

m=zeros(16,16);

m(1:5,1:5) = E;
m(6:8,1:5) = CHAT*E;

m(6:8,6:8) = B;

%c
c = zeros(16,8);

c(1:5,1:5) = S;
c(6:8,1:5) = CHAT*S;
c(9:16,1:8) = eye(8);


%% D. Measurement equations


%set up d1 and d2

d1=zeros(3,16);
d2=zeros(3,16);

%d1
%Does not include standard devations of measurement error. These are now in
%r
d1(1,6)=1;d1(1,3)=1;d1(1,5)=delta;
d1(2,7)=1;d1(2,1)=1;
d1(3,8)=1;d1(3,1)=1+alpha;d1(3,5)=1;


d2(1,6)=-1;
d2(2,7)=-1;
d2(3,8)=-1;

%make measurement errors not permanent

% d2(1,14)=-sqrt(1.2304);
% d2(2,15)=-sqrt(.4862);
% d2(3,16)=-sqrt(.3208);

r=zeros(3,8);
%now has standard deviations of the measurement errors in r, rather than d1

r(1,6) = sigmay;
r(2,7) = sigmapi;
r(3,8) = sigmai;

%define number states, number observables. number shocks, how far to get
%the steady state KF

nos=16;
nobs=3;
nshk=8;
nt=1000;

prec=zeros(nos,nos,nt);
krec=zeros(nos,nobs,nt);
frec=zeros(nobs,nobs,nt);
capp=d1*m+d2;
%lam=d1*c+r;
lam=d1*c+r;
%compute filter using Nimark

prec(:,:,1)=eye(nos);


i=2;
while i<=nt;
    frec(:,:,i)=capp*prec(:,:,i-1)*capp'+lam*lam';
    
   
    krec(:,:,i)=(m*prec(:,:,i-1)*capp'+c*c'*d1'+c*r')*pinv(frec(:,:,i));
    
    prec(:,:,i)=m*prec(:,:,i-1)*m'+c*c'-krec(:,:,i)*(capp*prec(:,:,i-1)*capp'+lam*lam')*krec(:,:,i)';
    
    i=i+1;
end;
%define gain etc at end. This should be SS
k=zeros(nos,nobs);
k=krec(:,:,nt);
p=zeros(nos,nos);
p=prec(:,:,nt);

 pbar=m*p*m'+c*c'-k*(capp*p*capp'+lam*lam')*k';
 perr=pbar-p;
 disp('ssf p error')
 disp(perr);
 disp(trace(perr));

f=zeros(nobs,nobs);
f=frec(:,:,nt);

a=m*(eye(nos)-p*capp'*inv(f)*capp);
b=m*p*capp'*inv(f);
%now smooth using Kurz
nt=zeros(nos,nos,1000);
finv=inv(f);
nt1=capp'*finv*capp;
nt(:,:,1)=nt1;
j=2;
lt=m-k*(d1*m+d2);

disp('d1, d2 and M matrices')
disp(d1);
disp(d2);
disp(m);

%number of baclward steps
nrep=3000;
while j<=nrep;

nt(:,:,j)=nt1+lt'*nt(:,:,j-1)*lt;

j=j+1;
end;
ntend=nt(:,:,nrep);
psm=p-p*ntend*p';
disp('k filter gain');
disp(k);
disp('pt and pT for shocks');
disp('filtered');

disp(p(9:16,9:16));
disp('smoothed');
 disp(psm(9:16,9:16));
 





