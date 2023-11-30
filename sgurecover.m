%schmidt-grohe uribe model
%use values given bu uribe  
%write up was sguuribe.tex
%states will be dx zm dxm z dxr ytil pitil itil +8 shocks
%will be using the Nimark version of SSF to reduce states. Done by using
%the D1 and D2 matrices
%8 shocks
%state equations
m=zeros(16,16);
alpha=0.0;
delta=0;
m(1,1)=.2426;m(2,2)=.3298;m(3,3)=.2619;m(4,4)=.4254;m(5,5)=.311;
 m(6,6)=.2627;m(6,7)=.0187;m(6,8)=-.5031;m(6,1)=-.0956*.2426;m(6,3)=-.2603*.2619;m(6,4)=.4254;m(6,5)=-.0051*.311;
 m(7,7)=.3292;m(7,6)=.3129;m(7,8)=-.117;m(7,1)=-.4892*.2426;m(7,3)=.5632*.2619;m(7,4)=.8727*.4254;m(7,5)=.3651*.311;
 m(8,8)=.5048;m(8,6)=.2268;m(8,7)=.0977;m(8,1)=-1.3964*.2426;m(8,2)=.3298;m(8,3)=-.0309*.2619;m(8,4)=.2579*.4254;m(8,5)=-.2184*.311;
c=zeros(13,5);
c(1,1)=.4824;
c(2,2)=.6250;
c(3,3)=1.3624;
c(4,4)=1.0913;
c(5,5)=.4723;
c(6,1)=-.0956*.4824;
c(6,3)=-.2603*1.3624;
c(6,4)=1.0913;
c(6,5)=-.0051*.4723;
c(7,1)=-.4892*.4824;
c(7,3)=.5632*1.3624;
c(7,4)=.8727*1.0913;
c(7,5)=.3651*.4723;
c(8,1)=-1.3964*.4824;
c(8,2)=.6250;
c(8,3)=-.0309*1.3624;
c(8,4)=.2579*1.0913;
c(8,5)=-.2184*.4723;
c(9,1)=1;
c(10,2)=1;
c(11,3)=1;
c(12,4)=1;
c(13,5)=1;
c(14,6)=1;
c(15,7)=1;
c(16,8)=1;
%ons equations
% disp('c');
% disp(c);
delta=2;
%delta=8.3292;
%delta=4;
%need to map them into states
%formulas are below
%states will be dxm zm dx z dxr ytil pitil itil plus 8 shocks
%set up d1 and d2
d1=zeros(3,16);
d2=zeros(3,16);
d1(1,6)=1;d1(1,3)=1;d1(1,5)=delta;d1(1,14)=sqrt(1.2304);
d1(2,7)=1;d1(2,1)=1;d1(2,15)=sqrt(.4862);
d1(3,8)=1;d1(3,1)=1+alpha;d1(3,5)=1;d1(3,16)=sqrt(.3208);
d2(1,6)=-1;
d2(2,7)=-1;
d2(3,8)=-1;
%make measurement errors not permanent
% d2(1,14)=-sqrt(1.2304);
% d2(2,15)=-sqrt(.4862);
% d2(3,16)=-sqrt(.3208);
r=zeros(3,8);
 
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
 disp(p);
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
 
