

%MR recover
% 
% 
%now uses MEAN*****
alp1=1.48; %1.53;
alp2=-0.54;%-.54;
by=-.33;%-.32;
b1=.64;%.62;
alpr=-0.06;%-.05;


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
d1(1,2)=-alp1; %-1.53;
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
%d2(3,3)=-(b1*.1);




r=zeros(3,7);
 

%define number states, number observables. number shocks, how far to get
%the steady state KF

nos=15;
nobs=3;
nshk=7;
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
disp(p(9:15,9:15));
disp('smoothed');
 disp(psm(9:15,9:15));
 





