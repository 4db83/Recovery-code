clear; clc;

%Laubach Williams 2003 R Ec States paper. Baseline model
% 

% -----------------------------------------------------------------------------------------
% HLW Adrian's version with extra state variables (From Table 3 Buncic(2022) hlw.2e)
% -----------------------------------------------------------------------------------------
% a1 = 1.53991114;  a2 = -0.59855575; 
% ar = -0.06786964; by =  0.07859265; 
% sig_pi    = 0.78620285;
% sig_ytild = 0.33378693;
% sig_ystar = 0.57390968;
% sig_g     = 0.03073865;
% sig_z     = 0.17417263;

al=1.51;
a2=-.57;    % db: this should be -a2 not a1
by=.043;
cc=1.068;
ar=-.098;


s1=.387;    % sig_ytild = 0.33378693;
s2=.731;    % sig_pi    = 0.78620285;
s3=.323;    % sig_z     = 0.17417263;
s4=.605;    % sig_ystar = 0.57390968;
s5=.102;    % sig_g     = 0.03073865; need to be annualized, 0.102 is annualized already

%states will be y* y*(-1) g r* r*(-1) ep1 eps2 eps3 eps4 eps5 
%will be using the Nimark version of SSF to reduce states. Done by using
%the D1 and D2 matrices


%state equations

m=zeros(10,10);


m(1,1)=1;
m(1,3)=1;
m(3,3)=1;
m(4,4)=1;
m(2,1)=1;
m(5,4)=1;




c=zeros(10,5);
c(1,4)=s4;

c(3,5)=s5;


c(4,5)=cc*s5;
c(4,3)=s3;


c(6,1)=1;
c(7,2)=1;
c(8,3)=1;
c(9,4)=1;
c(10,5)=1;



%ons equations


% disp('c');
% disp(c);

%observables are psi1, psi2

%need to map them into states
%formulas are below



%set up d1 and d2

d1=zeros(2,10);
d2=zeros(2,10);

%curent shocks
%states will be y* y*(-1) g r* r*(-1) ep1 eps2 eps3 eps4 eps5 

d1(1,1)=1;
d1(1,2)=-al;
d1(1,6)=s1;
%lag 1 shocks
d1(1,5)=-ar/2;    
%lag 2 shocks 
d2(1,2)=-a2;      
d2(1,5)=-ar/2;    


%2nd obs equation

%current shocks 

d1(2,7)=s2;
%lag 1 shocks 

  
d1(2,2)=by;      

r=zeros(2,5);
 

%define number states, number observables. number shocks, how far to get
%the steady state KF

nos=10;
nobs=2;
nshk=5;
nt=3000;

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
    
   
    krec(:,:,i)=(m*prec(:,:,i-1)*capp'+c*c'*d1'+c*r')*inv(frec(:,:,i));
    
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
nrep=250;
vv=zeros(1,250);
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
disp(p(6:10,6:10));
disp('smoothed');
 disp(psm(6:10,6:10));
 
 ph=m-k*capp;
 disp('ph ');
 disp(ph);
 
 disp('k');
 disp(k);



