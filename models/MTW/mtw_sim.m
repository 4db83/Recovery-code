rng('default');
rng(1);
a=zeros(2,2);
a(1,2)=-.05;
a(2,2)=.95;
%work out bn shock weights from the var
g=inv(eye(2)-a);
disp(g);
%from var91) in observ
a(1,1)=-.043;
a(1,2)=-.049;
a(2,1)=-.112;
a(2,2)=.945;
g=inv(eye(2)-a);
disp(g);
%set cov matrix of var shocks
sigma=zeros(2,2);
sigma(1,1)=.1125;
sigma(1,2)=.1;
sigma(2,1)=.1;
sigma(2,2)=.1;
mu=zeros(1,2);
nn=50000;
%generate var shocks
v=mvnrnd(mu,sigma,nn);
%generate measurement error shocks
vms=mvnrnd(mu,.05*sigma,nn);
%compute dr dx and drobs
dr=zeros(nn,1);
dx=zeros(nn,1);
drobs=zeros(nn,1);
dxobs=zeros(nn,1);
i=2;
while i<=nn;
dr(i)=-.05*dx(i-1)+v(i,1);
dx(i)=.95*dx(i-1)+v(i,2);
drobs(i)=dr(i)+vms(i,1)-vms(i-1,1);
dxobs(i)=dx(i)+vms(i,2)-vms(i-1,2);
i=i+1;
end;
z=horzcat(dr,drobs,dx,dxobs,v(:,1:2),vms(:,1:2));
save('C:\Users\timr\Dropbox\Neutral real rate\2023\Code summary\MTW\mtwsim.txt','z','-ascii','-double');