var dy ddcpi di dxm zm dx dxr ytil pitil itil drstar z;
varexo eps1 eps2 eps3 eps4 eps5 eps6 eps7 eps8;

parameters rho1 rho2 rho3 rho4 rho5 s1 s2 s3 s4 s5 s6 s7 s8 sy spi si alpha delta b11 b12 b13 b21 b22 b23 b31 b32 b33 c11 c13 c15 c21 c23 c24 c25 c31 c33 c34 c35;

rho1=.2426;
rho2=.3298;
rho3=.2619;
rho4=.4254;
rho5=.3110 ;
b11=.2627;
b22=.3292;
b33=.5048;

b12=0.0187;
b13=-.5031;
b21=.3129;
b23=-.1170;
b31=.2268;
b32=-.0977;


c11= -.0956;
c13=-.2603;
 c15=-.0051;
 c21=-.4892;
 c23=.5632;
 c24 =.8727;
c25 =.3651;
c31 =1.3964;
c33=-.0309;
c34 =.2579;
c35=-.2184;


alpha=0;

//delta=8.3292;
delta=2;

s1=.4824  ;
s2=.6250 ;


s3=1.3624  ;  
s4=1.0913 ;  
s5=.4723    ;
s6  =sqrt(1.2304)  ;
s7 =sqrt(.4862)   ;
s8 =sqrt(.3208)  ;




model(linear);
dy=ytil-ytil(-1)+dx+delta*dxr+s6*eps6;

ddcpi=pitil-pitil(-1)+ dxm+s7*eps7;

di=itil-itil(-1)+(1+alpha)*dxr+s8*eps8;

dxm=rho1*dxm(-1)+s1*eps1;
zm=rho2*zm(-1)+s2*eps2;
dx=rho3*dx(-1)+s3*eps3;
z=rho4*z(-1)+s4*eps4;
dxr=rho5*dxr(-1)+s5*eps5;
ytil=b11*ytil(-1)+b12*pitil(-1)+b13*itil(-1)+c11*dxm+c13*dx+c15*dxr+z;
pitil=b21*ytil(-1)+b22*pitil(-1)+b23*itil(-1)+c21*dxm+c23*dx+c25*dxr+c24*z;
itil=b31*ytil(-1)+b32*pitil(-1)+b33*itil(-1)+c31*dxm+zm+c33*dx+c34*z+c35*dxr;

drstar=dxr+alpha*dxm;
end;



initval;
dxm =0;
zm =0;
dx=0;
 dxr=0;
 ytil=0;
 pitil=0;
 itil=0;
 drstar=0;
  z=0;
dy=0;
di=0;
ddcpi=0;

end;

shocks;
var eps1= 1;
var eps2= 1;

var eps3= 1;

var eps4= 1;

var eps5= 1;

var eps6= 1;

var eps7= 1;

var eps8= 1;


end;


stoch_simul(irf=40) dy ddcpi di drstar;

