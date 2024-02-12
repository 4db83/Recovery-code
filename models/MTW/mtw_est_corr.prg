'mtw_est.prg

'A. Create Workfile

wfcreate mtw u 1 50000

'B. Load data

import(type=txt) "C:\Users\timr\Dropbox\Neutral real rate\2023\Code summary\MTW\mtwsim.txt" ftype=ascii rectype=crlf skip=0 fieldtype=delimited delim=space colhead=0 eoltype=pad badfield=NA @freq U @id @date @destid @date @smpl @all

'C. rename the series
rename series01 dr
rename series02 drobs
rename series03 dx
rename series04 dxobs
rename series05 v1
rename series06 v2
rename series07 vms1
rename series08 vms2

'D. Estimate a VAR(1)
var varout.ls(noconst) 1 1 drobs dxobs
varout.results
varout.makeresids v1hat v2hat

matrix A = varout.@coefmat
A = @transpose(A)

matrix B = @identity(2) - A
matrix Binv = @inverse(B)

'E. Construct the true BN
series perm = v1 - v2
scalar permstd = @stdev(perm)

series perm_norm = perm/@stdev(perm)

'F. Construct the preliminary BN estimate
series permhat = Binv(1,1)*v1hat+Binv(1,2)*v2hat
scalar permhatstd = @stdev(permhat)

series permhat_norm = permhat/permhatstd

'G. Correlogram
freeze(correl_out) permhat.correl

'H. Fit an ARMA(1,2) to permhat
equation eq_arma.ls permhat ar(1) ma(1) ma(2)
eq_arma.makeresids omegahat
vector bb = eq_arma.@coefs

'I. Robust perm component
scalar permrobustcoef = (1+ bb(2) + bb(3))/(1-bb(1))
series permrobust = permrobustcoef*omegahat
scalar permrobuststd = @stdev(permrobust)
scalar omeghatstd = @stdev(omegahat)

series permrobust_norm = permrobust/permrobuststd

'J. Recovery R2
equation eq_recovery.ls perm permrobust
equation eq_recovery_prelim.ls perm permhat

'Normalised versions (for correl)

equation eq_recovery_norm.ls perm_norm c permrobust_norm
equation eq_recorvery_prelim_norm.ls perm_norm c permhat_norm 


'K. Omegahat equations
'K.1. Full
series dvms1 = d(vms1)
series dvms2 = d(vms2)
equation eq_omegahat_full.ls omegahat v1(0to-10) v2(0to-10) dvms1(0to-10) vms2(0to-10)

equation eq_omegahat_full_alt.ls omegahat v1(0to-10) v2(0to-10) vms1(0to-10) vms2(0to-10)

equation eq_omegahat_full_alt_d.ls omegahat v1(0to-10) v2(0to-10) dvms1(0to-10) dvms2(0to-10)

'K.2 No me

equation eq_omegahat_nome.ls omegahat v1(0to-10) v2(0to-10) 

'K.3 predict true BN
equation eq_true_BN.ls perm vms1(0to-10) vms2(0to-10)

'L. Workfile save
wfsave mtw_out_corr


