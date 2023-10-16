function [PtT_out, Ptt_out] = Kurz_steadystate_P(D1, D2, R, A, C, TT)
% TT = default number of iterations, optional
dim_X = size(D1,2);  % rows X(t)

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
Ptt = zeros(dim_X,dim_X,TT); 
% initialize state vector MSE 
P00 = eye(dim_X);
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
NT = zeros(dim_X);
Nt = NT;
LL = (A-Kt*G);
FF = Ft;
GinvFG = G'/FF*G;
% storage matrix if needed
Ntt = zeros(dim_X,dim_X,TT); 
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
