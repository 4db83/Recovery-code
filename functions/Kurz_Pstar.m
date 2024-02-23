function Pstar = Kurz_Pstar(D1, D2, R, A, C, eps0)
% TT = default number of iterations, optional
dim_X = size(D1,2);  % rows X(t)

% set default number of iterations if not supplied
if nargin < 6 || isempty(eps0); eps0 = 1e-15; end

% how far back to go to check convergence
T0  = 1e2;
% pre-compute some quantitities
Lam = (D1*C+R);
G   = (D1*A+D2);
CCT = C*C';
% initialize state vector MSE 
P00 = eye(dim_X);
% -----------------------------------------------------------------------------------------
% FILTERING
Ptt = P00;
% WHILE LOOP FORWARD RECURSIONS
norm_dPtt = 1; 
while norm_dPtt > eps0
  Ft  =  G*Ptt*G' + Lam*Lam';
  Kt  = (A*Ptt*G' + CCT*D1' + C*R') / Ft;
  Pt1 =  A*Ptt*A' + CCT - Kt*Ft*Kt';
  % CONVERGENCE CHECKING
  norm_dPtt = norm(Pt1 - Ptt);
  Ptt = Pt1;
end
% CONVERGENCE CHECK: compute the difference between Last entry and Last-T0 entry
Ptt_conv = norm_dPtt;
fprintf('Convergence of Steady-State P(t|t): %d\n', norm_dPtt)
% -----------------------------------------------------------------------------------------
% SMOOTHING 
NT = zeros(dim_X);
Nt = NT;
LL = (A-Kt*G);
FF = Ft;
GinvFG = G'/FF*G;
% WHILE LOOP BACKWARD RECURSIONS
norm_dNt = 1; 
while norm_dNt > eps0
  % THIS USES THE STEADY-STATE Kt AND Ft VALUES FROM ABOVE
  % Nt = G'/Ft*G + (A-Kt*G)'*Nt*(A-Kt*G);
  Nt1 = GinvFG + LL'*Nt*LL;
  % CONVERGENCE CHECKING
  norm_dNt = norm(Nt1 - Nt);
  Nt  = Nt1;
end
% CONVERGENCE CHECK: compute the difference between first entry and 10th entry
fprintf('Convergence of Steady-State N(t):   %d\n', norm_dNt)

% P(t|T) = smoothed steady-state MSE
PtT = Ptt - Ptt*Nt*Ptt;

% RETURN THE FOLLOWING
Pstar.tt = Ptt;
Pstar.tT = PtT;
Pstar.norm_dPtt = norm_dPtt;
Pstar.norm_dNt  = norm_dNt;
end











































