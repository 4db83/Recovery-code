function Pstar = Kurz_Pstar(D1, D2, R, A, Q, P00, eps0)
% function Pstar = Kurz_Pstar(D1, D2, R, A, Q, P00, eps0)
% --------------------------------------------------------------------------------------------------
% My Notation for Kurz State-Space Form (SSF): (Kurz notation: Q --> C).
% --------------------------------------------------------------------------------------------------
%   Observed: Z(t) = D1*X(t) + D2*X(t-1) + R*ε(t),    X(t) = latent States
%   State:    X(t) =  A*X(t-1)           + Q*ε(t),    ε(t) ~ MN(0,I)
% --------------------------------------------------------------------------------------------------
% Kurz_Pstar can be called without any data on the observed variable Z(t) in the SSF. 
% In Kurz_Filter, I have added a call to Kurz_Pstar to not have to call it afterwards again if this
% function is already called once. In Kurz_Filter, the P*(t|t) based estimate of a(t|t) is also
% computed so that the weigths can be ready of directly.
% --------------------------------------------------------------------------------------------------

dim_X = size(D1,2);  % rows X(t)

% INITIALIZE STATE VECTOR MSE 
% P00 = eye(dim_X);

% set default number of iterations if not supplied
if nargin < 6 || isempty(P00);  P00   = eye(dim_X); end
if nargin < 7 || isempty(eps0); eps0  = 1e-12;      end

% pre-compute some quantities
Lam = (D1*Q + R);
G   = (D1*A + D2);
QQ  = Q*Q';
QR  = Q*R';

% --------------------------------------------------------------------------------------------------
% FILTERING STEADY-STATE P(t|t)
% --------------------------------------------------------------------------------------------------
Pt_1t_1 = P00; % initialization
% WHILE LOOP FORWARD RECURSIONS (see p.43 [2.1] & [2.2] in Kurz (2018))
norm_dPt = 1; 
while norm_dPt > eps0
  Ft =  G*Pt_1t_1*G' + Lam*Lam';
  Kt = (A*Pt_1t_1*G' + QQ*D1' + QR) / Ft;
  Ptt_1 = A*Pt_1t_1*A' + QQ ;       % P(t|t-1) = ϕ*P(t-1|t-1)*ϕ' + QQ
  Ptt   = Ptt_1 - Kt*Ft*Kt';        % P(t|t)   = P(t|t-1) - K(t)*F(t)*K(t)'
  % CONVERGENCE CHECKING
  norm_dPt = norm(Ptt - Pt_1t_1);
  Pt_1t_1 = Ptt;                    % resetting: P(t-1|t-1) = P(t|t)
end
% CONVERGENCE CHECK
if norm_dPt <= eps0
  fprintf('Convergence to Steady-State P*(t|t): %d\n', norm_dPt)
  Pstar.tt_1  = Ptt_1;    % P*(t|t-1)
  Pstar.tt    = Ptt;      % P*(t|t)
else
  error('   Error: P(t|t) did not converge to steady-state P*(t|t) via recursions')
end

% --------------------------------------------------------------------------------------------------
% SMOOTHING STEADY-STATE P(t|t)
% --------------------------------------------------------------------------------------------------
NT = zeros(dim_X);
Nt = NT;
FF = G*Pstar.tt*G' + Lam*Lam';
KK = (A*Pstar.tt*G' + QQ*D1' + QR) / FF;
LL = (A - KK*G);
GinvFG = G'/FF*G;

% WHILE LOOP BACKWARD RECURSIONS (see p.44 [4.13] in Kurz (2018))
norm_dNt = 1; 
while norm_dNt > eps0
  % THIS USES THE STEADY-STATE Kt AND Ft VALUES FROM ABOVE
  % Nt = Gt'/Ft*Gt + (A-Kt*Gt)'*Nt*(A-Kt*Gt);
  Nt1 = GinvFG + LL'*Nt*LL;
  % CONVERGENCE CHECKING
  norm_dNt = norm(Nt1 - Nt);
  Nt  = Nt1;
end
% CONVERGENCE CHECK
fprintf('Convergence of Steady-State N(t):    %d\n', norm_dNt)
% P*(t|T) = smoothed steady-state MSE
Pstar.tT = Pstar.tt - Pstar.tt*Nt*Pstar.tt;


% RETURN THE FOLLOWING
% Pstar.tt_1  = Pstar_tt_1;
% Pstar.tt    = Pstar_tt;
% Pstar.tT    = Pstar_tT;
Pstar.norm_dPtt = norm_dPt;
Pstar.norm_dNt  = norm_dNt;
% also return FF/KK/LL if needed later
Pstar.FF = FF;
Pstar.KK = KK;
Pstar.LL = LL;




end












































