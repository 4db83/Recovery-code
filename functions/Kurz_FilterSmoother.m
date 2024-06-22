function [negLogLike, Kurz_FS] = Kurz_FilterSmoother(Z, D1, D2, R, A, Q, a00, P00)
% function [negLogLike, resStruct] = Kurz_FilterSmoother(Z, D1, D2, R, A, Q, a00, P00)
% NOW ADDS THE SMOOTHER AND FILTER OUTPUT TO ONE STRUCTURE. 
% --------------------------------------------------------------------------------------------------
% My Notation for Kurz State-Space Form (SSF): (Kurz notation: Q --> C).
% --------------------------------------------------------------------------------------------------
%   Observed: Z(t) = D1*X(t) + D2*X(t-1) + R*ε(t),    X(t) = latent States
%   State:    X(t) =  A*X(t-1)           + Q*ε(t),    ε(t) ~ MN(0,I)
% --------------------------------------------------------------------------------------------------
% MODIFIEDFILTER Nimark's (2015) modified Kalman filter for SSMwLS
% Purpose
%        The function computes Nimark's (2015) modified Kalman filter for
%        State Space Models with Lagged State (SSMwLS) in the measurement
%        equation.
%
%
% Usage
%        For SSMwLS
%           [negLogLike, resStruct] = modifiedFilter(Z, D1, D2, A, C, R)
%        For a classical SSM, i.e., without lagged state
%           [negLogLike, resStruct] = modifiedFilter(Z, D1, zeros(size(D1)), A, C, R)
%
%
% Model Equation
%       Measurement Equation
%           Z_t = D_1 X_t + D_2 X_t-1 + R u_t
%       State Equation
%           X_t = A X_t-1 + C u_t
%
% Inputs
%       Z  = (nObs x dimObs) vector of observables
%       D1 = (dimObs x dimState) matrix from the measurement equation
%       D2 = (dimObs x dimState) matrix from the measurement equation
%       A  = (dimState x dimState) matrix from the state equation
%       C  = (dimState x dimDisturbance) matrix from the state equation
%       R  = (dimObs x dimDisturbance) matrix from the measurement equation 
%
%
% Outputs
%      negLogLike = The negative log-likelihood
%      resStruct  = A structure containing
%                       Z_tilde Errors
%                       att     Filtered states
%                       Ptt     Filtered variances
%                       Pt1t    One-step ahead predictors of the variances
%                       Finv    Second term of the Kalman gain
%                       K       Kalman gain
%                       U       First term of the Kalman gain
%
%
% References
%      Nimark, K. P. 2015. "A low dimensional Kalman filter for systems
%         with lagged states in the measurement equation". Economics
%         Letters 127: 10-13.
%
%
%
% Author: Malte S. Kurz

% check and extract dimensions
[dimObs, dimState] = Kurz_checkDims_SSM(D1, D2, A, Q, R);
assert(size(Z,2) == dimObs)
T = size(Z,1);

D_tilde = (D1*A + D2);
QQ = Q*Q';

% intialize struct for the results
resStruct         = struct();
resStruct.Z_tilde = nan(T, dimObs);
resStruct.att     = nan(T, dimState);
resStruct.Ptt     = nan(dimState, dimState, T);
resStruct.Pt1t    = nan(dimState, dimState, T);
resStruct.Finv    = nan(dimObs,   dimObs, T);
resStruct.K       = nan(dimState, dimObs, T);
resStruct.U       = nan(dimState, dimObs, T);

% initialize filter
att = a00;
Ptt = P00;
negLogLike = 0;

% pre-compute eye and log(2*pi)/2
II = eye(dimObs);
log_2_pi_2 = 1.837877066409345483560659472811235279722 / 2;

for t = 1:T
    Z_tilde = Z(t, :)' - D_tilde*att;
    U = A*Ptt*D_tilde' + QQ*D1' + Q*R';
    F = D_tilde*Ptt*D_tilde' + (D1*Q + R)*(D1*Q + R)';
    
    Finv = II / F;
    % Finv = eye(size(F)) / F;
    % Finv = pinv(F);
    K = U*Finv;
    
    att = A*att + K*Z_tilde;
    Ptt = A*Ptt*A' + QQ - K*F*K';
    % P(t+1|t)
    Pt1t = A*Ptt*A' + QQ;
    
    resStruct.Z_tilde(t,:)   = Z_tilde;
    resStruct.att(t,:)       = att;
    resStruct.Ptt(:,:,t)     = Ptt;
    resStruct.Pt1t(:,:,t)    = Pt1t;
    resStruct.Finv(:,:,t)    = Finv;
    resStruct.K(:,:,t)       = K;
    resStruct.U(:,:,t)       = U;
    
    % negLogLike =  negLogLike + dimObs*log(2*pi)/2 + 0.5* (log(det(F)) + Z_tilde' * Finv * Z_tilde);
    negLogLike =  negLogLike + dimObs*log_2_pi_2 + 0.5*(log(det(F)) + Z_tilde'*Finv*Z_tilde);
end
% also return some input matrices
resStruct.a00 = a00;
resStruct.P00 = P00;
resStruct.D1  = D1;
resStruct.D2  = D2;
resStruct.R   = R;
resStruct.A   = A;
resStruct.Q   = Q;
resStruct.Z   = Z;
resStruct.LogLike = -negLogLike;
% lst(negLogLike)

% call to Pstar function to get steady-state P
Pstar = Kurz_Pstar(D1, D2, R, A, Q, Ptt);
resStruct.Pstar = Pstar;
% Pstar.FF; 
% Pstar.KK;

% compute / double check
FF  = D_tilde*Ptt*D_tilde' + (D1*Q + R)*(D1*Q + R)'; 
KK  = ( A*Ptt*D_tilde' + QQ*D1' + Q*R' ) / FF ; 
PSI = (A - KK*D_tilde); 
% --------------------------------------------------------------------------------------------------
% COMPUTE FILTERED STATE ESTIMATES USING STEADY-STATE PP = P*(t|t) KK = PP*M'/( M*PP*M' + RR ); 
% --------------------------------------------------------------------------------------------------
att_Pstar = zeros(dimState,T);
% W = [];
W = nan(dimState, dimObs, T);
for t = 1:T
	% FORECAST ERROR AND ITS MSE
  if t == 1
    att_Pstar(:,t) = PSI*a00              + KK*Z(t, :)';
  else
    att_Pstar(:,t) = PSI*att_Pstar(:,t-1) + KK*Z(t, :)';
  end
  W(:,:,t) = PSI^t*KK;
end

resStruct.PSI = PSI;
resStruct.KK  = KK;
resStruct.FF  = FF;
resStruct.WW  = W;
% also return the Pstar output
resStruct.att_Pstar = att_Pstar';


% now run the Smoother as well and cat results to resStr  
Kurz_FS = Kurz_Smoother( D1, D2, R, A, Q, resStruct); 


% resStruct = catstruct(resStruct, Kurz_SM);



end





































