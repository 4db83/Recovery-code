function [negLogLike, resStruct] = Kurz_Filter(Z, D1, D2, R, Phi, Q, a00, P00)
% function [negLogLike, resStruct] = Kurz_Filter(Z, D1, D2, R, Phi, Q, a00, P00)
% --------------------------------------------------------------------------------------------------
% My Notation for Kurz State-Space Form (SSF): (Kurz notation: Phi --> A, Q --> Q).
% --------------------------------------------------------------------------------------------------
%   Observed: Z(t) = D1*X(t)  + D2*X(t-1) + Rε(t)
%   State:    X(t) =  ϕ*X(t-1)            + Qε(t), where   Var(ε(t)) = I. 
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
[dimObs, dimState] = Kurz_checkDims_SSM(D1, D2, Phi, Q, R);
assert(size(Z,2) == dimObs)
nObs = size(Z,1);

D_tilde = (D1*Phi + D2);
QQT = Q*Q';

% intialize struct for the results
resStruct         = struct();
resStruct.Z_tilde = nan(nObs, dimObs);
resStruct.att     = nan(nObs, dimState);
resStruct.Ptt     = nan(dimState, dimState, nObs);
resStruct.Pt1t    = nan(dimState, dimState, nObs);
resStruct.Finv    = nan(dimObs, dimObs, nObs);
resStruct.K       = nan(dimState, dimObs, nObs);
resStruct.U       = nan(dimState, dimObs, nObs);

% initialize filter
att = a00;
Ptt = P00;

negLogLike = 0;

% pre-compute eye and log(2*pi)/2
II = eye(dimObs);
log_2_pi_2 = 1.837877066409345483560659472811235279722 / 2;

for iObs = 1:nObs
    
    Z_tilde = Z(iObs, :)' - D_tilde*att;
    U = Phi * Ptt * D_tilde' + QQT * D1' + Q * R';
    F = D_tilde * Ptt * D_tilde' + (D1 * Q + R) * (D1 * Q + R)';
    
    Finv = II / F;
    % Finv = eye(size(F)) / F;
    % Finv = pinv(F);
    K = U * Finv;
    
    att = Phi * att + K * Z_tilde;
    Ptt = Phi * Ptt * Phi' + QQT - K*F*K';
    % P(t+1|t)
    Pt1t = Phi * Ptt * Phi' + QQT;
    
    resStruct.Z_tilde(iObs,:)   = Z_tilde;
    resStruct.att(iObs,:)       = att;
    resStruct.Ptt(:,:,iObs)     = Ptt;
    resStruct.Pt1t(:,:,iObs)    = Pt1t;
    resStruct.Finv(:,:,iObs)    = Finv;
    resStruct.K(:,:,iObs)       = K;
    resStruct.U(:,:,iObs)       = U;
    
    % negLogLike =  negLogLike + dimObs*log(2*pi)/2 + 0.5* (log(det(F)) + Z_tilde' * Finv * Z_tilde);
    negLogLike =  negLogLike + dimObs*log_2_pi_2 + 0.5* (log(det(F)) + Z_tilde' * Finv * Z_tilde);
    
end
% also return some input matrices
resStruct.a00 = a00;
resStruct.P00 = P00;
resStruct.D1  = D1;
resStruct.D2  = D2;
resStruct.R   = R;
resStruct.Phi = Phi;
resStruct.Q   = Q;
% lst(negLogLike)

end





































