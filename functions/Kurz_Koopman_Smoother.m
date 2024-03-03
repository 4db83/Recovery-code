function resStruct = Kurz_Koopman_Smoother(D1, D2, R, Phi, Q, Kurz_KF)
% function resStruct = Kurz_Koopman_Smoother(D1, D2, R, Phi, Q, Kurz_KF)
% --------------------------------------------------------------------------------------------------
% My Notation for Kurz State-Space Form (SSF): (Kurz notation: Phi --> A, Q --> Q).
% --------------------------------------------------------------------------------------------------
%   Observed: Z(t) = D1*X(t)  + D2*X(t-1) + Rε(t)
%   State:    X(t) =  ϕ*X(t-1)            + Qε(t), where   Var(ε(t)) = I. 
% --------------------------------------------------------------------------------------------------
%MODIFIEDKOOPMANSMOOTHER Modified Koopman (1993) smoother for SSMwLS 
% Purpose
%        The function computes the modifed Koopman (1993) smoother for 
%        State Space Models with Lagged State (SSMwLS) in the measurement
%        equation, which is derived in Kurz (2018, Eq. (4.14)-(4.15)). The
%        smoother complements the modified Kalman filter of Nimark (2015).
%
%
% Usage
%           resStruct = modifiedKoopmanSmoother(D1, D2, A, C, R, Z_tilde, Finv, K)
%
%
% Model Equation
%       Measurement Equation
%           Z_t = D_1 X_t + D_2 X_t-1 + R u_t
%       State Equation
%           X_t = A X_t-1 + C u_t
%
% Inputs
%       D1       = (dimObs x dimState) matrix from the measurement equation
%       D2       = (dimObs x dimState) matrix from the measurement equation
%       A        = (dimState x dimState) matrix from the state equation
%       C        = (dimState x dimDisturbance) matrix from the state equation
%       R        = (dimObs x dimDisturbance) matrix from the measurement equation 
%       Z_tilde  = Errors -- Output of modifiedFilter()
%       Finv     = Second term of the Kalman gain -- Output of modifiedFilter()
%       K        = Kalman gain -- Output of modifiedFilter()
%
%
% Outputs
%      resStruct  = A structure containing
%                       a_t_T   Smooth states
%
%
% References
%      Koopman, S. J. 1993. "Disturbance Smoother for State Space Models".
%         Biometrika 80 (1): 117-126.
%      Kurz, M. S. 2018. "A note on low-dimensional Kalman smoother for
%         systems with lagged states in the measurement equation".
%      Nimark, K. P. 2015. "A low dimensional Kalman filter for systems
%         with lagged states in the measurement equation". Economics
%         Letters 127: 10-13.
%
%
% Author: Malte S. Kurz

% get the output from Kurz_KF function call
Z_tilde = Kurz_KF.Z_tilde; 
Finv    = Kurz_KF.Finv; 
K       = Kurz_KF.K; 
att     = Kurz_KF.att; 
Ptt     = Kurz_KF.Ptt;

% check and extract dimensions
[dimObs, dimState, dimDisturbance] = Kurz_checkDims_SSM(D1, D2, Phi, Q, R);
assert(size(Z_tilde,2) == dimObs)
nObs = size(Z_tilde,1);


D_tilde = (D1*Phi +D2);
D1CR    = (D1 * Q + R);

% intialize struct for the results
resStruct = struct();
resStruct.atT = nan(nObs, dimState);

% smooth disturbances
u_t_T = nan(nObs, dimDisturbance);

% initialize the smoother
r = zeros(dimState, 1);

% disturbance smoother (backward recursion)
for iObs = nObs:-1:1
    Finv_t     = Finv(:,:, iObs);
    Z_tilde_t  = Z_tilde(iObs,:)';
    K_t        = K(:,:, iObs);
    
    M = Q - K_t *D1CR;
    u_t_T(iObs, :) = D1CR' * Finv_t * Z_tilde_t + M' * r;
    
    if iObs == nObs
        L = 0;
    else
        L = Phi - K_t * D_tilde;
    end
    r = D_tilde' * Finv_t * Z_tilde_t + L' * r;
    
end

% initialization
% [a_0_0, P_0_0] = initializeSSM(A, C, dimState);
% initialize filter
% [a_t_t, P_t_t] = initializeSSM(A, C, dimState);
a00 = zeros(dimState, 1);
P00 = 1*eye(dimState, dimState);

atT = a00 + P00 * r;

% state smoother (forward recursion)
for iObs = 1:nObs
    atT = Phi * atT + Q * u_t_T(iObs, :)';
    resStruct.atT(iObs, :) = atT;
end

resStruct.att = att;
resStruct.Ptt = Ptt;

end
























