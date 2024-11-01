function resStruct = Kurz_DeJongKohnAnsley_Smoother(D1, D2, A, Kurz_KF)
% function resStruct = Kurz_DeJongKohnAnsley_Smoother(D1, D2, A, Kurz_KF)
% --------------------------------------------------------------------------------------------------
% My Notation for Kurz State-Space Form (SSF): (Kurz notation: Q --> C).
% --------------------------------------------------------------------------------------------------
%   Observed: Z(t) = D1*X(t) + D2*X(t-1) + R*ε(t),    X(t) = latent States
%   State:    X(t) =  A*X(t-1)           + Q*ε(t),    ε(t) ~ MN(0,I)
% --------------------------------------------------------------------------------------------------
% function resStruct = Kurz_DeJongKohnAnsley_Smoother(D1, D2, A, Z_tilde, Finv, K, att, Ptt)
% MODIFIEDDEJONGKOHNANSLEYSMOOTHER Modified de Jong (1988, 1989) and Kohn and Ansley (1989) smoother for SSMwLS 
% Purpose
%         The function computes the modifed de Jong (1988, 1989) and Kohn and Ansley (1989) 
%         smoother for State Space Models with Lagged State (SSMwLS) in the measurement 
%         equation, which is derived in Kurz (2018, Eq. (4.11)). The smoother complements the 
%         modified Kalman filter of Nimark (2015). The variance of the smooth states is also 
%         computed (Eq. (4.13) in Kurz (2018)).
%
%
% Usage
%           resStruct = modifiedDeJongKohnAnsleySmoother(D1, D2, A, Z_tilde, Finv, K, atT, PtT)
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
%       Z_tilde  = Errors -- Output of modifiedFilter()
%       Finv     = Second term of the Kalman gain -- Output of modifiedFilter()
%       K        = Kalman gain -- Output of modifiedFilter()
%       atT      = Filtered states -- Output of modifiedFilter()
%       PtT      = Filtered variances -- Output of modifiedFilter()
%
%
% Outputs
%      resStruct  = A structure containing
%                   - atT Smooth states
%                   - PtT Smooth variances
%
%
% References
%      de Jong, P. 1988. "A Cross-Validation Filter for Time Series Models". Biometrika 75 (3): 594-600.
%      de Jong, P. 1989. "Smoothing and Interpolation with the State-Space Model". Journal of the American Statistical Association 84 (408): 1085-1088.
%      Kohn, R., and C. F. Ansley. 1989. "A Fast Algorithm for Signal Extraction, Influence and Cross-Validation in State Space Models". Biometrika 76 (1): 65-79.
%      Kurz, M. S. 2018. "A note on low-dimensional Kalman smoother for systems with lagged states in the measurement equation".
%      Nimark, K. P. 2015. "A low dimensional Kalman filter for systems with lagged states in the measurement equation". Economics Letters 127: 10-13.
%
% Author: Malte S. Kurz

% get the output from Kurz_KF function call
Z_tilde = Kurz_KF.Z_tilde; 
Finv    = Kurz_KF.Finv; 
K       = Kurz_KF.K; 
att     = Kurz_KF.att; 
Ptt     = Kurz_KF.Ptt;

% check and extract dimensions
[dimObs, dimState] = Kurz_checkDims_SSM(D1, D2, A);
assert(size(Z_tilde,2) == dimObs)
nObs = size(Z_tilde,1);


D_tilde = (D1*A + D2);

% intialize struct for the results
resStruct = struct();
resStruct.atT = nan(nObs, dimState);
resStruct.PtT = nan(dimState, dimState, nObs);

% initialize the smoother
r = zeros(dimState, 1);
N = zeros(dimState, dimState);

% for iObs = nObs
for iObs = nObs:-1:1
    Finv_t     = Finv(:,:, iObs);
    Z_tilde_t  = Z_tilde(iObs,:)';
    K_t        = K(:,:, iObs);
    a_filtered = att(iObs,:)';
    P_filtered = Ptt(:,:, iObs);
    
    resStruct.atT(iObs,:)   = a_filtered + P_filtered * r;
    resStruct.PtT(:,:,iObs) = P_filtered - P_filtered * N * P_filtered;
    
    if iObs == nObs
        L = 0;
    else
        L = A - K_t * D_tilde;
    end
    r = D_tilde' * Finv_t * Z_tilde_t + L' * r;
    N = D_tilde' * Finv_t * D_tilde + L' * N * L;
end

% uncomment to also return the Filtered series
% resStruct.att = att;
% resStruct.Ptt = Ptt;


end
