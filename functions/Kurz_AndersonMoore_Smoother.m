function resStruct = Kurz_AndersonMoore_Smoother(D1, D2, A, Z_tilde, Finv, K, att, Ptt)
%MODIFIEDANDERSONMOORESMOOTHER Modified Anderson and Moore (1979) smoother for SSMwLS 
% Purpose
%        The function computes the modifed Anderson and Moore (1979)
%        smoother for State Space Models with Lagged State (SSMwLS) in the
%        measurement equation, which is derived in Kurz (2018, Eq. (4.3)).
%        The smoother complements the modified Kalman filter of Nimark
%        (2015). The variance of the smooth states is also computed (Eq. 
%        (4.6) in Kurz (2018)).
%
%
% Usage
%           resStruct = modifiedAndersonMooreSmoother(D1, D2, A, Z_tilde, Finv, K, atT, PtT)
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
%       Z_tilde  = Errors -- Output of modifiedFilter()
%       Finv     = Second term of the Kalman gain -- Output of modifiedFilter()
%       K        = Kalman gain -- Output of modifiedFilter()
%       atT    = Filtered states -- Output of modifiedFilter()
%       PtT    = Filtered variances -- Output of modifiedFilter()
%
%
% Outputs
%      resStruct  = A structure containing
%                       atT   Smooth states
%                       PtT   Smooth variances
%
%
% References
%      Anderson, B. D. O., and J. B. Moore. 1979. Optimal filtering. Ed.
%         by T. Kailath. Prentice Hall information and system sciences
%         series. Englewood Cliffs, N.J.: Prentice-Hall.
%      Kurz, M. S. 2018. "A note on low-dimensional Kalman smoother for
%         systems with lagged states in the measurement equation".
%      Nimark, K. P. 2015. "A low dimensional Kalman filter for systems
%         with lagged states in the measurement equation". Economics
%         Letters 127: 10-13.
%
%
% Author: Malte S. Kurz



% check and extract dimensions
[dimObs, dimState] = Kurz_checkDims_SSM(D1, D2, A);
assert(size(Z_tilde,2) == dimObs)
nObs = size(Z_tilde,1);


D_tilde = (D1*A +D2);


% intialize struct for the results
resStruct = struct();
resStruct.atT = nan(nObs, dimState);
resStruct.PtT = nan(dimState, dimState, nObs);

resStruct.atT(nObs, :)     = att(nObs,:);
resStruct.PtT(:, :, nObs)  = Ptt(:,:, nObs);

for iObs = nObs-1:-1:1
    a_filtered = att(iObs,:)';
    P_filtered = Ptt(:,:, iObs);
    
    a_filtered_tp1    = att(iObs+1,:)';
    P_filtered_tp1    = Ptt(:,:, iObs+1);
    Finv_tp1          = Finv(:,:, iObs+1);
    K_tp1             = K(:,:, iObs+1);
    Z_tilde_tp1       = Z_tilde(iObs+1,:)';
    
    % one-step ahead smoother
    xx = P_filtered * D_tilde';
    atTp1 = a_filtered + xx * Finv_tp1 * Z_tilde_tp1;
    U_t_tp1 = P_filtered - xx * Finv_tp1 * xx';
    
    % components for J
    % P_filtered_tp1_Inv = eye(size(P_filtered_tp1)) / P_filtered_tp1;
    P_filtered_tp1_Inv = pinv(P_filtered_tp1);
    
    xx = A * P_filtered - K_tp1 * D_tilde * P_filtered;
    PtTp1_given_tp1 = xx';
    
    J = PtTp1_given_tp1 * P_filtered_tp1_Inv;
    
    resStruct.atT(iObs, :) = atTp1 + ...
        J * (resStruct.atT(iObs + 1, :)' - a_filtered_tp1);
    resStruct.PtT(:,:,iObs) = U_t_tp1 + ...
        J * (resStruct.PtT(:,:,iObs + 1) - P_filtered_tp1) * J';
    
    
end

end
