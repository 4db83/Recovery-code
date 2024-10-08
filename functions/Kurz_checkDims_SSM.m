function [dimObs, dimState, dimDisturbance] = Kurz_checkDims_SSM(D1, D2, A, C, R)
% function [dimObs, dimState, dimDisturbance] = checkDimsModifiedSSM(D1, D2, A, C, R)
% --------------------------------------------------------------------------------------------------
% Kurz State-Space Form (SSF):
% --------------------------------------------------------------------------------------------------
%   Observed: Z(t) = D1*X(t)  + D2*X(t-1) + Rε(t)
%   State:    X(t) = A*X(t-1)             + Cε(t), where   Var(ε(t)) = I.
% --------------------------------------------------------------------------------------------------
%CHECKDIMSSSMwLS Check and extract dimensions of SSMwLS

% check dimensions of observables and states
assert(isequal(size(D1), size(D2)))
assert(size(A,1) == size(A,2))

% extract dimensions
[dimObs, dimState] = size(D1);

if exist('C', 'var') && exist('R', 'var')
if not(isempty(C)) && not(isempty(R))
    % check and extract disturbance dimensions
    dimDisturbance = size(C,2);
    assert(size(C,2) == size(R,2))
elseif not(isempty(C)) && isempty(R)
    % extract disturbance dimensions
    dimDisturbance = size(C,2);
elseif not(isempty(R)) && isempty(C)
    % extract disturbance dimensions
    dimDisturbance = size(R,2);
else
    % nothing to do
end
end

end

% function [dimObs, dimState, dimDisturbance] = Kurz_checkDims_SSM(D1, D2, Phi, Q, R)
% % ------------------------------------------------------------------------------------------------
% % Kurz State-Space Form (SSF):
% % ------------------------------------------------------------------------------------------------
% %   Observed: Z(t) = D1*X(t)  + D2*X(t-1) + Rε(t)
% %   State:    X(t) = A*X(t-1)             + Cε(t), where   Var(ε(t)) = I.
% % ------------------------------------------------------------------------------------------------
% %CHECKDIMSSSMwLS Check and extract dimensions of SSMwLS
% 
% % check dimensions of observables and states
% assert(isequal(size(D1), size(D2)))
% assert(size(Phi,1) == size(Phi,2))
% 
% % extract dimensions
% [dimObs, dimState] = size(D1);
% 
% if exist('C', 'var') && exist('R', 'var')
% if not(isempty(Q)) && not(isempty(R))
%     % check and extract disturbance dimensions
%     dimDisturbance = size(Q,2);
%     assert(size(Q,2) == size(R,2))
% elseif not(isempty(Q)) && isempty(R)
%     % extract disturbance dimensions
%     dimDisturbance = size(Q,2);
% elseif not(isempty(R)) && isempty(Q)
%     % extract disturbance dimensions
%     dimDisturbance = size(R,2);
% else
%     % nothing to do
% end
% end
% 
% end