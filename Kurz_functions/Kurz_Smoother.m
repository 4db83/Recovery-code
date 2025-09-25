function KFS_out = Kurz_Smoother(D1, D2, R, A, Q, Kurz_KF, Smoother_type)
% function KFS_out = Kurz_Smoother(D1, D2, A, Kurz_KF, Smoother_type)
% --------------------------------------------------------------------------------------------------
% My Notation for Kurz State-Space Form (SSF): (Kurz notation: Q --> C).
% --------------------------------------------------------------------------------------------------
%   Observed: Z(t) = D1*X(t) + D2*X(t-1) + R*ε(t),    X(t) = latent States
%   State:    X(t) =  A*X(t-1)           + Q*ε(t),    ε(t) ~ MN(0,I)
% --------------------------------------------------------------------------------------------------

SetDefaultValue(8, 'Smoother_type', 1)
% DeJongKohnAnsley is default, NO INV, NO INITVALS FOR STATES

% set which smoother to use
if Smoother_type == 1
  KFS_out = Kurz_DeJongKohnAnsley_Smoother(D1, D2, A, Kurz_KF); 
elseif Smoother_type == 2 
  KFS_out = Kurz_AndersonMoore_Smoother(D1, D2, A, Kurz_KF);
elseif Smoother_type == 3 % does not return PtT
  KFS_out = Kurz_Koopman_Smoother(D1, D2, R, A, Q, Kurz_KF);
else 
  error(' Only Choices 1, 2, 3 apply ')
end


% Merge structures using struct function
KFS_out = catstruct(KFS_out, Kurz_KF);
