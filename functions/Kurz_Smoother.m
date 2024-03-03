function KFS_out = Kurz_Smoother(D1, D2, R, Phi, Q, Kurz_KF, Smoother_type)
% function KFS_out = Kurz_Smoother(D1, D2, Phi, Kurz_KF, Smoother_type)
% --------------------------------------------------------------------------------------------------
% My Notation for Kurz State-Space Form (SSF): (Kurz notation: Phi --> A, Q --> Q).
% --------------------------------------------------------------------------------------------------
%   Observed: Z(t) = D1*X(t)  + D2*X(t-1) + Rε(t)
%   State:    X(t) =  ϕ*X(t-1)            + Qε(t), where   Var(ε(t)) = I. 
% --------------------------------------------------------------------------------------------------

SetDefaultValue(8, 'Smoother_type', 1)
% DeJongKohnAnsley is default, NO INV, NO INITVALS FOR STATES

% set which smoother to use
if Smoother_type == 1
  KFS_out = Kurz_DeJongKohnAnsley_Smoother(D1, D2, Phi, Kurz_KF); 
elseif Smoother_type == 2 
  KFS_out = Kurz_AndersonMoore_Smoother(D1, D2, Phi, Kurz_KF);
elseif Smoother_type == 3 % does not return PtT
  KFS_out = Kurz_Koopman_Smoother(D1, D2, R, Phi, Q, Kurz_KF);
else 
  error(' Only Choices 1, 2, 3 apply ')
end

% add Filter output
KFS_out.KF = Kurz_KF;