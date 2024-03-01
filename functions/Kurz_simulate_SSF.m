function [Zt,Xt,Ut] = Kurz_simulate_SSF(D1, D2, R, Phi, Q, T, X00, BurnIn)
% function [Zt,Xt,Ut] = Kurz_simulate_SSF(D1, D2, R, Phi, Q, T, X00, BurnIn)
% --------------------------------------------------------------------------------------------------
% My Notation for Kurz State-Space Form (SSF): (Kurz notation: Phi --> A, Q --> Q).
% --------------------------------------------------------------------------------------------------
%   Observed: Z(t) = D1*X(t)  + D2*X(t-1) + Rε(t)
%   State:    X(t) =  ϕ*X(t-1)            + Qε(t), where   Var(ε(t)) = I. 
% --------------------------------------------------------------------------------------------------
%     - Z (KyxT) column vector of observed/measurment data
%     - T Number of draws to generate
%     - X00 initial condition of the state vector (default is 0) 
%     - BurnIn number of burn-in observations to reduce the impact of initialisation (default is 100)
% --------------------------------------------------------------------------------------------------

  dim_Z = size(D1,1);
  dim_X = size(Phi ,1);
  dim_R = size(R ,2);

  %Set Default BurnIn = 100, if not supplied
  if (nargin < 8 || isempty(BurnIn)); BurnIn = 100; end
  
  TT = T + BurnIn;
  % Simulate from Kurz SSF
  Xt = zeros(dim_X, TT);  % states initialized at 0;
  Zt = zeros(dim_Z, TT);  % observed/measurement variable
  Ut = randn(dim_R, TT);  % random draw from U(t)
  
  % initialize state vector at X00 state initial value (a00 in main file notation)
  % if supplied, otherwise, default to all zero values
  if ~(nargin < 7 || isempty(X00)); Xt(:,1) = X00; end

  for t = 2:TT
      Xt(:,t) = Phi * Xt(:,t-1) +  Q * Ut(:,t);
      Zt(:,t) =  D1 * Xt(:,t)   + D2 * Xt(:,t-1) + R * Ut(:,t);
  end
  % return (TxK) matrices of Z(t) and X(t), dropping initial values X(1) and 
  Zt = Zt(:,BurnIn+1:end)'; Xt = Xt(:,BurnIn+1:end)';
  % also return the U error terms if needed to be used somewhere else
  Ut = Ut(:,BurnIn+1:end)';
end