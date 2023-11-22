function [Zt,Xt,Ut] = Kurz_simulate_SSF(D1, D2, R, A, C, T, X00, BurnIn)
  % function [Zt,Xt,Ut] = Kurz_simulate_SSF(D1, D2, R, A, C, dim_Z, dim_X, dim_R, TT)
  dim_Z = size(D1,1);
  dim_X = size(A ,1);
  dim_R = size(R ,2);

  %Set Default BurnIn = 100, if not supplied
  if (nargin < 8 || isempty(BurnIn)); BurnIn = 100; end
  
  TT = T + BurnIn;
  % Simulate from Kurz SSF
  Xt = zeros(dim_X, TT);  % states initialized at 0;
  Zt = zeros(dim_Z, TT);  % observed/measurement variable
  Ut = randn(dim_R, TT);  % random draw from U(t)
  
  % initialze state vector at X00 state initial value (a00 in main file notation)
  % if supplied, otherwise, default to all zero values
  if ~(nargin < 7 || isempty(X00)); Xt(:,1) = X00; end

  for t = 2:TT
      Xt(:,t) =  A * Xt(:,t-1) + C  * Ut(:,t);
      Zt(:,t) = D1 * Xt(:,t)   + D2 * Xt(:,t-1) + R * Ut(:,t);
  end
  % return (TxK) matrices of Z(t) and X(t), dropping initial values X(1) and 
  Zt = Zt(:,BurnIn+1:end)'; Xt = Xt(:,BurnIn+1:end)';
  % also return the U error terms if needed to be used somewhere else
  Ut = Ut(:,BurnIn+1:end)';
end