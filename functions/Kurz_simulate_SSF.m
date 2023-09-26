function [Zt,Xt] = Kurz_simulate_SSF(D1,D2,R,A,C,dimZ,dimX,dimR,TT)
% Simulate from Kurz SSF
  Xt = zeros(dimX,TT);  % states initialized at 0;
  Zt = zeros(dimZ,TT); 
  Ut = randn(dimR,TT);  % random draw from U(t)
  
  for t = 2:TT
      Xt(:,t) =  A * Xt(:,t-1) + C  * Ut(:,t);
      Zt(:,t) = D1 * Xt(:,t)   + D2 * Xt(:,t-1) + R * Ut(:,t);
  end
  % return (TxK) matrices of Z(t) and X(t), dropping initial values X(1) and 
  Zt = Zt'; Xt = Xt';
end