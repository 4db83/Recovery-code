function rho = corr_theory(STDs, STDs_ET, PHIs)
% this function computes the theoretical correlation between the model implied shocks or states and
% the Kalman Filtered or Smoothed ones.
% INPUTS: STDs      = (kx1) vector of standard deviations from the theoretical model
%         STDs_ET   = (kx1) vector of the Kalman Smoothed (or filtered) stdev computed from the
%                     simulated SSF.
%         PHIs      = diag of Pstar.("P*(t|T)") or Pstar.("P*(t|t)").
% --------------------------------------------------------------------------------------------------

rho = 0.5.*(STDs_ET.^2 + STDs.^2 - PHIs )./(STDs_ET.*STDs);


end