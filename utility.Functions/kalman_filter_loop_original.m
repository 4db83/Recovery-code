% ---------------------------------------------------------------------------------------------------	
% MAIN KF LOOP FUNCTION (ONLY THIS PART NEEDS TO BE IN .MEX)
% ---------------------------------------------------------------------------------------------------	
% 		Observed:	y(t)		 = D(t) + M*alpha(t)			+ u(t);		Var[u(t)] = H.
% 		State:		alpha(t) = C(t) + Phi*alpha(t-1)	+ S*n(t);	Var[n(t)] = Q.
% ---------------------------------------------------------------------------------------------------	
function sumLL = kalman_filter_loop_original(y, D, M, H, C, Phi, SQS, att_1, Ptt_1, T, Ky) 
	log_2_pi = 1.837877066409345483560659472811235279722;
	sumLL = 0;
	for t = 1:T
		% FORECAST ERROR AND ITS MSE
		uhat	= y(:,t) - ( D(:,t) + M*att_1 );
		Ft		= M*Ptt_1*M' + H;
		% PRECOMPUTE PM_invF = P(t|t-1)*M'*inv(F(t))
		P_MT_invF = Ptt_1*M'/Ft;
		% UPDATING/FILTERING a(t|t) and P(t|t)	
		att		= att_1 + P_MT_invF*uhat;
		Ptt		= Ptt_1 - P_MT_invF*M*Ptt_1;
		% PREDICTION/FORECASTING a(t|t-1) and P(t|t-1)
		att_1	= Phi*att + C(:,t);
		Ptt_1	= Phi*Ptt*Phi' + SQS;
		% LOG-LIKELIHOOD 
		LL_t	= -0.5*( Ky*log_2_pi + log(det(Ft)) + uhat'/Ft*uhat );
		sumLL	= sumLL + LL_t; % return sum(LogLike) only
	end
% end of function
end



