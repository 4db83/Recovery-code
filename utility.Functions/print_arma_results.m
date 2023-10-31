function [] = print_arma_results(armaout,HAC_,diagout,addon)

% CHECK IF FIRST INPUT IS A STRUCTURE THAT CONTAINS ALL THE FIELDS
% if nargin < 2
% 	diagout = armaout.diagnostics;
% 	addon		= armaout.addon;
% 	armaout = armaout.armaout;
% end

SetDefaultValue(2 ,'HAC_'		, 1);
SetDefaultValue(3 ,'diagout', armaout.diagnostics);
SetDefaultValue(4 ,'addon'	, armaout.diagnostics);


if HAC_ == 1
	outp = [armaout.pars armaout.se_HAC armaout.tstat_HAC armaout.pval_HAC];
else
	outp = [armaout.pars armaout.se armaout.tstat armaout.pval];
end
	
p0 = max(diagout.P);
q0 = max(diagout.Q);

% number of X regressors included
nX = diagout.nX;

% p = 1:p0
p = diagout.P';
% q = 1:q0
q = diagout.Q';

if HAC_ == 1
	in.cnames = strvcat('Estimate', 'HAC std.err', 'HAC t-stat', 'HAC p-value');
else 
	in.cnames = strvcat('Estimate', 'Std.err', 't-stat', 'p-value');
end
	
in.rnames = strvcat('Parameters','Constant', [repmat('AR(',length(p),1) ...
	num2str(p') repmat(')',length(p),1)], [repmat('MA(',length(q),1) num2str(q') repmat(')',length(q),1)]);

if diagout.C == 0
	in.rnames = strvcat('Parameters', [repmat('AR(',length(p),1) num2str(p') repmat(')',length(p),1)], ...
		[repmat('MA(',length(q),1) num2str(q') repmat(')',length(q),1)]);
end

if nX > 0
	in.rnames = strvcat(in.rnames, [repmat('X',nX,1) num2str((1:nX)')]);
end

in.fmt    = strvcat('%16.6f');
% in.width	= 100;

if isempty(q0)
	QQ0 = 0;
else 
	QQ0 = q0;
end

if isempty(p0)
	PP0 = 0;
else 
	PP0 = p0;
end

% fprintf('-------------------------------------------------------------------------------\n');
fprintf('===============================================================================\n');
if isfield(addon, 'dates_out')
	fprintf(' ARMA(%d,%d) Model estimation results.', [PP0 QQ0]);
	fprintf(' Full Sample:  ')
	fprintf([datestr(addon.dates_out(1),'yyyy:qq') '-' datestr(addon.dates_out(end),'yyyy:qq')])
	fprintf(' [T = %d].\n', diagout.T)
	fprintf(' Effective Sample without Pre-sample values:			 ')
	fprintf([datestr(addon.dates_out(1+addon.Max_pq),'yyyy:qq') '-' datestr(addon.dates_out(end),'yyyy:qq')])
	fprintf(' [T = %d].\n', diagout.adjT)
else
	fprintf(' ARMA(%d,%d) Model estimation results. \n', [PP0 QQ0]);
end

fprintf('===============================================================================\n');
		myprint(outp,in,[],[],[],79);
fprintf('-------------------------------------------------------------------------------\n');

% print some extra results 
% ev1s = strjust(num2str(struct2mat(addon.eviews1)',in.fmt),'right');
% ev2s = strjust(num2str(struct2mat(addon.eviews2)',in.fmt),'right');

ev1s = strjust(num2str(struct2mat(armaout.addon.eviews1)',in.fmt),'right');
ev2s = strjust(num2str(struct2mat(armaout.addon.eviews2)',in.fmt),'right');

nam1 = strvcat(	'R-squared',...
				'Rbar-squared',...
				'SE of regression',...
				'Sum of Squared Errors',...
				'Log-likelihood',...	 		
				'F-statistic',...		 	
				'Nobs (T)');

nam2 = 		strvcat('|     Mean(y)',...
					'|     Std. Deviation(y)',...
					'|     AIC',...
					'|     BIC',...
					'|     HQC',...
					'|     Pr(F-statistic)',...
					'|     DW-statistic');


for ii = 1:size(nam1,1)								
	fprintf(' %s :  %s	  '	, nam1(ii,:), ev1s(ii,:))								
	fprintf('%s :  %s \n'		, nam2(ii,:), ev2s(ii,:))								
end
fprintf('===============================================================================\n');

% Compute AR lag polinomial roots
p_max = max(diagout.P);
if isempty(p_max); p_max = 0; end

q_max = max(diagout.Q);
if isempty(q_max); q_max = 0; end

NP  = max(p_max,q_max);

% lagPols = cell(NP+1,4);
% fill with nans rather than zeros.
lagPols = num2cell(nan(NP+1,4));

if p_max~=0
	lagPols(1,1) = {'AR(p)'};
	lagPols(1,2) = {'Modulus'};
	lagPols(2:p_max+1,1)= num2cell(diagout.aL_roots);
	lagPols(2:p_max+1,2)= num2cell(abs(diagout.aL_roots));
end
if q_max~=0
	lagPols(1,3) = {'MA(q)'};
	lagPols(1,4) = {'Modulus'};
	lagPols(2:q_max+1,3)= num2cell(diagout.bL_roots);
	lagPols(2:q_max+1,4)= num2cell(abs(diagout.bL_roots));
end


fprintf(' NOTE: max(p,q) = %d Observations are ''lost'' due to setting pre-sample values \n', armaout.addon.Max_pq);
fprintf(' for y at y(1:p), for e(1:q) at 0. Sample size including these is T = %d.\n', diagout.T)

fprintf('-------------------------------------------------------------------------------\n');
fprintf(' Unconditional mean of the model: mu = %2.4f.	\n', armaout.mu)
fprintf('-------------------------------------------------------------------------------\n');
fprintf(' Lag-polynomial roots (if Modulus larger than 1 -> stationary/Invertible)\n');
fprintf('-------------------------------------------------------------------------------\n');
% fprintf('		   AR(p)			     MA(q)\n');
% fprintf('-------------------------------------------------------------------------------\n');
disp(lagPols)
fprintf(' ------------------------------------------------------------------------------\n');
fprintf('  Lippi and Reichlin (1992) persistance measure PSI(1) = %2.4f where PSI(L)\n', armaout.addon.psi_1);
fprintf('  is equal to MA(infinity), that is: PSI(L) = THETA(L)/PHI(L), where PHI(L)\n');
fprintf('  and THETA(L) are the AR(p) and MA(q) lag polynomials.				\n');
fprintf('-------------------------------------------------------------------------------\n');















%% Compute AR lag polinomial roots
% p_max = max(diagout.P);
% q_max = max(diagout.Q);
% 
% if p_max ~= 0;
% %	aL = [1 -armaout.pars(1+diagout.C:diagout.C+max(p_max))'];
% 	fprintf(' AR(p) Lag-polynomial roots (if bigger than 1 -> stationary)\n');
% %	aL_roots = roots(aL);
% 	% fprtint(diagout.aL_roots')
% 	fprintf('  %s \n', mat2str(diagout.aL_roots,4))
% end;
% fprintf('-------------------------------------------------------------------------------\n');
% if q_max ~= 0;
% %	bL = [1 armaout.pars((1+diagout.C+p_max):end)'];
% 	fprintf(' MA(q) Lag-polynomial roots (if bigger than 1 -> invertible)\n');
% %	bL_roots = roots(bL);
% 	% disp(diagout.bL_roots')
% 	fprintf('  %s \n', mat2str(diagout.bL_roots,4))
% end;
% 

% aL = [1 -armaout.pars(1+diagout.C:end-q0)'];
% bL = [1 armaout.pars(2:end-q0)'];
% for ii = 1:NP;
% 	if ii < p_max+1;
% 		%lagPols{ii,1} = num2str(diagout.aL_roots(ii),'right');
% 		lagPols{ii,1} = strjust(mat2str(diagout.aL_roots(ii),5),'center');
% 	end;
% 	if ii < q_max+1;
% 		% lagPols{ii,2} = mat2str(diagout.bL_roots(ii),6);
% 		lagPols{ii,2} = strjust(mat2str(diagout.bL_roots(ii),5),'center');
% 	end
% end;
