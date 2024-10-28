function [] = print_fullols(fullols_output,variable_names,prnt_hac,prnt2xls)
% F: Prints full OLS results in Eviews Style to screen
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % clc;
% % %delete in
% % fullols_output = fullols(Y(1+h:end),X_us(1:end-h,:),[],[],[],1);

% width of line printing
WDTH = 133;

ols_out = fullols_output;
KK			= ols_out.k - ~ols_out.no_const;

default_var_names = cellstr([repmat('X',KK,1) num2str((1:KK)') ]);
SetDefaultValue(2, 'variable_names', default_var_names)
SetDefaultValue(3, 'prnt_hac'			,1)
SetDefaultValue(4, 'prnt2xls'			,0)

if prnt_hac == 1
	outp = [ols_out.bhat ols_out.hac_se ols_out.hac_tstat ols_out.hac_pval];
	Fstat				= ols_out.hac_Fstat;
	pval_Fstat	= ols_out.hac_pval_Fstat;
else
	outp = [ols_out.bhat ols_out.se ols_out.tstat ols_out.pval];
	Fstat				= ols_out.Fstat;
	pval_Fstat	= ols_out.pval_Fstat;
end

Nvn = size(variable_names,2);
if Nvn > 1
	variable_names = variable_names';
end

if length(variable_names) ~= KK
  variable_names = default_var_names;
end

if prnt_hac == 1
	in.cnames = {'Estimate'; 'stderr(HAC)'; 't-stat(HAC)'; 'p-value'};
else
	in.cnames = {'Estimate'; 'stderr'; 't-stat'; 'p-value'};
end

in.rnames = append(' ', [{'Variable';'Constant'}; variable_names]);
in.width	= WDTH;

if ols_out.no_const == 1
 	in.rnames = append(' ', ['Variable'; variable_names]);
end

% if no regressors only a constant is in the regression, remove variable_names.
if KK == 0 
	in.rnames = [{'Variable';'Constant'}];
end

in.fmt    = strvcat('%26.6f');
% fprintf('===============================================================================\n');
sep(WDTH,'=')
		myprint(outp,in,[],[],[],WDTH);
sep(WDTH,'-')
% fprintf('-------------------------------------------------------------------------------\n');

% all regressor including constant (it if exists)
K_all = ols_out.k; 
% standard error or regression MLE (no DF adjustment) sqrt(SSE/T).
stderr_MLE = sqrt(ols_out.sse/ols_out.N);

% ols_out.INCLUDE_PW

%%
% --------------------------------------------------------------------------------------------------
extra_info1 = [	ols_out.R2;
								ols_out.R2_bar;
								sqrt(ols_out.sig2);
								ols_out.sse;
								ols_out.LL;
								Fstat;
								pval_Fstat;
								ols_out.N;
								stderr_MLE;
								ols_out.INCLUDE_PW];

extra_info2 = [	KK;
								ols_out.k	
								mean(ols_out.y);
								std(ols_out.y,1); 
								ols_out.aic; 
								ols_out.aicc; 
								ols_out.bic;  
								ols_out.hq;
								ols_out.DW;
								ols_out.L];


% print some extra results 
ev1s = strjust(num2str((extra_info1),in.fmt),'right');
ev2s = strjust(num2str((extra_info2),in.fmt),'right');

% nam1 = strvcat(	'R-squared'			,...
% 				'Rbar-squared'		,...
% 				'SE of regression'	,...
% 				'Sum Squared Errors',...
% 				'Log-likelihood'	,...
% 				'F-statistic'		,...		 	
% 				'Pr(F-statistic)'	,...
% 				'Number of Regressors' , ...
% 				'Plus Const. if exists');

nam1 = strvcat(	'R-squared'			,...
								'Rbar-squared'					,...
								'SE of regression'			,...
								'Sum Squared Errors'		,...
								'Log-likelihood'				,...
								'F-statistic'						,...		 	
								'Pr(F-statistic)'				,...
								'No. of observations'		,...
								'Std.err.MLE (div by T)',...
								'Include Pre-whitening');

nam2 = strvcat(	'|   No. of Regressors'	, ...
								'|   Plus Const.(if exist)', ...
								'|   Mean(y)'						,...
								'|   Stdev(y)'					,...
								'|   AIC'								,...
								'|   AICc'							,...
								'|   BIC'								,...
								'|   HQIC'							,...
								'|   DW-stat.'					,...
								'|   HAC Trunct.Lag.'				);

for ii = 1:size(nam2,1)								
	fprintf(' %s :		  %s', nam1(ii,:), ev1s(ii,:))								
	fprintf(' %s :	   %s \n', nam2(ii,:), ev2s(ii,:))								
end
sep(WDTH,'=')


% --------------------------------------------------------------------------------------------------
%% print to excel as well
% xls_filename = 'ols_output.xls';
[Nr,Nc] = size(outp);

if prnt2xls ~= 0
	% check if it is a string, for filename decision
	if isstr(prnt2xls)
		xls_filename = prnt2xls;
	else
		xls_filename = 'ols_output.xls';
	end
	
	% main tstats etc
	% --------------------------------------------------------------------------------------------------
	xlswrite(xls_filename,in.rnames	,'Sheet1','A1')
	xlswrite(xls_filename,in.cnames','Sheet1','B1')
	xlswrite(xls_filename,outp			,'Sheet1','B2')

	% additional info
	% --------------------------------------------------------------------------------------------------
	Rbeg = num2str(Nr+2);

	xlswrite(xls_filename,cellstr(nam1)	,'Sheet1',['A' Rbeg])
	xlswrite(xls_filename,extra_info1		,'Sheet1',['B' Rbeg])

	xlswrite(xls_filename,cellstr(nam2(:,5:end))	,'Sheet1',['D' Rbeg])
	xlswrite(xls_filename,extra_info2				,'Sheet1',['E' Rbeg])
end





















% EOF



