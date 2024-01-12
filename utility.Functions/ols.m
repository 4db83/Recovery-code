function olsout = ols(y,x,no_const,Xnames,L,INCLUDE_PW,print_results,Recession_indicator)
% function olsout = ols(y,x,no_const,Xnames,L,INCLUDE_PW,print_results,Recession_indicator)
% olsout = ols(y,X,L,INCLUDE_PW,Recession_indicator,no_const);
% wrapper function to fullols to print results to screen.

SetDefaultValue(3,'no_const',0)						% default is to include a constant
SetDefaultValue(4,'Xnames',[]);           % cellstring of regressor names
SetDefaultValue(5,'L',[]);								% set default value for truncation lag to L = 
SetDefaultValue(6,'INCLUDE_PW',[]);				% set to 1 to include pre-whitening in LRV computation
SetDefaultValue(7,'print_results',1)			% print results to screen is default
SetDefaultValue(8,'Recession_indicator',[]);

ytt = isa(y,"timetable");

if ytt 
  dates = y.Properties.RowTimes;
  y     = y.Variables; 
end

xtt = isa(x,"timetable"); 
if xtt
  if isempty(Xnames) 
    Xnames  = x.Properties.VariableNames;
  else 
    Xnames = Xnames;
  end
  dates   = x.Properties.RowTimes;
  x       = x.Variables; 
end

II = anynans(y,x);
% olsout = fullols(y,x,L,INCLUDE_PW,Recession_indicator,no_const);
olsout = fullols(y(~II),x(~II,:),no_const,L,INCLUDE_PW,Recession_indicator);

% add the dates vector
if ytt || xtt
  Dates = dates(~II);
  olsout.dates = Dates;
end

% sqrt(diag(olsout.HAC_VCV))
% sqrt(diag(olsout.HAC_VCV_pw_ARMA11))
% sqrt(diag(olsout.HAC_VCV_pw_AR1))

if ~isempty(Xnames)
  if ~iscell(Xnames); Xnames = cellstr(Xnames); end
end

% print sample period if it exists. 
if ytt || xtt
  fprintf( ['     Estimation period:  ' date_string(Dates(1),'m') ' to ' date_string(Dates(end),'m') '   Sample Size N = ' num2str(olsout.N) '\n' ])
end

if ischar(print_results)
	print_fullols(olsout, Xnames, [], print_results);
elseif print_results==1
	print_fullols(olsout, Xnames);
% 	print_fullols(olsout, Xnames, [], print_results);
end


