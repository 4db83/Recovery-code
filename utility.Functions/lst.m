function []= lst(data_without_dates,decimls,unts)
% -------------------------------------------------------------------------------------
% function	lst() prints in nice format numerical values out. 
% call as:	lst(datamatrix_without_dates).
% or as		lst(datamatrix_without_dates,decimals,units).
% =====================================================================================
% db (March, 2014).
% -------------------------------------------------------------------------------------

% set default value for truncation lag to L = 
SetDefaultValue(2,'decimls',4);
SetDefaultValue(3,'unts',4);

prnt_format = ['%' num2str(unts) '.' num2str(decimls) 'f'];

if ischar(decimls)
  prnt_format = decimls;
end

x = data_without_dates;

if isa(x,'table')
  x = table2timetable(x);
end

if isstruct(x)
	fldnams   = fieldnames(x);
	fld_data	= fldnams(find(strcmp(fldnams,'data')));
	fld_dates	= fldnams(find(strcmp(fldnams,'dates')));
	% call to lstd with dates function
	lstd([x.(char(fld_dates)) x.(char(fld_data))],[],prnt_format);
elseif isa(x,'timetable')
  fldnams = ['Date' x.Properties.VariableNames];
  xout    = [datestr(x.Time) repmat('    ',size(x,1),1) num2str(x.Variables,[prnt_format '  '])];
%   disp(fldnams);
%   fprintf('%s \n',fldnams)
  fprintf([repmat('%s             ', 1, size(fldnams,2)) '\n'], fldnams{:})
  disp(xout);
else
	if_no_dates = 1;
	if if_no_dates
		% is0 = x==0
		xout = [num2str(x(:,1:end),[prnt_format '  '])];
	else
		xout = [datestr(x(:,1), date_format) repmat('    ',size(x,1),1) num2str(x(:,2:end),[prnt_format '  '])];
  end
	% return the required string
	strng_data		= xout;
	% varargout{1}	= strng_data;
	disp(strng_data);
end

	% % if nargin < 2
	% % 	if_no_dates = 0;
	% % 	prnt_format = '%6.6f';
	% % end;
	% % 
	% % if nargin < 3
	% % 	prnt_format = '%6.6f';
	% % end;
