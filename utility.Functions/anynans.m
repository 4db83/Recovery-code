function [Inan, varargout] = anynans(X,Y)
% find any rows with NANs in regressor matrix X and also in Y if inputed and return indicator Inan
% useful for workign with time series data that has been lagged by mlag or differenced by delta.
% can be used to trim out the require entries for OLS etc. 
% then return rows with anynans.

% find missing values in X 
if strcmp(class(X),'fints'); X = fts2mat(X); end
if istimetable(X); X = X.Variables; end

Ixall = isnan(X);   % for all individually columns
Ix = any(Ixall,2);  % for all jointly

if nargin > 1
  % find missing values in Y if supplied
  if strcmp(class(Y),'fints'); Y = fts2mat(Y); end
  if istimetable(Y); Y = Y.Variables; end

  Iyall = isnan(Y);   % for all individually columns
  Iy = any(Iyall,2);  % for all jointly
  % now join them for the output
  Inan = Ix|Iy;
else 
  Inan = Ix;
end

if nargout > 1
	varargout{1} = Ixall;
  if nargin > 1
    varargout{2} = Iyall;
  end
end
