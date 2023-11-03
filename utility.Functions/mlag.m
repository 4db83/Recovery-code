function Xout = mlag(Xin,nlag,return_Xt,intitValue)
% PURPOSE: generates a matrix of n lags from a matrix (or vector)
%          containing a set of vectors (For use in var routines)
%---------------------------------------------------
% USAGE:     xlag = mlag(x,nlag)
%       or: xlag1 = mlag(x), which defaults to 1-lag
% where:		 x = a matrix (or vector), nobs x nvar
%					nlag = # of contiguous lags for each vector in x (default = 1 lag)
%                if == 0, returns X(t)
% return_Xt		 = returns also the x(t) vector/matrix in the first column(s), default set to NO (default = 0)
%		intitValue = (optional) scalar value to feed initial missing values (default = 0)
%---------------------------------------------------
% RETURNS:
%        xlag = a matrix of lags (nobs x nvar*nlag)
%        x1(t-1), x1(t-2), ... x1(t-nlag), x2(t-1), ... x2(t-nlag) ...
% --------------------------------------------------
% COMMENT: CALLS THE lag() FUNCTION
%---------------------------------------------------
% orgininally by lesage, but modded by db.

SetDefaultValue(2,'nlag',1);
SetDefaultValue(3,'return_Xt',0);
SetDefaultValue(4,'intitValue',NaN);

% check if input is a TimeTable, if so, return a timetable as well.
isTT = isa(Xin,'timetable');

if isTT
  xDates = Xin.Date;
  var_names = getvarnames(Xin);
  X = Xin.Variables;
  % make lagged names 
  N_names = length(var_names);
  tlag_names = cellstr(make_table_names('(t-', 1:1:nlag, ')'));
  Tlag_names = repmat(tlag_names',N_names,1);
  tmp_ = [];
  % loop through nlags for variable name output creation
  for jj = 1:nlag
    tmp_ = [tmp_; cellstr([char(var_names) char(Tlag_names(:,jj))]) ]; 
  end
  % if nlag > 0
  %   lag_names = strrep( tmp_,' ','.'); --> not needed, only to remove spaces in file anme
  % end
  lag_names = tmp_;
else
  X = Xin;
end

if return_Xt
	Xlag = X;
else
	Xlag = [];
end

for i = 1:nlag
  Xlag = [Xlag lag_F(X, i, intitValue)];
end 

if nlag == 0
  Xlag = X;
  if isTT
    lag_names = var_names;
  end
end

if isTT
  if return_Xt % add the time t names
    lag_names = [var_names; lag_names];
  end
  Xout = make_timetable(Xlag, lag_names, xDates);
else
  Xout = Xlag;
end


end % EOF 


function z = lag_F(x,n,v)
% PURPOSE: creates a matrix or vector of lagged values
% -------------------------------------------------------
% USAGE: z = lag(x,n,v)
% where: x = input matrix or vector, (nobs x k)
%        n = order of lag
%        v = (optional) initial values (default=NaN)
% e.g.
%     z = lag(x) creates a matrix (or vector) of x, lagged 1 observations
%         with initial value being set to zero. so needs to be trimmed to
%         get rid of initial 0.
%     z = lag(x,n) creates a matrix (or vector) of x, lagged n observations
%     z = lag(x,n,v) creates a matrix (or vector) of x, lagged n observations,
%         with initial values taking a value v.
% ------------------------------------------------------
% RETURNS: z = matrix (or vector) of lags (nobs x k)
% ------------------------------------------------------
% NOTES: if n <= 0, z = [] is returned. While you may find this
%        preverse, it is sometimes useful.
%-------------------------------------------------------
% SEE ALSO: mlag() 
%-------------------------------------------------------

% written by:
% James P. LeSage, Dept of Economics
% University of Toledo
% 2801 W. Bancroft St,
% Toledo, OH 43606
% jpl@jpl.econ.utoledo.edu

switch(nargin)
  case 1
    n = 1; v = NaN;
    zt = ones(n,size(x,2))*v;
    z = [ zt; trimr_F(x,0,n)];
  case 2
    v = NaN;
    if n < 1
    z = [];
    else
    zt = ones(n,size(x,2))*v;
    z = [ zt; trimr_F(x,0,n)];
    end
  case 3
    if n < 1
    z = [];
    else
    zt = ones(n,size(x,2))*v;
    z = [ zt; trimr_F(x,0,n)];
    end
  otherwise
    error('lag: wrong # of input arguments');
  end
end


function z = trimr_F(x,n1,n2)
% PURPOSE: return a matrix (or vector) x stripped of the specified rows.
% -----------------------------------------------------
% USAGE: z = trimr_F(x,n1,n2)
% where: x = input matrix (or vector) (n x k)
%       n1 = first n1 rows to strip
%       n2 = last  n2 rows to strip
% NOTE: modeled after Gauss trimr_F function
% -----------------------------------------------------
% RETURNS: z = x(n1+1:n-n2,:)
% -----------------------------------------------------

% written by:
% James P. LeSage, Dept of Economics
% University of Toledo
% 2801 W. Bancroft St,
% Toledo, OH 43606
% jpl@jpl.econ.utoledo.edu

  [n junk] = size(x);
  if (n1+n2) >= n; 
     error('Attempting to trim too much in trimr_F');
  end;
  h1 = n1+1;   
  h2 = n-n2;
  z = x(h1:h2,:);
end  