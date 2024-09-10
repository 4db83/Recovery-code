function z = lag(x,n,v)
% NOTE: there is a name clash. lag is also available with the timetable functions in matlab, but
% leave the name the same and let matlab automatically switch when input is timetable or double
% PURPOSE: creates a matrix or vector of lagged values:
% -------------------------------------------------------
% NOTE: works with timetable objects, but NOT with table.
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
   n = 1; 
   v = NaN;
   zt = ones(n,cols(x))*v;
   z = [ zt; trimr_F(x,0,n)];

case 2
   v = NaN;
    if n < 0
      z = [];
   % if n = 0, ie, no lag, just return the series itself
   elseif n == 0
      z = x;
   else 
      zt = ones(n,cols(x))*v;
      z = [ zt; trimr_F(x,0,n)];
   end

case 3
   if n < 0
      z = [];
   % if n = 0, ie, no lag, just return the series itself
   elseif n == 0
      z = x;
   else 
      zt = ones(n,cols(x))*v;
      z = [ zt; trimr_F(x,0,n)];
   end

otherwise
error('lag: wrong # of input arguments');
end
end




% calls the function
function z = trimr_F(x,n1,n2)
  % PURPOSE: return a matrix (or vector) x stripped of the specified rows.
  % -----------------------------------------------------
  % USAGE: z = trimr(x,n1,n2)
  % where: x = input matrix (or vector) (n x k)
  %       n1 = first n1 rows to strip
  %       n2 = last  n2 rows to strip
  % NOTE: modeled after Gauss trimr function
  % -----------------------------------------------------
  % RETURNS: z = x(n1+1:n-n2,:)
  % -----------------------------------------------------
  
  % written by:
  % James P. LeSage, Dept of Economics
  % University of Toledo
  % 2801 W. Bancroft St,
  % Toledo, OH 43606
  % jpl@jpl.econ.utoledo.edu
  
  % if only trimming from the front
  if nargin < 3; n2 = 0; end
  
  [n junk] = size(x);
  if (n1+n2) >= n 
     error('Attempting to trim too much in trimr');
  end
  h1 = n1+1;   
  h2 = n-n2;
  z = x(h1:h2,:);
end  
