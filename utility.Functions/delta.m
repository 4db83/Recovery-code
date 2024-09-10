function dyout = delta(y,k,remove_nans)
% Function: Computes the kth difference/change in y(t): 
%                       Δᵏy(t) = (1-L)ᵏy(t)
% -------------------------------------------------------------------------------------------------
% NOTE: This is different form the function delta_long or long_diff which 
%       compute the long-difference, that is:
%                       Δₖ(y(t)) = (1-Lᵏ)y(t) = y(t) - y(t-k)
% -------------------------------------------------------------------------------------------------
% DESCRIPTION:
%     Computes the kth difference/change in y(t): Δᵏy(t) = (1-L)ᵏy(t)
% -------------------------------------------------------------------------------------------------
% USAGE:	dy = delta(y,k).
% -------------------------------------------------------------------------------------------------
% INPUT:
%     y = dependent variable  (Txn).                                     
%    	k = order of change to be computed.al for the
%  
% OUTPUT:
%	   dy = (Txn) vector of kth changes in y(t), with NaNs in
%         first k lagged postions, so needs to be trimmed.
% -------------------------------------------------------------------------------------------------
%   		Created by Daniel Buncic on 27.10.2023.
% 			Modified on:                27.10.2023.
% -------------------------------------------------------------------------------------------------

SetDefaultValue(2,'k',1);
SetDefaultValue(3,'remove_nans',0);

isTT = isa(y,'timetable');

% if isa(class(y),'fints')
if isTT 
	% get fieldnames
	names0	= y.Properties.VariableNames;
	% get dates first
	dates0	= y.Properties.RowTimes;
	% datamat 
	matx0		= y.Variables;
	[~,cy]  = size(matx0);
	dy_tmp  = [NaN(k,cy); diff(matx0,k)];
  dy      = make_timetable(dy_tmp,names0,dates0);
  switch k
    case 1;   dyout = prefix_varnames(['Δ' ],dy);
    case 2;   dyout = prefix_varnames(['Δ²'],dy);
    case 3;   dyout = prefix_varnames(['Δ³'],dy);
    case 4;   dyout = prefix_varnames(['Δ⁴'],dy);
    case 5;   dyout = prefix_varnames(['Δ⁵'],dy);
    otherwise
      dyout = prefix_varnames(['Δ' num2str(k)],dy);
  end
else
  % return actual series y if k == 0
  if k == 0 
    dyout  = y;	
  else
      % if k ~= 0
	  [~,cy] = size(y);
	  dyout  = [NaN(k,cy); diff(y,k)];	
  end
end

if remove_nans  
  dyout = removenans(dyout);
end

