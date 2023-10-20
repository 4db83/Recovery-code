function dyout = delta(y,k,remove_nans)
% function dy = delta(varargin)
% Function: Computes the kth change in y(t) as: y(t)-y(t-k)
%______________________________________________________________
%
% DESCRIPTION:
%
%	Computes the kth change in y(t) as: y(t)-y(t-k)
%______________________________________________________________
%
% USAGE:	dy = delta(y,k).
%______________________________________________________________
%
% INPUT:
%		
%     y = dependent variable  (Txn).                                     
%    	k = order of change to be computed.al for the
%  
% OUTPUT:
%
%	   dy = (Txn) vector of kth changes in y(t), with NaNs in
%         first k lagged postions, so needs to be trimmed.
%	
% NOTES:-------------------------------------------------------
%		NONE.
%______________________________________________________________
%   
%   		Created by Daniel Buncic on 29/8/2005.
% 			Modified on: 29/8/2005.

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
	[~,cy] = size(matx0);
	dy_tmp  = [NaN(k,cy); trimr_F(matx0,k,0)-trimr_F(matx0,0,k)];
	% dy			= fints(dates0,dy_tmp,names0(end-cy+1:end));
  dy      = make_timetable(dy_tmp,names0,dates0);
  switch k
    case 1;   dyout = prefix_varnames(['Δ' ],dy);
    % case 2;   dyout = prefix_varnames(['Δ₂'],dy);
    % case 3;   dyout = prefix_varnames(['Δ³'],dy);
    % case 4;   dyout = prefix_varnames(['Δ⁴'],dy);
    % case 5;   dyout = prefix_varnames(['Δ⁵'],dy);
    otherwise
      dyout = prefix_varnames(['Δ' num2str(k)],dy);
  end
else
	[~,cy] = size(y);
	dyout  = [NaN(k,cy); trimr_F(y,k,0) - trimr_F(y,0,k)];	
end

if remove_nans  
dyout = removenans(dyout);
end

% % [cax,args,nargs] = axescheck(varargin{:});
% % y = args{1};
% % [ry,cy] = size(y);
% % 
% % if nargs == 1;
% %   dy  = [NaN(1,cy); trimr_F(y,1,0)-trimr_F(y,0,1)];
% % else
% %   k   = args{2};
% %   dy  = [NaN(k,cy); trimr_F(y,k,0)-trimr_F(y,0,k)];
% % end






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

  n = size(x);
  if (n1+n2) >= n 
     error('Attempting to trim too much in trimr_F');
  end
  h1 = n1+1;   
  h2 = n-n2;
  z = x(h1:h2,:);
  