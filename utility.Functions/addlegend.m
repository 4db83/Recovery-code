function varargout = addlegend(varargin)%
% QUICKLY ADD LEGEND TO PLOT WITH STANDARD SETTINGS
% CALL AS:	addlegend(LG,legnames,1,16,'Latex');
%   addlegend(input_subselection_of_lines, CellofLegnames, Position either 1 3 5 7, Fontsize, Latex 1 Tex 2)
% or simply without the input selection lines
% 	addlegend(CellofLegnames, Position either 1 3 5 7, Fontsize, 'Latex' or 'Tex')
%
% MINIMUM call: addlegend(CellofLegnames) with default location and fonts size being TopRigh, 12 with latex font.
% ---------------------------------------------------------------------------------------------------
% NOTE: input_subselection_of_lines selects from a plot with many lines the ones that are
% supposed to get a legend entry. 
% BUFFER may need to be adjuste at times to move the legend box up or down
% varargin{1}
% iscell(varargin{2})
% function h_leg = addlegend(legNames,Anchor,FN,LATEX)
% or function h_leg = addlegend(PNT,legNames,Anchor,FN,LATEX)
% ---------------------------------------------------------------------------------------------------

% SetDefaultValue(2,'TEXX', 'Latex')

% to stop writing Latex font warnings to screen
warning off

FNS_0 = get(gca,'FontSize');
% FNS0	= FNS_0 - 1;
FNS0	= FNS_0 + 1;

Length_vargin = nargin;

% check if first entry is a cell of legend names, if not set to 1
if ischar(varargin{1}); varargin{1} = cellstr(varargin{1}); end
input_subselection_of_lines = ~iscell(varargin{1});

% if Length_vargin
% default is Latex
% INTerPret = 'Tex';
INTerPret = 'Latex';

if input_subselection_of_lines
	
	legName = varargin{2};
  if ischar(legName); legName = cellstr(legName); end

	if Length_vargin == 5 %|| isempty(varargin{4})
		INTerPret = varargin{5};
		if isempty(varargin{4}) 
			FN	= FNS0; 
		else
			FN	= varargin{4};
		end
		% check that Anchor is non-empty, if empty replcace with Default 3
		if isempty(varargin{3}) 
			Anchor	= 3; 
		else
			Anchor	= varargin{3};
		end
	end
	
	if Length_vargin == 4 %|| isempty(varargin{4})
		FN	= varargin{4};
		% check that Anchor is non-empty, if empty replcace with Default 3
		if isempty(varargin{3}) 
			Anchor	= 3; 
		else
			Anchor	= varargin{3};
		end
	end
	
	% Anchor
	if Length_vargin == 3 %|| isempty(varargin{3})
 		FN			= FNS0;
		Anchor	= varargin{3};
	end

	if Length_vargin == 2 %|| isempty(varargin{3})
 		Anchor		= 3;
		FN				= FNS0;
		INTerPret = 'Latex';
	end

% % replace _ with .
% for ii = 1:length(legName)
%   legName = strrep(legName,'_','.');
% %   if any(ismember(legName{ii}, '_'))
% %     legName{ii} = strrep(legName(ii),'_','.');
% %   end
%   if any(ismember(char(legName{ii}),'!'))
%     legName{ii} = strrep(legName{ii},'!','_');
%   end
% end

	legendflex(varargin{1}, ... 
						 legName, ...
						'fontsize'	, FN, ...
						'anchor'		, Anchor*ones(1,2), ... % 3=top Right, 2=top Center, 1=top Left, 5=bottom Right, 7=bottom Left, 
						'buffer'		,	[0 0], ... 
						'Interpreter', INTerPret);
else
	
	
legName = varargin{1};

	if Length_vargin == 4 %|| isempty(varargin{4})
% 		FN	= varargin{4};
% 		% check that Anchor is non-empty, if empty replcace with Default 3
% 		if isempty(varargin{3}) 
% 			Anchor	= 3; 
% 		else
% 			Anchor	= varargin{2};
% 		end
		if isempty(varargin{2}); Anchor	= 3; else; Anchor	= varargin{2}; end
		if isempty(varargin{3}); FN	= FNS0; else; FN	= varargin{3}; end
		if isempty(varargin{4}); INTerPret = 'Latex'; else; INTerPret	= varargin{4}; end
	end
	
	if Length_vargin == 3
		if isempty(varargin{2}); Anchor	= 3; else; Anchor	= varargin{2}; end
		if isempty(varargin{3}); FN	= FNS0; else; FN	= varargin{3}; end
		INTerPret = 'Latex';
	end

	% Anchor
	if Length_vargin == 2 %|| isempty(varargin{3})
		if isempty(varargin{2}); Anchor	= 3; else; Anchor	= varargin{2}; end
		FN				= FNS0;
		INTerPret = 'Latex';
	end

	if Length_vargin == 1 %|| isempty(varargin{3})
 		Anchor		= 3;
		FN				= FNS0;
		INTerPret = 'Latex';
	end

%   legName = strrep(legName,'_','.');

	legendflex(legName, ...
						'fontsize'	, FN, ...
						'anchor'		, Anchor*ones(1,2), ... % 3=top Right, 2=top Center, 1=top Left, 5=bottom Right, 7=bottom Left, 
						'buffer'		,	[0 0], ... 
						'Interpreter',INTerPret);
end

warning on

% varible output arguments
if nargout > 0
	varargout{1} = h_leg;
end


					
% SetDefaultValue(2,'Anchor', 3)
% SetDefaultValue(3,'FN', 12)
% SetDefaultValue(3,'PD', 12)

% % 	'fontsize'	,	11, ... 
% % 	'padding'		,	[5 4 5], ... 
% % 	'anchor'		,	{'nw','nw'}, ...
% % 	'buffer'		,	[0 0], ... 
% % 	'xscale'		,	.5);

% legendflex({'T-Bill(3M)';'Fed-Funds rate (LW03)';'Interest rate (HLW)'})
% 
% legendflex({'$\Delta y_t$';
% 						'MUE($\hat\lambda$)';
% 						'MUE($\hat\lambda^{\mathrm{Up}}$)'
% 						'MPLE';
% 						'MMLE';
% 						'Clark UC0';
% 						'Clark UC' ;}, 'fontsize', Fns-1, 'anchor',3*ones(1,2) ,'Interpreter','Latex')

% % 
% % if nargin == 0
% % 	fnames = get(gca,'Legend');
% % end
% % 
% % vargin
% % if nargin > 0
% % 	%SetDefaultValue('fontsize', 11)
% % 	nargin
% % end
% % 
% % h_leg = legendflex(fnames.String, ... 
% % 	'fontsize'	,	11, ... 
% % 	'padding'		,	[5 4 5], ... 
% % 	'anchor'		,	{'nw','nw'}, ...
% % 	'buffer'		,	[0 0], ... 
% % 	'xscale'		,	.5);
% % 
% % 
% % 
