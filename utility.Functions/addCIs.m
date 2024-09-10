function [varargout] = addCIs(CI, colr, alpha) % ,color,edge,add,transparency)
  % Adds Confidence Interval (CI) shading to plot
  % --------------------------------------------------------------------------------------------------
  % colr optional: 
  % 'r' - Red,    % 'g' - Green,  % 'b' - Blue,  % 'c' - Cyan,  % 'm' - Magenta
  % 'y' - Yellow, % 'k' - Black,  % 'w' - White
  % --------------------------------------------------------------------------------------------------
  
  alpha_default = .15;
  Xgrid = (1:length(CI));
  
  upper_Y = CI(:,2)';
  lower_Y = CI(:,1)';
  
  if (nargin < 2) || isempty(colr)
    colr = .8*[.2 .4 .8]; % default color is grey;
    alpha = alpha_default;
  end
  
  if (nargin < 3) || isempty(alpha)
    alpha = alpha_default;
  end
  
  filled	= [ upper_Y, fliplr(lower_Y) ];
  XX = [Xgrid, fliplr(Xgrid)];
  
  fillhandle = fill(XX, filled, colr, 'EdgeColor','none', 'FaceAlpha', alpha);
  
  if nargout > 0
	  varargout{1} = fillhandle;
  end	
end

function cout = CLRs(ii)
color_base = [
      0         0.4470    0.7410; % Blue
      0.8500    0.3250    0.0980; % Red
      0.9290    0.6940    0.1250; % Yellow
      0.4940    0.1840    0.5560; % Purple
      0.4660    0.6740    0.1880; % Green
      0.3010    0.7450    0.9330; % Cyan
      0.6350    0.0780    0.1840  % Dark Red
  ];

  cout = color_base(ii,:);
end
