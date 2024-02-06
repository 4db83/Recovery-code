function [] = add2yaxislabel(rightAlign)

if (nargin<1)||isempty(rightAlign); rightAlign = 0; end

% function create a second y-axis label.
ylim_					= get(gca,'YLim');
ylim_ticks 		= get(gca,'YTick');
ylim_lables_0 = get(gca,'YTickLabel');
ylim_lables 	= strrep(ylim_lables_0,' ','');
% if rightAlign second axis lables
if rightAlign == 1
  ylim_lables = num2str(ylim_ticks');
end

% ylim_lables
% set second axis here
ax2 = gca();
yyaxis(ax2,'right');
ylim(ylim_);
set(ax2,'YColor',[0 0 0],'YTickLabel',ylim_lables,'YTick',ylim_ticks);


