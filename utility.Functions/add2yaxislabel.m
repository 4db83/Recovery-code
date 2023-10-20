function [] = add2yaxislabel
% function create a second y-axis label.
ylim_					= get(gca,'YLim');
ylim_ticks 		= get(gca,'YTick');
ylim_lables_0 = get(gca,'YTickLabel');
ylim_lables 	= strrep(ylim_lables_0,' ','');
% set second axis here
ax2 = gca();
yyaxis(ax2,'right');
ylim(ylim_);
set(ax2,'YColor',[0 0 0],'YTickLabel',ylim_lables,'YTick',ylim_ticks);







% ax = get(gca);
% abc0 = ax.YAxis.Exponent;
% ylim_	= get(gca,'YLim');
% ylim_ticks = get(gca,'YTick');
% ylim_lables_0 = get(gca,'YTickLabel');
% % ylim_lables = strrep(ylim_lables_0,' ','')
% ylim_lables = ylim_lables_0

% % if abc0 ~= 0
% yyaxis(gca,'right');
% % ylim(ylim_);
% % % abc0
% % set(gca,'YColor',[0 0 0])%,'YTickLabel',ylim_lables);
% % ax2 = gca();
% % abc = ax2.YAxis(2).Exponent
% % ax2.YAxis(2).Exponent = abc0;
% % disp('hell0')
% % else
% % ax2 = gca()
% % % make the zere entry alligned to the right
% % % ylim_lables(find(ylim_ticks==0)) = {' aa'}
% % ylim_lables = strrep(ylim_lables_0,' ','')
% % ylim_lables = ylim_lables_0
% % % set(gca,'YColor',[0 0 0],'YTickLabel',ylim_lables);
% % set(gca,'YColor',[0 0 0],'YTickLabel',ylim_lables,'YTick',ylim_ticks);
% % ax2.YTick(2) = ylim_ticks
% % end

% set(gca,'YColor',[0 0 0],'YTickLabel',ylim_lables,'YTick',ylim_ticks);


% ax.YAxis.Exponent 
% ax.YAxis.Exponent = 2
% set(gca,'YColor',[0 0 0])
% set(gca,'YColor',[0 0 0],'YTickLabel',ylim_lables,'YTick',ylim_ticks);


% function [] = add2yaxislabel
% function create a second y-axis label.


% % make the zere entry alligned to the right
% ylim_lables(find(ylim_ticks==0)) = {' aa'}
% 
% set(gca,'YColor',[0 0 0],'YTickLabel',ylim_lables,'YTick',ylim_ticks);
% 
% a = 1


% ax = gca;
% % ax.YRuler.TickLabelFormat = '%.1f';
% ax.YRuler.Exponent
% setoutsideTicks; add2yaxislabel 
% yyaxis(gca,'right');
% ax2 = gca;
% % set(gca,'YColor',[0 0 0],'YTickLabels',flipud(ylim_lables_0));
% set(ax2,'YColor',[0 0 0]);% 'YTickLabel',flipud(ylim_lables),'YTick',flipud(ylim_ticks));
% ax2.YTickLabels = 
