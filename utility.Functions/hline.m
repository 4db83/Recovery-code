function hhh=hline(y,in1,in2,LW,CLR,LOS,LFNTS)
% function h=hline(y, linetype, label,linewidth,color,label_offset,label_fnts)
% 
% Draws a horizontal line on the current axes at the location specified by 'y'.  Optional arguments are
% 'linetype' (default is 'r:') and 'label', which applies a text label to the graph near the line.  The
% label appears in the same color as the line.
%
% The line is held on the current axes, and after plotting the line, the function returns the axes to
% its prior hold state.
%
% The HandleVisibility property of the line object is set to "off", so not only does it not appear on
% legends, but it is not findable by using findobj.  Specifying an output argument causes the function to
% return a handle to the line, so it can be manipulated or deleted.  Also, the HandleVisibility can be 
% overridden by setting the root's ShowHiddenHandles property to on.
%
% h = hline(42,'g','The Answer')
%
% returns a handle to a green horizontal line on the current axes at y=42, and creates a text object on
% the current axes, close to the line, which reads "The Answer".
%
%
% draws three lines with the appropriate labels and colors.
% 
% By Brandon Kuczenski for Kensington Labs.
% brandon_kuczenski@kensingtonlabs.com
% 8 November 2001
% aLW = get(gca,'LineWidth')
% SetDefaultValue(4, 'LW'			, 1.25);
SetDefaultValue(4, 'LW'			, get(gca,'LineWidth'));
SetDefaultValue(5, 'CLR'		, 0.15*ones(1,3) );
SetDefaultValue(6, 'LOS'		, [.2 .2]);
SetDefaultValue(7, 'LFNTS'	, get(gca,'FontSize'));

if length(y)>1  % vector input
    for I=1:length(y)
        switch nargin
        case 1
            linetype='r:';
            label='';
        case 2
            if ~iscell(in1)
                in1={in1};
            end
            if I>length(in1)
                linetype=in1{end};
            else
                linetype=in1{I};
            end
            label='';
        case 3
            if ~iscell(in1)
                in1={in1};
            end
            if ~iscell(in2)
                in2={in2};
            end
            if I>length(in1)
                linetype=in1{end};
            else
                linetype=in1{I};
            end
            if I>length(in2)
                label=in2{end};
            else
                label=in2{I};
						end
				case 4
            if ~iscell(in1)
                in1={in1};
            end
            if ~iscell(in2)
                in2={in2};
            end
            if I>length(in1)
                linetype=in1{end};
            else
                linetype=in1{I};
            end
            if I>length(in2)
                label=in2{end};
            else
                label=in2{I};
						end
				case 5
            if ~iscell(in1)
                in1={in1};
            end
            if ~iscell(in2)
                in2={in2};
            end
            if I>length(in1)
                linetype=in1{end};
            else
                linetype=in1{I};
            end
            if I>length(in2)
                label=in2{end};
            else
                label=in2{I};
						end
				case 6
            if ~iscell(in1)
                in1={in1};
            end
            if ~iscell(in2)
                in2={in2};
            end
            if I>length(in1)
                linetype=in1{end};
            else
                linetype=in1{I};
            end
            if I>length(in2)
                label=in2{end};
            else
                label=in2{I};
            end
        end
        h(I)=hline(y(I),linetype,label);
    end
else
    switch nargin
    case 1
%        linetype='r:';
        linetype='k-';
        label='';
    case 2
        linetype=in1;
        label='';
    case 3
        linetype=in1;
        label=in2;
		case 4
        linetype=in1;
        label=in2;
		case 5
        linetype=in1;
        label=in2;
		case 6
        linetype=in1;
        label=in2;
    end
  
    g=ishold(gca);
    hold on

		vertical_offset		= LOS(1);
		horizontal_offset = LOS(2);

    x=get(gca,'xlim');
    if isempty(linetype);linetype = '-'; end

    h=plot(x,[y y],linetype,'Linewidth',LW,'Color',CLR);

    if ~isempty(label)
        yy=get(gca,'ylim');
        yrange=yy(2)-yy(1);
        yunit=(y-yy(1))/yrange;
        if yunit<0.2
            text(x(1) + horizontal_offset*(x(2)-x(1)),y + vertical_offset*yrange,label,'color',CLR,'Fontsize',LFNTS);
        else
            text(x(1) + horizontal_offset*(x(2)-x(1)),y - vertical_offset*yrange,label,'color',CLR,'Fontsize',LFNTS);
        end
    end

    if g==0
    hold off
    end
    set(h,'tag','hline','handlevisibility','off') % this last part is so that it doesn't show up on legends
end % else

if nargout
    hhh=h;
end






























