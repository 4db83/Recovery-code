function varargout = addgrid(LineWidth,Alpha,LineStyle,GridLayer)
% return only grid axis handle

box on;
grid on;
LW0 = 1.11;
LW0 = 4/5;
% set(hg,'LineWidth') = 4/3;
Alpha0 = 1/3;
LS0 = ':';
SetDefaultValue(1,'LineWidth' , LW0)
SetDefaultValue(2,'Alpha'     , Alpha0)
SetDefaultValue(3,'LineStyle' , LS0)
SetDefaultValue(4,'GridLayer' , 0)
set(gca,'LineWidth',LineWidth)
hg = gca;


% set(gca,'GridLineStyle',':','GridAlpha',1/3)
% hg.GridColor = [0 .5 .5];
hg.GridLineStyle  = LineStyle;
% check if property exists for older versions of matlab
if isprop(hg,"GridLineWidth"); hg.GridLineWidth = LineWidth; end
hg.GridAlpha      = Alpha;
hg.Layer          = 'bottom';
if GridLayer == 1; hg.Layer = 'top'; end

for k = 1:nargout
  varargout{k} = hg;
end

tickshrink(.5)