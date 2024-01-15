function varargout = addgrid(LineWidth,Alpha,LineStyle)
% return only grid axis handle

box on;
grid on;
LW0 = 1.11;
% set(hg,'LineWidth') = 4/3;
Alpha0 = 1/3;
LS0 = ':';
SetDefaultValue(3,'Alpha'     , Alpha0)
SetDefaultValue(2,'LineStyle' , LS0)
SetDefaultValue(1,'LineWidth' , LW0)
set(gca,'LineWidth',LineWidth)
hg = gca;


% set(gca,'GridLineStyle',':','GridAlpha',1/3)
% hg.GridColor = [0 .5 .5];
hg.GridLineStyle  = LineStyle;
hg.GridLineWidth  = LineWidth;
hg.GridAlpha      = Alpha;
hg.Layer          = 'top';

for k = 1:nargout
  varargout{k} = hg;
end

