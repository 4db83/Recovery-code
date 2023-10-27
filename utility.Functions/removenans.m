function [Xno_nan, I] = removenans(x)
% F: removes nan ROWS.
% See also addnans(xin, beg_x, end_x) to pad a series with nans at top and bottom for plotting

I = anynans(x);

if sum(I)==0
	Xno_nan = x;
else
	Xno_nan = x(~I,:);
end

