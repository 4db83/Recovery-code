function [Xno_nan, I] = removenans(x)
% F: removes nan ROWS.

I = anynans(x);

if sum(I)==0
	Xno_nan = x;
else
	Xno_nan = x(~I,:);
end

