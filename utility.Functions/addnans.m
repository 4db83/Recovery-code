function xout = addnans(xin, beg_x, end_x)
% add nans to front and/or end of a vector/matrix xin for plotting purposes.
% See also removenans
SetDefaultValue(2, 'beg_x', 0) % default value = nan.
SetDefaultValue(3, 'end_x', 0) % default value = nan.

[T,K] = size(xin);

xout = [nan(beg_x,K); xin; nan(end_x,K)];


