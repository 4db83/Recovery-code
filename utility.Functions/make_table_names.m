function names_out = make_table_names(prefix, numbers, suffix, frmt)
% function to easliy make cell array of indexed names. 
% call as: aout = make_table_names('ridge_z_', 1:12, '%0.4g').
% --------------------------------------------------------------------------------------------------

nS = length(numbers);
if nargin < 3 || isempty(suffix)
  suffix = '';
end
nX = length(suffix);

if nS == 1
  d = floor(log10(numbers))+1;
  name_end = num2str([1:numbers]', ['%0' num2str(d) 'd']);
  if nX == 1; suffix = repmat(suffix,numbers,1); end
  names_out = [repmat(prefix, numbers, 1) name_end suffix];
else
  names_out = cell(nS,1);
  for ii = 1:nS
    z_i = numbers(ii);
    % frmt = '%0.5g';
    % frmt = []
    if nargin < 4 || isempty(frmt)
      names_out(ii) = cellstr([prefix strrep(num2str(z_i),'.','') suffix]);
    else
      names_out(ii) = cellstr([prefix strrep(num2str(z_i,frmt),'.','') suffix]);
    end
  end
end


end