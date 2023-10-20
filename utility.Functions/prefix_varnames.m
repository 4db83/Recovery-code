function TT_out = prefix_varnames(string_name_beg, TT_in, string_name_end, prefix_rownames)
% ------------------------------------------------------------------------------------------
% FUNCTION: Add a prefix 'string_name' to each variable in timetable TT_in and to end if
% string_name_beg is past through the function as well.
% ------------------------------------------------------------------------------------------
%         Call as: TT_out = prefix_varnames(string_name, TT_in)
% this is useful when wanting to merge to timetable with the same variable name structure. 
% db: Stockholm 30.05.2023
% ------------------------------------------------------------------------------------------
% SEE ALSO:
%   - TT_out = remove_prefix(hml2r, TT_in), where hml2r = how_many_letters_2_remove
%   - TT_out = make_timetable(MatX, Date, VarNames)
%   - TT_out = prefix_varnames(string_name, TT_in)
%   - TT_out = as_timetable(FredData, Varnames)
%   - aout = getFredData()
% ------------------------------------------------------------------------------------------
% nargin
  if nargin < 4
    prefix_rownames = 0;
  end 

  TT_varnmaes   = TT_in.Properties.VariableNames';
  n             = size(TT_varnmaes,1);
  % adding to the front of name
  
  % adding to the back of name
  if nargin < 3
    string_name_end = '';
  end
  % strng_char_beg = repmat(string_name_beg,n,1);
  % strng_char_end = repmat(string_name_end,n,1);
  
  strng_char_beg = string_name_beg;
  strng_char_end = string_name_end;
 
  % char(TT_varnmaes)

  TT_varnmaes   = strip(TT_varnmaes);
  new_varnames  = cell(n,1);
  for ii = 1:n
    new_varnames{ii}  = [strng_char_beg TT_varnmaes{ii} strng_char_end];
  end
  
  % new_varnames}  = cellstr([ strng_char_beg (TT_varnmaes) strng_char_end]);
  % now create the new table
  % new_varnames = cell(n,1)
  % for ii = 1:length(TT_varnmaes)
  %   new_varnames{ii}  = cellstr([ strng_char_beg char(TT_varnmaes{ii}) strng_char_end]);
  % end
  
  if prefix_rownames == 1
    Nrows = size(TT_in.Properties.RowNames,1);
    new_rownames = [repmat(string_name_beg, Nrows, 1)  ...
                    char(TT_in.Properties.RowNames)    ...
                    repmat(string_name_end, Nrows, 1) ]; 
  
    TT_in.Properties.RowNames = cellstr(new_rownames);
  else
    TT_in.Properties.VariableNames = new_varnames;
  end
  % output the table
  TT_out        = TT_in;
end




































