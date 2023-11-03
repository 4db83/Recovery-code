function vnames = getvarnames(time_table_in)
isTT = isa(time_table_in,'timetable');
isTL = isa(time_table_in,'table');
isST = isa(time_table_in,'struct');

if isTT
  vnames = time_table_in.Properties.VariableNames';
elseif isST
  vnames = fieldnames(time_table_in);
end