function aout = print_sample_period(datevector_in, date_format, add_sep, date_text)
% set default dateformate
SetDefaultValue(2,'date_format'  , 'd');
SetDefaultValue(3,'add_sep'      , [1 1]);
SetDefaultValue(4,'date_text'    , '       Period: ');

if length(date_format) == 1
  switch lower(date_format)
    case 'd'
      date_format = 'dd.mm.yyyy';
    case 'm'
      date_format = 'mmm-yyyy';
    case 'q'
      date_format = 'yyyy:QQ';
  end
else
  date_format = lower(date_format);
end

% check if input is timetable or datevectors
if istimetable(datevector_in)
  dates_full = datevector_in.Properties.RowTimes;
elseif istable(datevector_in)
  dates_full = datevector_in.Properties.UserData;
elseif isa(datevector_in,'timerange')
  dates_struct = struct(datevector_in);
else
  dates_full = datevector_in;
end 

if numel(add_sep) == 1;
  add_sep = add_sep*ones(1,2);
end

if isnumeric(date_text)
  date_text = ['      Tis = ' num2str(date_text) ',  Time Period: '];
end

if isa(datevector_in,'timerange')
  date0 = datestr(dates_struct.first, date_format);
  dateT = datestr(dates_struct.last, date_format);
else
  date0 = datestr( dates_full(1),   date_format);
  dateT = datestr( dates_full(end), date_format);
end

if nargout == 1
  aout = [ date0 ' to ' dateT];
  % aout = [ 'Sample: ' datestr(dates_full(1), date_format) ' to ' datestr(dates_full(end), date_format)];
else 
  % just print to screen if no output arguments
  if add_sep(1) == 1; sep(180); end
  fprintf([ date_text ' %s to']   , date0);
  fprintf(' %s\n'                 , dateT);
  if add_sep(2) == 1; sep(180); end
end

end