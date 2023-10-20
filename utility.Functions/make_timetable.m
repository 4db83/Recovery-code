function ttout = make_timetable(MatX, VarNames, Dates, Dates_String_format)
% --------------------------------------------------------------------------------------------------
% Function: to make TimeTable from MatX,Date,VarNames. 
%          CALL AS: TT = make_timetable(MatX, VarNames, Dates, Dates_String_format)
%                        
% OLD make_timetable(MatX, Date, VarNames)
% --------------------------------------------------------------------------------------------------
% this is to compliment Matlab's array2timetable function, but converts automatically the annoying
% default to Time for the 'datetime' column and makes the Dates input more flexible.
% db: stockholm 30.05.2023.
% --------------------------------------------------------------------------------------------------
% SEE ALSO:
%   - TT_out = remove_prefix(hml2r, TT_in), where hml2r = how_many_letters_2_remove
%   - TT_out = make_timetable(MatX, Date, VarNames)
%   - TT_out = prefix_varnames(string_name, TT_in)
%   - TT_out = as_timetable(FredData, Varnames)
%   - aout = getFredData()
% --------------------------------------------------------------------------------------------------

% note the MM for month. --> see also datetime for the formatting
SetDefaultValue(4,'Dates_String_format','yyyy.MM.dd');

isDN = isa(Dates,'double');
isST = isa(Dates,'char');
isDT = isa(Dates,'datetime');


if (~isDN)*(~isDT)*(~isST)
  error('ERROR: Date vector must be a Datenum, Char/String, or Datetime object')
end

if isDN
  Dates = datetime( Dates, 'ConvertFrom','datenum');
elseif isST
  Dates = datetime( num2str(Dates), 'InputFormat', Dates_String_format);
else
  Dates = Dates;
end



ttout = array2timetable(MatX, ...
                        'RowTimes', Dates, ...
                        'VariableNames', VarNames ...
                        );

% change Time Column to Date
ttout.Properties.DimensionNames{1} = 'Date';  

end

