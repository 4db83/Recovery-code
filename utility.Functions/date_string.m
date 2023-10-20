function string_out = date_string(xx,frmat)
% simple date_string conversion function to daily, monthly, qarterly print
frmat = lower(frmat);
switch frmat
case 'd'
  string_out = datestr(xx,'dd.mm.yyyy');
case 'm'
  string_out = datestr(xx,'yyyy:mm');
case 'q'
  string_out = datestr(xx,'yyyy:qq');
end