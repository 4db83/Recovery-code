function [] = printstrng2xls(xlsfilename,strng,Cellin,Sheetin)
% F: print string to xls file.
% --------------------------------------------------------------------------------------------------
% Call as: printstrng2xls(xlsfilename,strng) as minimum call.
% or
% as:		printstrng2xls(xlsfilename,strng,Cellin,Sheetin).
% for a full call.
% DB.
% --------------------------------------------------------------------------------------------------

if nargin < 4; Sheetin = 'Sheet1';	end;
if nargin < 3; Cellin	 = 'A1';			end;

% xlswrite(xlsfilename,cellstr(prntstrng2xls(strng)),Sheetin,Cellin);
xlswrite(xlsfilename,cellstr(strng),Sheetin,Cellin);

