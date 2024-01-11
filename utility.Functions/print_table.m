function [] = print_table(y_in, FMT, no_index, table_header, xlsout, counter_exists )
% FUNCTION: print table or timetable to screen / or excel
% 
% USAGE:      print_table(y_in, FMT, no_index, table_header, xlsout )
% -----------------------------------------------------------------------------------------
% FMT = '%10.4f';
% FMT = {'%2.6f'};
% FMT  = cellstr(repmat('%10.4f',7,1))'; FMT(1) = {'%5.0f'};
% FMT = 14;
% FMT = [0 9*ones(1,6)]
% these are the possible format values. 
% ----------------------------------------------------------------------------------------- 
SetDefaultValue(2, 'FMT', 4)
SetDefaultValue(3, 'no_index', 0)
SetDefaultValue(4, 'table_header', [])
SetDefaultValue(5, 'xlsout', 0)
SetDefaultValue(6, 'counter_exists', []) % set ot zero to force not using the counter check in column 1
% --------------------------------------------------------------------------------------------------

% PRINT TO SCREEN
% check if table or timetable
isT   = istable(y_in);
isTT  = istimetable(y_in);

if ~(isT|isTT); error(' ERROR: Input must be a table or a timetable'); end

ColNames_tmp  = y_in.Properties.VariableNames(1:end);
if isTT % if TT, get dates
  RowNames_TT   = y_in.Properties.RowTimes;
else
  RowNames_TT   = y_in.Properties.RowNames;
end
y             = table2array(y_in);

% if isT
%   % if isempty(RowNames_TT)
%   % RowNames_tmp  = eval(['y_in.' char(y_in.Properties.VariableNames(1))]);
%   RowNames_tmp  = y_in.Properties.RowNames;
%   isempty(RowNames_tmp)
%   y					    = table2array(y_in(:,1:end));
%   ColNames_     = ColNames_tmp(2:end);
% 
%   % if isa(RowNames_tmp,'double')
%   %   RowNames_tmp  = y_in.Properties.RowNames;
%   % end
%   % else
%   %   RowNames_tmp  = y_in.Properties.RowNames;
%   % end
% elseif isTT
%   y					    = table2array(y_in(:,1:end));
%   ColNames_     = ColNames_tmp(1:end);
%   RowNames_tmp  = y_in.Properties.RowTimes;
% end 

mmy = max(max(abs(y)));
Nc  = ceil(log10(mmy));

% das ist kcakck

% --------------------------------------------------------------------------------------------------
% MAKE HERE THE TAILORED FORMATTING ON THE LEFT OF THE DECIMAL FOR THE INDIVIDUAL COMLUMNS
kk = 1.4; % scaling factor

frmat_dble = ceil(log10( nanmax(abs(ceil(y))) + 1 )) ;
frmat_dble(isnan(frmat_dble)) = 0;  % set formatting of nan values to 0
frmat_dble = ceil(frmat_dble*kk + FMT*kk) ; 
tmp1 = strcat('%', num2str(frmat_dble'), '.', num2str(FMT), 'f');
% tmp1 = strcat('%', num2str(11*ones(size(tmp1,1),1)), '.', num2str(FMT), 'f')
% [{'%5.0f'} cellstr(repmat(frmat_dble, length(ColNames_), 1))']
% --------------------------------------------------------------------------------------------------

[Nr, Nc] = size(y);
% add an index to the table
nn = (1:Nr)';

% Force counter_exits to 0
if isempty(counter_exists)
  % check if counter exists already
  counter_exists = any(mod(y(:,1),1)) == 0;
else
  counter_exists = counter_exists;
end

if counter_exists
  X = [y];
else
  if no_index
    X = y;
  else % add counter if it does not exist
  X = [nn y];
  end
end

if Nr == length(RowNames_TT)
  RowNames_tmp = RowNames_TT;
end
if Nc == length(ColNames_tmp)
  ColNames_ = ColNames_tmp;
end

% if no_index
%   X = y;
% end

% CHECK FMT FORMAT: these all return a cellstructure, must be converted to char later on
if isa(FMT, 'char')
  FMT0 = FMT;
  % add the counter formattinig
  FMT1 = [{'%5.0f'} cellstr(repmat(FMT0, length(ColNames_)-counter_exists, 1))'];
  
elseif isa(FMT, 'cell')
  disp('cell')
  if length(FMT) > 1
    disp('>1')
    FMT0 = FMT;
    % add the counter formattinig
    FMT1 = [{'%5.0f'} FMT0];
  else
    FMT0 = char(FMT);
    FMT1 = [{'%5.0f'} cellstr(repmat(FMT0, length(ColNames_) -counter_exists, 1))'];    
  end
elseif isa(FMT, 'double') % if double 
  if length(FMT) == 1     % if scalar input
    FMT1 = [{'%5.0f'} cellstr(tmp1(1+counter_exists:end,:))'];
  else 
    FMT = FMT';
    L   = length(FMT);
    s   = FMT<6;
    N5  = 7*ones(L,1);
    b   = FMT>=5 & FMT<=8;
    N10 = 5*ones(L,1);
    vb  = FMT>8;
    N18 = 14*ones(L,1);
    JN  = N5.*s + N10.*b + N18.*vb;
    fmtout0   = [ repmat('%',L,1) num2str(JN) repmat('.',L,1) num2str(FMT) repmat('f',L,1)];
    FMT0 = cellstr(fmtout0);
    % add the counter formattinig
    FMT1 = [{'%5.0f'}; FMT0];
  end
end
% FMT1

% ADD ROW COUNTER COLUMN 'N' TO LIST OF COLNAMES
if counter_exists
  if ~strcmp({'N'},ColNames_(1))
    ColNames_ = [{'N'} ColNames_];
  else
    ColNames_ = ColNames_;
  end
else
  if ~ no_index
    ColNames_ = [{'N'} ColNames_];
  else
    ColNames_ = ColNames_;
  end
end

if isTT
  RowNames_ = char(['Date'; cellstr(datestr(datenum(y_in.Properties.RowTimes)))]);
else
  % RowNames_ = ['Model'; RowNames_tmp];
  RowNames_ = [' '; RowNames_tmp];
end

TTin = [];
TTin.cnames = char(ColNames_);
TTin.rnames = char(RowNames_);
if length(ColNames_) == length(FMT1)
  TTin.fmt = char(FMT1(1:end));
else
  TTin.fmt = char(FMT1(2:end));
end

header_empty = isempty(table_header);
date_empty   = isempty(y_in.Properties.UserData);

% make empty space for oos period.
oos = [];       % if empty
if ~date_empty  % else fill with oos period.
  oos = print_sample_period(y_in.Properties.UserData, 'm',0);
end

if table_header == 0
  mprint(X, TTin)
else
  if header_empty
    sep(180)
    fprintf([oos '\n'] );
    sep(180)
    mprint(X, TTin)
  else
    sep(180)
    fprintf([table_header oos '\n'] );
    sep(180)
    mprint(X, TTin)
  end
end

% sep(180)

% IF PRINTING TO EXCEL
if (xlsout~=0)
	if ischar(xlsout)
		xlsout_name = [xlsout '.xlsx'];
	else
		% xlsout_name = 'matlab_output.xlsx';
    xlsout_name = [inputname(1) '.xlsx'];
	end

  % NOW WRITE TO XLS.
  if ~header_empty
    printstrng2xls(xlsout_name, [table_header oos]  , 'B1');
 	  xlswrite(xlsout_name, ColNames_,	'Sheet1', 'B2');
	  xlswrite(xlsout_name, RowNames_,  'Sheet1', 'A2');
	  xlswrite(xlsout_name, X,					'Sheet1', 'B3');
  else  
	  xlswrite(xlsout_name, ColNames_,	'Sheet1', 'B1');
	  xlswrite(xlsout_name, RowNames_,  'Sheet1', 'A1');
	  xlswrite(xlsout_name, X,					'Sheet1', 'B2');
  end
end
end


























% --------------------------------------------------------------------------------------------------
% HELPER FUNCTIONS
% --------------------------------------------------------------------------------------------------
function mprint(y,info)
% PURPOSE: print an (nobs x nvar) matrix in formatted form
%---------------------------------------------------
% USAGE:     mprint(x,info) 
% where: x         = (nobs x nvar) matrix (or vector) to be printed
%        info      = a structure containing printing options
%        info.begr = beginning row to print,    (default = 1)
%        info.endr = ending row to print,       (default = nobs)
%        info.begc = beginning column to print, (default = 1
%        info.endc = ending column to print,    (default = nvar)        
%        info.cnames = an (nvar x 1) string vector of names for columns (optional)
%                      e.g. info.cnames = strvcat('col1','col2');
%                      (default = no column headings)
%        info.rnames = an (nobs+1 x 1) string vector of names for rows (optional)
%                      e.g. info.rnames = strvcat('Rows','row1','row2');
%                      (default = no row labels)
%        info.fmt    = a format string, e.g., '%12.6f' or '%12d' (default = %10.4f)
%                      or an (nvar x 1) string containing formats
%                      e.g., info.fmt=strvcat('%12.6f','%12.2f','%12d'); for nvar = 3
%        info.fid    = file-id for printing results to a file
%                      (defaults to the MATLAB command window)
%                      e.g. fid = fopen('file.out','w'); 
%        info.rflag  = 1 for row #'s printed, 0 for no row #'s (default = 0) 
%        info.width  = # of columns before wrapping occurs (default = 80)                                                  
%---------------------------------------------------
% e.g.   in.cnames = strvcat('col1','col2');
%        in.rnames = strvcat('rowlabel','row1','row2');
%        mprint(y,in), prints entire matrix, column and row headings
%        in2.endc = 3; in2.cnames = strvcat('col1','col2','col3');
%    or: mprint(y,in2), prints 3 columns of the matrix, just column headings 
%    or: mprint(y), prints entire matrix, no column headings or row labels 
% NOTES: - defaults are used for info-elements not specified
%        - default wrapping occurs at 80 columns, which varies depending on the
%          format you use, e.g. %10.2f will wrap after 8 columns
%---------------------------------------------------
% SEE ALSO: tsprint, mprint_d, lprint
%---------------------------------------------------

% written by:
% James P. LeSage, Dept of Economics
% University of Toledo
% 2801 W. Bancroft St,
% Toledo, OH 43606
% jpl@jpl.econ.utoledo.edu
% header line default length and style
% cwidth  = 120;
  cwidth = 180;
  lsep    = '-';
  LL_ = [repmat(lsep,1,cwidth)];
  
  % setup defaults
  fid = 1; rflag = 0; cflag = 0; rnum = 0; nfmts = 1; 
  [nobs, nvars] = size(y);
  begr = 1; endr = nobs; begc = 1; endc = nvars; 
  fmt = '%10.4f';
  if nargin == 1
  % rely on defaults
  elseif nargin == 2
    if ~isstruct(info)
      error('mprint: you must supply the options as a structure variable'); 
    end;
  fields = fieldnames(info);
  nf = length(fields);
  for i=1:nf
      if strcmp(fields{i},'fmt')
          fmts = info.fmt; 
    [nfmts, ~] = size(fmts);
    if nfmts == nvars
     fmt = fmts;
    elseif nfmts == 1
     fmt = fmts;
    else
     error('mprint: wrong # of formats in string -- need nvar');
    end
      elseif strcmp(fields{i},'fid')
          fid = info.fid;
      elseif strcmp(fields{i},'begc')
          begc = info.begc;
      elseif strcmp(fields{i},'begr')
          begr = info.begr;
      elseif strcmp(fields{i},'endc')
          endc = info.endc;
      elseif strcmp(fields{i},'endr');
          endr = info.endr;
      elseif strcmp(fields{i},'width');
        cwidth = info.width;
      elseif strcmp(fields{i},'cnames');
          cnames = info.cnames;
          cflag = 1;
      elseif strcmp(fields{i},'rnames');
          rnames = info.rnames;
          rflag = 1;
      elseif strcmp(fields{i},'rflag');
          rnum = info.rflag;
      end;
  end;
  
  else
  error('Wrong # of arguments to mprint');
     
  end; % end of if-elseif input checking
  
  
  % see if the user supplied row names and set rnum
  % correct her mistake if she did this
  if rflag == 1
  rnum = 0;
  end;
  
  % parse formats
  if nfmts == 1
     f1 = strtok(fmt,'%');
     f2 = strtok(f1,'.'); 
      if strcmp(f1,f2)
       f2 = strtok(f2,'d');
       dflag = 1;
       fflag = 0;
      else
       tmp1 = strtok(fmt,'f');
       % ----------------------------------------------------------------------------------------- 
       % tmp2 = strtok(tmp1,'.');
       % ----------------------------------------------------------------------------------------- 
       [tmp2,a] = strtok(tmp1,'.');
       tmp1 = tmp1(2:length(tmp1));
       tmp2 = tmp2(2:length(tmp2));
       opoint = num2str(str2num(tmp1) - str2num(tmp2));
       % ----------------------------------------------------------------------------------------- 
       % decimal = opoint(1,length(opoint));
       % ----------------------------------------------------------------------------------------- 
       % I HAVE CHANGED THIS TO MAKE THE DECIMLA PRINTING NOT STOP AT 9 DIGITS AS IN ORIGINAL
       % MPRINT.M
       % ----------------------------------------------------------------------------------------- 
       decimal = a(2:length(a));
       f2 = strtok(f2,'f');
       fflag = 1;
       dflag = 0;
      end;
     f2 = str2num(f2);
     nwide = floor(cwidth/f2); % 80 columns divided by format
     nvar = endc-begc+1;
     nsets = ceil(nvar/nwide);
  else %  wrapping in this case is based on widest format in the list
  nwidev = zeros(nfmts,1);
  nsetsv = zeros(nfmts,1);
  f2v = zeros(nfmts,1);
  dflagv = zeros(nfmts,1);
  fflagv = zeros(nfmts,1);
  decimalv = zeros(nfmts,1);
     for ii=1:nfmts;
     f1 = strtok(fmt(ii,:),'%');
     f2 = strtok(f1,'.');
      if strcmp(f1,f2)
       f2 = strtok(f2,'d');
       dflagv(ii,1) = 1;
       fflagv(ii,1) = 0;     
      else
       tmp1 = strtok(fmt(ii,:),'f');
       % ----------------------------------------------------------------------------------------- 
       % tmp2 = strtok(tmp1,'.');
       % ----------------------------------------------------------------------------------------- 
       [tmp2,a] = strtok(tmp1,'.');
       tmp1 = tmp1(2:length(tmp1));
       tmp2 = tmp2(2:length(tmp2));
       opoint = num2str(str2num(tmp1) - str2num(tmp2));
       % ----------------------------------------------------------------------------------------- 
       decimalv(ii,1) = opoint(1,length(opoint));
       % ----------------------------------------------------------------------------------------- 
       % decimalv(ii,1) = a(2:length(a));
       f2 = strtok(f2,'f');
       fflagv(ii,1) = 1;
       dflagv(ii,1) = 0;     
      end;
     f2v(ii,1) = str2num(f2);
     nwidev(ii,1) = floor(cwidth/f2v(ii,1)); % cwidth columns divided by format
     nvar = endc-begc+1;
     nsetsv(ii,1) = ceil(nvar/nwidev(ii,1));   
  end;
  nsets = min(nsetsv); 
  nwide = max(nwidev);
  end; 
  
  % if we have row and column labels
  % adjust variable labels and column heading strings
  % to match the width of the printing format
  
  if rnum == 1
  dstr = 'Obs#';
  end;
  
  if cflag == 1 % we have column headings
   [vsize nsize] = size(cnames); % error check cnames argument
   if vsize ~= nvars; error('Wrong # cnames in mprint'); end;    
   if nfmts == 1 % case of only 1 format string
    nmax = max(f2,nsize); % build format strings 
                          % based on widest format              
    sfmt = ['%', num2str(nmax)];
    sfmt = [sfmt,'s ']; 
    ffmt = ['%', num2str(nmax)];
     if dflag == 1
     ffmt = [ffmt,'d '];
     elseif fflag == 1
     ffmt = [ffmt,'.'];
     ffmt = [ffmt,decimal];
     ffmt = [ffmt,'f '];
     end;
   else % we have multiple format strings, process each
   sfmtv = []; fmtv = [];
    for ii=1:nfmts % find and parse multiple formats
    nmax = max(f2v(ii,:),nsize); % build format strings 
                          % based on widest format              
    sfmtv{ii} = ['%', num2str(nmax)];
    sfmtv{ii} = [sfmtv{ii},'s ']; 
    ffmtv{ii} = ['%', num2str(nmax)];
     if dflagv(ii,1) == 1
     ffmtv{ii} = [ffmtv{ii},'d '];
     elseif fflagv(ii,1) == 1
     ffmtv{ii} = [ffmtv{ii},'.'];
     ffmtv{ii} = [ffmtv{ii},decimalv(ii,1)];    
     ffmtv{ii} = [ffmtv{ii},'f '];
     end;
    end; % end of for ii loop
   end; % end of if-else
  elseif cflag == 0 % we have no column headings
   if nfmts == 1 % case of only 1 format string
    nmax = f2; % augment format string with a space (the hard way) 
    ffmt = ['%', num2str(nmax)];
     if dflag == 1
     ffmt = [ffmt,'d '];
     elseif fflag == 1
     ffmt = [ffmt,'.'];
     ffmt = [ffmt,decimal];
     ffmt = [ffmt,'f '];
     end;
   else % we have multiple format strings, process each
   sfmtv = []; fmtv = [];
    for ii=1:nfmts % find and parse multiple formats
    nmax = f2v(ii,:); % augment format strings with a space 
    ffmtv{ii} = ['%', num2str(nmax)];
     if dflagv(ii,1) == 1
     ffmtv{ii} = [ffmtv{ii},'d '];
     elseif fflagv(ii,1) == 1
     ffmtv{ii} = [ffmtv{ii},'.'];
     ffmtv{ii} = [ffmtv{ii},decimalv(ii,1)];    
     ffmtv{ii} = [ffmtv{ii},'f '];
     end;
    end; % end of for ii loop
   end; % end of if-else    
  end; % end of if-elseif cflag == 0,1
     
  if rflag == 1 % we have row labels
   [vsize nsize] = size(rnames); % error check cnames argument
   if vsize ~= nobs+1; error('Wrong # rnames in mprint'); end;  
   rfmt = ['%', num2str(nsize)]; 
   rfmt = [rfmt,'s ']; 
  end; % end of if rflag == 1
  
  if (rflag == 0 & cflag == 0)
      ffmt = fmt;
  end;
  
  % print matrix
  for j=1:nsets;
   if nfmts == 1 % print row header and column headers
   if rnum == 1;fprintf(fid,'%5s ',dstr);     
       elseif rflag == 1    
    fprintf(fid,rfmt,rnames(1,:));
       end;  
       if cflag == 1
      for i = (j-1)*nwide+begc:j*nwide+begc-1
    if i <= endc
  % find version #; 
    %[version,junk] = version; vers = str2num(version);
     %if vers == 5.2
     fprintf(fid,sfmt,strjust(cnames(i,:),'right'));
     %else
     %fprintf(fid,sfmt,strjust(cnames(i,:)));
     %end;
    end;
   end;
       end;
    % fprintf(fid,'\n');
  % fprintf(fid,'\n------------------------------------------------------------\n');
  fprintf(fid,['\n' LL_ '\n']);
   else % we have multiple formats
   if rnum == 1;fprintf(fid,'%5s ',dstr);     
      elseif rflag == 1   
   fprintf(fid,rfmt,rnames(1,:));
      end;
      if cflag == 1
     for i = (j-1)*nwide+begc:j*nwide+begc-1
    if i <= endc
  % find version #; 
    %[version,junk] = version; vers = str2num(version);
     %if vers == 5.2
     fprintf(fid,sfmtv{i},strjust(cnames(i,:),'right'));
     %else
     %fprintf(fid,sfmtv{i},strjust(cnames(i,:)));
     %end;
    end;
     end;
      end;
   % fprintf(fid,'\n');
    % fprintf(fid,'\n------------------------------------------------------------\n');
   fprintf(fid,['\n' LL_ '\n']);
   end; % end of if-else nfmts
   
   for k = begr:endr; % print row labels and numbers in matrix
    if rnum == 1; fprintf(fid,'%5d ',k);
          elseif rflag == 1        
    fprintf(fid,rfmt,rnames(k+1,:));
          end;
    for l = (j-1)*nwide+begc:j*nwide+begc-1
     if l <= endc
      if nfmts == 1
      fprintf(fid,ffmt,y(k,l));
      else
      fprintf(fid,ffmtv{l},y(k,l));
      end;
     end;
    end; % end of for l
    fprintf(fid,'\n');
   end; % end of for k
  % fprintf(fid,'\n');
  % fprintf(fid,'\n------------------------------------------------------------\n');
   % fprintf(fid,[LL_ '\n']);
  end; % end of for j
end




















