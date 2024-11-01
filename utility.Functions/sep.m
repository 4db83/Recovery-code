function [] = sep(LN,symbl,NL,FileID)
%F: prints line separtions.
% FileID

width = 181;

if nargin < 3 || isempty(NL)
  NL = '';
else
  if NL == 1
    NL = '\n';
  else
    NL = '';
  end
end

if nargin == 1
if LN == 0
	width = 150;
	fprintf([NL repmat('-',1,width) '\n']);	
end
end

if nargin < 1
	fprintf([NL repmat('-',1,width) '\n']);	
elseif nargin < 2 || isempty(symbl)
	if ischar(LN)
		fprintf([NL repmat(LN,1,width) '\n']);
	else
		fprintf([NL repmat('-',1,LN) '\n']);	
	end
elseif nargin < 4 || isempty(FileID)
	fprintf([NL repmat(symbl,1,LN) '\n']);	
else 
	fprintf(FileID,[NL repmat(symbl,1,LN) '\n']);	
end


% if nargin < 3
% 	fprintf([repmat(symbl,1,LN) '\n']);	
% else
% 	fprintf(FileID,[repmat(symbl,1,LN) '\n']);
% 	disp('ehll')
% end

% % 	fprintf(FileID,[repmat(symbl,1,LN) '\n']);
% % if nargin < 3
% % 	if isempty(symbl)
% % 		symbl = '-';
% % 	end;
% % 	fprintf([repmat(symbl,1,LN) '\n']);
% % elseif nargin < 2
% % 	if ischar(LN)
% % 		fprintf([repmat(LN,1,70) '\n']);
% % 	else
% % 		fprintf([repmat('-',1,LN) '\n']);
% % 	end
% % else
% % % elseif nargin == 0
% % % 	symbl = '-';
% % % 	LN = 70;
% % % 	fprintf(FileID,[repmat(symbl,1,LN) '\n']);
% % % end
% % 	fprintf(FileID,[repmat(symbl,1,LN) '\n']);
% % end