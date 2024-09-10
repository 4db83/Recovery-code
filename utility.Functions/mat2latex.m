function [] = mat2latex(M, D, rowNames, colNames, ColSpace)
% D = number of decimals
SetDefaultValue(2, 'D', 4)
SetDefaultValue(3, 'rowNames', [])
SetDefaultValue(4, 'colNames', [])
SetDefaultValue(5, 'ColSpace', 12)

% rowNames = {'Row1';'Row2';'Row3'};
if ~isempty(rowNames)
  rowNames = strcat(rowNames,' &');
end

[~,NrowNames] = size(char(rowNames));
if NrowNames == 0; NrowNames = ColSpace-1; end

% Get matrix dimensions
[nR,nC] = size(M);
% Create first line
s = sprintf('  \\begin{bmatrix}\n  '); 
% ColNames = [{'empty    '}; strcat({'Column '},num2str((1:nC)'))]
if isempty(colNames)
  print_colNames = strcat({'Col('},num2str((1:nC)'),')');
else
  print_colNames = colNames;
end

[~,NcolNames] = size(char(print_colNames));

% first column name
firstName = repmat(' ', 1 , NrowNames);

% extra_space = repmat(' ', 1, ColSpace - NcolNames - 1);
extra_space = repmat(' ', 1, ColSpace - NcolNames );
% size(extra_space)

% make column names first
for C = 1:nC
  % s = sprintf(['%s &     '  print_colNames{C} ], s); 
  if isempty(rowNames)
    if C == 1; 
      s = sprintf(['%s      '   print_colNames{C} ], s);  
    else
      s = sprintf(['%s  & '     print_colNames{C} ], s); 
    end
  else 
    if C == 1; s = sprintf(['%s  ' firstName ], s); end
    s = sprintf(['%s &  '  extra_space  print_colNames{C} ], s); 
  end
end
s = sprintf('%s     \\\\ \n  ', s);
% s

% Add matrix content
for R = 1:nR
    for C = 1:nC
      if C == 1
        if ~isempty(rowNames)
          s = sprintf(['%s    ' rowNames{R} ' % ' num2str(ColSpace) '.' num2str(D) 'f'], s, M(R, C)); % print 3 decimal places, align to 6 characters
        else 
          s = sprintf(['%s % ' num2str(ColSpace) '.' num2str(D) 'f'], s, M(R, C)); % print 3 decimal places, align to 6 characters
        end
      else
        s = sprintf(['%s % ' num2str(ColSpace) '.' num2str(D) 'f'], s, M(R, C)); % print 3 decimal places, align to 6 characters
      end
      if C < nC
        s = sprintf('%s &', s);
      end
    end

    if R < nR
      s = sprintf('%s       \\\\ ', s);
    end
    s = sprintf('%s\n  ', s);
end
% Add last line
s = sprintf('%s\\end{bmatrix}', s);
% Print the result
disp(s);
