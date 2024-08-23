function [] = mat2latex(M,D)
% D = number of decimals
if nargin < 2
  D = 4;
end

% Get matrix dimensions
[nR,nC] = size(M);
% Create first line
s = sprintf('  \\begin{bmatrix}\n  ');
% Add matrix content
for k = 1:nR
    for l = 1:nC
        s = sprintf(['%s % 10.' num2str(D) 'f'], s, M(k, l)); % print 3 decimal places, align to 6 characters
        if l < nC
            s = sprintf('%s &', s);
        end
    end
    if k < nR
        s = sprintf('%s \\\\ ', s);
    end
    s = sprintf('%s\n  ', s);
end
% Add last line
s = sprintf('%s\\end{bmatrix}', s);
% Print the result
disp(s);
