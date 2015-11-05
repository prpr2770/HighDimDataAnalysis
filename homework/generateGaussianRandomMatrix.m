function GR = generateGaussianRandomMatrix(n)
% Q26.b) Random Matrices:: Gaussian 

% create diagonal first
dia = 0 + 2*randn(1,n);
diaM = diag(dia);

% creating off-diagonal parts
one = ones(n);
oneU = ~tril(one);
oneM = oneU + oneU';
bodyM = randn(n) .* oneM;

% combining all parts
GR = diaM + bodyM;

% % % Assumption: The diagonal and off-diagonal entries are sampled from
% % % std.normal.
% % GR = randn(n);
end