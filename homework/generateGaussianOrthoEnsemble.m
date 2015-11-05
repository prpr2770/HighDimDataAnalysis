function GOE = generateGaussianOrthoEnsemble(n)
% Q22.b) Random Wigner Matrices:: Gaussian Orthogonal Ensembles(GOE)

% create diagonal first
dia = 0 + 2*randn(1,n);
diaM = diag(dia);

% creating off-diagonal parts
one = ones(n);
oneU = ~tril(one);
body = randn(n) .* oneU;
bodyM = body + body';

% combining all parts

GOE = diaM + bodyM;

end