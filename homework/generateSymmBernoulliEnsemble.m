function SBE = generateSymmBernoulliEnsemble(n,p)
% Q22.a) Random Wigner Matrices:: Symmetric Bernoulli Ensemble with p=0.5
% create diagonal first
dia = rand(1,n);
diaPos = dia > p;
diaNeg = dia <= p;
diaBin = diaPos - diaNeg;
diaM = diag(diaBin);

% creating off-diagonal parts
one = ones(n);
oneU = ~tril(one);

body = rand(n);
bodyPos = (body > p);
bodyNeg = (body <= p);

bodyU = oneU.*(bodyPos - bodyNeg); % extract upper-diagonal
bodyM = bodyU + bodyU';

% combining all parts

SBE = diaM + bodyM;
end
