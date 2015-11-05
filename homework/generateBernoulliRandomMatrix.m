function BR = generateBernoulliRandomMatrix(n,p)
% Q26.a) Random Bernoulli Matrices

body = rand(n);
bodyPos = (body > p);
bodyNeg = (body <= p);

BR = (bodyPos - bodyNeg); % extract upper-diagonal

end
