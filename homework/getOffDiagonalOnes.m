function M = getOffDiagonalOnes(n)
A = ones(n);
upTri = ~tril(A);
M = upTri + upTri';
end