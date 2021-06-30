function [i,j] = triind2sub(siz,ind)
%TRIIND2SUB Convert triangular indices to matrix subscripts
%   This function is for use with pdist which calculates distances between
%   pairs of points in an m-by-n matrix.  Results are returned as a row
%   vector in order
%       (2,1), (3,1), ..., (m,1), (3,2), ..., (m,2), ..., (m,m–1))
%   triind2sub takes the size of the input matrix to pdist and an array
%   of indices of the distances returned.  The result is two arrays of the
%   same dimensions of the indices containing the column and row
%   subscripts.  Since the results of pdist represent a strictly lower
%   triangular matrix, the column subscripts in i will be greater than the
%   corresponding row subscripts in j
%
%   Example:
%      square = [0 0; 0 1; 1 1; 1 0];
%      % find the distances between the corners of the unit square
%      dists = pdist(square);
%      % report the pairs of points with maximum separation
%      [col, row] = triind2sub(size(square), find(dists == max(dists)))
%      % report the pairs of points with minimum separation
%      [col, row] = triind2sub(size(square), find(dists == min(dists)))
%
%   The approach above can avoid the use of the squareform function while
%   using less space and avoiding duplicate results.
%
%   See also PDIST, SQUAREFORM, SUB2TRIIND.

% Copyright 2017 James Ashton

n = siz(1); % number of columns given to pdist
ntri = n * (n - 1) / 2; % number of results returned by pdist
ind = (ntri + 1) - ind; % reverse the order
j = ceil((sqrt(ind * 8 + 1) - 1) / 2); % inverse of ind = j*(j+1)/2
i = n + 1 + j .* (j - 1) ./ 2 - ind;
j = n - j;
end