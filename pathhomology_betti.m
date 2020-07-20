function betti = pathhomology_betti(bdry,varargin)
% Computes Betti numbers for use in pathhomology.m.
% For further documentation, see pathhomology.m.
% 
% Input:
%     bdry: cell of matrices representing the chain boundary operator
% Output:
%     betti: row vector of Betti numbers
% 
% Copyright (c) 2020, BAE Systems. All rights reserved.
%
% This work is licensed under the Creative Commons
% Attribution-NonCommercial-ShareAlike 4.0 International License. To view a
% copy of this license, visit
% http://creativecommons.org/licenses/by-nc-sa/4.0/ or send a letter to
% Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.

global pmax

%% Compute dimensions of boundary map images and kernels
dim_im = zeros(1,pmax+1);
dim_ker = zeros(1,pmax+1);
for p = 0:pmax
    % Kill zero rows before computing rank (for speedup)
    B = bdry{p+1}(any(bdry{p+1},2),:);
    dim_im(p+1) = rank(full(B));
    % Use the rank-nullity theorem
    dim_ker(p+1) = size(bdry{p+1},2)-dim_im(p+1);
end

%% Compute Betti numbers
betti = dim_ker(1:pmax)-dim_im(2:end);

end