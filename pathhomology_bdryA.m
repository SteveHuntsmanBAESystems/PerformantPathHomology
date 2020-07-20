function bdryA = pathhomology_bdryA(allPaths,varargin)
% Computes the unrestricted boundary operator for use in pathhomology.m.
% For further documentation, see pathhomology.m.
% 
% Input:
%     allPaths: cell of matrices, each having rows representing allowed
%         paths (entries being vertex indices)
% Output:
%     bdryA: cell of matrices, representing unrestricted boundary operator
%         (with domain as allowed paths and codomain as elementary paths)
% 
% Copyright (c) 2020, BAE Systems. All rights reserved.
%
% This work is licensed under the Creative Commons
% Attribution-NonCommercial-ShareAlike 4.0 International License. To view a
% copy of this license, visit
% http://creativecommons.org/licenses/by-nc-sa/4.0/ or send a letter to
% Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.

global pmax n
% Helper method to parse a path as a base-n integer
lexValue = @pathhomology_lex;
% Boundary map
bdryA = cell(1,pmax+1);

for p = 0:pmax
    %% Initialize modified boundary map
    % bdryA has domain allowed p-paths and codomain elementary p-1 paths
    % (i.e., formal sums)
    bdryA{p+1} = spalloc(n^p,size(allPaths{p+1},1),...
        (p+1)*size(allPaths{p+1},1));
    %% Populate columns of the modified boundary map
    for i1 = 1:size(allPaths{p+1},1)
        curPath = allPaths{p+1}(i1,:);
        for i2 = 0:p
            % Kill i2th vertex
            temp = curPath([1:i2,(i2+2):(p+1)]);
            % Elementary p-1 paths are ordered lexicographically.
            bdryA{p+1}(lexValue(temp,p-1),i1) = (-1)^i2;
        end
    end
end
end