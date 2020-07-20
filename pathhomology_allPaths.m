function [allPaths,indA,leaves] = pathhomology_allPaths(D,varargin)
% Computes allowed paths for use in pathhomology.m.
% For further documentation, see pathhomology.m.
% 
% Input:
%     D: digraph
% Output:
%     allPaths: cell of matrices, each having rows representing allowed
%         paths (entries being vertex indices)
%     indA: cell of row vectors, each representing the inclusion of allowed
%         paths into elementary paths.
%     leaves: logical row vector representing the vertices pruned as leaves
% 
% Copyright (c) 2020, BAE Systems. All rights reserved.
%
% This work is licensed under the Creative Commons
% Attribution-NonCommercial-ShareAlike 4.0 International License. To view a
% copy of this license, visit
% http://creativecommons.org/licenses/by-nc-sa/4.0/ or send a letter to
% Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.

global pmax n
A = adjacency(D);

%% Recursively find leaves
leaves = false(1,n);
if ~ismember('keepLeaves',varargin)
    kill = 1;
    while any(kill)
        % restrA is the pruned version of the graph
        restrA = (A.*~leaves).*~leaves';
        % Total (in+out) degrees of each vertex
        degree = sum(restrA,1)+sum(restrA,2)';
        % Newly identified leaves
        kill = (degree==1)&~leaves;
        % When two leaves identified at once have an edge, they are the last
        % two vertices in their connected component, so they should not both be
        % pruned.
        kill = kill&~any(restrA(kill,:),1);
        leaves = leaves|kill;
    end
end
nonleaves = find(~leaves);

%% Generate allowed paths
% All allowed paths (double entendre)
allPaths = cell(1,pmax+1);
for p = 0:pmax
    % For length-0 paths, there are no previous paths to build from
    if p == 0
        allPaths{p+1} = nonleaves';
        continue
    end
    % Preallocate space for more paths than necessary
    allPaths{p+1} = NaN((n-1)*size(allPaths{p},1),p+1);
    i2 = 1;
    for i1 = 1:size(allPaths{p},1)
        % Build allowed paths from the shorter paths
        for v = nonleaves
            if ~A(allPaths{p}(i1,p),v)
                continue
            end
            allPaths{p+1}(i2,:) = [allPaths{p}(i1,:) v];
            i2 = i2 + 1;
        end
    end
    % Trim extra preallocated space
    allPaths{p+1} = allPaths{p+1}(1:i2-1,:);
end

%% Indices of allowed p-paths in set of elementary p-paths
indA = cell(1,pmax+1);
% Helper method to parse a path as a base-n integer
lexValue = @pathhomology_lex;
for p=0:pmax
    if ~isempty(allPaths{p+1})
        indA{p+1} = lexValue(allPaths{p+1},p);
    end
end