function y = pathhomology_leafPaths(D,y,leaves,varargin)
% Reincorporates pruned leaves into homology for use in pathhomology.m.
% For further documentation, see pathhomology.m.
% 
% Input:
%     D: digraph
%     y: struct representing homology
%     leaves: logical row vector describing which vertices were pruned
% Output:
%     y: struct representing homology
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
% Helper method to parse a path as a base-n integer
lexValue = @pathhomology_lex;
% Unpack struct input
allPaths = y.allPaths;
bdryA = y.bdryA;
omega = y.omega;
bdry = y.bdry;

%% Get proper indices
% Indexed inclusion of old allowed paths into all allowed paths
inclOld = cell(1,pmax+1);
% Indexed inclusion of new allowed paths into all allowed paths
inclNew = cell(1,pmax+1);
% Source (shorter path) for each new allowed path
sourceNew = cell(1,pmax+1);
% New allowed paths
newPaths = cell(1,pmax+1);
for p = 0:pmax
    if p == 0
        % Each new vertex is a new path
        newPaths{p+1} = find(leaves)';
        % Source is the empty path
        sourceNew{p+1} = ones(1,length(newPaths{p+1}));
    else
        % Preallocate space for more paths than necessary
        maxLen = (n-1)*size(allPaths{p},1);
        newPaths{p+1} = NaN(maxLen,p+1);
        sourceNew{p+1} = NaN(1,maxLen);
        i2 = 1;
        % Construct paths ending in new vertices
        for i1 = 1:size(allPaths{p},1)
            currpath = allPaths{p}(i1,:);
            for v = find(leaves)
                if ~A(currpath(end),v)
                    continue
                end
                newPaths{p+1}(i2,:) = [currpath v];
                sourceNew{p+1}(i2) = i1;
                i2 = i2+1;
            end
        end
        % Construct new paths ending in old vertices
        for i1 = 1:size(newPaths{p},1)
            currpath = newPaths{p}(i1,:);
            for v = find(~leaves)
                if ~A(currpath(end),v)
                    continue
                end
                newPaths{p+1}(i2,:) = [currpath v];
                sourceNew{p+1}(i2) = inclNew{p}(i1);
                i2 = i2+1;
            end
        end
        % Trim excess preallocated space
        newPaths{p+1} = newPaths{p+1}(1:i2-1,:);
        sourceNew{p+1} = sourceNew{p+1}(1:i2-1);
    end
    % Set up indexed inclusions
    inclOld{p+1} = 1:size(allPaths{p+1},1);
    inclNew{p+1} = (1:size(newPaths{p+1},1))+size(allPaths{p+1},1);
    % allPaths in lexicographical order (could be optimized as a sorted
    % merge but MATLAB's quicksort is good enough)
    [allPaths{p+1},indices] = sortrows([allPaths{p+1};newPaths{p+1}]);
    [~,indices] = sort(indices);
    inclOld{p+1} = indices(inclOld{p+1});
    inclNew{p+1} = indices(inclNew{p+1});
end

%% Update data
for p = 0:pmax
    %% Insert columns into bdryA
    % Add empty columns (i.e., preallocate a larger sparse matrix)
    bdryTemp = spalloc(size(bdryA{p+1},1),size(allPaths{p+1},1),...
        size(allPaths{p+1},1)*(p+1));
    bdryTemp(:,inclOld{p+1}) = bdryA{p+1};
    bdryA{p+1} = bdryTemp;
    % Populate columns
    for i1 = 1:length(inclNew{p+1})
        if p == 0
            bdryA{p+1}(1,inclNew{p+1}(i1)) = 1;
            continue
        end
        % Source of new path (i.e., index in previous bdryA)
        s = sourceNew{p+1}(i1);
        % Index of new path
        ind = inclNew{p+1}(i1);
        % Copy boundary operator from source
        for j = find(bdryA{p}(:,s))
            % Obtain new lexicographical index from old one
            target = (j-1)*n+allPaths{p+1}(ind,end);
            % Construct boundary operator
            bdryA{p+1}(target,ind) = bdryA{p}(j,s);
        end
        % This last entry isn't in the source's boundary
        bdryA{p+1}(lexValue(allPaths{p}(s,:),p-1),ind) = (-1)^p;
    end
    
    %% Insert empty rows into omega
    % The paths involving leaves are never invariant, so we don't worry
    % about them
    omegaTemp = zeros(size(allPaths{p+1},1),size(omega{p+1},2));
    omegaTemp(inclOld{p+1},:) = omega{p+1};
    omega{p+1} = omegaTemp;
    
    %% Add pruned edges into omega (low dimension only)
    if p <= 1
        omegaNew = zeros(size(omega{p+1},1),length(inclNew{p+1}));
        omegaNew(inclNew{p+1},:) = eye(length(inclNew{p+1}));
        % This leaves the indexing a disorganized mess but it's easy and
        % technically correct
        omega{p+1} = [omega{p+1} omegaNew];
    end
    
    %% Fix bdry in low dimensions
    % Add blank rows
    if ismember(p,[1 2])
        bdryTemp = zeros(size(omega{p},2),size(bdry{p+1},2));
        bdryTemp(1:size(bdry{p+1},1),:) = bdry{p+1};
        bdry{p+1} = bdryTemp;
    end
    % Add columns (since omega has more vectors)
    if p <= 1
        oldIn = size(bdry{p+1},2);
        % This is a hack that only works because bdryA{p+1}'s codomain
        % (elementary paths) is equal to omega{p}'s codomain (allowed
        % paths) in low dimension; otherwise we'd have to use indA
        bdryNew = bdryA{p+1}*omega{p+1}(:,oldIn+1:end);
        if p > 0
            bdryNew = omega{p}\bdryNew;
        end
        bdry{p+1} = [bdry{p+1} bdryNew];
    end
end

% Save output values into struct
y.allPaths = allPaths;
y.bdryA = bdryA;
y.omega = omega;
y.bdry = bdry;
end