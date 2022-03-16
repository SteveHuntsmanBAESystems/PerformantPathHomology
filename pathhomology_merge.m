function y = pathhomology_merge(yi,bins,varargin)
% Merges multiple homologies for use in pathhomology.m.
% For further documentation, see pathhomology.m.
% 
% Input:
%     yi: cell containing structs, representing homologies
%     bins: row vector describing which component each vertex is in
% Output:
%     y: struct containing data accumulated from input homologies
% 
% Copyright (c) 2020, BAE Systems. All rights reserved.
%
% This work is licensed under the Creative Commons
% Attribution-NonCommercial-ShareAlike 4.0 International License. To view a
% copy of this license, visit
% http://creativecommons.org/licenses/by-nc-sa/4.0/ or send a letter to
% Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.

global pmax n
C = length(yi);
y = struct();
% Helper method to parse a path as a base-n integer
lexValue = @pathhomology_lex;

%% Find proper indices for inclusion
% Indexed inclusion for elementary paths (in lexicographic order)
inclLex = cell(C,pmax);
for i = 1:C
    component = find(bins==i);
    % For p == 0, inclusion is trivial
    inclLex{i,1} = component';
    for p = 1:pmax
        % Build inclusion in recursive construction (similar to a fractal)
        suffixInd = repmat(inclLex{i,p},length(component),1);
        prefixInd = repelem((component'-1)*n^p,length(inclLex{i,p}));
        inclLex{i,p+1} = prefixInd+suffixInd;
    end
end
y.allPaths = cell(1,pmax+1);
% Indexed inclusion for allowed paths
inclAll = cell(C,pmax+1);
for p = 0:pmax
    % Accumulate all paths into one list
    y.allPaths{p+1} = NaN(0,p+1);
    for i = 1:C
        % Keep track of indices for each component
        inclAll{i,p+1} = (1:size(yi{i}.allPaths{p+1},1))+...
            size(y.allPaths{p+1},1);
        newPaths = arrayfun(@(x) inclLex{i,1}(x),yi{i}.allPaths{p+1});
        y.allPaths{p+1} = [y.allPaths{p+1};newPaths];
    end
    % Sort list and get indexed inclusion
    [y.allPaths{p+1}, indices] = sortrows(y.allPaths{p+1});
    indices(indices) = 1:length(indices);   % bug fix 15.03.2022 from MY
    for i = 1:C
        inclAll{i,p+1} = indices(inclAll{i,p+1});
    end
end

%% Merge data from disconnected components
y.bdryA = cell(1,pmax+1);
y.omega = cell(1,pmax+1);
y.bdry = cell(1,pmax+1);
for p = 0:pmax
    % initialize accumulating arrays
    y.bdryA{p+1} = sparse(n^p,size(y.allPaths{p+1},1));
    y.omega{p+1} = zeros(size(y.allPaths{p+1},1),0);
    y.bdry{p+1} = [];
    for i = 1:C
        % Accumulate data (allowed paths are in lex order, invariant paths
        % are ordered by component)
        if p >= 1
            y.bdryA{p+1}(inclLex{i,p},inclAll{i,p+1}) = yi{i}.bdryA{p+1};
        end
        omegaNew = zeros(size(y.allPaths{p+1},1),size(yi{i}.omega{p+1},2));
        omegaNew(inclAll{i,p+1},:) = yi{i}.omega{p+1};
        y.omega{p+1} = [y.omega{p+1} omegaNew];
        if p == 0
            y.bdry{p+1} = [y.bdry{p+1} yi{i}.bdry{p+1}];
        else
            y.bdry{p+1} = blkdiag(y.bdry{p+1},yi{i}.bdry{p+1});
        end
    end
    if p == 0
        % Exception necessary because no inclLex{i,0}
        y.bdryA{p+1} = sparse(ones(1,size(y.allPaths{p+1},1)));
    end
    % Cast to symbolic if requested (there's probably a better way to do
    % this)
    if ismember('symbolic',varargin)
        y.omega{p+1} = sym(y.omega{p+1});
        y.bdry{p+1} = sym(y.bdry{p+1});
    end
end

y.betti = zeros(1,pmax);
for i = 1:C
    y.betti = y.betti + yi{i}.betti;
end
% Definition of Betti 0
y.betti(1) = C-1;

if ismember('representatives',varargin)
    y.hom = cell(1,pmax);
    for p = 0:pmax-1
        % Initialize accumulator
        y.hom{p+1} = zeros(size(y.allPaths{p+1},1),0);
        for i = 1:C
            % Accumulate representatives ordered by component
            homNew = zeros(size(y.allPaths{p+1},1),size(yi{i}.hom{p+1},2));
            homNew(inclAll{i,p+1},:) = yi{i}.hom{p+1};
            y.hom{p+1} = [y.hom{p+1} homNew];
        end
        % Generate representatives for dimension 0
        if p == 0
            for i = 2:C
                homNew = zeros(size(y.allPaths{p+1},1),1);
                homNew(find(bins==1,1)) = 1;
                homNew(find(bins==i,1)) = -1;
                y.hom{p+1} = [y.hom{p+1} homNew];
            end
        end
        % Cast to symbolic if requested (there's probably a better way to do
        % this)
        if ismember('symbolic',varargin)
            y.hom{p+1} = sym(y.hom{p+1});
        end
    end
end

end
