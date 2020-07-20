function y = pathhomology(D,pMax,varargin)
% Computes the path homology of a digraph.
% For background, details, etc. see
%     [PPHDN] Chowdhury, S. and Memoli, F. "Persistent path homology of
%     directed networks." SODA (2018).
%     SIAM version:	https://doi.org/10.1137/1.9781611975031.75
%     ACM mirror:	https://dl.acm.org/citation.cfm?id=3175345
%     arXiv:        https://arxiv.org/abs/1701.00565
% or https://arxiv.org/abs/1207.2834 and other papers on the subject by
% Grigor'yan, Yong, Muranov, Yau, and collaborators. 
%
% The reduced and non-reduced Betti numbers differ only in dimension 0. The
% non-reduced zeroth Betti number is just the number of weakly connected
% components of D, which is easily obtained via the command
% conncomp(D,'Type','weak'). We thank Samir Chowdhury for this observation
% (see also 3.6 of https://arxiv.org/abs/1207.2834).
%
% The use of non-regular versus regular path homology is motivated by
% example 3.14 of https://arxiv.org/abs/1207.2834 in light of potential
% applications where 2-cycles seem worthy of "hole" status, but note that
% our preference is precisely opposite to the paper's.
%
% Homology representatives (if requested) are computed using MATLAB's
% rational basis for null space. Therefore, they are somewhat sparse and
% rational, but are not necessarily as sparse as possible, nor are they
% integral.
% 
% WARNING: symbolic computation doesn't seem to play nice with
% representatives
% 
% Input:
%     D: digraph object
%     pMax: highest dimension of homology (longest paths) to compute
% Input flags:
%     'keepLeaves' prevents pruning of leaves during computation
%     'removeLeaves' cancels reintroduction of leaves (i.e. faster
%         computation for less meaningful indexing); has no effect when
%         'keepLeaves' is specified
%     'representatives' adds homology representatives (.hom field) to the
%         output
%     'suppressRecursion' prevents breakup into connected components during
%         computation
%     'symbolic' specifies use of symbolic computation wherever applicable
%         (to avoid floating-point error)
% Output struct y has the following fields:
%     allPaths: cell of matrices, each having rows representing allowed
%         paths (entries being vertex indices)
%     bdryA: cell of sparse matrices representing the unrestricted boundary
%         operator (with domain as allowed paths and codomain as elementary
%         paths)
%     omega: cell of matrices representing the invariant chain complex
%         (each matrix being the inclusion from the invariant space to the
%         space of allowed paths)
%     bdry: cell of matrices representing the chain boundary operator (with
%         domain and codomain in the invariant chain complex)
%     betti: row vector of Betti numbers
%     hom [if requested]: cell of matrices representing the homology groups
%         (each matrix being the inclusion from the homology group to the
%         space of allowed paths)
% 
% Copyright (c) 2020, BAE Systems. All rights reserved.
%
% This work is licensed under the Creative Commons
% Attribution-NonCommercial-ShareAlike 4.0 International License. To view a
% copy of this license, visit
% http://creativecommons.org/licenses/by-nc-sa/4.0/ or send a letter to
% Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.

y = struct();
% All subroutines use pmax and n so they're global
global pmax n
pmax = pMax;
n = numnodes(D);

%% Break into connected components for faster computation
if ~ismember('suppressRecursion',varargin)
    % Find connected components
    [bins,binsizes] = conncomp(D,'type','weak');
    C = length(binsizes);
    % Compute each component's homology
    if C >= 2
        yi = cell(1,C);
        for i = 1:C
            Di = rmnode(D,find(bins~=i));
            yi{i} = pathhomology(Di,pmax,varargin{:});
        end
        % n gets changed in each pathhomology call
        n = size(D.Nodes,1);
        % Merge homologies from connected components
        y = pathhomology_merge(yi,bins,varargin{:});
        return
    end
end

%% Parse node names, remove loops, normalize edge weights
pathhomology_parse(D,varargin{:});

%% Find all allowed paths (excluding leaves)
[y.allPaths,indA,leaves] = pathhomology_allPaths(D,varargin{:});

%% Construct boundary map
y.bdryA = pathhomology_bdryA(y.allPaths,varargin{:});

%% Get boundary-invariant paths and restricted boundary map
[y.omega,y.bdry] = pathhomology_omega(y.bdryA,indA,varargin{:});

%% Reintroduce leaves
if ~(ismember('removeLeaves',varargin)||ismember('keepLeaves',varargin))
    y = pathhomology_leafPaths(D,y,leaves,varargin{:});
end

%% Compute Betti numbers (and homology representatives, if requested)
if ismember('representatives',varargin)
    [y.hom,y.betti] = pathhomology_homology(y.bdry,y.omega,varargin{:});
else
    y.betti = pathhomology_betti(y.bdry,varargin{:});
end

end