function [omega,bdry] = pathhomology_omega(bdryA,indA,varargin)
% Constructs the chain (path) complex for use in pathhomology.m.
% For further documentation, see pathhomology.m.
% 
% Input:
%     bdryA: cell of matrices, representing the unrestricted boundary
%         operator
%     indA: cell of row vectors, representing the inclusion of allowed
%         paths into elementary paths
% Output:
%     omega: cell of matrices, each representing the invariant spaces (each
%         matrix is the inclusion from invariant space to allowed paths)
%     bdry: cell of matrices, representing the boundary operator
% 
% Copyright (c) 2020, BAE Systems. All rights reserved.
%
% This work is licensed under the Creative Commons
% Attribution-NonCommercial-ShareAlike 4.0 International License. To view a
% copy of this license, visit
% http://creativecommons.org/licenses/by-nc-sa/4.0/ or send a letter to
% Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.

global pmax n
% Input flags for symbolic computation and/or homology representatives
symbolic = ismember('symbolic',varargin);
reps = ismember('representatives',varargin);

%% Construct the chain complex
% Space of invariant p-paths (represented as an inclusion)
omega = cell(1,pmax+1);
% Boundary map restricted to invariant space
bdry = cell(1,pmax+1);
for p = 0:pmax
    if isempty(bdryA{p+1})
        % When there are no allowed paths, there are no invariant paths
        omega{p+1} = zeros(size(bdryA{p+1},2),0);
        % Boundary operator with trivial domain
        if p == 0
            bdry{p+1} = zeros(1,0);
        else
            bdry{p+1} = zeros(size(omega{p},2),0);
        end
        continue
    end
    %% Find invariant space
    % The invariant space is the kernel of the boundary operator projected
    % onto disallowed paths.
    if p == 0
        bdryTemp = sparse(0,size(bdryA{p+1},2));
    else
        bdryTemp = bdryA{p+1}(setdiff(1:n^p,indA{p}),:);
    end
    % Kill zero rows before computing kernel to save LOTS of time
    bdryTemp = bdryTemp(any(bdryTemp,2),:);
    if symbolic
        % Symbolic computation, if requested.
        bdryTemp = sym(bdryTemp);
        omega{p+1} = null(bdryTemp);
    else
        bdryTemp = full(bdryTemp);
        if reps
            % Representatives are cleaner with a rational basis.
            omega{p+1} = null(bdryTemp,'r');
        else
            omega{p+1} = null(bdryTemp);
        end
    end
    % Explicitly specify invariant space if it's empty (MATLAB
    % sometimes misbehaves with empty null spaces).
    if ~numel(omega{p+1})
        if symbolic
            omega{p+1}=sym(zeros(size(bdryA{p+1},2),0));
        else
            omega{p+1}=zeros(size(bdryA{p+1},2),0);
        end
    end
    
    %% Restrict boundary operator
    % Restrict domain
    bdryTemp = bdryA{p+1}*omega{p+1};
    % Restrict codomain
    if p == 0
        bdryTemp = full(bdryTemp);
    else
        bdryTemp = omega{p}\bdryTemp(indA{p},:);
    end
    bdry{p+1} = bdryTemp;
end
end