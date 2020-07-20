function [hom,betti] = pathhomology_homology(bdry,omega,varargin)
% Computes homology representatives for use in pathhomology.m.
% For further documentation, see pathhomology.m.
% 
% Input:
%     bdry: cell of matrices, representing the chain boundary operator
%     omega: cell of matrices, representing the chain complex
% Output:
%     hom: cell of matrices, representing the homology groups
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
% Input flag for symbolic computation
symbolic = ismember('symbolic',varargin);

%% Compute homology representatives
% Compute the kernel of each boundary map
ker = cell(1,pmax+1);
for p = 0:pmax
    if symbolic
        ker{p+1} = null(bdry{p+1});
        % MATLAB dislikes empty symbolic null spaces.
        if isempty(ker{p+1})
            ker{p+1} = sym(null(double(bdry{p+1}),'r'));
        end
    else
        ker{p+1} = null(bdry{p+1},'r');
    end
end
% Compute homology representatives
hom = cell(1,pmax);
for p = 0:pmax-1
    if ~size(omega{p+1},2)
        % No kernel to quotient
        hom{p+1} = zeros(size(omega{p+1},1),0);
        if symbolic
            hom{p+1}=sym(hom{p+1});
        end
        continue
    end
    if ~size(omega{p+2},2)
        % No image to quotient by
        hom{p+1} = omega{p+1}*ker{p+1};
        if symbolic
            hom{p+1}=sym(hom{p+1});
        end
        continue
    end
    % Quotient computed as cokernel of boundary operator restricted
    % into the kernel (the cokernel is computed as the null space of
    % the transpose).
    if symbolic
        hom{p+1} = null(bdry{p+2}'*ker{p+1});
        if isempty(hom{p+1})
            % MATLAB dislikes empty symbolic null spaces.
            hom{p+1} = sym(null(double(bdry{p+2}'*ker{p+1}),'r'));
        end
    else
        hom{p+1} = null(bdry{p+2}'*ker{p+1},'r');
    end
    % Convert from kernel basis to invariant basis
    hom{p+1} = ker{p+1}*hom{p+1};
    % Convert from invariant basis to allowed basis
    hom{p+1} = omega{p+1}*hom{p+1};
end

%% Compute Betti numbers
betti = NaN(1,pmax);
for p = 0:p
    betti(p+1) = size(hom{p+1},2);
end

end