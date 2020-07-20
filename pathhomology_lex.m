function i = pathhomology_lex(paths, p)
% Generates lexicographical index for use in pathhomology.m.
% For further documentation, see pathhomology.m.
% 
% Input:
%     paths: matrix whose rows are paths (each entry being a node index)
%     p: length of paths
% Output:
%     i: column vector whose entries are integers (obtained by treating the
%         paths as base n integers)
% 
% Copyright (c) 2020, BAE Systems. All rights reserved.
%
% This work is licensed under the Creative Commons
% Attribution-NonCommercial-ShareAlike 4.0 International License. To view a
% copy of this license, visit
% http://creativecommons.org/licenses/by-nc-sa/4.0/ or send a letter to
% Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.

global n
i = 1+(paths-1)*(n.^(p:-1:0)');
end