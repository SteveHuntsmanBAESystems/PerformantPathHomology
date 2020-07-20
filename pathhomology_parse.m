function pathhomology_parse(D,varargin)
% Preprocesses a digraph for use in pathhomology.m.
% For further documentation, see pathhomology.m.
% 
% Performs some cosmetic alterations to digraph D, and removes any loops
% that exist.
% 
% Copyright (c) 2020, BAE Systems. All rights reserved.
%
% This work is licensed under the Creative Commons
% Attribution-NonCommercial-ShareAlike 4.0 International License. To view a
% copy of this license, visit
% http://creativecommons.org/licenses/by-nc-sa/4.0/ or send a letter to
% Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.

global n

% If D.Nodes.Name does not exist, make a simple lex ordered one
if ~any(cellfun(@numel,strfind(fieldnames(D.Nodes),'Name')))
    D.Nodes.Name = cellstr(dec2base((1:n)',10,ceil(log10(n))));
end
% Enforce lex order for convenience
D = reordernodes(D,sort(D.Nodes.Name));
% Warn if there are any loops and then remove them
loops = findedge(D,1:n,1:n);
if any(loops)
    warning('removing loops!');
    D = rmedge(D,1:n,1:n);
end
% Set any weights to unity
if any(cellfun(@numel,strfind(fieldnames(D.Edges),'Weight')))
    D.Edges.Weight = ones(size(D.Edges,1),1);
end

end
