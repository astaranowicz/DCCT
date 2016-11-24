function string = strjoin(items, adjoiner)
% Concatenate a cell array of strings.
%
% Input arguments:
% items:
%    a vector of items to join
% adjoiner:
%    string separating each neighboring element
%
% See also: cell2mat

% Copyright 2008-2011 Levente Hunyadi

error(nargchk(1, 2, nargin, 'struct'));
assert(isempty(items) || isvector(items), ...
    'string:strjoin:DimensionMismatch', ...
    'Only vectors are supported for string concatenation but argument has dimensions [%s].', int2str(size(items)));
if nargin > 1 && ~isempty(adjoiner)
    validateattributes(adjoiner, {'char'}, {'row'});
else
    adjoiner = '';
end

if iscellstr(items)
    p = items;
elseif isnumeric(items)
    if all(items == floor(items))
        p = arrayfun(@int2str, items, 'UniformOutput', false);
    else
        p = arrayfun(@num2str, items, 'UniformOutput', false);
    end
else
    error('string:strjoin:ArgumentTypeMismatch', ...
        'Type "%s" is not supported for string concatenation.', class(items));
end

% arrange substrings into cell array of strings
parts = cell(1,2*numel(items)-1);  % must be row vector
parts(1:2:end) = p;
parts(2:2:end) = {adjoiner};

% concatenate substrings preserving spaces
string = cell2mat(parts);