function [V,W,C,BW,E] = spm_histN(X,varargin)
% _________________________________________________________________________
%
% Compute the (joint) histogram of a (multidimensional) dataset
%
% FORMAT [values,count] = spm_histN(X,nbins,...)
% FORMAT [values,count] = spm_histN(X,centres,...)
%
% MANDATORY
% ---------
% X        - NxP matrix of observed values
% 
% OPTIONAL
% --------
% nbins    - 1x1 or 1xP number of bins [64]
%   or
% centres  - Mx1 ordered bin centres (or 1xP cell of Mpx1 bin centres)
%
% KEYWORD
% -------
% Weights  - Nx1 Observations weights [1]
% KeepZero - Keep bins with zero observations [true]
% Missing  - Keep rows with missing data [false]
%            Additional bins are created for missing values.
% Reshape  - Reshape output so that their lattice is M1xM2x... [false]
% Smooth   - FWHM of the smoothing kernel (in bins) [0]
% Verbose  - Verbosity level [0]
%
% OUTPUT
% ------
% values  - MxP matrix of multidimensional values (or M1x...xP if Reshape)
% count   - Mx1 vector of counts                  (or M1x...x1 if Reshape)
% centres - 1xP cell of Mpx1 bin centres
% widths  - 1xP vector (or 1xP cell of Mpx1 vectors) bin widths
%
% (M can be smaller that the specified number of bins if KeepZero = false)
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging
    
% -------------------------------------------------------------------------
% Parse inputs
p = inputParser;
p.FunctionName = 'spm_histN';
p.addRequired('X',                @isnumeric);
p.addOptional('B',         64,    @(arg) isnumeric(arg) || iscell(arg));
p.addParameter('Weights',  1,     @(arg) isscalar(arg)  || (numel(arg) == size(X,1)));
p.addParameter('KeepZero', true,  @isscalar);
p.addParameter('Missing',  false, @isscalar);
p.addParameter('Reshape',  false, @isscalar);
p.addParameter('Smooth',   0,     @isnumeric);
p.addParameter('Verbose',  0,     @isscalar);
p.parse(X, varargin{:});
B  = p.Results.B;
Wi = p.Results.Weights;
Wi = Wi(:);

if p.Results.Reshape && ~p.Results.KeepZero
    error(['To reshape, all values must be kept. ' ...
           'Use (''KeepZero'',true) with (''Reshape'',true).'])
end

% -------------------------------------------------------------------------
% Discard missing values
if ~p.Results.Missing
    missing = any(isnan(X),2);
    X       = X(~missing,:);
    if ~isscalar(Wi)
        Wi  = Wi(~missing);
    end
end

% -------------------------------------------------------------------------
% Compute bin centres / edges
P      = size(X,2); % Number of channels
N      = size(X,1); % Number of observations
minval = min(X, [], 1, 'omitnan'); % Min value / channel
maxval = max(X, [], 1, 'omitnan'); % Max value / channel
if ~iscell(B) && size(B,1) == 1
% Number of bins provided
    E = B;
    if numel(B) < P
        E = spm_padarray(E, [0, P-numel(B)], 'replicate', 'post');
    end
    BW = (maxval - minval)./B;
    E  = num2cell(E);
else
% Bin centres provided
    if ~iscell(B)
        if size(B,2) == 1
            B = repmat(B(:),1,P);
        end
        B = num2cell(B, 1);
    end
        
    E  = cell(1,P);
    BW = cell(1,P);
    for c=1:P
        E{c}  = (B{c}(2:end) + B{c}(1:end-1))/2;
        E{c}  = [minval(c); E{c}; maxval(c)]';
        BW{c} = E{c}(2:end) - E{c}(1:end-1);
    end
end
clear B

% -------------------------------------------------------------------------
% Discretize data
I  = cell(1,P);      % Bin index for each observation
V  = cell(1,P);      % Bin edges
dim = zeros(1,P);
hasnan = zeros(1,P,'logical');
for c=1:P
    [I{c},V{c}]       = discretize(X(:,c),E{c}); % /!\ in a toolbox?
    I{c}              = single(I{c});
    I{c}(isnan(I{c})) = numel(V{c});
    V{c}              = (V{c}(2:end) + V{c}(1:end-1))/2; % Edge to centre
    V{c}              = V{c}(:);
    hasnan(c)         = any(isnan(X(:,c)));
    dim(c)            = numel(V{c}) + hasnan(c);
    if hasnan(c)
        V{c}(end+1)   = NaN;
    end
end
clear E X

% -------------------------------------------------------------------------
% Count
if numel(dim) == 1
    linI = [I{:}];
else
    linI = sub2ind(dim, I{:}); % Convert P-dimensional index to (unique) linear index
end
linI = uint64(linI);
clear I
if isscalar(Wi)
    Wi = Wi * ones(size(linI));
end
Wi = double(Wi);
W  = spm_hist(linI, Wi, prod(dim)); % TODO: push new spm_hist
W  = W(:);
clear linI Wi

% -------------------------------------------------------------------------
% Create matrix of multivariate bin centres
C  = V;
V  = mycombvec(V{:});

% -------------------------------------------------------------------------
% Smooth
if p.Results.Smooth
    W = reshape(W, dim);
    lim = ceil(4/2.355*p.Results.Smooth);
    ker = spm_smoothkern(p.Results.Smooth, -lim:lim, 0);
    ker = ker(ker~=0);
    for c=1:P
        if hasnan(c)
            W1        = W;
            subs      = cell(1,P);
            [subs{:}] = deal(':');
            subs{c}   = 1:size(W,c)-1;
            W = subsref(W1, struct('type', '()', 'subs', {subs}));
        end
        W = convn(W, reshape(ker, [ones(1,c-1) numel(ker) 1]), 'same');
        if hasnan(c)
            [W1,W] = deal(W,W1);
            W = subsasgn(W, struct('type', '()', 'subs', {subs}), W1);
            clear W1
        end
        
    end
    W = W(:);
end

% -------------------------------------------------------------------------
% Reshape
if p.Results.Reshape
    W = reshape(W, dim);
    V = reshape(V, [dim P]);
end

% -------------------------------------------------------------------------
% Remove empty bins
if ~p.Results.KeepZero
    empty = W == 0;
    W     = W(~empty);
    V     = V(~empty,:);
end

function x = mycombvec(varargin)
% Create all possible combinations of input values.
%
% FORMAT x = mycombvec(y1,y2,...,yN)
% y1,y2,...,yN  must be column vectors
% x             is an MxN matrix, where M is the number of combinations
%
% /!\ It does not have the exact same behaviour as Matlab's combvec.
%     It is equivalent to: `x = combvec(y1(:)',y2(:)',...)'`
x = varargin{1}(:);
varargin = varargin(2:end);
while ~isempty(varargin)
    y  = varargin{1}(:);
    Nx = size(x,1);
    Ny = size(y,1);
    x  = repmat(x,  [Ny 1]);
    y  = repmat(y.', [Nx 1]);
    x  = [x y(:)];
    varargin = varargin(2:end);
end