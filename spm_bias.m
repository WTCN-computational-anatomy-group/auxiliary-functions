function varargout = spm_bias(X, varargin)
%__________________________________________________________________________
%
% Estimate a spatial bias field by fitting a GMM to the data.
%
% FORMAT [B,X,...] = spm_bias(X,...)
% 
% MANDATORY
% ---------
% X - [LAT]xP Multichannel image or volume
%
% KEYWORD
% -------
% VoxelSize  - 1x[D] Vector of voxel sizes [1] (mm/vox)
% FWHM       - 1x[D] Vector of FWHM to choose the number of bases [60] (mm)
% NbBases    - 1x[D] Number of bases along each dimension [0=use FWHM]
% RegParam   - 1x3   Regularisation parameters [0 0 10] (abs/mem/ben)
% GMM        - Cell of GMM options
% IterMax    - Maximum number of EM iterations [1000]
% Tolerance  - Convergence criterion (~ lower bound gain) [1e-4]
% BinWidth   - 1x[P] Bin width (histogram mode: add bits of variance) [0]
% InputDim   - Input space dimension [0=try to guess]
% Verbose    - Verbosity level: [0]=quiet, 1=write, 2=plot, 3=plot more
%
% OUTPUT
% ------
% B    - [LAT]xP   bias field
% X    - [LAT]xP   corrected input volume
% ...  - Same output as spm_gmm, in the same order > help spm_gmm
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

% -------------------------------------------------------------------------
% Parse inputs
p = inputParser;
p.FunctionName = 'spm_bias';
p.addRequired('X',                      @isnumeric);
p.addParameter('VoxelSize',  1,         @isnumeric);
p.addParameter('FWHM',       60,        @isnumeric);
p.addParameter('NbBases',    0,         @isnumeric);
p.addParameter('RegParam',   1E7,       @isnumeric);
p.addParameter('GMM',        {},        @iscell);
p.addParameter('IterMax',    1000,      @(X) isscalar(X) && isnumeric(X));
p.addParameter('SubIterMax', 1000,      @(X) isscalar(X) && isnumeric(X));
p.addParameter('LineSearch', 6,         @(X) isscalar(X) && isnumeric(X));
p.addParameter('Tolerance',  1e-4,      @(X) isscalar(X) && isnumeric(X));
p.addParameter('BinWidth',   0,         @isnumeric);
p.addParameter('JointOptim', true,      @(X) isscalar(X) && islogical(X));
p.addParameter('InputDim',   0,         @(X) isscalar(X) && isnumeric(X));
p.addParameter('Verbose',    0,         @(X) isscalar(X) && (isnumeric(X) || islogical(X)));
p.parse(X, varargin{:});
vs         = p.Results.VoxelSize;
fwhm       = p.Results.FWHM;
nbbases    = p.Results.NbBases;
RegParam   = p.Results.RegParam;
GMMopt     = p.Results.GMM;
P          = p.Results.InputDim;
BinWidth   = p.Results.BinWidth;
IterMax    = p.Results.IterMax;
SubIterMax = p.Results.SubIterMax;
Verbose    = p.Results.Verbose;
Tolerance  = p.Results.Tolerance;
LineSearch = p.Results.LineSearch;
JointOptim = p.Results.JointOptim;

if ~isempty(GMMopt) && ~ischar(GMMopt{1})
    K = GMMopt{1};
    GMMopt = GMMopt(2:end);
else
    K = 2;
end

% -------------------------------------------------------------------------
% Prepare input
[P,latX] = dimFromObservations(P, X);
dim      = numel(latX);
C        = spm_gmm_lib('obs2code', reshape(X, [], P));     % Code image
L        = unique(C);                                      % List of codes

% -------------------------------------------------------------------------
% Prepare bases
if sum(nbbases) == 0
    nbbases = spm_bias_lib('fwhm2nbcomp', latX, vs, fwhm);
end
[Bx,By,Bz] = spm_bias_lib('dcbasis', latX, nbbases);
switch dim
    case 2, bases = {Bx,By};
    case 3, bases = {Bx,By,Bz};
end
ICO   = spm_bias_lib('regulariser', 'bending', latX, nbbases, vs);
coeff = zeros([nbbases P], 'like', X);
field = spm_bias_lib('reconstruct', bases, coeff, 'mult');

% -------------------------------------------------------------------------
% GMM prior
MU0 = zeros(P,K);
BX  = reshape(field.*X, [], P);
minval = min(BX, [], 'omitnan');
maxval = max(BX, [], 'omitnan');
clear BX
for p=1:P
    tmp = linspace(minval(p), maxval(p), K+1);
    MU0(p,:) = (tmp(1:end-1) + tmp(2:end))/2;
end
A0 = diag(((maxval - minval)./(2*K)).^(-2));
A0 = repmat(A0, [1 1 K]);
n0 = 10 * ones(1,K);
b0 = 10 * ones(1,K);
V0 = bsxfun(@rdivide, A0, reshape(n0, [1 1 K]));
a0 = 0 * ones(1,K);

mean    = {MU0, b0};
prec    = {V0, n0};
cluster = {mean prec};
a       = a0;
PI      = repmat(1/K, [1 K]);
Z       = [];

% -------------------------------------------------------------------------
% EM loop
lb  = struct('sum', [], 'last', NaN, 'X', NaN, 'XB', 0, 'B', 0, 'Z', NaN, 'P', NaN, 'MU', NaN, 'A', NaN);
ThisRegParam = RegParam(1);
RegParam     = RegParam(2:end);
for em=1:IterMax
    
    
    % ---------------------------------------------------------------------
    % Save previous lower bound value
    obj0 = lb.last;
    
    % ---------------------------------------------------------------------
    % Optimise GMM
    [Z,cluster,prop,lb] = spm_gmm_loop(reshape(field.*X, [], P), ...
        cluster, {'Prop', PI, 'Dir', a}, ...
        'PropPrior',      a0, ...
        'LowerBound',     lb, ...
        'Missing',        true, ...
        'MissingCode',    {C,L}, ...
        'IterMax',        IterMax, ...
        'Tolerance',      Tolerance, ...
        'BinUncertainty', (bsxfun(@times, reshape(field, [], P), reshape(BinWidth, 1, [])).^2)/12, ...
        'Verbose',        Verbose, ...
        GMMopt{:});
%         'GaussPrior',     {MU0, b0, V0, n0}, ...
%         'SubIterMax',     SubIterMax, ...
%         'SubTolerance',   SubTolerance, ...
    MU = cluster.MU;
    b  = cluster.b;
    A  = cluster.A;
    V  = cluster.V;
    n  = cluster.n;
    PI = prop.Prop;
    a  = prop.Dir;
    if sum(b) > 0, mean = {MU b};
    else,          mean = {MU};
    end
    if sum(n) > 0, prec = {V n};
    else,          prec = {A};
    end
    cluster = {mean prec};
    
    % ---------------------------------------------------------------------
    % Optimise Bias Field
    [field,coeff,lb,ok] = spm_bias_loop(X, Z, cluster, bases, ...
        'Coefficients',   coeff, ...
        'BiasField',      field, ...
        'LowerBound',     lb, ...
        'RegPrecision',   ICO, ...
        'RegParam',       ThisRegParam, ...
        'MissingCode',    {C,L}, ...
        'IterMax',        IterMax, ...
        'Tolerance',      Tolerance, ...
        'LineSearch',     LineSearch, ...
        'BinWidth',       BinWidth, ...
        'JointOptim',     JointOptim, ...
        'Verbose',        Verbose);
    if ~ok
        break;
    end
    
    % ---------------------------------------------------------------------
    % Check convergence
    obj  = lb.last;
    den  = max(lb.sum(:))-min(lb.sum(:));
    gain = check_convergence(obj, obj0, den, em, Verbose);
    if gain < Tolerance
        if ~isempty(RegParam)
            PrevRegParam = ThisRegParam;
            ThisRegParam = RegParam(1);
            RegParam     = RegParam(2:end);
            lb.B(end+1)  = lb.B(end) * ThisRegParam / PrevRegParam;
            lb.X(end+1)  = lb.X(end);
            lb.XB(end+1) = lb.XB(end);
            lb.Z(end+1)  = lb.Z(end);
            lb.P(end+1)  = lb.P(end);
            lb.MU(end+1) = lb.MU(end);
            lb.A(end+1)  = lb.A(end);
        else
            break
        end
    end
end

% -------------------------------------------------------------------------
% Prepare output
Z = reshape(Z, [latX K]);
varargout{1}  = field;
varargout{2}  = field.*X;
varargout{3}  = Z;
varargout{4}  = MU;
varargout{5}  = A;
varargout{6}  = PI;
varargout{7}  = b;
varargout{8}  = V;
varargout{9}  = n;
varargout{10} = a;

% =========================================================================
function gain = check_convergence(obj, obj0, den, it, verbose)
% FORMAT gain = check_convergence(obj, obj0, den, it, verbose)

gain = abs(obj - obj0)/den;
if verbose >= 1
    if     isnan(obj - obj0),           incr = '';
    elseif obj > obj0,                  incr = '(+)';
    elseif obj < obj0,                  incr = '(-)';
    else,                               incr = '(=)';
    end
    fprintf('%-5s | %4d | lb = %-10.6g | gain = %-10.4g | %3s\n', 'b+g', it, obj, gain, incr);
end

% =========================================================================
function [P,latX] = dimFromObservations(P, X)
dimX = size(X);
if P == 0
    if numel(dimX) == 2
        if dimX(1) == 1
            % row-vector case
            P = dimX(1);
        else
            % matrix case
            P = dimX(2);
        end
    else
        % N-array case
        P = dimX(end);
    end
end
if P == 1
    latX = dimX;
else
    latX = dimX(1:end-1);
end
