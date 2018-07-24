function varargout = spm_bias_lib(varargin)
%__________________________________________________________________________
%
% Library of functions for Bias correction
%
% Bias correction is performed by optimising a GMM fit to the data.
%
%--------------------------------------------------------------------------
% Basis functions
% ---------------
%
% basis = spm_bias_lib('dcbasis',     lattice, nb_components)
% prec  = spm_bias_lib('regulariser', mode, lattice, nb_components)
% field = spm_bias_lib('reconstruct', basis, coefficients, ['mult']/'add')
% [TODO] coeff = spm_bias_lib('rescale',     coeff, centre)
%
%--------------------------------------------------------------------------
% Optimisation
% ------------
%
% [g,H] = spm_bias_lib('derivatives', p, obs, basis, resp, cluster, codes)
% ll    = spm_bias_lib('objective',   obs, resp, bias, mean, prec, codes, binvar)
% [TODO] ll    = spm_bias_lib('prior',       coeff, precision)
% 
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

if nargin == 0
    help spm_bias_lib
    error('Not enough argument. Type ''help spm_bias_lib'' for help.');
end
id = varargin{1};
varargin = varargin(2:end);
switch lower(id)
    case 'dcbasis'
        [varargout{1:nargout}] = dcbasis(varargin{:});   
    case 'reconstruct'
        [varargout{1:nargout}] = reconstruct(varargin{:});  
    case 'regulariser'
        [varargout{1:nargout}] = regulariser(varargin{:});   
    case 'derivatives'
        [varargout{1:nargout}] = derivatives(varargin{:});  
    case 'objective'
        [varargout{1:nargout}] = objective(varargin{:});           
    otherwise
        help spm_bias_lib
        error('Unknown function %s. Type ''help spm_bias_lib'' for help.', id)
end

% =========================================================================
function varargout = dcbasis(lattice, varargin)
% FORMAT [Bx,By,Bz,...] = dcbasis(lattice, nb_component)
% FORMAT [Bx,By,Bz,...] = dcbasis(lattice, voxel_size, fwhm)
%
% lattice      - Dimensions of the lattice [dx dy ...]
% nb_component - Number of basis functions along each dimension [nx ny ...]
% voxel_size   - Voxel size of the lattice [vx vy ...]
% fwhm         - Full-width half-max of the highest frequency basis (mm)
%
% If voxel_size and fwhm are provided, the number of components is chosen
% so that the full-width half-max of the highest frequency basis function 
% is smaller than fwhm. The bias field cannot model effects whose spatial
% frequency is higher than this value.
%
% If only one value is provided for nb_component, voxel_size or fwhm, the
% same value is used along all dimensions.

ndim = numel(lattice);

% -------------------------------------------------------------------------
% Compute number of components per direction
switch numel(varargin)
    case 1
        nb_component = reshape(varargin{1}, 1, []);
        if numel(nb_component) < ndim
            nb_component = padarray(nb_component, [0 ndim-numel(nb_component)], 'replicate', 'post');
        end
    case 2
        vs = reshape(varargin{1}, 1, []);
        if numel(vs) < ndim
            vs = padarray(vs, [0 ndim-numel(vs)], 'replicate', 'post');
        end
        fwhm = reshape(varargin{2}, 1, []);
        if numel(fwhm) < ndim
            fwhm = padarray(fwhm, [0 ndim-numel(fwhm)], 'replicate', 'post');
        end
        nb_component = ceil(2 * vs .* lattice ./ fwhm);
        nb_component = max(nb_component, 1);
end

% -------------------------------------------------------------------------
% Compute each basis
varargout = cell(1,min(ndim, nargout));
for d=1:min(ndim, nargout)
    varargout{d} = spm_dctmtx(lattice(d),nb_component(d));
end

% =========================================================================
function L = regulariser(mode, lattice, varargin)
% FORMAT L = regulariser(param, lattice, nb_component)
% FORMAT L = regulariser(param, lattice, voxel_size, fwhm)
%
% FORMAT L = regulariser(mode,  lattice, nb_component)
% FORMAT L = regulariser(mode,  lattice, voxel_size, fwhm)
%
% param        - Parameters for absolute, membrane and bending energies
% mode         - Name of a single energy ('absolute'/'membrane'/'bending')
% lattice      - Dimensions of the lattice [dx dy ...]
% nb_component - Number of basis functions along each dimension [nx ny ...]
% voxel_size   - Voxel size of the lattice [vx vy ...]
% fwhm         - Full-width half-max of the highest frequency basis (mm)
%
% If numerical parameters are provided, a weighted combination of the  
% three types of regularisation is returned.
% If an energy name is provided, the matrix that allows to compute it is
% returned (without weighting: the regularisation parameter should be 
% multiplied with this matrix)
%
% If voxel_size and fwhm are provided, the number of components is chosen
% so that the full-width half-max of the highest frequency basis function 
% is smaller than fwhm. The bias field cannot model effects whose spatial
% frequency is higher than this value.
%
% If only one value is provided for nb_component, voxel_size or fwhm, the
% same value is used along all dimensions.

% -------------------------------------------------------------------------
% Special case: mixture of regularisers
if ~ischar(mode)
    param = mode;
    L = 0;
    for i=1:numel(param)
        if param(i) ~=0
            switch i
                case 1
                    mode = 'absolute';
                case 2
                    mode = 'membrane';
                case 3
                    mode = 'bending';
            end
            L = L + param(i) * regulariser(mode, lattice, varargin{:});
        end
    end
    return
end


% -------------------------------------------------------------------------
% Compute number of components per direction
ndim = numel(lattice);
switch numel(varargin)
    case 1
        nb_component = reshape(varargin{1}, 1, []);
        if numel(nb_component) < ndim
            nb_component = padarray(nb_component, [0 ndim-numel(nb_component)], 'replicate', 'post');
        end
    case 2
        vs = reshape(varargin{1}, 1, []);
        if numel(vs) < ndim
            vs = padarray(vs, [0 ndim-numel(vs)], 'replicate', 'post');
        end
        fwhm = reshape(varargin{2}, 1, []);
        if numel(fwhm) < ndim
            fwhm = padarray(fwhm, [0 ndim-numel(fwhm)], 'replicate', 'post');
        end
        nb_component = ceil(2 * vs .* lattice ./ fwhm);
        nb_component = max(nb_component, 1);
end

% -------------------------------------------------------------------------
% Mode-specific options
switch lower(mode)
    case 'absolute'
        maxdiff = 0;
    case 'membrane'
        maxdiff = 1;
    case 'bending'
        maxdiff = 2;
    otherwise
        error('Unknown mode %s, should be ''absolute'', ''membrane'' or ''bending''.', mode);
end

% -------------------------------------------------------------------------
% Compute each basis + square it
basis = cell(ndim, maxdiff);
for d=1:ndim
    for diff=0:maxdiff
        switch maxdiff
            case 0
                basis{d,diff+1} = spm_dctmtx(lattice(d),nb_component(d));
            case 1
                basis{d,diff+1} = spm_dctmtx(lattice(d),nb_component(d),'diff');
            case 2
                basis{d,diff+1} = spm_dctmtx(lattice(d),nb_component(d),'diff2');
        end
        basis{d,diff+1} = basis{d,diff+1}.' * basis{d,diff+1};
    end
end

% -------------------------------------------------------------------------
% Compute precision matrix
switch lower(mode)
    case 'absolute'
        L = 1;
        for d=1:ndim
            L = spm_krutil(basis{d,1}, L);
        end
    case 'membrane'
        L = 0;
        for dd=1:ndim               % Which dimension to differentiate
            L1 = 1;
            for d=1:ndim            % Kronecker loop
                if d == dd
                    L1 = spm_krutil(basis{d,2}, L1);
                else
                    L1 = spm_krutil(basis{d,1}, L1);
                end
            end
            L = L + L1;
        end
    case 'bending'
        L = 0;
        for dd1=1:ndim              % First dimension to differentiate
            L1 = 1;
            for d=1:ndim            % Kronecker loop
                if d == dd1
                    L1 = spm_krutil(basis{d,3}, L1);
                else
                    L1 = spm_krutil(basis{d,1}, L1);
                end
            end
            L = L + L1;
            for dd2=dd1+1:ndim      % Second dimension to differentiate
                L1 = 1;
                for d=1:ndim        % Kronecker loop
                    if d == dd1 || d == dd2
                        L1 = spm_krutil(basis{d,2}, L1);
                    else
                        L1 = spm_krutil(basis{d,1}, L1);
                    end
                end
                L = L + 2 * L1;
            end
        end
end

% =========================================================================
function field = reconstruct(basis, coeff, mode)
% FORMAT field = spm_bias_lib('reconstruct', {Bx,By,...}, coefficients)

if nargin < 3
    mode = 'mult';
end

% -------------------------------------------------------------------------
% Get number of components per basis
ndim = numel(basis);
ncomp = zeros(1,ndim);
lat   = zeros(1,ndim);
for i=1:ndim
    lat(i)   = size(basis{i}, 1);
    ncomp(i) = size(basis{i}, 2);
end

% -------------------------------------------------------------------------
% Coefficients provided
if ~isempty(coeff)
    field = coeff;
    for d=1:ndim
        field    = reshape(field, ncomp(1), []); % Coeffs in matrix form
        field    = basis{d} * field;             % Basis x Coeff
        ncomp(1) = size(field, 1);               % Update size (nbcoeffs -> nbvoxels)
        field    = reshape(field, ncomp);        % Coeffs in ND-array form
        ncomp    = circshift(ncomp, -1);         % Shift dimensions
        field    = shiftdim(field, 1);           % Shift dimensions
    end
    switch lower(mode)
        case 'add'
        case 'mult'
            field = single(exp(field));
        otherwise
            error('Unknown bias mode %s. Should be ''add'' or ''mult''.', mode)
    end
    
% -------------------------------------------------------------------------
% No coefficients provided
else
    switch lower(mode)
        case 'add'
            field  = zeros(lat, 'single');
        case 'mult'
            field  = ones(lat, 'single');
        otherwise
            error('Unknown bias mode %s. Should be ''add'' or ''mult''.', mode)
    end
end

% =========================================================================
function [g,H] = derivatives(p, X, B, Z, cluster, codes, binvar)
% FORMAT [g,H] = derivatives(p, obs, basis, resp, cluster, codes, binvar)
%
% p        -     Channel to process. If empty: compute joint gradient/Hessian
% obs      - NxP Bias corrected image
% basis    -     Smooth basis functions {Bx,By,...}
% resp     - NxK Cluster responsibilites
% cluster  -     Either {MU,A} or {MU,V,n} -> Gaussian mixture parameters
% codes    - Nx1 List of cdes encoding missing configuations: C or {C,L}
% binvar   - NxP Bias-modulated binning uncertainty
%
% g - (J1xJ2x...)xP                Gradient w.r.t. bias coefficients
% H - (J1xJ2x...)xPx(J1xJ2x...)xP  Hessian  w.r.t. bias coefficients
%
% Compute gradient and Hessian of the conditional term


MU = [];
A  = [];
V  = [];
n  = [];
C  = [];
L  = [];
if nargin < 7
    binvar = 0;
end


%--------------------------------------------------------------------------
% Read input arguments
if ~iscell(cluster)
    MU = cluster;
else
    if numel(cluster) >= 1
        MU = cluster{1};
        if numel(cluster) >= 2
            A = cluster{2};
            if numel(cluster) >= 3
                V = A;
                n = cluster{3};
            end
        end
    end
end
if nargin >= 6
    if ~iscell(codes)
        C = codes;
    else
        if numel(codes) >= 1
            C = codes{1};
            if numel(codes) >= 2
                L = codes{2};
            end
        end
    end
    if isempty(L)
        L = unique(C);
    end
end

% -------------------------------------------------------------------------
% Dimensions
P  = size(MU,1);
K  = size(MU,2);
if isempty(L), L = 2^P - 1;    end % None missing
ndim  = numel(B);
ncomp = zeros(1,ndim);
lat   = zeros(1,ndim);
for i=1:ndim
    lat(i)   = size(B{i}, 1);
    ncomp(i) = size(B{i}, 2);
end

% -------------------------------------------------------------------------
% Initialise arrays to store statistics for gradient and Hessian
if isempty(p)
    g = zeros(N,1);    % <- 0.5 * Spp * x_p^2 - x_p * [Sp*(mu-0.5*x)] 
    H = zeros(N,1);    % <- 1.5 * Spp * x_p^2 - x_p * [Sp*(mu-0.5*x)]
else
    g = zeros(N,P);    % [p]   <- 0.5 * Spp * x_p^2 - x_p * [Sp*(mu-0.5*x)]
    H = zeros(N,P,P);  % [p,p] <- 1.5 * Spp * x_p^2 - x_p * [Sp*(mu-0.5*x)]
                       % [p,q] <- 0.5 * Spq * x_p * x_q
end

% -------------------------------------------------------------------------
% For each combination of missing voxels
for i=1:numel(L)
    
    % ---------------------------------------------------------------------
    % Get mask of missing modalities (with this particular code)
    c        = L(i);
    observed = code2bin(c, P);
    missing  = ~observed;
    if ~isempty(p) && missing(p), continue; end
    if isempty(C), msk = ones(N, 1, 'logical');
    else,          msk = (C == c);
    end
    Pm      = sum(missing);
    Po      = P-Pm;
    Nc      = sum(msk);
    if Nc == 0, continue; end
    
    % ---------------------------------------------------------------------
    % Convert channel indices to observed indices
    if isempty(p)
        list_p = 1:Po;
    else
        list_p = 1:P;
        list_p = list_p(observed);
        [~,list_p] = find(list_p == p);
    end
    
    X1 = X(msk,observed);
    
    for k=1:K
        
        Z1 = Z(msk,k);
        
        % -----------------------------------------------------------------
        % Compute expected precision (see GMM+missing data)
        if sum(n) > 0
            Ao = V(observed,observed,k) - V(observed,missing,k)*(V(missing,missing,k)\V(missing,observed,k));
            Ao = (n(k)-Pm) * Ao;
        else
            Ao = A(observed,observed,k) - A(observed,missing,k)*(A(missing,missing,k)\A(missing,observed,k));
        end
        MUo = MU(observed,k);
        
        % -----------------------------------------------------------------
        % Compute statistics
        sk1 = zeros(Nc,numel(list_p));
        sk2 = zeros(Nc,numel(list_p),numel(list_p));
        for q=list_p
            tmp1 = Ao(q,q) * X1(:,q).^2;
            tmp2 = X1(:,q) .* (bsxfun(@minus, MUo.', 0.5*X1) * Ao(q,:).');
            sk1(:,q)   = (0.5 * tmp1 - tmp2);
            sk2(:,q,q) = (1.5 * tmp1 - tmp2);
            if numel(binvar) > 1
                sk2(:,q,q) = sk2(:,q,q) + 2 * binvar(msk,q);
            end
            clear tmp1 tmp2
            if numel(list_p) > 1
                for qq=q+1:Po
                    sk2(:,q,qq) = 0.5 * X1(:,q) .* X1(:,qq);
                end
            end
        end
        sk1 = bsxfun(@times, sk1, Z1);
        sk2 = bsxfun(@times, sk2, Z1);
        
        % -----------------------------------------------------------------
        % Accumulate
        if isempty(p)
            g(msk,observed)          = g(msk,observed)          + sk1;
            H(msk,observed,observed) = H(msk,observed,pbserved) + sk2;
        else
            g(msk) = g(msk) + sk1;
            H(msk) = H(msk) + sk2;
        end
        clear sk1 sk2
        
    end
    
    % ---------------------------------------------------------------------
    % Normalisation term
    g(msk,:) = g(msk,:) - 1;
    
end

% -------------------------------------------------------------------------
% Multiply with basis functions
if isempty(p)
    % ---------------------------------------------------------------------
    % Gradient
    dimG = [lat P];
    g    = reshape(g, dimG);
    for d=1:ndim
        g       = reshape(g, dimG(1), []);  % Stats in matrix form
        g       = B{d}.' * g;               % Basis x Stats
        dimG(1) = size(g, 1);               % Update size (nbvoxels -> nbcoeffs)
        g       = reshape(g, dimG);         % Coeffs in ND-array form
        dimG    = circshift(dimG, -1);      % Shift dimensions
        g       = shiftdim(g, 1);           % Shift dimensions
    end
    if P > 1
        g       = shiftdim(g, 1);
    end
    g = reshape(g, [ncomp P]);
    % ---------------------------------------------------------------------
    % Hessian
    BB = ones(1, 'like', B{1});
    for d=1:ndim
        BB = spm_krutil(B{d}, BB);
    end
    BB = reshape(BB, [lat ncomp]);
    H  = bsxfun(@times, BB, reshape(H, [lat 1 P P]));
    BB = reshape(BB, [], prod(ncomp));
    H  = BB' * reshape(H, [], prod(ncomp)*P*P); clear BB
    H  = reshape(H, [prod(ncomp) prod(ncomp) P P]);
    H  = reshape(permute(H, [1 3 2 4]), [ncomp P ncomp P]);
else
    % ---------------------------------------------------------------------
    % Gradient
    dimG = lat;
    g    = reshape(g, dimG);
    for d=1:ndim
        g       = reshape(g, dimG(1), []);  % Stats in matrix form
        g       = B{d}.' * g;               % Basis x Stats
        dimG(1) = size(g, 1);               % Update size (nbvoxels -> nbcoeffs)
        g       = reshape(g, dimG);         % Coeffs in ND-array form
        dimG    = circshift(dimG, -1);      % Shift dimensions
        g       = shiftdim(g, 1);           % Shift dimensions
    end
    g = reshape(g, ncomp);
    % ---------------------------------------------------------------------
    % Hessian
    BB = ones(1, 'like', B{1});
    for d=1:ndim
        BB = spm_krutil(B{d}, BB);
    end
    BB = reshape(BB, [lat ncomp]);
    H  = bsxfun(@times, BB, reshape(H, lat));
    BB = reshape(BB, [], prod(ncomp));
    H  = BB' * reshape(H, [], prod(ncomp)); clear BB
    H  = reshape(H, [ncomp ncomp]);
end


% =========================================================================
function lb = objective(X, Z, B, mean, prec, codes, binvar)
% FORMAT lb = spm_bias_lib('objective', obs, resp, bias, mean, prec, codes, binvar)
%
% MANDATORY
% ---------
% obs    - NxP      observations (non-corrected) 
% resp   - NxK      responsibilities
% bias   - NxP      bias field (non-exponentiated)
% mean   - PxK      GMM mean:      MU or {MU,b} 
% prec   - PxPxK    GMM precision: A  or {V,n}
%
% OPTIONAL
% --------
% codes  - Nx1   image of missing codes (and code list): C or {C,L}
% binvar - 1xP   binning uncertainties
%
% Compute the conditional data term of the objective function:
% sum_n { sum_k log zk * N( Bx | MUk, Ak ) + log |B] }

if nargin < 7, binvar = 0; end

C  = [];
L  = [];

%--------------------------------------------------------------------------
% Read input arguments
if nargin >= 6
    if ~iscell(codes)
        C = codes;
    else
        if numel(codes) >= 1
            C = codes{1};
            if numel(codes) >= 2
                L = codes{2};
            end
        end
    end 
    if isempty(L)
        L = unique(C);
    end
end

%--------------------------------------------------------------------------
% Normalisation term
lb = sum(B(:));
B  = exp(B);

%--------------------------------------------------------------------------
% Compute GMM likelihood from bias-corrected data
X = X .* B;
[lSS0,lSS1,lSS2] = spm_gmm_lib('SuffStat', 'base', X, Z, codes);
if sum(binvar) > 0
    SS2b = spm_gmm_lib('SuffStat', 'bin', binvar, Z, 1, B, codes);
end
lb = lb + spm_gmm_lib('MarginalSum', lSS0, lSS1, lSS2, mean, prec, L, SS2b);
