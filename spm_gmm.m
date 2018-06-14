function [Z,MU,A,PI,b,V,n,a,X] = spm_gmm(X, K, varargin)
% _________________________________________________________________________
%
% Fit a [Bayesian] Gaussian mixture model to observed [weighted] data.
%
% FORMAT [Z,MU,A,PI,(b,V,n),(a),(X)] = spm_gmm(X,K,...)
% 
% MANDATORY
% ---------
% X - NxP matrix of observed values
% K - Number of cluster
% 
% OPTIONAL
% --------
% W  - Nx1 Vector of weights associated with each observation [1]
%
% KEYWORD
% -------
% PropPrior  - 1x[K] vector of Dirichlet priors [0=ML]
%                or NxK matrix of fixed observation-wise proportions.
% GaussPrior - {MU(PxK),b(1x[K]),V(PxPx[K]),n(1x[K])} Gauss-Wishart prior
%              [{}=non-Bayesian GMM]
% Prune      - Threshold on proportions to prune uninformative clusters
%              [0=no pruning]
% Missing    - Infer missing data [true]
% Start      - Starting method: METHOD or {METHOD, PRECISION} with
%    METHOD    = ['linspace'],'kmeans','prior','sample','uniform',MU(PxK)
%    PRECISION = A([PxP]x[K]) [default: diag(a) with a = (range/(2K))^(-2)]
% KMeans     - Cell of KMeans options [{}].
% IterMax    - Maximum number of EM iterations [1000]
% Tolerance  - Convergence criterion (~ lower bound gain) [1e-5]
% Verbose    - Verbosity level: [0]=quiet, 1=write, 2=plot
%
% OUTPUT
% ------
% Z    - NxK   cluster reponsibility
% MU   - PxK   means                                    [E[MU] if Bayesian]
% A    - PxPxK precision matrices                        [E[A] if Bayesian]
% PI   - 1xK   proportions                              [E[PI] if Bayesian]
% b    - 1xK   posterior mean degrees of freedom              [if Bayesian]
% V    - PxPxK posterior scale matrix                         [if Bayesian]
% n    - 1xK   posterior precision degrees of freedom         [if Bayesian]
% a    - 1xK   posterior Dirichlet                            [if Bayesian]
% X    - NxK   obs and inferred values                         [if Missing]
%
% (Note: MAP labels can be obtained with [~,L] = max(Z,[],2)
% _________________________________________________________________________
%
% Use a learned mixture to segment an image.
%
% FORMAT Z = spm_gmm(X,{MU,A},{PI},...)     > Classical
% FORMAT Z = spm_gmm(X,{MU,b,V,n},{a},...)  > Bayesian
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

% TODO
% - include weights when computing lower bound
% - code the "apply" mode
% - decide what default GaussPrior should be:
%   > Bayesian with small DF is more robust than non-Bayesian
% - Better way to include Prop prior for entirely missing voxels?
%   (when all(isnan(X),2))

% -------------------------------------------------------------------------
% Parse inputs
p = inputParser;
p.FunctionName = 'spm_gmm';
p.addRequired('X', @isnumeric);
p.addRequired('K', @isnumeric);
p.addOptional('W',           1,             @isnumeric);
p.addParameter('PropPrior',  0,             @isnumeric);
p.addParameter('GaussPrior', {},            @iscell);
p.addParameter('Prune',      0,             @isscalar);
p.addParameter('Missing',    true,          @isscalar);
p.addParameter('Start',      'kmeans',      @(X) ischar(X) || isnumeric(X));
p.addParameter('KMeans',     {},            @iscell);
p.addParameter('IterMax',    1000,          @isnumeric);
p.addParameter('Tolerance',  1e-5,          @isnumeric);
p.addParameter('Verbose',    0,             @isscalar);
p.parse(X, K, varargin{:});
W         = p.Results.W;
Start     = p.Results.Start;
KMeans    = p.Results.KMeans;
PropPrior = p.Results.PropPrior;

KMeans = [KMeans {'Missing', p.Results.Missing}];
if ~iscell(Start)
    Start = {Start};
end

% -------------------------------------------------------------------------
% Special case: Apply model
% > Here, we use a learned GMM to segment an image
if numel(K) > 1
    Z = gmm_apply(X, K, p.Results.Missing);
    return
end

% -------------------------------------------------------------------------
% Guess clusters from provided initial values
if ~ischar(Start{1})
    K = size(Start{1}, 2);
end

% -------------------------------------------------------------------------
% Proportion / Dirichlet
if size(PropPrior, 1) > 1
    PI              = PropPrior;
    logPI           = log(max(PI,eps));
    voxelwise_prior = true;
    a0              = zeros(1,K);
else
    a0              = PropPrior(:)';
    if numel(a0) < K
        a0          = padarray(a0, [0 K - numel(a0)], 'replicate', 'post');
    end
    voxelwise_prior = false;
    PI              = [];
end

% -------------------------------------------------------------------------
% Vector case
row_vector = size(X,1) == 1 && numel(size(X)) == 2;
if row_vector
    X = X';
end
dim = size(X);
X   = reshape(X, [], dim(end));
X   = double(X);
P   = size(X,2);

% -------------------------------------------------------------------------
% Prepare weights and "missing code"
if numel(W) == 1
    W = W * ones(size(X,1),1);
end
W = double(W(:));
if p.Results.Missing % Deal with missing data
    code      = obs2code(X);
    code_list = unique(code);
else % Discard rows with missing values
    N0        = size(X,1);
    missing   = any(isnan(X),2);
    W         = W(~missing);
    X         = X(~missing,:);
    code      = double2int((2^P-1) * ones(size(X,1),1));
    code_list = 2^P-1;
    if voxelwise_prior
        PI0   = PI(missing,:);
        PI    = PI(~missing,:);
    end
end

% -------------------------------------------------------------------------
% Initialise Gauss-Wishart prior
pr = initialise_prior(p.Results.GaussPrior, K, P);
    
% -------------------------------------------------------------------------
% Initialise mixture
[~, MU, A, PIstart,logPIstart] = start(Start, X, W, K, a0, pr, KMeans);
if isempty(PI)
    PI    = PIstart;
    logPI = logPIstart;
end
clear PIstart

% -------------------------------------------------------------------------
% Default prior mean/precision if needed
if isempty(pr.MU) && ~isempty(pr.b)
    pr.MU = MU;
end
if isempty(pr.V) && ~isempty(pr.n)
    pr.V = bsxfun(@rdivide, A, reshape(pr.n, [1 1 K]));
end

% -------------------------------------------------------------------------
% Initialise posterior
if isempty(pr.b) && isempty(pr.n)
    b      = [];
    n      = [];
    V      = [];
    const  = constant_term(MU,[],A,[]);
    const0 = constant_term(MU,[],A,[],code_list);
elseif isempty(pr.b)
    b      = [];
    n      = pr.n;
    V      = bsxfun(@rdivide,A,reshape(n,[1 1 K]));
    const  = constant_term(MU,[],V,n);
    const0 = constant_term(MU,[],V,n,code_list);
elseif isempty(pr.n)
    b      = pr.b;
    n      = [];
    V      = [];
    const  = constant_term(MU,b,A,[]);
    const0 = constant_term(MU,b,A,[],code_list);
else
    b      = pr.b;
    n      = pr.n;
    V      = bsxfun(@rdivide,A,reshape(n,[1 1 K]));
    const  = constant_term(MU,b,V,n);
    const0 = constant_term(MU,b,V,n,code_list);
end
          
% EM loop
lb  = struct('sum', [], 'X', [], 'Z', [], 'P', [], 'MU', [], 'A', []);
ll0 = nan;
for em=1:p.Results.IterMax
    
    if em == 1
        % Compute responsibilities from non-missing values only
        Z = responsibility(X, MU, A, logPI, code, code_list, const0);
        Z(all(isnan(X),2),:) = repmat([0.9 repmat(0.1/(K-1), [1 K-1])], [sum(all(isnan(X),2)) 1]);
    else
        % Infer missing data
        X = infer_missing(X, Z, MU, A, code, code_list);

        % Compute responsibilities
        [Z,lb.X(end+1),lb.Z(end+1)] = responsibility(X, MU, A, logPI, code, code_list, const);
    end
    
    % Compute sufficient statistics
    [SS0,SS1,SS2] = suffstat(X, W, Z, A, code, code_list);
    
    % Update Proportions
    if em > 1
        lb.P(end+1) = kl_dirichlet(a, a0);
    end
    if ~voxelwise_prior
        [PI,logPI,a] = update_proportions(SS0, a0);
    end
    
    % Update GMM
    if em > 1
        if isempty(V)
            [lb.MU(end+1),lb.A(end+1)] = kl_gauss_wishart(MU,b,A,n,pr.MU,pr.b,pr.V,pr.n);
        else
            [lb.MU(end+1),lb.A(end+1)] = kl_gauss_wishart(MU,b,V,n,pr.MU,pr.b,pr.V,pr.n);
        end
    end
    [MU,b,V,n,A] = update_gmm(SS0,SS1,SS2,pr);
    if isempty(V),  const = constant_term(MU,b,A,[]);
    else,           const = constant_term(MU,b,V,n);   end
    
    % Check convergence
    if em > 1
        lb.sum(end+1) = lb.X(end) + lb.Z(end) + lb.P(end) + lb.MU(end) + lb.A(end);
        if p.Results.Verbose >= 2
            subplot(2, 3, 1);
            plot(lb.sum)
            title('Lower Bound')
            subplot(2, 3, 2);
            plot(lb.X)
            title('Observations (E[ln p(X)] - E[ln q(H)])')
            subplot(2, 3, 3);
            plot(lb.Z)
            title('Responsibilities (KL)')
            subplot(2, 3, 4);
            plot(lb.P)
            title('Proportions (KL)')
            subplot(2, 3, 5);
            plot(lb.MU)
            title('Means (KL)')
            subplot(2, 3, 6);
            plot(lb.A)
            title('Precisions (KL)')
            drawnow
        end
        gain = abs((ll0 - lb.sum(end))/(max(lb.sum(:))-min(lb.sum(:))));
        if p.Results.Verbose
            fprintf('%3d | lb = %6g | gain = %6g\n', em, lb.sum(end), gain);
        end
        if gain < p.Results.Tolerance
            break;
        end
        ll0 = lb.sum(end);
    end
    
end

% -------------------------------------------------------------------------
% Prune clusters
if ~voxelwise_prior && p.Results.Prune > 0
    kept = PI >= p.Results.Prune;
    Z    = Z(:,kept);
    MU   = MU(:,kept);
    A    = A(:,:,kept);
    PI   = PI(kept);
    if ~isempty(pr)
        b = b(kept)';
        n = n(kept)';
        V = V(:,:,kept);
    end
end

% -------------------------------------------------------------------------
% Replace discarded missing values
if ~p.Results.Missing
    Z1            = Z;
    Z             = zeros([N0 size(Z1,2)], 'like', Z1);
    Z(~missing,:) = Z1; clear Z1;
    if voxelwise_prior
        PI1            = PI;
        PI             = zeros(N0, K);
        PI(~missing,:) = PI1; clear PI1
        PI(missing,:)  = PI0;
    end
end

% =========================================================================
function X = infer_missing(X, Z, MU, A, code, code_list)
% FORMAT X = infer_missing(X, Z, MU, A, code, code_list)
% X         - NxP   observations
% Z         - NxK   responsibilities
% MU        - PxK   (expected) means
% A         - PxPxK (expected) precision matrices
% code      - "Missing value" code image
% code_list - List of existing codes
%
% Compute the mean expected value of missing voxels.

for i=1:numel(code_list)
    c       = code_list(i);
    msk     = code == c;
    missing = ~code2bin(c, size(X,2));
    if sum(msk(:)) == 0 || sum(missing) == 0, continue; end
    
    K = size(Z,2);
    X(msk,missing) = 0;
    for k=1:K
        X1k = 0;
        X1k = X1k + MU(missing,k).' * A(missing,missing,k);
        X1k = X1k + (MU(~missing,k).' - X(msk,~missing)) * A(~missing,missing,k);
        X1k = X1k / A(missing,missing,k);
        X(msk,missing) = X(msk,missing) + bsxfun(@times, X1k, Z(msk,k));
    end
end


% =========================================================================
function [logpX, logqX] = gmm_elogp(X, MU, A, code, code_list, const)
% FORMAT [logpX, logqX] = gmm_elogp(X, MU, A, code, code_list, const)
% 
% X         - NxP observed values
% MU        - PxK (expected) means
% A         - PxPxK (expected) precision matrices
% code      - Image of "missing" codes
% code_list - List of existing codes
% const     - [M]xK constant terms. If M > 1, discard missing values
%
% logpX     - NxK (expected) log-likelihood of belonging to each class
% logqX     - 1x1 Posterior part of the lower bound: sum E[log q(h)]

if size(const, 1) == 1
    const = repmat(const, [numel(code_list) 1]);
    discard_missing = false;
else
    discard_missing = true;
end

N  = size(X,1);
K  = size(MU,2);
logpX  = zeros([N K]);
logqX  = 0;

% -------------------------------------------------------------------------
% For each combination of missing voxels
for i=1:numel(code_list)
    c       = code_list(i);
    msk     = code == c;
    missing = ~code2bin(c, size(X,2));
    if sum(msk(:)) == 0, continue; end
    
    logpX(msk,:) = repmat(const(i,:), [sum(msk), 1]);
    X1 = X(msk,:)';
    for k=1:K
        % Compute expected log-likelihood (eq (23) of GMM-MD paper)
        l = zeros([1 sum(msk)]);
        % 1) quadratic in observed values
        part1 = A(~missing,~missing,k) * (X1(~missing,:) - 2*MU(~missing,k));
        part1 = dot(part1, X1(~missing,:), 1);
        l = l - 0.5 * part1; clear part1
        if ~discard_missing && sum(missing) > 0
            % 2) quadratic in missing values
            part2 = A(missing,missing,k) * (X1(missing,:) - 2*MU(missing,k));
            part2 = dot(part2, X1(missing,:), 1);
            part2 = part2 + sum(missing);
            l = l - 0.5 * part2; clear part2
            % 3) quadratic in observed x missing values
            part3 = MU(missing,k)' * A(missing,~missing,k) * X1(~missing,:);
            l = l + part3; clear part3
            part4 = A(missing,~missing,k) * (X1(~missing,:) - MU(~missing,k));
            part4 = dot(part4, X1(missing,:), 1);
            l = l - part4; clear part4
        end
        % Reshape as a column vector
        logpX(msk,k) = logpX(msk,k) + l'; clear l
        
        if  ~discard_missing && sum(missing) > 0
            % E[ln q(H)] (posterior ~ missing values) 
            logqX = logqX + 0.5*sum(msk)*( ...
                - spm_matcomp('LogDet',A(missing,missing,k)) ...
                - numel(missing)*(1+log(2*pi)) );
        end
    end
end

% =========================================================================
function [Z,llX,klZ] = responsibility(X, MU, A, logPI, code, code_list, const)
% FORMAT [Z,ll] = responsibility(X, MU, A, logPI, code, code_list, const)
%
% X          - NxK   observed (+ inferred) values
% MU         - PxK   (Expected) means
% A          - PxPxK (Expected) precision matrices
% logPI      - 1xK   (Expected) log proportions
% code       - Nx1   Code image of missing values
% code_list  -       List of exisinting codes (saves time)
% const      - 1xK   Constant term (across voxels) of E[ln p(Xn|MUk,Ak)]
%
% Z          - NxK   (Posterior) responsibilities
% ll         -       Lower bound =   E[Z]E[ln p(X|MU,A)] - E[ln q(H)]
%                                  + E[ln p(Z|PI)]       - E[ln q(Z)]

klZ = 0;
% log-likelihood part: E[log p(X|MUk,Ak)]

[logpX, llX] = gmm_elogp(X, MU, A, code, code_list, const);

% prior part: E[log PI]
Z = logpX + logPI;

% Exponentiate and normalise
Z = Z - max(Z, [], 2);
Z = exp(Z);
Z = bsxfun(@rdivide, Z, sum(Z, 2));

% E[ln p(Z|PI)] (prior ~ responsibilities)
klZ = klZ + sum(sum(Z .* logPI));

% -E[ln q(Z)] (posterior ~ responsibilities))
klZ = klZ - sum(sum(Z .* log(max(Z,eps))));

% E[Z]E[ln p(X|A,MU)] (prior ~ intensities) 
llX = -llX + sum(sum(logpX .* Z));

% =========================================================================
function [SS0,SS1,SS2] = suffstat(X, W, Z, A, code, code_list)
% FORMAT [SS0,SS1,SS2] = suffstat(X, W, Z, code, code_list);
%
% X    - NxP Observed + Inferred values
% W    - Nx1 Observation weights
% Z    - NxK Responsibilities
% A    - PxPxK precision matrices (useful for uncertainties)
% code - Nx1 Missing values "code"
% code_list - List of codes present (saves a tiny bit of time if provided)
%
% SS0 - 1xK   0th order suff stat (sum of resp)
% SS1 - PxK   1st order suff stat (weighted sum of intensities)
% SS2 - PxPxK 2nd order suff stat (weighted sum of squared intensities)
%
% Compute sufficient statistics up to 2nd order, taking into account
% inferred values and their uncertainty.
%
% Note that, because we use the code image, no uncertainties are added 
% if missing values were discarded ("~Missing" case)

if nargin < 6
    code_list = unique(code);
end

N = size(X,1);
P = size(X,2);
K = size(A,3);
Z   = bsxfun(@times, Z, W); % Multiply resp with observation count
SS0 = sum(Z, 1);
SS1 = sum(bsxfun(@times, X, reshape(Z, [N 1 K])), 1, 'omitnan');
SS1 = reshape(SS1, [P K]);
SS2 = zeros(P,P,K);
% Add quadratic terms ~ observed/inferred
for i=1:P
    SS2(i,i,:) = reshape(sum(bsxfun(@times, Z, X(:,i).^2),1,'omitnan'), [1 1 K]);
    for j=i+1:P
        SS2(i,j,:) = reshape(sum(bsxfun(@times, Z, X(:,i).*X(:,j)),1,'omitnan'), [1 1 K]);
        SS2(j,i,:) = SS2(i,j,:);
    end
end
% Add uncertainty ~ inferred values
for k=1:K
    A(:,:,k) = spm_matcomp('Inv', A(:,:,k));
end
for i=1:numel(code_list)
    c       = code_list(i);
    missing = ~code2bin(c, size(X,2));
    if sum(missing) > 0
        msk = (code == c) & (~any(isnan(X),2));
        if sum(msk) > 0
            Zc  = Z(msk,:);
            sumA = bsxfun(@times, ...
                reshape(Zc, [size(Zc,1) 1 1 K]), ...
                reshape(A(missing,missing,:), [1 sum(missing) sum(missing) K]));
            sumA = reshape(sum(sumA, 1), [sum(missing) sum(missing) K]);
            SS2(missing,missing,:) = SS2(missing,missing,:) + sumA;
        end
    end
end

% =========================================================================
function const = constant_term(MU,b,V,n,code_list)
% FORMAT const = constant_term(MU,b,(V|A),n)
% MU - (Expected) mean
% b  - Mean df (if empty or 0 -> no Bayesian prior)
% V  - Scale matrix     (if not n empty or 0) 
% A  - Precision matrix (if n empty or 0)
% n  - Precision df (if empty or 0 -> no Bayesian prior)
%
% Compute the constant term (w.r.t. voxels) of each Gaussian 
% (expected) log-likelihood.

P = size(MU,1);
K = size(MU,2);
    
if nargin == 5
% Discard missing values
    const = zeros(numel(code_list), K);
    for i=1:numel(code_list)
        c = code_list(i);
        missing = ~code2bin(c, P);
        for k=1:K
            const(i,k) = - 0.5 * numel(~missing) * log(2*pi);
            if ~isempty(n) && n(k) > 0
                const(i,k) = const(i,k) + 0.5 * spm_prob('W','ELogDet',V(~missing,~missing,k),n(k)) ...
                                    - 0.5 * n(k) * MU(~missing,k)' * V(~missing,~missing,k) * MU(~missing,k);
            else
                const(i,k) = const(i,k) + 0.5 * spm_matcomp('LogDet',V(~missing,~missing,k)) ...
                                    - 0.5 * MU(~missing,k)' * V(~missing,~missing,k) * MU(~missing,k);
            end
            if ~isempty(b) && b(k) > 0
                const(i,k) = const(i,k) - 0.5 * numel(~missing) / b(k);
            end
        end
        
    end
    
else
% With missing values
    const = zeros(1,K);
    for k=1:K
        const(k) = - 0.5 * P * log(2*pi);
        if ~isempty(n) && n(k) > 0
            const(k) = const(k) + 0.5 * spm_prob('W','ELogDet',V(:,:,k),n(k)) ...
                                - 0.5 * n(k) * MU(:,k)' * V(:,:,k) * MU(:,k);
        else
            const(k) = const(k) + 0.5 * spm_matcomp('LogDet',V(:,:,k)) ...
                                - 0.5 * MU(:,k)' * V(:,:,k) * MU(:,k);
        end
        if ~isempty(b) && b(k) > 0
            const(k) = const(k) - 0.5 * P / b(k);
        end
    end
end

% =========================================================================
function [MU,b,V,n,A] = update_gmm(SS0,SS1,SS2,pr)
% FORMAT [MU,b,V,n,A] = update_gmm_bayes(SS0,SS1,SS2,pr)
% SS0 - 0th order sufficient statistics (sum Z_i)
% SS1 - 1st order sufficient statistics (sum Z_i * X_i)
% SS2 - 2nd order sufficient statistics (sum Z_i * (X_i * X_i'))
% pr  - Structure of prior Gauss-Wishart parameters.
%
% Compute posterior GMM parameters from suff stats.

K  = numel(SS0);

% -------------------------------------------------------------------------
% Mean
if sum(pr.b) == 0
    % ---------------------------------------------------------------------
    % Without prior
    b = [];
    for k=1:K
        SS2(:,:,k) = SS2(:,:,k) - SS1(:,k) * SS1(:,k).';
    end
    MU  = bsxfun(@rdivide, SS1, SS0);
else
    % ---------------------------------------------------------------------
    % With prior
    b  = pr.b + SS0;
    MU = bsxfun(@rdivide, SS1 + bsxfun(@times,pr.b,pr.MU), b);
    for k=1:K
        SS2(:,:,k) = SS2(:,:,k) + pr.b(k) * pr.MU(:,k) * pr.MU(:,k).' ...
                                -    b(k) *    MU(:,k) *    MU(:,k).';
    end
end

% -------------------------------------------------------------------------
% Scale/Precision
if sum(pr.n) == 0
    % ---------------------------------------------------------------------
    % Without prior
    n   = [];
    SS2 = bsxfun(@rdivide, SS2, reshape(SS0, [1 1 K]));
else
    % ---------------------------------------------------------------------
    % With prior
    n = pr.n + SS0;
    for k=1:K
        SS2(:,:,k) = SS2(:,:,k) + spm_matcomp('Inv',pr.V(:,:,k));
    end
end
V = SS2;
for k=1:K
    V(:,:,k) = spm_matcomp('Inv', V(:,:,k));
end
if sum(n) > 0
    A = bsxfun(@times, V, reshape(n, [1 1 K]));
else
    A = V;
    V = [];
end


% =========================================================================
function [klMU,klA] = kl_gauss_wishart(MU,b,V,n,MU0,b0,V0,n0)

P = size(MU,1);
K = size(MU,2);
LogDetA = zeros(1,K);
if sum(n) > 0
    A = bsxfun(@times, V, reshape(n, [1 1 K]));
    for k=1:K
        LogDetA(k) = spm_prob('W','ELogDet',V(:,:,k),n(k));
    end
else
    A = V;
    for k=1:K
        LogDetA(k) = spm_matcomp('LogDet',A(:,:,k));
    end
end

% Lower bound
klMU = 0;
klA  = 0;
for k=1:K
    % + prior
    if sum(b0) > 0
        % prior
        klMU = klMU - P*log(2*pi) ...
                    + P*log(b0(k)) ...
                    + LogDetA(k) ...
                    - b0(k)*(MU(:,k)-MU0(:,k)).'*A(:,:,k)*(MU(:,k)-MU0(:,k)) ...
                    - P*b0(k)/b(k);
        % posterior
        klMU = klMU + P*log(2*pi) ...
                    - P*log(b(k)) ...
                    - LogDetA(k) ...
                    + P;
    end
    if sum(n0) > 0
        klA = klA - spm_prob('W', 'kl', V(:,:,k), n(k), V0(:,:,k), n0(k));
    end
end
klMU = 0.5 * klMU;
klA  = 0.5 * klA;


% =========================================================================
function [PI,logPI,a,ll] = update_proportions(SS0, a0)
% FORMAT [PI,logPI,a,ll] = update_proportions(SS0, a0)
%
% SS0 - 1xK 0th order sufficient statistics (sum of responsibilities)
% a0  - 1xK Dirichlet prior (can be 0)
%
% PI    - 1xK Cluster proportion posterior expected value
% logPI - 1xK ln(PI) or E[ln(PI)] (if Bayesian)
% a     - Dirichlet posterior (if Bayesian)
% ll    - 1x1 Lower bound: E[ln p(PI|a)] - E[ln q(PI)]
%
% Bayesian or ML update of cluster proportions.

K = numel(a0);
a = a0 + SS0;
if sum(a0(:))
% Bayesian
    % expected values
    logPI = psi(a) - psi(sum(a));
    PI    = a ./ sum(a(:));
else
% Maximum Likelihood
    ll    = 0;
    a     = max(a, eps);
    PI    = a ./ sum(a(:));
    logPI = log(PI);
end

% =========================================================================
function klP = kl_dirichlet(a, a0)

klP = 0;
K   = numel(a);
if sum(a0) > 0
    % prior
    klP = gammaln(sum(a0)) - sum(gammaln(a0));
    klP = klP + sum((a0-1) .* (psi(a) - K*psi(sum(a))));
    % posterior
    klP = klP - gammaln(sum(a)) - sum(gammaln(a));
    klP = klP - sum((a-1) .* (psi(a) - K*psi(sum(a))));
end

% =========================================================================
function [Z,MU,A,PI,logPI] = start(method, X, W, K, a0, pr, kmeans)
% FORMAT C = start(method, X, W, K, a0, pr, kmeans)
%
% method - Method to use to select starting centroids
%               'kmeans', 'sample', 'uniform' or provided matrix
% X      - Vector of NxP observed values
% W      - Vector of Nx1 weights
% K      - Number of clusters
% a0     - 1xK dirichlet priors (or empty)
% pr     - Structure of Gauss-Wishart prior parameters
% kmeans - Options for spm_kmeans
%
% MU - PxK   means
% A  - PxPxK precision matrices
% PI - 1xK   proportions
%
% Compute starting estimates

if ~iscell(method)
    method = {method};
end

switch method{1}
    
    case 'kmeans'
    % Use kmeans to produce a first clustering of the data
        P = size(X,2);
        N = size(X,1);
        % Get labelling and centroids from K-means
        [L,MU] = spm_kmeans(X, K, W, kmeans{:});
        MU = MU';
        % Convert to "responsibility"
        Z = zeros(N,K);
        for k=1:K
            Z(:,k) = (L == k) + eps;
        end
        Z(all(isnan(X),2),:) = 1/K;
        clear L
        Z = bsxfun(@rdivide, Z, sum(Z,2));
        % Compute proportion from labelling
        PI = zeros(1,K);
        for k=1:K
            PI(k) = sum(Z(:,k)) + a0(k);
        end
        if sum(a0) > 0
            logPI = psi(PI) - psi(sum(PI));
            PI    = PI ./ sum(PI);
        else
            PI    = max(PI,eps);
            PI    = PI ./ sum(PI);
            logPI = log(PI);
        end
        if numel(method) == 1
            % Compute precision from intra-class sample variance
            % 1) SuffStat2 = sum_i { Wi * (Xi-MU)*(Xi-MU)' }
            X2 = X - reshape(MU, [1 P K]);
            A  = zeros(N,P,P,K);
            for c1=1:P
                A(:,c1,c1,:) = reshape(X2(:,c1,:).^2, [N 1 1 K]);
                for c2=c1+1:P
                    A(:,c1,c2,:) = reshape(X2(:,c1,:) .* X2(:,c2,:), [N 1 1 K]);
                    A(:,c2,c1,:) = A(:,c1,c2,:);
                end
            end
            clear X2
            A = sum(bsxfun(@times, A, reshape(bsxfun(@times, Z, W), [N 1 1 K])), 1, 'omitnan');
            A = reshape(A, [P P K]);
            % 2) SuffStat0 = sum_i { Wi * Present*Present' }
            X   = ~isnan(X);
            SUM = zeros(N,P,P);
            for c1=1:P
                SUM(:,c1,c1) = X(:,c1,:).^2;
                for c2=c1+1:P
                    SUM(:,c1,c2) = X(:,c1) .* X(:,c2);
                    SUM(:,c2,c1) = SUM(:,c1,c2);
                end
            end
            SUM = sum(bsxfun(@times, SUM, reshape(bsxfun(@times, Z, W), [N 1 1 K])), 1, 'omitnan');
            SUM = reshape(SUM, [P P K]);
            % 3) Normalise and invert
            A = A ./ SUM; clear SUM X
            for k=1:K
                A(:,:,k) = A(:,:,k) + A(:,:,k)'; % Ensure symmetric
                A(:,:,k) = spm_matcomp('Inv', A(:,:,k));
            end
        end
        
    case 'prior'
    % Use prior expected value
        if isempty(pr.MU) || isempty(pr.V) || isempty(pr.n)
            error('To initialise from prior, a full prior must be provided')
        end
        MU = pr.MU;
        A  = bsxfun(@times, pr.V, reshape(pr.n, [1 1 K]));
        
    case 'sample'
    % Sample uniform
    % All centroids are selected at random from the observed values.
    % They are all unique (to avoid starting with several identical
    % centroids)
        X   = X(~(any(W==0,2)|any(isnan(X),2)),:); % Remove all rows w/ NaNs or w/o obs
        i   = randperm(size(X,1));
        i   = i(1:K);
        MU = X(i,:)';
        
    case 'uniform'
    % Range uniform
    % All centroids are selected at random from the continuous range of 
    % observed values.
    % They are all unique (to avoid starting with several identical
    % centroids)
        minval = min(X, [], 1, 'omitnan');
        maxval = max(X, [], 1, 'omitnan');
        MU = rand(K,size(X,2));
        MU = bsxfun(@times, MU, maxval - minval) + minval;
        MU = MU';
        
    case 'linspace'
    % Range uniform
    % All centroids are selected at random from the continuous range of 
    % observed values.
    % They are all unique (to avoid starting with several identical
    % centroids)
        minval = min(X, [], 1, 'omitnan');
        maxval = max(X, [], 1, 'omitnan');
        MU = rand(size(X,2),K);
        for p=1:size(X,2)
            ticks = linspace(minval(p), maxval(p), 2*K+1);
            MU(p,:) = ticks(2:2:end-1);
        end
        
    otherwise
    % Provided
        MU = method{1};
end

% Precision matrix
% - If nothing is provided, we split the range in 2K pieces and set the 
%   standard deviation (sigma) to that value. This way, half of each 
%   Gaussian covers (1/K)th of the range of values.
% - If scalars,are provided, we use them as precisions in diagonal
%   precision matrices.
% - If matrices are provided, we just copy them.
if numel(method) >= 2 || ~any(strcmpi(method{1}, {'kmeans','prior'}))
    if numel(method) > 1
        A = method{2};
    else
        minval = min(X, [], 1, 'omitnan');
        maxval = max(X, [], 1, 'omitnan');
        A = diag(((maxval - minval)./(2*K)).^(-2));
    end
    if size(A,3) < K
        A = padarray(A, [0 0 K-size(A,3)], 'replicate', 'post');
    end
    if size(A,1) == 1
        A = bsxfun(@times, A, eye(size(X,2)));
    end
end
% Proportion
% - If a Dirichlet prior is provided, we use its expected value
% - Else, we use 1/K
if ~any(strcmpi(method{1}, {'kmeans'}))
    if sum(a0) > 0
        PI    = a0;
        logPI = psi(PI) - psi(sum(PI));
        PI    = PI ./ sum(PI);
    else
        PI    = 1/K * ones(1,K);
        PI    = PI ./ sum(PI);
        logPI = log(PI);
    end
end
% Responsibility
% - We use the expected proportion.
if ~any(strcmpi(method{1}, {'kmeans'}))
    Z = repmat(PI(:)', [size(X,1) 1]);
    allnan = all(isnan(X), 2);
    Z(allnan,:) = repmat([0.9 repmat(0.1/(K-1), [1 K-1])], [sum(allnan) 1]);
end

% =========================================================================
function pr = initialise_prior(pr0, K,P)
% FORMAT pr = initialise_prior(pr0)
%
% pr0 - Cell of paramters {MU, b, V, n}
% K   - Number of classes
% P   - Number of channels
% pr  - Structure with fields MU, b, V, n
%
% Each input parameter can be missing or empty, in which case:
%   isempty(MU) ->  isempty(b)  -> no prior over MU (ML optimisation)
%               -> ~isempty(b)  -> 'Start' method
%   isempty(b)  -> ~isempty(MU) -> eps (most uninformative prior)
%   isempty(V)  ->  isempty(n)  -> no prior over V (ML optimisation)
%               -> ~isempty(n)  -> eye(P)
%   isempty(n)  -> ~isempty(V)  -> P (most uninformative prior)

pr = struct;
pr.MU = [];
pr.b  = [];
pr.V  = [];
pr.n  = [];
if isempty(pr0)
    return
end

% -------------------------------------------------------------------------
% Mean
pr.MU = pr0{1};
if ~isempty(pr.MU)
    P = size(pr.MU, 1);
    K = size(pr.MU, 2);
end
% -------------------------------------------------------------------------
% Mean df
if numel(pr0) > 1
    pr.b = pr0{2};
end
% -------------------------------------------------------------------------
% Mean df [default]
if ~isempty(pr.MU) && isempty(pr.b)
    pr.b = eps;
end
if ~isempty(pr.b)
    pr.b = pr.b(:)';
    if numel(pr.b) < K
        pr.b = padarray(pr.b, [0 K - numel(pr.b)], 'replicate', 'post');
    end
end
% -------------------------------------------------------------------------
% Scale
if numel(pr0) > 2
    pr.V = pr0{3};
end
% -------------------------------------------------------------------------
% Scale df
if numel(pr0) > 3
    pr.n = pr0{4};
end
% -------------------------------------------------------------------------
% Scale df [default]
if ~isempty(pr.V) && isempty(pr.n)
    pr.n = P;
end
if ~isempty(pr.n)
    pr.n = pr.n(:)';
    if numel(pr.n) < K
        pr.n = padarray(pr.n, [0 K - numel(pr.n)], 'replicate', 'post');
    end
end
% -------------------------------------------------------------------------
% Scale [default]
if ~isempty(pr.V)
    if size(pr.V,3) < K
        pr.V = padarray(pr.V, [0 0 K - size(pr.V,3) ], 'replicate', 'post');
    end
    if size(pr.V,1) == 1
        pr.V = bsxfun(@times, pr.V, eye(P));
    end
end

% =========================================================================
function code = obs2code(X)
% FORMAT code = obs2code(X)
%
% Compute a "missing code" image for the input observation matrix.

code = double2int(sum(bsxfun(@times, ~isnan(X), 2.^(0:size(X,2)-1)), 2));

% =========================================================================
function bin = code2bin(code, C)
    bin = dec2bin(code,C) == '1';
    bin = bin(end:-1:1);

% =========================================================================
function L = double2int(L)
% FORMAT L = double2int(L)
%
% Find the best suited integer type to convert L, based on min and max
% values

minval = min(L(:));
maxval = max(L(:));
type   = range2int(maxval,minval);
func   = str2func(type);
L      = func(L);

% =========================================================================
function type = range2int(maxval,minval)
% FORMAT type = range2int(maxval,minval)
%
% Find the best suited integer type to store integer values in the range
% [minval,maxval]

if nargin < 2
    minval = 0;
end

type     = 'int';
unsigned = minval >= 0;
if unsigned
    type   = ['u' type];
    minval = 0;
else
    minval = numel(dec2base(-minval,2));
end
maxval = numel(dec2base(maxval,2));
nbits  = max(minval,maxval);
if unsigned
    nbits = nbits + 1;
end
if nbits <= 8
    type = [type '8'];
elseif nbits <= 16
    type = [type '16'];
elseif nbits <= 32
    type = [type '32'];
elseif nbits <= 64
    type = [type '54'];
else
    type = 'double';
end