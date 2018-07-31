function [Z,cluster,prop,lb] = spm_gmm_loop(obs, cluster, prop, varargin)
%__________________________________________________________________________
%
% Fit a [Bayesian] Gaussian mixture model to observed [weighted] data.
%
% This function is the core of the fitting process. However, it needs all
% inputs to be well formatted and initialised and is, thus, not the usual
% entry point. To fit a GMM without having to bother with these issues, use
% spm_gmm instead.
%
% FORMAT [resp,cluster,prop,lb] = spm_gmm_loop(obs,cluster,prop,...)
%
% MANDATORY
% ---------
%
% obs <- X, {X}, or {X,W}
%   X - NxP observations
%   W - Nx1 weights [1]
%
% cluster <- {MU,A}, {{MU,b},A}, {MU,{V,n}}, or {{MU,b},{V,n}}
%   MU - PxK   means
%   b  - 1xK   mean d.f. [0=ML]
%   A  - PxPxK precision matrices
%   V  - PxPxK scale matrices
%   n  - 1xK   precision d.f. [0=ML]
%
% prop <- LogPi or {('LogProp', LogPi), ('Prop', Pi), ('Dir', a)}
%   LogPi - NxK Fixed voxel-wise log-proportions
%           1xK Pre-computed log(Pi) or E[log(Pi)]
%   Pi    - NxK Fixed voxel-wise proportions
%           1xK Pre-computed Pi or E[Pi]
%   a     - 1xK Posterior concentration parameter (Dirichlet)
%
% KEYWORD
% -------
%
% LowerBound     - Pre-computed lower bound structure with fields:
%                   sum, last, X, Z, P, MU, A
% Resp           - NxK Pre-computed responsibilities
% GaussPrior     - {MU0,b0,V0,n0} [{}=ML]
% PropPrior      - a0 [0=ML]
% Missing        - Infer missing data [true]
% MissingCode    - C, {C}, or {C,L} [recompute]
%   C - Nx1 Image of missing code
%   L - List of unique codes
% IterMax        - Max number of EM iterations [1024]
% Tolerance      - Gain tolerance to stop the EM algorithm [1e-4]
% SubIterMax     - Max number of sub-EM iterations (Missing == true) [1024]
% SubTolerance   - Sub-EM gain tolerance (Missing == true) [1e-4]
% BinUncertainty - 1xP Binning uncertainty
%                  NxP Bias-modulated binning uncertainty
% Verbose        - Verbosity level: [0]=quiet, 1=write, 2=plot, 3=plot more,
%                                   4=plot more more
% dm             - Original image dimensions (2d or 3d), necessary when
%                  Verbose=4 [[]]
% Template       - [NxK] Voxel-vise probabalistic template [[]]
% 
% OUTPUT
% ------
% 
% Z       - NxK responsibilities
% cluster - Structure with fields: MU, b, A, V, n
% prop    - Structure with fields: LogProp, Prop, Dir
% lb      - Structure with fields: sum, last, X, Z, P, MU, A
% 
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

lb0 = struct('sum', [], 'last', NaN, ...
             'X', [], 'Z', [], 'P', [], 'MU', [], 'A', []);

% -------------------------------------------------------------------------
% Parse inputs
p = inputParser;
p.FunctionName = 'spm_gmm_loop';
p.addParameter('LowerBound',     lb0,   @isstruct);
p.addParameter('Resp',           [],    @isnumeric);
p.addParameter('GaussPrior',     {},    @iscell);
p.addParameter('PropPrior',      0,     @isnumeric);
p.addParameter('Missing',        true,  @islogical);
p.addParameter('MissingCode',    {},    @(X) isnumeric(X) || iscell(X));
p.addParameter('IterMax',        1024,  @(X) isscalar(X) && isnumeric(X));
p.addParameter('Tolerance',      1e-4,  @(X) isscalar(X) && isnumeric(X));
p.addParameter('SubIterMax',     1024,  @(X) isscalar(X) && isnumeric(X));
p.addParameter('SubTolerance',   1e-4,  @(X) isscalar(X) && isnumeric(X));
p.addParameter('BinUncertainty', 0,     @isnumeric);
p.addParameter('Verbose',        0,     @(X) isscalar(X) && (isnumeric(X) || islogical(X)));
p.addParameter('dm',             [],    @isnumeric);
p.addParameter('Template',       [],    @isnumeric);
p.parse(varargin{:});
lb = p.Results.LowerBound;
Z  = p.Results.Resp;
E  = p.Results.BinUncertainty;
a0 = p.Results.PropPrior;
C  = p.Results.MissingCode;
GaussPrior   = p.Results.GaussPrior;
Missing      = p.Results.Missing;
IterMax      = p.Results.IterMax;
Tolerance    = p.Results.Tolerance;
SubIterMax   = p.Results.SubIterMax;
SubTolerance = p.Results.SubTolerance;
Verbose      = p.Results.Verbose;
dm           = p.Results.dm;
Template     = p.Results.Template;

% -------------------------------------------------------------------------
% Unfold inputs
W   = 1;  % Observations weight
b   = 0;  % Mean degrees of freedom (posterior)
n   = 0;  % Precision degrees of freedom (posterior)
V   = []; % Scale matrix (posterior)
b0  = 0;  % Mean degrees of freedom (prior)
n0  = 0;  % Precision degrees of freedom (prior)
MU0 = []; % Mean (prior)
V0  = []; % Scale matrix (prior)
C   = []; % Image of missing codes
L   = []; % List of unique codes
logPI = []; % [expected] log-proportions
PI    = []; % [expected] proportions
a     = []; % Dirichlet posterior

if ~iscell(obs)
    X = obs;
else
    X = obs{1};
    if numel(obs) >= 2
        W = obs{2};
    end
end
if iscell(C)
    codes = C;
    C = codes{1};
    if numel(codes) >= 2
        L = codes{2};
    else
        L = unique(C);
    end
    clear codes
end
if Missing && isempty(C) && any(any(isnan(X)))
    C = spm_gmm_lib('obs2code', X);
    L = unique(C);
end
if ~iscell(cluster) || numel(cluster) < 2
    error('[spm_gmm_loop] At least one mean and one precision matrix are needed.');
else
    if ~iscell(cluster{1})
        MU = cluster{1};
    else
        MU = cluster{1}{1};
        if numel(cluster{1}) >= 2
            b = cluster{1}{2};
        end
    end
    if ~iscell(cluster{2})
        A = cluster{2};
    else
        A = cluster{2}{1};
        if numel(cluster{2}) >= 2
            n = cluster{2}{2};
            if sum(n) > 0
                V = A;
                A = bsxfun(@times, V, reshape(n, 1, 1, []));
            end
        end
    end
end

% Prepare proportions
if ~iscell(prop)
    logPI = prop;
else
    i = 1;
    while i < numel(prop)
        switch lower(prop{i})
            case 'logprop'
                i = i + 1;
                logPI = prop{i};
            case 'prop'
                i = i + 1;
                PI = prop{i};
            case 'dir'
                i = i + 1;
                a = prop{i};
            otherwise
                logPI = prop{i};
        end
        i = i + 1;
    end
end
if numel(GaussPrior) >= 1
    MU0 = GaussPrior{1};
    if numel(GaussPrior) >= 2
        b0 = GaussPrior{2};
        if numel(GaussPrior) >= 3
            V0 = GaussPrior{3};
            if numel(GaussPrior) >= 4
                n0 = GaussPrior{4};
            end
        end
    end
end

mean = {MU,b};
prec = {V,n};

if ~isempty(Template)
    % Compute logPI by combining Template and [1xK] proportions in PI
    logPI = bsxfun(@times,Template,PI);    
    logPI = log(bsxfun(@times,logPI,1./sum(logPI,2)));
end

% -------------------------------------------------------------------------
% Compute log-prop if needed
if isempty(logPI)
    if sum(a) > 0
        if isempty(PI)
            PI = a ./ sum(a);
        end
        logPI = psi(a) - psi(sum(a));
    elseif ~isempty(PI)
        logPI = log(bsxfun(@rdivide, PI+eps, sum(PI+eps, 2)));
    else
        error('At least one of Prop, LogProp or Dir must be provided.');
    end
end

% -------------------------------------------------------------------------
% Choose how to start
if isempty(Z)
    if Missing, const = spm_gmm_lib('Const', mean, prec, L);
    else,       const = spm_gmm_lib('Const', mean, prec);
    end
    logpX = spm_gmm_lib('Marginal', X, [{MU} prec], const, {C,L}, E);
    Z     = spm_gmm_lib('Responsibility', logpX, logPI);
    clear logpX
end

% -------------------------------------------------------------------------
% EM loop
K = size(logPI,2); % Number of classes
for em=1:IterMax
    
    % ---------------------------------------------------------------------
    % Compute sufficient statistics (bin uncertainty part)
    if sum(E) > 0
    	SS2b = spm_gmm_lib('SuffStat', 'bin', E, Z, W, 1, {C,L});
    else
        SS2b = 0;
    end
    
    if Missing
    % ---------------------------------------------------------------------
    % sub-EM algorithm to update Mean/Precision with missing data
    % . Responsibilities (E[z]) are kept fixed
    % . Missing values (E[z*h], E[z*hh']) are updated
    % . Cluster parameters (MU,b,A/V,n) are updated
    
        % -----------------------------------------------------------------
        % Compute fast sufficient statistics:
        % > sum{E[z]}, sum{E[z]*g}, sum{E[z]*gg'}
        %   for each configuration of missing data
        [lSS0,lSS1,lSS2] = spm_gmm_lib('SuffStat', 'base', X, Z, W, {C,L});
        
        LB = NaN(1,SubIterMax);
        for i=1:SubIterMax
            % -------------------------------------------------------------
            % Infer missing suffstat
            % sum{E[z]}, sum{E[z*x]}, sum{E[z*xx']}
            [SS0,SS1,SS2] = spm_gmm_lib('SuffStat', 'infer', lSS0, lSS1, lSS2, {MU,A}, L);
            SS2 = SS2 + SS2b;

            % -------------------------------------------------------------
            % Update GMM
            [MU,A1,b,V1,n] = spm_gmm_lib('UpdateClusters', ...
                                       SS0, SS1, SS2, {MU0,b0,V0,n0});
            for k=1:size(MU,2)
                [~,cholp] = chol(A1(:,:,k));
                if cholp == 0
                    A(:,:,k) = A1(:,:,k);
                    if sum(n) > 0
                        V(:,:,k) = V1(:,:,k);
                    end
                end
            end
            mean = {MU,b};
            if ~sum(n), prec = {A};
            else,       prec = {V,n};   end
            
            % -------------------------------------------------------------
            % Marginal / Objective function
            [L1,L2]    = spm_gmm_lib('KL', 'GaussWishart', {MU,b}, prec, {MU0,b0}, {V0,n0});
            [L3,const] = spm_gmm_lib('MarginalSum', lSS0, lSS1, lSS2, mean, prec, L, SS2b);
            LB(i+1)    = L1+L2+L3;
            subgain    = abs(LB(i+1)-LB(i))/(max(LB(2:i+1))-min(LB(2:i+1)));
            if Verbose > 0
                fprintf('sub | %4d | lb = %6.10f | gain = %6.10f\n', i, LB(i+1), subgain);
            end
            if subgain < SubTolerance
                break
            end
        end
        
    else
    % ---------------------------------------------------------------------
    % Classical M-step
    
        % -----------------------------------------------------------------
        % Compute sufficient statistics
        [SS0,SS1,SS2] = spm_gmm_lib('SuffStat', X, Z, W);
        SS2 = SS2 + SS2b;
        
        % -------------------------------------------------------------
        % Update GMM
        [MU,A1,b,V1,n] = spm_gmm_lib('UpdateClusters', ...
                                   SS0, SS1, SS2, {MU0,b0,V0,n0});
        for k=1:size(MU,2)
            [~,cholp] = chol(A1(:,:,k));
            if cholp == 0
                A(:,:,k) = A1(:,:,k);
                if sum(n) > 0
                    V(:,:,k) = V1(:,:,k);
                end
            end
        end
        mean = {MU,b};
        if ~sum(n), prec =  {A};
        else,       prec = {V,n};   end
        const = spm_gmm_lib('const', mean, prec, L);
        
    end
                 
    % ---------------------------------------------------------------------
    % Update Proportions
    if ~isempty(Template)
        % Update proportions when a template is given
        logPI = bsxfun(@times,Template,PI);    
        logPI = log(bsxfun(@times,logPI,1./sum(logPI,2)));    
        
        mgm = 1./(Template*PI');
        mgm = mgm'*Template;                
        
        PI = zeros(1,K);
        for k=1:K                     
            PI(k) = (SS0(k) + 1)/(mgm(k) + K);
        end
        PI = PI/sum(PI);
    elseif size(PI,1) == 1
        [PI,logPI,a] = spm_gmm_lib('UpdateProportions', SS0, a0);
    end        

    % ---------------------------------------------------------------------
    % Plot GMM
    if p.Results.Verbose >= 3
        spm_gmm_lib('Plot', 'GMM', {X,W}, {MU,A}, PI);
    end
   
    % ---------------------------------------------------------------------
    % Marginal / Objective function
    logpX = spm_gmm_lib('Marginal', X, [{MU} prec], const, {C,L}, E);
    
    % ---------------------------------------------------------------------
    % Compute lower bound
    lb.P(end+1) = spm_gmm_lib('KL', 'Dirichlet', a, a0);
    lb.Z(end+1) = spm_gmm_lib('KL', 'Categorical', Z, W, logPI);
    if ~Missing
        [lb.MU(end+1),lb.A(end+1)] = spm_gmm_lib('KL', 'GaussWishart', ...
            {MU,b}, prec, {MU0,b0}, {V0,n0});
        lb.X(end+1) = sum(sum(bsxfun(@times, logpX, bsxfun(@times, Z, W)),'omitnan'),'omitnan');
    else
        lb.MU(end+1) = L1;
        lb.A(end+1)  = L2;
        lb.X(end+1)  = L3;
    end

    % ---------------------------------------------------------------------
    % Check convergence
    [lb,gain] = check_convergence(lb, em, Verbose);
    if gain < Tolerance
        break;
    end
    
    % ---------------------------------------------------------------------
    % Compute responsibilities
    Z = spm_gmm_lib('Responsibility', logpX, logPI);
    clear logpX
end

% ---------------------------------------------------------------------
% Plot ML of responsibilities and template
if p.Results.Verbose >= 4 && ~isempty(dm)
    spm_gmm_lib('Plot', 'ml', dm, Z, Template);
end

% -------------------------------------------------------------------------
% Format output
cluster = struct('MU', MU, 'b', b, 'A', A, 'V', V, 'n', n);
prop    = struct('LogProp', logPI, 'Prop', PI, 'Dir', a);

% =========================================================================
function [lb,gain] = check_convergence(lb, em, verbose)
% FORMAT [lb,gain] = check_convergence(lb, em, verbose)
% lb      - Lower bound structure with fields X, Z, P, MU, A, sum, last
% em      - EM iteration
% verbose - Verbosity level (>= 0)
%
% Compute lower bound (by summing its parts) and its gain
% + print info

lb.sum(end+1) = lb.X(end) + lb.Z(end) + lb.P(end) + lb.MU(end) + lb.A(end);
if verbose >= 2
    spm_gmm_lib('plot', 'LB', lb)
end
gain = abs((lb.last - lb.sum(end))/(max(lb.sum(:))-min(lb.sum(:))));
if verbose >= 1
    if     numel(lb.sum) < 2,           incr = '';
    elseif lb.sum(end) > lb.sum(end-1), incr = '(+)';
    elseif lb.sum(end) < lb.sum(end-1), incr = '(-)';
    else,                               incr = '(=)';
    end
    fprintf('%3d | lb = %10.6g |�gain = %10.4g | %3s\n', em, lb.sum(end), gain, incr);
end
lb.last = lb.sum(end);