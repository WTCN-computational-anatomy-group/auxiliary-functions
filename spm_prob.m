function varargout = spm_prob(varargin)
%__________________________________________________________________________
% Collection of tools for probability distributions (PDF, KL-div, ...):
%   [n/normal/gaussian]
%   [g/gamma]
%   [ig/inverse-gamma]
%   [w/wishart]
%   [iw/inverse-wishart]
%   [ng/normal-gamma]
%   [nw/normal-wishart]
%
% FORMAT out = spm_prob('Normal', ...)
% FORMAT out = spm_prob('Gamma', ...)
%
% FORMAT help spm_prob>function
% Returns the help file of the selected function.
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

    if nargin == 0
        help spm_prob
        error('Not enough argument. Type ''help spm_prob'' for help.');
    end
    id = varargin{1};
    varargin = varargin(2:end);
    switch lower(id)
        case {'normal', 'n', 'gaussian'}
            [varargout{1:nargout}] = normal(varargin{:});
        case {'gamma', 'g'}
            [varargout{1:nargout}] = gamma(varargin{:});
        otherwise
            help spm_prob
            error('Unknown function %s. Type ''help spm_prob'' for help.', id)
    end
end

% =========================================================================
%   NORMAL
% =========================================================================

function varargout = normal(varargin)
%__________________________________________________________________________
% Characteristic functions of the (uni/multivariate) Normal distribution:
%   [pdf]       Probability density function
%   [ll/logpdf] Log-probability density function
%   [kl]        Kullback-Leibler divergence
%   [up/update] Conjugate (or Bayesian) update
%
% FORMAT pdf = spm_prob('Normal', 'pdf',    x, mu, sigma)
% FORMAT pdf = spm_prob('Normal', 'pdf',    x, mu, lambda, 'precision')
% FORMAT ll  = spm_prob('Normal', 'logpdf', x, mu, sigma)
% FORMAT ll  = spm_prob('Normal', 'logpdf', x, mu, lambda, 'precision')
%   >> (Log) Probability density function.
%
% FORMAT pdf = spm_prob('Normal', 'pdf',    x, mu, n, sigma)
% FORMAT pdf = spm_prob('Normal', 'pdf',    x, mu, n, lambda, 'precision')
% FORMAT ll  = spm_prob('Normal', 'logpdf', x, mu, n, sigma)
% FORMAT ll  = spm_prob('Normal', 'logpdf', x, mu, n, lambda, 'precision')
%   >> Reparameterisation of the (log)-PDF when the Normal distribution
%      is used as a distribution for a Normal mean with known variance.
%
% FORMAT kl  = spm_prob('Normal', 'kl', mu1, sigma1,  mu0, sigma0)
% FORMAT kl  = spm_prob('Normal', 'kl', mu1, lambda1, mu0, lambda0, 'precision')
%   >> Kullback-Leibler divergence from N0 to N1 = KL(N1||N0)
%
% FORMAT kl  = spm_prob('Normal', 'kl', mu1, n1, mu0, n0, sigma)
% FORMAT kl  = spm_prob('Normal', 'kl', mu1, n1, mu0, n0, lambda, 'precision')
%   >> Reparameterisation of the KL-divergence when the Normal distribution
%      is used as a conjugate prior for a Normal mean with known variance.
%
% FORMAT [mu1, n1] = spm_prob('Normal', 'update', mu, n, mu0, n0)
%   >> Posterior parameters of the Normal distribution over a Normal mean
%      computed from prior parameters and the sample mean.
%
% FORMAT help spm_prob>function
%   >> Returns the help file of the selected function.
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging
    if nargin == 0
        help spm_prob>normal
        error('Not enough argument. Type ''help spm_prob>normal'' for help.');
    end
    id = varargin{1};
    varargin = varargin(2:end);
    switch lower(id)
        case 'pdf'
            [varargout{1:nargout}] = normal_pdf(varargin{:});
        case {'logpdf', 'll'}
            [varargout{1:nargout}] = normal_logpdf(varargin{:});
        case 'kl'
            [varargout{1:nargout}] = normal_kl(varargin{:});
        case {'up', 'update'}
            [varargout{1:nargout}] = normal_up(varargin{:});
        otherwise
            help spm_prob>normal
            error('Unknown function %s. Type ''help spm_prob>normal'' for help.', id)
    end
end

% -------------------------------------------------------------------------

function pdf = normal_pdf(x, mu, varargin)
% FORMAT pdf = spm_prob('Normal', 'pdf',    x, mu,    sigma)
% FORMAT pdf = spm_prob('Normal', 'pdf',    x, mu,    lambda, 'precision')
% FORMAT pdf = spm_prob('Normal', 'pdf',    x, mu, n, sigma)
% FORMAT pdf = spm_prob('Normal', 'pdf',    x, mu, n, lambda, 'precision')

    % Check if we are in the reparameterised case
    if nargin == 5
        if startsWith(varargin{3}, 'p', 'IgnoreCase', true);
            pdf = normal_pdf(x, mu, varargin{2}*varargin{1}, 'precision');
        else
            pdf = normal_pdf(x, mu, varargin{2}/varargin{1});
        end
        return
    elseif nargin == 4 && ~ischar(varargin{2})
        pdf = normal_pdf(x, mu, varargin{2}/varargin{1});
        return
    end
    
    % Else, set default values
    if nargin < 4
        mode = 'covariance';
        if nargin < 2
            mu = zeros(size(x));
        end
        if nargin < 3
            varargin{1} = eye(numel(mu));
        end
    else
        mode = varargin{2};
    end
    precision = startsWith(mode, 'p', 'IgnoreCase', true);
    
    % Usual PDF
    if precision
        pdf = (det(varargin{1}/(2*pi)))^(0.5) * exp(-0.5*(x(:)-mu(:))'*varargin{1}*(x(:)-mu(:)));
    else
        pdf = (det(varargin{1}*2*pi))^(-0.5) * exp(-0.5*(x(:)-mu(:))'*(varargin{1}\(x(:)-mu(:))));
    end
    
end

% -------------------------------------------------------------------------

function pdf = normal_logpdf(x, mu, varargin)
% FORMAT ll  = spm_prob('Normal', 'logpdf', x, mu,    sigma)
% FORMAT ll  = spm_prob('Normal', 'logpdf', x, mu,    lambda, 'precision')
% FORMAT ll  = spm_prob('Normal', 'logpdf', x, mu, n, sigma)
% FORMAT ll  = spm_prob('Normal', 'logpdf', x, mu, n, lambda, 'precision')

    % Check if we are in the reparameterised case
    if nargin == 5
        if startsWith(varargin{3}, 'p', 'IgnoreCase', true);
            pdf = normal_pdf(x, mu, varargin{2}*varargin{1}, 'precision');
        else
            pdf = normal_pdf(x, mu, varargin{2}/varargin{1});
        end
        return
    elseif nargin == 4 && ~ischar(varargin{2})
        pdf = normal_pdf(x, mu, varargin{2}/varargin{1});
        return
    end
    
    % Else, set default values
    if nargin < 4
        mode = 'covariance';
        if nargin < 2
            mu = zeros(size(x));
        end
        if nargin < 3
            varargin{1} = eye(numel(mu));
        end
    else
        mode = varargin{2};
    end
    precision = startsWith(mode, 'p', 'IgnoreCase', true);
    K = size(varargin{1}, 1);
    
    if precision
        pdf = -0.5*( K*log(2*pi) - spm_matcomp('LogDet', varargin{1}) + (x(:)-mu(:))'*varargin{1}*(x(:)-mu(:)) );
    else
        pdf = -0.5*( K*log(2*pi) + spm_matcomp('LogDet', varargin{1}) + (x(:)-mu(:))'*(varargin{1}\(x(:)-mu(:))) );
    end
end

% -------------------------------------------------------------------------

function kl = normal_kl(mu1, par1, mu0, par0, varargin)
% FORMAT kl  = spm_prob('Normal', 'kl', mu1, sigma1,  mu0, sigma0)
% FORMAT kl  = spm_prob('Normal', 'kl', mu1, lambda1, mu0, lambda0, 'precision')
% FORMAT kl  = spm_prob('Normal', 'kl', mu1, n1,      mu0, n0,      sigma)
% FORMAT kl  = spm_prob('Normal', 'kl', mu1, n1,      mu0, n0,      lambda, 'precision')

    % Check if we are in the reparameterised case
    if nargin == 6
        K = size(varargin{1}, 1);
        if startsWith(varargin{2}, 'p', 'IgnoreCase', true);
            kl = 0.5 * ( K * (par0/par1 - log(par0/par1) - 1) ...
                         + par0 * (mu0(:)-mu1(:))'*varargin{1}*(mu0(:)-mu1(:)) );
        else
            kl = 0.5 * ( K * (par0/par1 - log(par0/par1) - 1) ...
                         + par0 * (mu0(:)-mu1(:))'*(varargin{1}\(mu0(:)-mu1(:))) );
        end
        return
    elseif nargin == 5 && ~ischar(varargin{1})
        K  = size(varargin{1}, 1);
        kl = 0.5 * ( K * (par0/par1 - log(par0/par1) - 1) ...
                     + par0 * (mu0(:)-mu1(:))'*(varargin{1}\(mu0(:)-mu1(:))) );
        return
    end
    
    % Else, set default values
    if nargin < 5
        mode = 'covariance';
    else
        mode = varargin{1};
    end
    precision = startsWith(mode, 'p', 'IgnoreCase', true);
    K = size(par1, 1);
    
    % Common KL-divergence
    if precision
        kl = 0.5 * ( trace(par1\par0) ...
                     - spm_matcomp('LogDet', par0) ...
                     + spm_matcomp('LogDet', par1) ...
                     - K ...
                     + (mu0(:)-mu1(:))'*par0*(mu0(:)-mu1(:)) );
    else
        kl = 0.5 * ( trace(par0\par1) ...
                     - spm_matcomp('LogDet', par1) ...
                     + spm_matcomp('LogDet', par0) ...
                     - K ...
                     + (mu0(:)-mu1(:))'*(par1\(mu0(:)-mu1(:))) );
    end
    
end

% -------------------------------------------------------------------------

function [mu1, n1] = normal_up(mu, n, mu0, n0)
% FORMAT [mu1, n1] = spm_prob('Normal', 'update', mu, n, mu0, n0)

    n1 = n + n0;
    mu1 = (n0 * mu0 + n * mu)/n1;

end

% =========================================================================
%   GAMMA
% =========================================================================

function varargout = gamma(varargin)
%__________________________________________________________________________
% Characteristic functions of the Gamma distribution:
%   [pdf]       Probability density function
%   [ll/logpdf] Log-probability density function
%   [kl]        Kullback-Leibler divergence
%   [up/update] Conjugate (or Bayesian) update
%
% FORMAT pdf = spm_prob('Gamma', 'pdf',    x, alpha, beta)
% FORMAT ll  = spm_prob('Gamma', 'logpdf', x, alpha, alpha)
%   >> (Log) Probability density function.
%
% FORMAT pdf = spm_prob('Gamma', 'pdf',    x, lambda, n, K)
% FORMAT ll  = spm_prob('Gamma', 'logpdf', x, lambda, n, K)
%   >> Reparameterisation of the (log)-PDF when the Gamma distribution is
%      used as a distribution for an univariate Normal precision with known 
%      mean (K=1) or for a multivariate Normal precision magnitude with 
%      known mean (K>1).
%
% FORMAT kl  = spm_prob('Gamma', 'kl', alpha1, beta1, alpha0, beta0)
%   >> Kullback-Leibler divergence from G0 to G1 = KL(G1||G0).
%
% FORMAT kl  = spm_prob('Gamma', 'kl', lambda1, n1, lambda0, n0, K)
%   >> Reparameterisation of the KL-divergence when the Gamma distribution
%      is used as a conjugate prior for an univariate Normal precision with  
%      known mean (K=1) or for a multivariate Normal precision magnitude  
%      with known mean (K>1).
%
% FORMAT [lam1, n1] = spm_prob('Gamma', 'up', lam, n,     lam0, n0)
% FORMAT [lam1, n1] = spm_prob('Gamma', 'up', s0, s1, s2, lam0, n0, 
%                                             (mu=0), (Lam=eye))
%   >> Posterior parameters of the Gamma distribution over an univariate 
%      Normal precision or a multivariate Normal precision magnitude 
%      with known mean.
%
% FORMAT help spm_prob>function
%   >> Returns the help file of the selected function.
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging
    if nargin == 0
        help spm_prob>gamma
        error('Not enough argument. Type ''help spm_prob>gamma'' for help.');
    end
    id = varargin{1};
    varargin = varargin(2:end);
    switch lower(id)
        case 'pdf'
            [varargout{1:nargout}] = gamma_pdf(varargin{:});
        case {'logpdf', 'll'}
            [varargout{1:nargout}] = gamma_logpdf(varargin{:});
        case 'kl'
            [varargout{1:nargout}] = gamma_kl(varargin{:});
        case {'up', 'update'}
            [varargout{1:nargout}] = gamma_up(varargin{:});
        otherwise
            help spm_prob>gamma
            error('Unknown function %s. Type ''help spm_prob>gamma'' for help.', id)
    end
end

% -------------------------------------------------------------------------

function pdf = gamma_pdf(x, varargin)
% FORMAT pdf = gamma_pdf(x, alpha,  beta)
% FORMAT pdf = gamma_pdf(x, lambda, n,    K)

    % Check if we are in the reparameterised case
    if nargin == 4
        pdf = gamma_pdf(x, 0.5*varargin{3}*varargin{2}, ...
                           0.5*varargin{3}*varargin{2}/varargin{1});
        return
    end
    
    % Usual pdf
    alpha = varargin{1};
    beta  = varargin{2};
    pdf   = beta^alpha * x^(alpha-1) .* exp(-beta*x) / gamma(alpha);
    
end

% -------------------------------------------------------------------------

function pdf = gamma_logpdf(x, varargin)
% FORMAT pdf = gamma_logpdf(x, alpha,  beta)
% FORMAT pdf = gamma_logpdf(x, lambda, n,    K)

    % Check if we are in the reparameterised case
    if nargin == 4
        pdf = gamma_logpdf(x, 0.5*varargin{3}*varargin{2}, ...
                              0.5*varargin{3}*varargin{2}/varargin{1});
        return
    end
    
    % Usual pdf
    alpha = varargin{1};
    beta  = varargin{2};
    pdf   = alpha*log(beta) + (alpha-1)*log(x) - beta*x - gammaln(alpha);
    
end

% -------------------------------------------------------------------------

function kl = gamma_kl(varargin)
% FORMAT pdf = gamma_pdf(alpha1,  beta1, alpha0,  beta0)
% FORMAT pdf = gamma_pdf(lambda1, n1,    lambda0, n0,    K)

    % Check if we are in the reparameterised case
    if nargin == 5
        kl = gamma_kl(0.5*varargin{5}*varargin{2}, ...
                      0.5*varargin{5}*varargin{2}/varargin{1}, ...
                      0.5*varargin{5}*varargin{4}, ...
                      0.5*varargin{5}*varargin{4}/varargin{3});
        return
    end
    
    % Usual KL
    alpha1 = varargin{1};
    beta1  = varargin{2};
    alpha0 = varargin{3};
    beta0  = varargin{4};
    kl     = - alpha0 * log(beta0/beta1)         ...
             + alpha1 * (beta0/beta1 - 1)        ...
             + gammaln(alpha0) - gammaln(alpha1) ...
             + (alpha1- alpha0) * psi(alpha1);
    
end

% -------------------------------------------------------------------------

function [lambda1, n1] = gamma_up(varargin)
% FORMAT [lambda1, n1] = gamma_up(lambda, n,     lambda0, n0)
% FORMAT [lambda1, n1] = gamma_up(ss0, ss1, ss2, lambda0, n0, (mu=0), (Lambda=eye))

    if nargin > 4
        % Sufficient statistics case
        ss0     = varargin{1};
        ss1     = varargin{2};
        ss2     = varargin{3};
        lambda0 = varargin{4};
        n0      = varargin{5};
        K  = size(ss2, 1);
        if nargin < 7
            Lambda = eye(K);
            if nargin < 6
                mu = zeros(size(ss1));
            else
                mu = varargin{6};
            end
        else
            Lambda = varargin{7};
        end
        n      = ss0;
        lambda = trace(ss2*Lambda) ...
                 - 2 * mu'*Lambda*ss1 ...
                 + ss0 * mu'*Lambda*mu;
        lambda = K*ss0/lambda;
        [lambda1, n1] = gamma_up(lambda, n, lambda0, n0);
    else
        % Average case
        lambda  = varargin{1};
        n       = varargin{2};
        lambda0 = varargin{3};
        n0      = varargin{4};
        n1      = n + n0;
        lambda1 = n1 / (n0/lambda0 + n/lambda);
    end
    
    
end