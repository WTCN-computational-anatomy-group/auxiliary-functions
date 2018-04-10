function varargout = spm_imbasics(varargin)
%__________________________________________________________________________
% Collection of tools for image calculation (gradient, suff stat, ...).
%
% FORMAT div = dive(Dx,Dy,Dz,vx)
% FORMAT [Dx,Dy,Dz] = grad(X,vx)
% FORMAT smooth_img_in_mem(img,fwhm) 
% FORMAT [mg,mn,vr] = fit_gmm2hist(h,x,K,verbose)
%
% FORMAT help spm_imbasics>function
% Returns the help file of the selected function.
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging
if nargin == 0
    help spm_imbasics
    error('Not enough argument. Type ''help spm_imcalc'' for help.');
end
id = varargin{1};
varargin = varargin(2:end);
switch lower(id)
    case 'dive'
        [varargout{1:nargout}] = dive(varargin{:});        
    case 'grad'
        [varargout{1:nargout}] = grad(varargin{:});
    case 'smooth_img_in_mem'
        [varargout{1:nargout}] = smooth_img_in_mem(varargin{:});          
    case 'fit_gmm2hist'
        [varargout{1:nargout}] = fit_gmm2hist(varargin{:});                 
    otherwise
        help spm_imcalc
        error('Unknown function %s. Type ''help spm_imcalc'' for help.', id)
end
%==========================================================================

%==========================================================================
function div = dive(Dx,Dy,Dz,vx)  
% Computes the divergence of an image
% FORMAT div = dive(Dx,Dy,Dz,vx) 
% [Dx,Dy,Dz] - Gradients in x-,y- and z-direction
% vx         - Voxel size
% div        - Divergence
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging
if nargin<4, vx = ones(3,1); end

if size(Dx,3) == 1
    Du = [-Dx(:,1), -diff(Dx(:,1:end-1),1,2), Dx(:,end-1)];
    Dv = [-Dy(1,:); -diff(Dy(1:end-1,:),1,1); Dy(end-1,:)];
    div = Du./vx(2) + Dv./vx(1);
else
    Du = cat(2, -Dx(:,1,:), -diff(Dx(:,1:end-1,:),1,2), Dx(:,end-1,:)); 
    Dv = cat(1, -Dy(1,:,:), -diff(Dy(1:end-1,:,:),1,1), Dy(end-1,:,:));
    Dw = cat(3, -Dz(:,:,1), -diff(Dz(:,:,1:end-1),1,3), Dz(:,:,end-1));
    div = Du./vx(2) + Dv./vx(1) + Dw./vx(3);
end
%==========================================================================

%==========================================================================
function [Dx,Dy,Dz] = grad(X,vx) 
% Calculate 2D or 3D gradient of an image
% FORMAT [Dx,Dy,Dz] = grad(X,vx)
% X          - Image
% vx         - voxel size
% [Dx,Dy,Dz] - Gradients in x-,y- and z-direction
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging   
if nargin<2, vx = ones(3,1); end

precision = get_type(X);

if size(X,3)==1
    Dx = [diff(X,1,2),zeros(size(X,1),1,precision)]./vx(2);
    Dy = [diff(X,1,1);zeros(1,size(X,2),precision)]./vx(1);
    Dz = 0;
else
    Dx = cat(2,diff(X,1,2),zeros(size(X,1),1,size(X,3),precision))./vx(2);
    Dy = cat(1,diff(X,1,1),zeros(1,size(X,2),size(X,3),precision))./vx(1);
    Dz = cat(3,diff(X,1,3),zeros(size(X,1),size(X,2),1,precision))./vx(3);  
end
%==========================================================================

%==========================================================================
function smooth_img_in_mem(img,fwhm) 
% Smooth an image with a Gaussian kernel
% FORMAT smooth_img_in_mem(img,fwhm) 
% img          - Image
% fwhm         - Full-width at half maximum
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging   
if nargin<2, fwhm = 10; end

if numel(fwhm)==1
    fwhm = fwhm*ones(1,3);
end

lim = ceil(2*fwhm);
x   = -lim(1):lim(1); x = spm_smoothkern(fwhm(1),x); x = x/sum(x);
y   = -lim(2):lim(2); y = spm_smoothkern(fwhm(2),y); y = y/sum(y);
z   = -lim(3):lim(3); z = spm_smoothkern(fwhm(3),z); z = z/sum(z);
i   = (length(x) - 1)/2;
j   = (length(y) - 1)/2;
k   = (length(z) - 1)/2;
spm_conv_vol(img,img,x,y,z,-[i j k]);
%==========================================================================

%==========================================================================
function [mg,mn,vr] = fit_gmm2hist(h,x,K,verbose)
% Fit a GMM to image histogram
% FORMAT [mg,mn,vr] = fit_gmm2hist(h,x,K,verbose)
% h - Histogram counts
% x - Intensity values
% K - Clusters
% verbose - Output level [0]
% mg - Mixing weights
% mn - Means
% vr - Variances
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging   
if nargin<4, verbose = 0; end

mg = ones(K,1)/K;
mn = linspace(min(x),max(x),K)'./K;
sd = ones(K,1)*(max(x) - min(x))./(K);  
vr = zeros(size(sd));

m0    = zeros(K,1);
m1    = zeros(K,1);
m2    = zeros(K,1);
ll(1) = -Inf;
for iter=1:10000
    p  = zeros(numel(x),K);
    for k=1:K
        % Product Rule
        % p(class=k, intensity | mg, nu, sig) = p(class=k|mg) p(intensity | nu, sig, class=k)
        p(:,k) = mg(k)*normpdf(x(:),mn(k),sd(k));
    end

    % Sum Rule
    % p(intensity | mg, nu, sig) = \sum_k p(class=k, intensity | mg, nu, sig)
    sp         = sum(p,2)+eps;
    ll(iter+1) = sum(log(sp).*h(:));
    if ll(iter+1) - ll(iter) < 1e-8*sum(h)
        if verbose==2
            figure(4001);
            md = mean(diff(x));
            plot(x(:),(h/sum(h))/md,'b-',x(:),sp,'r-'); hold on
            plot(x(:),p,'--');        
            set(gca,'TickLabelInterpreter','latex');  
            xlabel('Image intensity','Interpreter','latex')
            ylabel('Probability','Interpreter','latex')
            legend({'Empirical','Fit','Air','Tissue'},'Interpreter','latex');
            drawnow;
        end
        break; 
    end

    if verbose == 3
        figure(4001);
        subplot(121); plot(0:numel(ll)-2,ll(2:end))  
        md = mean(diff(x));
        subplot(122); plot(x(:),p,'--',x(:),h/sum(h)/md,'b.',x(:),sp,'r'); 
        drawnow
    end

    % Bayes Rule
    % p(class=k | intensity, mg, nu, sig) = p(class=k, intensity | mg, nu, sig) / p(intensity | mg, nu, sig)
    p = bsxfun(@rdivide,p,sp);

    % Compute moments from the histograms, weighted by the responsibilities (p).
    for k=1:K
        m0(k) = sum(p(:,k).*h(:));             % Number of voxels in class k
        m1(k) = sum(p(:,k).*h(:).*x(:));       % Sum of the intensities in class k
        m2(k) = sum(p(:,k).*h(:).*x(:).*x(:)); % Sum of squares of intensities in class k
    end
    mg = m0/sum(m0);
    for k=1:K
        mn(k) = m1(k)./m0(k);                                % Mean
        vr(k) = (m2(k)-m1(k)*m1(k)/m0(k)+1e-6)/(m0(k)+1e-6); % Variance
    end
    sd = sqrt(vr);
end
%==========================================================================

%==========================================================================
% HELPER FUNCTIONS
%==========================================================================

%==========================================================================
function out = get_type(var)
tmp = whos('var');
out = tmp.class;
%==========================================================================