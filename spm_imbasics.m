function varargout = spm_imbasics(varargin)
%__________________________________________________________________________
% Collection of tools for image calculation (gradient, suff stat, ...).
%
% FORMAT div = dive(Dx,Dy,Dz,vx)
% FORMAT [Dx,Dy,Dz] = grad(X,vx)
% FORMAT smooth_img_in_mem(img,fwhm) 
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
        otherwise
            help spm_imcalc
            error('Unknown function %s. Type ''help spm_imcalc'' for help.', id)
    end
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

    fwhm = fwhm*ones(1,3);

    lim = ceil(2*fwhm);
    x   = -lim(1):lim(1); x = spm_smoothkern(fwhm(1),x); x = x/sum(x);
    y   = -lim(2):lim(2); y = spm_smoothkern(fwhm(2),y); y = y/sum(y);
    z   = -lim(3):lim(3); z = spm_smoothkern(fwhm(3),z); z = z/sum(z);
    i   = (length(x) - 1)/2;
    j   = (length(y) - 1)/2;
    k   = (length(z) - 1)/2;
    spm_conv_vol(img,img,x,y,z,-[i j k]);
end
%==========================================================================

%==========================================================================
% HELPER FUNCTIONS
%==========================================================================

%==========================================================================
function out = get_type(var)
    tmp = whos('var');
    out = tmp.class;
end
%==========================================================================