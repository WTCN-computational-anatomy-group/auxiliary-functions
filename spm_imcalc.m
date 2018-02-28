function varargout = spm_imcalc(varargin)
%__________________________________________________________________________
% Collection of tools for image calculation (gradient, suff stat, ...).
%
% FORMAT [Dx,Dy,Dz] = grad(X,vx)
%
% FORMAT help spm_imcalc>function
% Returns the help file of the selected function.
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

    if nargin == 0
        help spm_imcalc
        error('Not enough argument. Type ''help spm_imcalc'' for help.');
    end
    id = varargin{1};
    varargin = varargin(2:end);
    switch lower(id)
        case 'grad'
            [varargout{1:nargout}] = grad(varargin{:});
        otherwise
            help spm_imcalc
            error('Unknown function %s. Type ''help spm_imcalc'' for help.', id)
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
% HELPER FUNCTIONS
%==========================================================================

%==========================================================================
function out = get_type(var)
    tmp = whos('var');
    out = tmp.class;
end
%==========================================================================