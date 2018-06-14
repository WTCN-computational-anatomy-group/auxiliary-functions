function varargout = spm_impreproc(varargin)
%__________________________________________________________________________
% Collection of tools for image pre-processing.
%
% FORMAT [Affine,bb] = spm_impreproc('atlas_crop',P,Affine,prefix,rem_neck)
% FORMAT spm_impreproc('nm_reorient',Vin,vx,type)
% FORMAT spm_impreproc('reset_origin',P)
% FORMAT R = spm_impreproc('rigid_align',P)
% FORMAT V = spm_impreproc('reg_and_reslice',V)
% FORMAT spm_impreproc('subvol',V,bb,prefix)
% FORMAT nfname = spm_impreproc('downsample_inplane',fname,vx1)
%
% FORMAT help spm_impreproc>function
% Returns the help file of the selected function.
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging
if nargin == 0
    help spm_impreproc
    error('Not enough argument. Type ''help spm_impreproc'' for help.');
end
id       = varargin{1};
varargin = varargin(2:end);
switch lower(id)
    case 'atlas_crop'
        [varargout{1:nargout}] = atlas_crop(varargin{:});  
    case 'nm_reorient'
        [varargout{1:nargout}] = nm_reorient(varargin{:});                   
    case 'coreg'
        [varargout{1:nargout}] = coreg(varargin{:});                 
    case 'reset_origin'
        [varargout{1:nargout}] = reset_origin(varargin{:});              
    case 'reslice'
        [varargout{1:nargout}] = reslice(varargin{:});                     
    case 'rigid_align'
        [varargout{1:nargout}] = rigid_align(varargin{:});                  
    case 'subvol'
        [varargout{1:nargout}] = subvol(varargin{:});    
    case 'downsample_inplane'
        [varargout{1:nargout}] = downsample_inplane(varargin{:});            
    otherwise
        help spm_impreproc
        error('Unknown function %s. Type ''help spm_impreproc'' for help.', id)
end
%==========================================================================

%==========================================================================
function [Affine,bb] = atlas_crop(P,prefix,rem_neck)
% Removes air outside of head
% FORMAT [Affine,bb] = atlas_crop(P,Affine,prefix,rem_neck)
% P        - Path to NIfTI file
% prefix   - File prefix (if empty -> overwrites) ['']
% rem_neck - Remove neck/spine [false]
% bb - Computed bounding box
%
% This function rigidly registers the SPM atlas to an image and then
% removes image data outside of the head.
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging
if nargin<2, prefix   = ''; end
if nargin<3, rem_neck = false; end

% Locate TPM.nii in SPM
pth_tpm = fullfile(spm('dir'),'tpm','TPM.nii,');

Vin  = spm_vol(P);
Vtpm = spm_vol(pth_tpm);

mat    = Vin.mat;
mattpm = Vtpm.mat;

tpm = spm_load_priors8(Vtpm);    

V = spm_vol(P);

M               = V(1).mat;
c               = (V(1).dim+1)/2;
V(1).mat(1:3,4) = -M(1:3,1:3)*c(:);
[Affine1,ll1]   = spm_maff8(V(1),8,(0+1)*16,tpm,[],'mni'); % Closer to rigid
Affine1         = Affine1*(V(1).mat/M);

% Run using the origin from the header
V(1).mat      = M;
[Affine2,ll2] = spm_maff8(V(1),8,(0+1)*16,tpm,[],'mni'); % Closer to rigid

% Pick the result with the best fit
if ll1>ll2, Affine  = Affine1; else Affine  = Affine2; end

Affine = spm_maff8(P,2,32,tpm,Affine,'mni');
Affine = spm_maff8(P,2,1,tpm,Affine,'mni');

% Voxel locations in TPM.nii
Ltpm1 = [120 72.2 37.3 1]'; Ltpm2 = [120 72.2 75.9 1]';
Rtpm1 = [3  72.2 37.3 1]'; Rtpm2 = [3  72.2 75.9 1]';

Stpm1 = [58.6 42.6 119 1]'; Stpm2 = [58.60 99.4 119 1]';
if rem_neck
    Itpm1 = [58.6 39.4 2.5   1]'; Itpm2 = [58.60 99.4 2.5  1]';    
else
    Itpm1 = [58.6 39.4 -200  1]'; Itpm2 = [58.60 99.4 -200 1]';
end
Atpm1 = [58.6 144 28.4 1]'; Atpm2 = [58.60 144 82.3 1]';
Ptpm1 = [58.6 3.5  28.4 1]'; Ptpm2 = [58.60 3.5 82.3 1]'; 

% Voxel locations in input
T  = mat\(Affine\mattpm);
L1 = T*Ltpm1; L2 = T*Ltpm2;
R1 = T*Rtpm1; R2 = T*Rtpm2;
U1 = T*Stpm1; U2 = T*Stpm2;
D1 = T*Itpm1; D2 = T*Itpm2;
A1 = T*Atpm1; A2 = T*Atpm2;
P1 = T*Ptpm1; P2 = T*Ptpm2;

% Bounding-box
bb = zeros(2,3);
for i=1:3
    X       = [L1(i) R1(i) U1(i) D1(i) A1(i) P1(i)...
               L2(i) R2(i) U2(i) D2(i) A2(i) P2(i)];
    bb(1,i) = max(X);
    bb(2,i) = min(X);
end

% Do cropping
spm_impreproc('subvol',Vin,bb,prefix);      
%==========================================================================

%==========================================================================
function nm_reorient(Vin,vx,type,prefix)
% Re-orient images
% FORMAT nm_reorient(Vin,vx,type)
% Vin  - SPM volume objects
% vx   - New voxel size
% type - Order of interpolation
% prefix - Prefix of file to write
%
% The function reslices the input images to a resolution of vx mm.
% Output images (with the prefix "pn_r") are written in the transverse
% orientation (using information from the ".mat" files).
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging
spm('defaults', 'FMRI');

if length(vx)<3
    vx=[vx vx vx];
end
if nargin<4, prefix = 'ro_'; end

% If no arguments, then prompt for images
%PP = spm_get([1 Inf],'*.img','Select files to reorient');

% Get information about the image volumes
VV = spm_vol(Vin);

for V=VV', % Loop over images

    % The corners of the current volume
    d = V.dim(1:3);
    c = [	1    1    1    1
        1    1    d(3) 1
        1    d(2) 1    1
        1    d(2) d(3) 1
        d(1) 1    1    1
        d(1) 1    d(3) 1
        d(1) d(2) 1    1
        d(1) d(2) d(3) 1]';

    % The corners of the volume in mm space
    tc = V.mat(1:3,1:4)*c;
    if spm_flip_analyze_images, tc(1,:) = -tc(1,:); end;

    % Max and min co-ordinates for determining a bounding-box
    mx = round(max(tc,[],2)');
    mn = round(min(tc,[],2)');

    % Translate so that minimum moves to [1,1,1]
    % This is the key bit for changing voxel sizes,
    % output orientations etc.
    mat = spm_matrix(mn)*diag([vx 1])*spm_matrix(-[1 1 1]);

    % Dimensions in mm
    dim = ceil((mat\[mx 1]')');

    % Output image based on information from the original
    VO               = V;

    % Create a filename for the output image (prefixed by 'r')
    [lpath,name,ext] = fileparts(V.fname);
    VO.fname         = fullfile(lpath,[prefix name ext]);

    % Dimensions of output image
    VO.dim(1:3)      = dim(1:3);

    % Voxel-to-world transform of output image
    if spm_flip_analyze_images, mat = diag([-1 1 1 1])*mat; end;
    VO.mat           = mat;

    % Initialise plot of how far reslicing has gone
    %spm_progress_bar('Init',dim(3),'reslicing...','planes completed');

    % Create .hdr and open output .img
    VO = spm_create_vol(VO);

    for i=1:dim(3), % Loop over slices of output image

        % Mapping from slice i of the output image,
        % to voxels of the input image
        M   = inv(spm_matrix([0 0 -i])*inv(VO.mat)*V.mat);

        % Extract this slice according to the mapping
        img = spm_slice_vol(V,M,dim(1:2),type);

        % Write this slice to output image
        spm_write_plane(VO,img,i);

        % Update the progress bar
        %spm_progress_bar('Set',i);

    end; % End loop over output slices

    % Get rid of the progress bar
    %spm_progress_bar('Clear');

end; % End loop over images
%==========================================================================

%==========================================================================
function reset_origin(P,orig)
% Reset origin of image
% FORMAT reset_origin(P)
% P = Path to NIfTI image
%
% OBS: Image will have the matrix in its header adjusted.
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging
if nargin<2, orig = []; end

V   = spm_vol(P);
M   = V.mat;
dim = V.dim;
vx  = sqrt(sum(M(1:3,1:3).^2));

if det(M(1:3,1:3))<0
    vx(1) = -vx(1); 
end;

if isempty(orig)
    orig = (dim(1:3)+1)/2;
end

off  = -vx.*orig;
M1   = [vx(1) 0      0      off(1)
           0      vx(2) 0      off(2)
           0      0      vx(3) off(3)
           0      0      0      1];

spm_get_space(P,M1);   
%==========================================================================

%==========================================================================
function R = rigid_align(P)
% Reposition an image by affine aligning to MNI space and Procrustes adjustment
% FORMAT rigid_align(P)
% P - name of NIfTI image
% R - Affine matrix
%
% OBS: Image will have the matrix in its header adjusted.
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

% Load tissue probability data
tpm = fullfile(spm('dir'),'tpm','TPM.nii,');
tpm = [repmat(tpm,[6 1]) num2str((1:6)')];
tpm = spm_load_priors8(tpm);

% Do the affine registration
V = spm_vol(P);

M               = V(1).mat;
c               = (V(1).dim+1)/2;
V(1).mat(1:3,4) = -M(1:3,1:3)*c(:);
[Affine1,ll1]   = spm_maff8(V(1),8,(0+1)*16,tpm,[],'rigid'); % Closer to rigid
Affine1         = Affine1*(V(1).mat/M);

% Run using the origin from the header
V(1).mat      = M;
[Affine2,ll2] = spm_maff8(V(1),8,(0+1)*16,tpm,[],'rigid'); % Closer to rigid

% Pick the result with the best fit
if ll1>ll2, Affine  = Affine1; else Affine  = Affine2; end

Affine = spm_maff8(P,2,32,tpm,Affine,'rigid'); % Heavily regularised
Affine = spm_maff8(P,2,1 ,tpm,Affine,'rigid'); % Lightly regularised

% Load header
Nii    = nifti(P);

% Generate mm coordinates of where deformations map from
x      = affind(rgrid(size(tpm.dat{1})),tpm.M);

% Generate mm coordinates of where deformation maps to
y1     = affind(x,inv(Affine));

% Weight the transform via GM+WM
weight = single(exp(tpm.dat{1})+exp(tpm.dat{2}));

% Weighted Procrustes analysis
[~,R]  = spm_get_closest_affine(x,y1,weight);

% Invert
% R      = inv(R);

% Write the new matrix to the header
Nii.mat = R\Nii.mat;
create(Nii);
%==========================================================================

%==========================================================================
function V = coreg(V)
% Co-register images
% FORMAT V = coreg(V)
% V - SPM volume object that can contain N different modalities (e.g. T1- 
% and T2-weighted MRIs.
%
% WARNING: This function overwrites orientation matrices!
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging
N = numel(V);
if N==1
    return;
end

% Get image with smallest voxel size and pick this image as reference
prod_vx = zeros(1,N);
for n=1:N
    vx         = spm_misc('vxsize',V(n).mat);
    prod_vx(n) = prod(vx);
end
[~,ref_ix] = min(prod_vx);

% Set options
matlabbatch{1}.spm.spatial.coreg.estimate.ref               = {V(ref_ix).fname};
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep      = [4 2];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol      = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm     = [7 7];

% Co-register
ixs       = 1:N;
source_ix = ixs(ixs~=ref_ix);
for n=source_ix
    matlabbatch{1}.spm.spatial.coreg.estimate.source = {V(n).fname};         

    spm_jobman('run',matlabbatch);
end
%==========================================================================

%==========================================================================
function V = reslice(V,ref_ix)
% Re-slice images
% FORMAT V = reslice(V)
% V - SPM volume object that can contain N different modalities (e.g. T1- 
% and T2-weighted MRIs.
% ref_ix - index of reference image in V
%
% Takes medical images of the same subject and re-slices the images to the
% same dimensions. If no reference index is given, the image with the largest field of view is chosen as
% reference for the re-slicing. First order interpolation is used not to
% introduce any negative values.
%
% WARNING: This function overwrites the input data!
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging
N = numel(V);
if N==1
    return;
end

if nargin<2
    % Get image with largest volume (for reslicing using this image as
    % reference)
    vol = zeros(N,3);
    for n=1:N
        vx       = spm_misc('vxsize',V(n).mat);
        vol(n,:) = vx.*V(n).dim;
    end
    vol        = prod(vol,2);
    [~,ref_ix] = max(vol);
end

% Set options
matlabbatch{1}.spm.spatial.coreg.write.ref             = {V(ref_ix).fname};
matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 1;
matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap   = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.write.roptions.mask   = 0;
matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'res_';        

% Re-slice
ixs       = 1:N;
source_ix = ixs(ixs~=ref_ix);
for n=source_ix
    matlabbatch{1}.spm.spatial.coreg.write.source = {V(n).fname};     

    output_list = spm_jobman('run',matlabbatch);

    delete(V(n).fname);
    V(n) = spm_vol(output_list{1}.rfiles{1});    
end
%==========================================================================

%==========================================================================
function VO = subvol(V,bb,prefix)
% Extract a subvolume
% FORMAT VO = subvol(V,bb,prefix)
% V      - SPM volume object
% bb     - bounding box
% prefix - file prefix (if empty -> overwrites)
% VO     - resized image
%
% Example:
%     V = spm_vol(spm_select(1,'image'));
%     subvol(V,[32 64 ; 1 64 ; 1 48]');
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging
if nargin<3, prefix = 'sv_'; end

bb      = round(bb);
bb      = sort(bb);
bb(1,:) = max(bb(1,:),[1 1 1]);
bb(2,:) = min(bb(2,:),V.dim(1:3));

VO            = V;
[pth,nam,ext] = fileparts(V.fname);
VO.fname      = fullfile(pth,[prefix nam ext]);
VO.dim(1:3)   = diff(bb)+1;
VO.mat        = V.mat*spm_matrix((bb(1,:)-1));

VO = spm_create_vol(VO);
for z=1:VO.dim(3),
    M   = V.mat\VO.mat*spm_matrix([0 0 z]);
    img = spm_slice_vol(V,M,VO.dim(1:2),0);
    VO  = spm_write_plane(VO,img,z);
end
%==========================================================================

%==========================================================================
function nfname = downsample_inplane(fname,vx1)
% Down-sample a NIfTI image in the high-resolution plane
% FORMAT nfname = downsample_inplane(fname,vx1)
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging  
if nargin<2, vx1 = 1; end

if numel(vx1)==1
    vx1 = [vx1(1) vx1(1) vx1(1)];
end

Nii  = nifti(fname);
mat0 = Nii.mat;             

% Get down-sampling factor
vx0       = spm_misc('vxsize',mat0);   
d         = ((vx0 < vx1).*vx0)./vx1;
d(d == 0) = 1;   

if sum(d)==3
    % Do not downsample if in-plane res Xhat equals in-plane res Y
    warning('do_dsinp::false')
    return
end

% Downsample
%--------------------------------------------------------------------------
D      = diag([d, 1]);          
mat_ds = mat0/D;
vx_ds  = spm_misc('vxsize',mat_ds);

X   = Nii.dat(:,:,:);     
dm0 = size(X);    

% fwhm = max(vx_ds./vx0 - 1,0.01);        
% spm_imbasics('smooth_img_in_mem',X,fwhm);                                                 
         
C               = spm_bsplinc(X,[0 0 0 0 0 0]); % Resample using 0th order b-splines (NN)
[x1,y1,z1]      = get_downsampling_grid(D,dm0);                  
X               = spm_bsplins(C,x1,y1,z1,[1 1 1 0 0 0]);
X(~isfinite(X)) = 0;

fname         = Nii.dat.fname;
[pth,nam,ext] = fileparts(fname);
nfname        = fullfile(pth,['ds_' nam ext]);

spm_misc('create_nii',nfname,X,mat_ds,Nii.dat.dtype,Nii.descrip);
%==========================================================================

%==========================================================================
% HELPER FUNCTIONS
%==========================================================================

%==========================================================================
function x = rgrid(d)
x = zeros([d(1:3) 3],'single');
[x1,x2] = ndgrid(single(1:d(1)),single(1:d(2)));
for i=1:d(3),
    x(:,:,i,1) = x1;
    x(:,:,i,2) = x2;
    x(:,:,i,3) = single(i);
end
%==========================================================================

%==========================================================================
function y1 = affind(y0,M)
y1 = zeros(size(y0),'single');
for d=1:3,
    y1(:,:,:,d) = y0(:,:,:,1)*M(d,1) + y0(:,:,:,2)*M(d,2) + y0(:,:,:,3)*M(d,3) + M(d,4);
end
%==========================================================================

%==========================================================================
function [x1,y1,z1] = get_downsampling_grid(M,dm)
T          = eye(4)/M;   
dm         = floor(M(1:3,1:3)*dm')';
[x0,y0,z0] = ndgrid(1:dm(1),...
                    1:dm(2),...
                    1:dm(3));

x1 = T(1,1)*x0 + T(1,2)*y0 + T(1,3)*z0 + T(1,4);
y1 = T(2,1)*x0 + T(2,2)*y0 + T(2,3)*z0 + T(2,4);
z1 = T(3,1)*x0 + T(3,2)*y0 + T(3,3)*z0 + T(3,4);  
%==========================================================================