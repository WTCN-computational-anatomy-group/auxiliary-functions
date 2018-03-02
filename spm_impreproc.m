function varargout = spm_impreproc(varargin)
%__________________________________________________________________________
% Collection of tools for image pre-processing.
%
% FORMAT Affine = atlas_crop(P,Affine,prefix,rem_neck)
% FORMAT V = crop(V,rem_neck)
% FORMAT nm_reorient(Vin,vx,type)
% FORMAT reset_origin(P)
% FORMAT rigid_align(P)
% FORMAT V = realign2mni(V)
% FORMAT V = reg_and_reslice(V)
% FORMAT subvol(V,bb,prefix)
% FORMAT skullstrip(V)
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
        case 'crop'
            [varargout{1:nargout}] = crop(varargin{:});   
        case 'nm_reorient'
            [varargout{1:nargout}] = nm_reorient(varargin{:});                       
        case 'realign2mni'
            [varargout{1:nargout}] = realign2mni(varargin{:});                   
        case 'reg_and_reslice'
            [varargout{1:nargout}] = reg_and_reslice(varargin{:});                 
        case 'reset_origin'
            [varargout{1:nargout}] = reset_origin(varargin{:});              
        case 'rigid_align'
            [varargout{1:nargout}] = rigid_align(varargin{:});             
        case 'skullstrip'
            [varargout{1:nargout}] = skullstrip(varargin{:});                  
        case 'subvol'
            [varargout{1:nargout}] = subvol(varargin{:});                
        otherwise
            help spm_impreproc
            error('Unknown function %s. Type ''help spm_impreproc'' for help.', id)
    end
end
%==========================================================================

%==========================================================================
function Affine = atlas_crop(P,Affine,prefix,rem_neck)
% Removes air outside of head
% FORMAT Affine = atlas_crop(P,Affine,prefix,rem_neck)
% P        - Path to NIfTI file
% Affine   - Affine matrix [[]]
% prefix   - File prefix (if empty -> overwrites) ['']
% rem_neck - Remove neck/spine [false]
%
% This function rigidly registers the SPM atlas to an image and then
% removes image data outside of the head.
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

    if nargin<2, Affine   = []; end
    if nargin<3, prefix   = ''; end
    if nargin<4, rem_neck = false; end

    % Locate TPM.nii in SPM
    pth_tpm = fullfile(spm('dir'),'tpm','TPM.nii,');

    Vin  = spm_vol(P);
    Vtpm = spm_vol(pth_tpm);

    mat    = Vin.mat;
    mattpm = Vtpm.mat;

    if isempty(Affine)
        tpm    = spm_load_priors8(Vtpm);        
        Affine = spm_maff_new(P,4,(0+1)*16,tpm,[],'mni');
    end

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
    for i=1:3
        X     = [L1(i) R1(i) U1(i) D1(i) A1(i) P1(i)...
                 L2(i) R2(i) U2(i) D2(i) A2(i) P2(i)];
        mx(i) = max(X);
        mn(i) = min(X);
    end

    % Do cropping
    spm_impreproc('subvol',Vin,[mn(1) mx(1);mn(2) mx(2);mn(3) mx(3)]',prefix);      
end
%==========================================================================

%==========================================================================
function V = crop(V,rem_neck)
% Removes air outside of head from N-channel data
% FORMAT V = crop(V,rem_neck)
% V        - SPM volume object
% rem_neck - Remove neck/spine [false]
%
% This function rigidly registers the SPM atlas to an image and then
% removes image data outside of the head.
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

    if nargin<2, rem_neck = false; end

    N      = numel(V);
    Affine = [];
    for n=1:N          
        if n==1
            Affine = spm_impreproc('atlas_crop',V(n).fname,[],'sv_',rem_neck);     
        else
            spm_impreproc('atlas_crop',V(n).fname,Affine,'sv_',rem_neck);                 
        end

        [pth,nam,ext] = fileparts(V(n).fname);
        delete(V(n).fname);

        nfname = fullfile(pth,['sv_' nam ext]);
        V(n)   = spm_vol(nfname);
    end
end
%==========================================================================

%==========================================================================
function nm_reorient(Vin,vx,type)
% Re-orient images
% FORMAT nm_reorient(Vin,vx,type)
% Vin  - SPM volume objects
% vx   - New voxel size
% type - Order of interpolation
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
        VO.fname         = fullfile(lpath,['ro_' name ext]);

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
end
%==========================================================================

%==========================================================================
function reset_origin(P)
% Reset origin of image
% FORMAT reset_origin(P)
% P = Path to NIfTI image
%
% OBS: Image will have the matrix in its header adjusted.
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

    V   = spm_vol(P);
    M   = V.mat;
    dim = V.dim;
    vx  = sqrt(sum(M(1:3,1:3).^2));

    if det(M(1:3,1:3))<0
        vx(1) = -vx(1); 
    end;

    orig = (dim(1:3)+1)/2;
    off  = -vx.*orig;
    M    = [vx(1) 0      0      off(1)
               0      vx(2) 0      off(2)
               0      0      vx(3) off(3)
               0      0      0      1];

    spm_get_space(P,M);    
end
%==========================================================================

%==========================================================================
function rigid_align(P)
% Reposition an image by affine aligning to MNI space and Procrustes adjustment
% FORMAT rigid_align(P)
% P - name of NIfTI image
%
% OBS: Image will have the matrix in its header adjusted.
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

    % Load tissue probability data
    tpm = fullfile(spm('dir'),'tpm','TPM.nii,');
    tpm = [repmat(tpm,[6 1]) num2str((1:6)')];
    tpm = spm_load_priors8(tpm);

    % Do the affine registration
    Affine = eye(4);
    Affine = spm_maff8(P,2,32,tpm,Affine,'mni'); % Heavily regularised
    Affine = spm_maff8(P,2,1 ,tpm,Affine,'mni'); % Lightly regularised

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
end
%==========================================================================

%==========================================================================
function V = realign2mni(V)
% Rigidly align subject's images to MNI space
% FORMAT V = realign2mni(V)
% V - SPM volume object 
%
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

    N = numel(V);
    
    for n=1:N              
        vx = vxsize(V(n).mat);

        spm_impreproc('nm_reorient',V(n).fname,vx,0);          

        [pth,nam,ext] = fileparts(V(n).fname);
        delete(V(n).fname);

        nfname = fullfile(pth,['ro_' nam ext]);
        V(n)   = spm_vol(nfname);

        spm_impreproc('reset_origin',V(n).fname);  

        spm_impreproc('rigid_align',V(n).fname);          
    end
end
%==========================================================================

%==========================================================================
function V = reg_and_reslice(V)
% Co-register and reslice images
% FORMAT V = reg_and_reslice(V)
% V - SPM volume object that can contain N different modalities (e.g. T1- 
% and T2-weighted MRIs.
%
% Takes medical images of the same subject and co-registers them and also
% re-slices the images to the same dimensions. The image with the largest
% field of view is chosen as reference for the re-slicing. First order
% interpolation is used not to introduce any negative values.
%
% WARNING: This function overwrites the input data!
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

    N = numel(V);
    if N==1
        return;
    end
    
    % Get image with largest volume and reslice using this image as reference
    vol = zeros(N,3);
    for n=1:N
        vx       = sqrt(sum(V(n).mat(1:3,1:3).^2));
        vol(n,:) = vx.*V(n).dim;
    end
    vol        = prod(vol,2);
    [~,ref_ix] = max(vol);

    % Set options
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep      = [4 2];
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol      = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm     = [7 7];
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp   = 1;
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap     = [0 0 0];
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask     = 0;
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix   = 'r';

    matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {V(ref_ix).fname};

    % Register and reslice
    ixs       = 1:N;
    source_ix = ixs(ixs~=ref_ix);
    for n=source_ix
        matlabbatch{1}.spm.spatial.coreg.estwrite.source = {V(n).fname};                        

        output_list = spm_jobman('run',matlabbatch);
        delete(V(n).fname);

        V(n) = spm_vol(output_list{1}.rfiles{1});            
    end
end
%==========================================================================

%==========================================================================
function skullstrip(V)
% Skull-strip subject images
% FORMAT V = skullstrip(V)
% V - SPM volume object that can contain N different modalities (e.g. T1- 
% and T2-weighted MRIs.
%
% WARNING: This function overwrites the input data!
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

    obj          = struct;
    obj.bb       = NaN(2,3);
    obj.bb       = [-90 -126 -72; 90 90 108];
    obj.vox      = 2;
    obj.cleanup  = 1;
    obj.mrf      = 2;
    obj.affreg   = 'mni';
    obj.reg      = [0 0.001 0.5 0.05 0.2]*0.1;
    obj.fwhm     = 1;
    obj.samp     = 4;
    obj.biasreg  = 0.001*(1/5);
    obj.biasfwhm = 60;

    tpmname   = fullfile(spm('dir'),'tpm','TPM.nii');
    obj.lkp   = 1:6;
    obj.tpm   = spm_load_priors8(tpmname);
    obj.image = V;

    % Initial affine registration.
    obj.Affine = spm_maff8(obj.image(1),6,obj.fwhm,obj.tpm,[],obj.affreg);

    % Run the actual segmentation
    res = spm_preproc8(obj);

    % Final iteration, so write out the required data.
    required        = false(max(obj.lkp),4);
    required(1:3,1) = true; % GM, WM, CSF
    spm_preproc_write8(res,required,false(1,2),false(1,2),obj.mrf,obj.cleanup,obj.bb,obj.vox);
    
    % Create mask
    pth   = fileparts(obj.image(1).fname);
    files = spm_select('FPList',pth,'^c.*\.nii$');
    V0    = spm_vol(files);
    K     = numel(V0);
    msk   = zeros(V0(1).dim,'single');
    for k=1:K
        msk = msk + V0(k).private.dat(:,:,:);
    end
    msk = msk>0;

    % Mask images
    N = numel(obj.image);
    for n=1:N
        img       = single(obj.image(n).private.dat(:,:,:));        
        img(~msk) = 0;      
        obj.image(n).private.dat(:,:,:) = img;
    end    
end
%==========================================================================

%==========================================================================
function subvol(V,bb,prefix)
% Extract a subvolume
% FORMAT subvol(V,bb,prefix)
% V      - SPM volume object
% bb     - bounding box
% prefix - file prefix (if empty -> overwrites)
%
% Example:
%     V = spm_vol(spm_select(1,'image'));
%     subvol(V,[32 64 ; 1 64 ; 1 48]');
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

    if nargin<3, prefix = ''; end

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
end
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
end
%==========================================================================

%==========================================================================
function y1 = affind(y0,M)
    y1 = zeros(size(y0),'single');
    for d=1:3,
        y1(:,:,:,d) = y0(:,:,:,1)*M(d,1) + y0(:,:,:,2)*M(d,2) + y0(:,:,:,3)*M(d,3) + M(d,4);
    end
end
%==========================================================================
