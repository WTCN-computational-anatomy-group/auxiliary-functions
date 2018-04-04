lat = [3 3 3];
vs = [0.5 0.4 0.3];
% vs = [1 1 1];
v = zeros([lat 3]);
v(1,1,1,1) = 1;

verbose = false;
threshold = 1e-7;

%% ABSOLUTE
prm = [1 0 0 eps 0];

% sparse matrix
L = spm_sparse('precision', 'diffeo', lat, vs, prm);
m1 = reshape(L*v(:), [lat 3]);

% convolution
m2 = spm_diffeo('vel2mom', single(v), [vs prm]);

if sum((double(m1(:)) - double(m2(:))).^2)/numel(m1) > threshold
    fprintf('Absolute vel2mem FAILED\n')
    if verbose
        m1(:,:,:,1)
        m2(:,:,:,1)
    end
else
    fprintf('Absolute vel2mem OK\n')
end


ker1 = spm_sparse('kernel', 'diffeo', prm, lat, vs);
ker2 = spm_diffeo('kernel', lat, [vs prm]);

if sum((double(ker1(:)) - double(ker2(:))).^2)/numel(ker1) > threshold
    fprintf('Absolute kernel FAILED\n')
    if verbose
        ker1
        ker2
    end
else
    fprintf('Absolute kernel OK\n')
end

%% MEMBRANE
prm = [0 1 0 eps 0];

% sparse matrix
L = spm_sparse('precision', 'diffeo', lat, vs, prm);
m1 = reshape(L*v(:), [lat 3]);

% convolution
m2 = spm_diffeo('vel2mom', single(v), [vs prm]);

if sum((double(m1(:)) - double(m2(:))).^2)/numel(m1) > threshold
    fprintf('Membrane vel2mem FAILED\n')
    if verbose
        m1(:,:,:,1)
        m2(:,:,:,1)
    end
else
    fprintf('Membrane vel2mem OK\n')
end


ker1 = spm_sparse('kernel', 'diffeo', prm, lat, vs);
ker2 = spm_diffeo('kernel', lat, [vs prm]);

if sum((double(ker1(:)) - double(ker2(:))).^2)/numel(ker1) > threshold
    fprintf('Membrane kernel FAILED\n')
    if verbose
        ker1
        ker2
    end
else
    fprintf('Membrane kernel OK\n')
end

%% BENDING
prm = [0 0 1 eps 0];

% sparse matrix
L = spm_sparse('precision', 'diffeo', lat, vs, prm);
m1 = reshape(L*v(:), [lat 3]);

% convolution
m2 = spm_diffeo('vel2mom', single(v), [vs prm]);

if sum((double(m1(:)) - double(m2(:))).^2)/numel(m1) > threshold
    fprintf('Bending vel2mem FAILED\n')
    if verbose
        m1(:,:,:,1)
        m2(:,:,:,1)
    end
else
    fprintf('Bending vel2mem OK\n')
end

ker1 = spm_sparse('kernel', 'diffeo', prm, lat, vs);
ker2 = spm_diffeo('kernel', lat, [vs prm]);

if sum((double(ker1(:)) - double(ker2(:))).^2)/numel(ker1) > threshold
    fprintf('Bending kernel FAILED\n')
    if verbose
        ker1
        ker2
    end
else
    fprintf('Bending kernel OK\n')
end

%% LINEAR-ELASTIC SYMJAC
prm = [0 0 0 1 0];

% sparse matrix
L = spm_sparse('precision', 'diffeo', lat, vs, prm);
m1 = reshape(L*v(:), [lat 3]);

% convolution
m2 = spm_diffeo('vel2mom', single(v), [vs prm]);

if sum((double(m1(:)) - double(m2(:))).^2)/numel(m1) > threshold
    fprintf('Linear-Elastic-SymJac vel2mem FAILED\n')
    if verbose
        m1(:,:,:,1)
        m2(:,:,:,1)
    end
else
    fprintf('Linear-Elastic-SymJac vel2mem OK\n')
end

ker1 = spm_sparse('kernel', 'diffeo', prm, lat, vs);
ker2 = spm_diffeo('kernel', lat, [vs prm]);

if sum((double(ker1(:)) - double(ker2(:))).^2)/numel(ker1) > threshold
    fprintf('Linear-Elastic-SymJac kernel FAILED\n')
    if verbose
        ker1
        ker2
    end
else
    fprintf('Linear-Elastic-SymJac kernel OK\n')
end

%% LINEAR-ELASTIC DIV
prm = [0 0 0 0 1];

% sparse matrix
L = spm_sparse('precision', 'diffeo', lat, vs, prm);
m1 = reshape(L*v(:), [lat 3]);

% convolution
m2 = spm_diffeo('vel2mom', single(v), [vs prm]);

if sum((double(m1(:)) - double(m2(:))).^2)/numel(m1) > threshold
    fprintf('Linear-Elastic-Div vel2mem FAILED\n')
    if verbose
        m1(:,:,:,1)
        m2(:,:,:,1)
    end
else
    fprintf('Linear-Elastic-Div vel2mem OK\n')
end

ker1 = spm_sparse('kernel', 'diffeo', prm, lat, vs);
ker2 = spm_diffeo('kernel', lat, [vs prm]);

if sum((double(ker1(:)) - double(ker2(:))).^2)/numel(ker1) > threshold
    fprintf('Linear-Elastic-Div kernel FAILED\n')
    if verbose
        ker1
        ker2
    end
else
    fprintf('Linear-Elastic-Div kernel OK\n')
end
