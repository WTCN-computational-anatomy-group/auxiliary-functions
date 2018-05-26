function varargout = spm_json_manager(varargin)
%__________________________________________________________________________
% Collection of functions for reading and organising data.
%
% FORMAT [dat,dict] = spm_json_manager('init_dat',dir_population,load_dat,dat)
% FORMAT model      = spm_json_manager('init_model',input,load_dat,dat)
% FORMAT spm_json_manager('modify_json_field',pth_json,field,val)
% FORMAT spm_json_manager('modify_pth_in_population',dir_population,field,npth)
% FORMAT spm_json_manager('make_pth_relative',input)
%
% FORMAT help spm_json_manager>function
% Returns the help file of the selected function.
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging
if nargin == 0
    help spm_json_manager
    error('Not enough argument. Type ''help spm_json_manager'' for help.');
end
id = varargin{1};
varargin = varargin(2:end);
switch lower(id) 
    case 'init_dat'
        [varargout{1:nargout}] = init_dat(varargin{:});   
    case 'init_model'
        [varargout{1:nargout}] = init_model(varargin{:});           
    case 'modify_json_field'
        [varargout{1:nargout}] = modify_json_field(varargin{:});             
    case 'modify_pth_in_population'
        [varargout{1:nargout}] = modify_pth_in_population(varargin{:});    
    case 'make_pth_relative'
        [varargout{1:nargout}] = make_pth_relative(varargin{:});              
    otherwise
        help spm_json_manager
        error('Unknown function %s. Type ''help spm_json_manager'' for help.', id)
end
%==========================================================================


%==========================================================================
function [dat,dict] = init_dat(input,load_dat,dat)
%__________________________________________________________________________
% Initialise a dat object from one or several JSON files.
% These JSON files can either be explicitely provided, or searched for in a
% directory:
%
% FORMAT [dat,dict] = spm_json_manager('init_dat',input,load_dat)
%
% input - Path to a JSON file or ot a directory with all the JSON files.
%         The input can also be a list (= cell) of paths.
% load_dat - Load (if exists) from population specific dat.mat files
%
% dat   - A hierarchical object (cell of structs) representing all input
%         files and metadata
% dict  - A dictionary mapping subject IDs to indices in dat.
%
%
% It is also possible to addup to an existing dat object:
%
% FORMAT [dat,dict] = spm_json_manager('init_dat', ..., dat)
%
% dat - An already initialised model structure.
%--------------------------------------------------------------------------
% JSON SYNTAX
% -----------
%
% JSON files should contain either a dictionary or a list of dictionaries.
% Each dictionary should follow the following syntax:
%
% mandatory
% ---------
% 'name':       Subject name
% 
% optional
% --------
% 'population': Populatipn name ('Healty'/'Lesion'/...). 
%               If provided, the unique subject id is <population>_<name>.
% 'modality':   Modality name (for imaging modality) ('CT'/'MRI'/...)
% + 'channel':  Channel name ('T1'/'T2'/...)
% + 'pth':      Path to image file (absolute or relative w.r.t. JSON file)
% 'rater':      Rater name (for manual segmentation)
% + 'pth':      Path to image file (absolute or relative w.r.t. JSON file)
% 'tissue':     Tissue name (for automated segmentation_map) ('GM'/'WM'/...)
% + 'type':     Segmentation type ('c'/'wc'/'rc'/...)
% + 'pth':      Path to image file (absolute or relative w.r.t. JSON file)
%
% Any other field name can be used to store additional metadata (age,
% weight, ...). They will be stored at the root of each subject.
%
% Be aware that all keys (name, population, modality, channel, rater, ...)
% are case sensitive! Consequently, 'T1' is not the same as 't1'.
%--------------------------------------------------------------------------
% OUTPUT OBJECT
% -------------
%
% The output object has the following form. Here, all fields are shown.
% However, in reality, fields only appear if the corresponding data was
% found in the JSON files.
%
% dat{s}.modality_map                         [map modality names to indices]
% dat{s}.modality{m}.name                                  [modality names]
% dat{s}.modality{m}.nii                      [nifti - single channel case]
% dat{s}.modality{m}.channel_map                [map channel names to indices]
% dat{s}.modality{m}.channel{c}.name                        [channel names]
% dat{s}.modality{m}.channel{c}.nii         [nifti - multiple channel case]
% dat{s}.rater_map                                [map rater names to indices]
% dat{s}.label{r}.name                                        [rater names]
% dat{s}.label{r}.nii                                               [nifti]
% dat{s}.segmentation_map                  [map segmentation types to indices]
% dat{s}.segmentation{t}.name                           [segmentation type]
% dat{s}.segmentation{t}.class{k}.name                        [tissue name]
% dat{s}.segmentation{t}.class{k}.nii                               [nifti]
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

if nargin<2
    load_dat = false;
end
if nargin<3
    dat = {};
end

% -------------------------------------------------------------------------
% Create a dictionary which will map subject names to dat indices
dict = containers.Map;
if ~isempty(dat)
    for s=1:numel(dat)
        dict(dat{s}.name) = s;
    end
end
pre_count = dict.Count;

% -------------------------------------------------------------------------
% Get all input json files
if ~iscell(input), input = {input}; end
json_files = [];
for i=1:numel(input)
    if numel(input{i}) >= 5 && strcmpi(input{i}(end-4:end), '.json')
        json_files = [json_files ; dir(input{i})];
    else
        json_files = [json_files ; rdir(input{i}, '*.json')];
    end
end
% clear input
J = numel(json_files);

% -------------------------------------------------------------------------
% Display the number of files read
base10 = floor(log10(J)) + 1;
str    = sprintf(['Initialising dat | %' num2str(base10) 'd of %' num2str(base10) 'd files read.'],0,J);
fprintf(1, ['%-' num2str(2*base10 + 50) 's'], str);
tic;

% -------------------------------------------------------------------------
% Used as a flag for saving the dat object to check if more than one population
% exists
odir_json   = '';
first_iter  = true;
loading_dat = false;

% -------------------------------------------------------------------------
% Loop over files
for j=1:J
    
    % ---------------------------------------------------------------------
    % Display the number of files read
    if ~mod(j,10)
        fprintf(1, repmat('\b',1,2*base10 + 50));
        str = sprintf(['Initialising dat | %' num2str(base10) 'd of %' num2str(base10) 'd files read.'],j,J);
        fprintf(1, ['%-' num2str(2*base10 + 50) 's'], str);
    end
    
    % ---------------------------------------------------------------------
    % Get path to JSON file
    pth_json = fullfile(json_files(j).folder,json_files(j).name);    
    
    if load_dat        
        % Loading from file (if exists)
        pth_dat = fullfile(json_files(j).folder,'dat.mat');
        if exist(pth_dat,'file')==2
            loading_dat = true;                        
            
            dir_json = json_files(j).folder;
            if strcmp(dir_json,odir_json)
                % This just increments the for-loop index
                continue
            end
            odir_json = dir_json;
            
            dat1 = load(pth_dat);
            dat1 = dat1.dat1;
            J1   = numel(dat1);                      
            
            for j1=1:J1
                dat{end + 1} = dat1{j1};
            end                              
            
            first_iter = true;
             
            continue
        end
    end
    
    if loading_dat
        error('Currently does not support both loading from dat.mat and reading from JSONs')
    end
    
    list_metadata = spm_jsonread(pth_json);  
    if ~iscell(list_metadata), list_metadata = {list_metadata}; end
    I = numel(list_metadata);
            
    % ---------------------------------------------------------------------
    % Loop over elements in the file
    % > In order to store multiple elements in a single JSON file, I
    %   support lists of dictionaries. Each element of the list is treated 
    %   as one json object.
    for i=1:I
        metadata = list_metadata{i};
        metadata = check_metadata_dat(metadata);
          
        % Get poulation and subject name
        name       = metadata.name;
        population = metadata.population;
                  
        if ~first_iter
            if ~strcmp(opopulation,population)
                % New population -> save dat object for previous population
                pth_dat = fullfile(json_files(j - 1).folder,'dat.mat');
                if exist(pth_dat,'file')==2
                    delete(pth_dat);
                end
                
                dat1 = get_population(dat,opopulation);                
                save(pth_dat,'dat1')
            end
        end
        opopulation = population;     
        
        % Create dictionary key
        if ~isempty(population),    key = [population '_' name];
        else,                       key = name;
        end

        if ~dict.isKey(key)
            % Subject not in dictionary -> add subject to dictionary
            dict(key)                 = dict.Count + 1;         
            dat{dict(key)}.name       = name;
            dat{dict(key)}.population = population;
        end
        s = dict(key);

        % -----------------------------------------------------------------
        % Process imaging data
        % -----------------------------------------------------------------
        if ~isempty(metadata.modality)
            
            modality = metadata.modality;
            channel  = metadata.channel;
            hospital = metadata.hospital;
            
            % -------------------------------------------------------------
            % Get path to image and read
            Nii = read_nifti(metadata.pth, json_files(j).folder, json_files(j).name, 'ro');
            if isempty(Nii)
                continue
            end
            
            % -------------------------------------------------------------
            % Create modality
            if ~isfield(dat{s},'modality')
                % No image data exists -> create image data fields
                dat{s}.modality_map         = containers.Map;
                dat{s}.modality             = {};
            end
            if ~dat{s}.modality_map.isKey(modality)
                m                           = numel(dat{s}.modality) + 1;
                dat{s}.modality{m}.name     = modality;
                dat{s}.modality_map(modality) = m;
                if ~isempty(hospital)
                    dat{s}.modality{m}.hospital = hospital;
                end
            end
            % -------------------------------------------------------------
            % Add modality (single channel)  
            m = dat{s}.modality_map(modality);
            if isempty(channel)
                if isfield(dat{s}.modality{m}, 'nii')
                    N = numel(dat{s}.modality{m}.nii);
                else
                    dat{s}.modality{m}.nii = nifti;
                    N = 0;
                end
                dat{s}.modality{m}.nii(N + 1)      = Nii;
                dat{s}.modality{m}.json(N + 1).pth = pth_json;
            else
                % ---------------------------------------------------------
                % Create channel
                if ~isfield(dat{s}.modality{m}, 'channel')
                    dat{s}.modality{m}.channel  = {};
                    dat{s}.modality{m}.channel_map = containers.Map;
                end
                if ~dat{s}.modality{m}.channel_map.isKey(channel)
                    c                                    = numel(dat{s}.modality{m}.channel) + 1;
                    dat{s}.modality{m}.channel{c}.name   = channel;
                    dat{s}.modality{m}.channel_map(channel) = c;
                end
                % ---------------------------------------------------------
                % Add channel
                c = dat{s}.modality{m}.channel_map(channel);
                if isfield(dat{s}.modality{m}.channel{c}, 'nii')
                    N = numel(dat{s}.modality{m}.channel{c}.nii);
                else
                    dat{s}.modality{m}.channel{c}.nii = nifti;
                    N = 0;
                end
                dat{s}.modality{m}.channel{c}.nii(N + 1)      = Nii;
                dat{s}.modality{m}.channel{c}.json(N + 1).pth = pth_json;
            end % < Channel
        end % < Modality

        % -----------------------------------------------------------------
        % Process label data
        % -----------------------------------------------------------------
        if ~isempty(metadata.rater)
            
            rater     = metadata.rater;
            protocol = metadata.protocol;
            
            % -------------------------------------------------------------
            % Get path to image and read
            Nii = read_nifti(metadata.pth, json_files(j).folder, json_files(j).name, 'ro');
            if isempty(Nii)
                continue
            end
            
            % -------------------------------------------------------------
            % Create label
            if ~isfield(dat{s},'label')
                % No image data exists -> create image data fields
                dat{s}.rater_map   = containers.Map;
                dat{s}.label    = {};
            end
            if ~dat{s}.rater_map.isKey(rater)
                r                       = numel(dat{s}.label) + 1;
                dat{s}.label{r}.name    = rater;
                dat{s}.rater_map(rater) = r;
            end
            % -------------------------------------------------------------
            % Add label
            r = dat{s}.rater_map(rater);
            if isfield(dat{s}.label{r}, 'nii')
                N = numel(dat{s}.label{r}.nii);
            else
                dat{s}.label{r}.nii = nifti;
                N = 0;
            end
            dat{s}.label{r}.nii(N+1)      = Nii;
            dat{s}.label{r}.json(N+1).pth = pth_json;
            dat{s}.label{r}.protocol     = protocol;
        end

        
        % -----------------------------------------------------------------
        % Process segmentation data
        % -----------------------------------------------------------------
        if ~isempty(metadata.tissue)
                
            type   = metadata.type;        
            tissue = metadata.tissue;   
            
            % -------------------------------------------------------------
            % Get path to image and read
            Nii = read_nifti(metadata.pth, json_files(j).folder, json_files(j).name, 'ro');
            if isempty(Nii)
                continue
            end

            % -------------------------------------------------------------
            % Create segmentation
            if ~isfield(dat{s},'segmentation')
                % No image data exists -> create image data fields
                dat{s}.segmentation_map = containers.Map;
                dat{s}.segmentation  = {};
            end
            if ~dat{s}.segmentation_map.isKey(type)
                t                             = numel(dat{s}.segmentation) + 1;
                dat{s}.segmentation{t}.name   = type;
                dat{s}.segmentation_map(type) = t;
            end
            % -------------------------------------------------------------
            % Create class
            t = dat{s}.segmentation_map(type);
            if ~isfield(dat{s}.segmentation{t}, 'class')
                dat{s}.segmentation{t}.class   = {};
                dat{s}.segmentaload_dattion{t}.class_map = containers.Map;
            end
            if ~dat{s}.segmentation{t}.class_map.isKey(tissue)
                c                                      = numel(dat{s}.segmentation{t}.class) + 1;
                dat{s}.segmentation{t}.class{c}.name   = tissue;
                dat{s}.segmentation{t}.class_map(tissue) = c;
            end
            % -------------------------------------------------------------
            % Add class
            if isfield(dat{s}.segmentation{t}.class{c}, 'nii')
                N = numel(dat{s}.segmentation{t}.class{c}.nii);
            else
                dat{s}.segmentation{t}.class{c}.nii = nifti;
                N = 0;
            end
            dat{s}.segmentation{t}.class{c}.nii(N+1)      = Nii;
            dat{s}.segmentation{t}.class{c}.json(N+1).pth = pth_json;
            
        end
    
        % -----------------------------------------------------------------
        % Append other meta data fields (if there are any)
        % -----------------------------------------------------------------
        fn = fieldnames(metadata);
        protected_fields = {'pth','modality','name','rater','channel','tissue','protocol','hospital','type'};
        for k=1:numel(fn)
            field_name = fn{k};

            if any(strcmpi(field_name, protected_fields))
                continue
            end

            dat{s}.(fn{k}) = metadata.(field_name);
        end            

        % Make sure fields are ordered alphabetically
        dat{s} = orderfields(dat{s});
        
    end % < Loop over elements (I)
    
    if ~first_iter || J==1
        if strcmp(opopulation,population) && j==J
            % Only one, or last population -> save dat object
            pth_dat = fullfile(json_files(j).folder,'dat.mat');
            if exist(pth_dat,'file')==2
                delete(pth_dat);
            end

            dat1 = get_population(dat,opopulation);                
            save(pth_dat,'dat1')
        end
    else
        first_iter = false;
    end
end % < Loop over files (J)

% -------------------------------------------------------------------------
% Display number of files read
fprintf(1, repmat('\b',1,2*base10 + 50));
str = sprintf(['Initialising dat | %' num2str(base10) 'd of %' num2str(base10) 'd files read.'],j,J);
fprintf(1, ['%-' num2str(2*base10 + 50) 's'], str);
fprintf('\n');

reswhos = whos('dat');
siz = reswhos.bytes;
unit = 'B';
if siz > 1024
    siz = ceil(siz/1024);
    unit = 'KB';
end
if siz > 1024
    siz = ceil(siz/1024);
    unit = 'MB';
end
if siz > 1024
    siz = ceil(siz/1024);
    unit = 'GB';
end

if numel(input) == 1
    fprintf('Initialising dat | Loaded %i subjects from %s in %0.1f seconds (%d%s).\n',numel(dat),input{1},toc,siz,unit);
else
    fprintf('Initialising dat | Loaded %i subjects in %0.1f seconds (%d%s).\n',numel(dat),toc,siz,unit);
end
%==========================================================================


%==========================================================================
function model = init_model(input, varargin)
%__________________________________________________________________________
% Initialise a model object from one or several JSON files.
% These JSON files can either be explicitely provided, or searched for in a
% directory.
%
% FORMAT model = spm_json_manager('init_model',input)
% input - Path to a JSON file or ot a directory with all the JSON files.
%         The input can also be a list (= cell) of paths.
%
% It is also possible to provide an already initialised 'dat' structure
% from which some information can be extracted (modality_map, hospital_map, ...)
%
% FORMAT model = spm_json_manager('init_model', ..., dat)
% dat - An already initialised dat structure from which we can get info.
%
% Finally, it is possible to addup to an existing model object. An input
% structure is detected as a model object if it has only one element
% (whereas dat objects have several elements).
%
% FORMAT model = spm_json_manager('init_model', ..., model)
% model - An already initialised model structure.
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

model = struct;
for i=1:numel(varargin)
    if isstruct(varargin{i}) && numel(varargin{i}) == 1
        model = varargin{i};
        varargin = [varargin{1:i-1} varargin{i+1:end}];
    end
end
if numel(varargin) == 1
    dat = varargin{1};
elseif numel(varargin) > 1
    warning('[spm_json_manager::init_model] Too many arguments. Don''t know what to do.')
end

% -------------------------------------------------------------------------
% Dictionary for modality_map / hospital_map
% -------------------------------------------------------------------------
if ~isfield(model, 'modality_map'), model.modality_map = containers.Map; end
if ~isfield(model, 'hospital_map'),  model.hospital_map  = containers.Map; end

% -------------------------------------------------------------------------
% Parse JSON files
% -------------------------------------------------------------------------

% JSON file can be used to provide already known elements of the model.
% It can be an existing template, an existing shape model (subspace,
% residual precision, latent precision...) or a GMM.

% -------------------------------------------------------------------------
% Get all input json files
if ~iscell(input), input = {input}; end
json_files = [];
for i=1:numel(input)
    if numel(input{i}) >= 5 && strcmpi(input{i}(end-4:end), '.json')
        json_files = [json_files ; dir(input{i})];
    else
        json_files = [json_files ; rdir(input{i}, '*.json')];
    end
end
clear input
J = numel(json_files);

% -------------------------------------------------------------------------
% Loop over files
for j=1:J
    pth_json      = fullfile(json_files(j).folder,json_files(j).name);
    list_metadata = spm_jsonread(pth_json);  
    if ~iscell(list_metadata), list_metadata = {list_metadata}; end
    I = numel(list_metadata);
    
    % ---------------------------------------------------------------------
    % Loop over elements in the file
    % > In order to store multiple elements in a single JSON file, I
    %   support lists of dictionaries. Each element of the list is treated 
    %   as one json object.
    for i=1:I
        metadata = list_metadata{i};
        metadata = check_metadata_model(metadata);
    
        switch lower(metadata.type)
            % -------------------------------------------------------------
            % Shape model
            case 'shape'
                model.shape = struct;
                if isfield(metadata, 'prm')
                    model.shape.prm = metadata.prm;
                end
                if isfield(metadata, 'boundary')
                    model.shape.boundary = metadata.boundary;
                end
                if isfield(metadata, 'subspace')
                    model.shape.subspace.nii = read_nifti(metadata.subspace, json_files(j).folder, json_files(j).name, 'rw');
                end
                if isfield(metadata, 'latent')
                    if ischar(metadata.latent)
                        model.shape.latent.A = read_mat(metadata.latent, json_files(j).folder, json_files(j).name);
                    else
                        model.shape.latent.A = metadata.latent;
                    end
                end
                if isfield(metadata, 'residual')
                    if ischar(metadata.residual)
                        model.shape.residual.lambda = read_mat(metadata.residual, json_files(j).folder, json_files(j).name);
                    else
                        model.shape.residual.lambda = metadata.residual;
                    end
                end
            % -------------------------------------------------------------
            % Template
            case 'template'
                model.template = struct;
                if isfield(metadata, 'prm')
                    model.template.prm = metadata.prm;
                end
                if isfield(metadata, 'boundary')
                    model.template.boundary = metadata.boundary;
                end
                if isfield(metadata, 'log')
                    model.template.log.nii = read_nifti(metadata.log, json_files(j).folder, json_files(j).name, 'rw');
                end
                if isfield(metadata, 'derivatives')
                    model.template.derivatives.nii = read_nifti(metadata.derivatives, json_files(j).folder, json_files(j).name, 'rw');
                end
            % -------------------------------------------------------------
            % Intensity priors
            case 'gmm'
                if ~isfield(metadata, 'gmm')
                    error('GMM json files must have a ''gmm'' field.');
                end
                if ~isfield(metadata, 'channel_map')
                    error('GMM json files must have a ''channel_map'' field.');
                end
                if ~isfield(model, 'modality'), model.modality = {}; end
                if ~model.modality_map.isKey(metadata.modality)
                    model.modality{end+1}.name          = metadata.modality;
                    model.modality_map(metadata.modality) = numel(model.modality);
                    model.modality{end}.channel_map        = metadata.channel_map;
                end
                m = model.modality_map(metadata.modality);
                if isfield(metadata, 'hospital')
                    if ~isfield(model.modality{m}, 'hospital')
                        model.modality{m}.hospital = {};
                    end
                    if ~model.hospital_map.isKey(metadata.hospital)
                        model.modality{m}.hospital{end+1}.name = metadata.hospital;
                        model.hospital_map(metadata.hospital)     = numel(model.modality{m}.hospital);
                    end
                    h = model.hospital_map(metadata.hospital);
                    model.modality{m}.hospital{h}.gmm = metadata.gmm;
                else
                    model.modality{m}.gmm = metadata.gmm;
                end
        end
    end
end

% -------------------------------------------------------------------------
% Extract info from dat structure
% -------------------------------------------------------------------------

% By looking to the dat structure, we can detect how many modality_map (and
% thus GMM) exist, as well as how many channel_map per GMM. It is important to
% look at the entire dat structure, because some individuals may have
% missing modality_map or missing channel_map.
% A dictionary is used to map channel names with indices in the GMM.

model.dim = 0;
if ~isempty(dat)
    for s=1:numel(dat)
        if isfield(dat{s}, 'modality')
            for m=1:numel(dat{s}.modality)
                if ~isfield(model, 'modality'), model.modality = {}; end
                % Create modality to store future GMM
                name = dat{s}.modality{m}.name;
                if ~model.modality_map.isKey(name)
                    model.modality{end+1}.name   = name;
                    model.modality_map(name)       = numel(model.modality);
                    model.modality{end}.channel_map = containers.Map;
                    model.modality{end}.gmm      =  {};
                end
                mm = model.modality_map(name);
                % Register channel names to map with channel number in GMM.
                if isfield(dat{s}.modality{m}, 'channel')
                    for c=1:numel(dat{s}.modality{m}.channel)
                        name = dat{s}.modality{m}.channel{c}.name;
                        if ~model.modality{mm}.channel_map.isKey(name)
                            model.modality{mm}.channel_map(name) = model.modality{mm}.channel_map.Count + 1;
                        end
                        dim = 3 - (size(dat{s}.modality{m}.channel{c}.nii, 3) == 1);
                        if model.dim && model.dim ~= dim
                            error('Input data with different dimensions (2D/3D)');
                        elseif ~model.dim
                            model.dim = dim;
                        end
                    end
                else
                    model.modality{mm}.channel_map(name) = 1;
                    dim = 3 - (size(dat{s}.modality{m}.nii, 3) == 1);
                    if model.dim && model.dim ~= dim
                        error('Input data with different dimensions (2D/3D)');
                    elseif ~model.dim
                        model.dim = dim;
                    end
                end
            end
        end
        if isfield(dat{s}, 'label')
            for r=1:numel(dat{s}.label)
                if ~isfield(model, 'rater_map'),     model.rater_map     = containers.Map; end
                if ~isfield(model, 'protocol_map'), model.protocol_map = containers.Map; end
                if ~isfield(model, 'rater'),      model.rater      = {}; end
                if ~isfield(model, 'protocol'),  model.protocol  = {}; end
                % Create rater (to store ???)
                name      = dat{s}.label{r}.name;
                protocol = dat{s}.label{r}.protocol;
                if ~model.rater_map.isKey(name)
                    model.rater{end+1}.name = name;
                    model.rater_map(name)      = numel(model.rater);
                end
                if ~model.protocol_map.isKey(protocol)
                    model.protocol{end+1}.name = protocol;
                    model.protocol_map(protocol) = numel(model.protocol);
                end
            end
        end
    end
end

%==========================================================================


%==========================================================================
function modify_json_field(pth_json,field,val)
% FORMAT spm_json_manager('modify_json_field',pth_json,field,val)
%
% pth_json - Path to JSON file.
% field    - Field to change.
% val      - Value to change to.
%
% Modifies a field in a JSON file.
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging
a         = spm_jsonread(pth_json);
a.(field) = val;
a         = orderfields(a);
spm_jsonwrite(pth_json,a);
%==========================================================================


%==========================================================================
function modify_pth_in_population(dir_population,field,npth)
% FORMAT spm_json_manager('modify_pth_in_population',dir_population,field,new_pth)
%
% dir_population - Path to a directory containing JSON files. 
%                  Each JSON file holds subject-specific meta data.
% field          - Field to change.
% new_pth        - New path
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

% Get all JSON files
json_files = dir(fullfile(dir_population,'*.json'));
J          = numel(json_files);

for s=1:J
    pth_json    = fullfile(dir_population,json_files(s).name);
        
    a           = spm_jsonread(pth_json);
    opth        = a.(field);
    ix1         = opth(1);
    [~,nam,ext] = fileparts(opth);
    if ~strcmp(ix1,filesep)
        ix1 = '';
    end
    nval        = fullfile([ix1 npth],[nam ext]);
    a.(field)   = nval;
    a           = orderfields(a);
    spm_jsonwrite(pth_json,a);
end
%==========================================================================

%==========================================================================
function make_pth_relative(input)
% FORMAT spm_json_manager('make_pth_relative',input)
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

% -------------------------------------------------------------------------
% Parse JSON files
% -------------------------------------------------------------------------

% JSON file can be used to provide already known elements of the model.
% It can be an existing template, an existing shape model (subspace,
% residual precision, latent precision...) or a GMM.

% -------------------------------------------------------------------------
% Get all input json files
if ~iscell(input), input = {input}; end
json_files = [];
for i=1:numel(input)
    if numel(input{i}) >= 5 && strcmpi(input{i}(end-4:end), '.json')
        json_files = [json_files ; dir(input{i})];
    else
        json_files = [json_files ; rdir(input{i}, '*.json')];
    end
end
clear input
J = numel(json_files);

% -------------------------------------------------------------------------
% Loop over files
for j=1:J
    pth_json      = fullfile(json_files(j).folder,json_files(j).name);
    list_metadata = spm_jsonread(pth_json); 
    is_list = iscell(list_metadata);
    if ~iscell(list_metadata)
        list_metadata = {list_metadata};
    end
    I = numel(list_metadata);
    
    fprintf('%s\n', pth_json);
    
    % ---------------------------------------------------------------------
    % Loop over elements in the file
    % > In order to store multiple elements in a single JSON file, I
    %   support lists of dictionaries. Each element of the list is treated 
    %   as one json object.
    for i=1:I
        metadata = list_metadata{i};
        
        if isfield(metadata, 'pth')
            [~, fname, ext] = fileparts(metadata.pth);
            new_pth = fullfile(json_files(j).folder, [fname ext]);
            if exist(new_pth, 'file')
                metadata.pth = [fname ext];
            end
        end
        
        list_metadata{i} = metadata;
    end
    
    if ~is_list
        list_metadata = list_metadata{1};
    end
    spm_jsonwrite(pth_json, list_metadata);
end
%==========================================================================



%==========================================================================
% HELPER FUNCTIONS
%==========================================================================


%==========================================================================
function files = rdir(pop_dir, pattern)
% FORMAT files = rdir(folder, pattern)
% folder  - Root folder where to start the search
% pattern - Pattern for file names (ex: *.json)
%
% Recursive file search.

if ~is_absolute(pop_dir)
    pop_dir = fullfile(pwd, pop_dir);
end
if isdir(pop_dir)
    files           = dir(fullfile(pop_dir, pattern));
    [files.folder]  = deal(pop_dir);
    dirs            = dir(pop_dir);
    [dirs.folder]   = deal(pop_dir);
    is_dir          = [dirs.isdir];
    dirs            = dirs(is_dir);
    for i=1:numel(dirs)
        if ~any(strcmpi(dirs(i).name, {'.','..'}))
            files = [files ; rdir(fullfile(dirs(i).folder, dirs(i).name), pattern)];
        end
    end
end
%==========================================================================


%==========================================================================
function test = is_absolute(path)
% FORMAT test = is_absolute(path)
% path  - A path
%
% Check of a path is absolute or relative.

if isempty(path)
    test = false;
elseif path(1) == filesep
    test = true;
elseif ~ispc
    test = false;
else
    res  = regexp(path, '^[A-Za-z]:\\', 'once');
    test = ~isempty(res);
end
%==========================================================================


%==========================================================================
function path = check_path(path, folder, json, possible_ext)

% If no path: try to find an image with the same name as the JSON file
if isempty(path)
    [~, fname, ~] = fileparts(json);
    for e=1:numel(possible_ext)
        path = fullfile(folder, [fname possible_ext{e}]);
        if exist(path, 'file')
            break
        end
        path = '';
    end
end

% If still no path: failure
if isempty(path)
    warning('[spm_json_manager] Missing path. Could not find the corresponding file.')
    return
end

% If not absolute: assume relative path w.r.t. JSON file location
if ~is_absolute(path)
    path = fullfile(folder, path);
end

% If file does not exist: failure
if ~exist(path, 'file')
    warning('[spm_json_manager] Missing path. Could not find the corresponding file.')
    path = '';
    return
end
%==========================================================================


%==========================================================================
function nii = read_nifti(path, folder, json, permission)

if nargin < 4
    permission = 'ro'; % Inputs should be read-only (for security)
end

path = check_path(path, folder, json, {'.nii', '.nii.gz', '.img', '.img.gz'});
if isempty(path)
    nii = [];
    return
end

% Read nifti file
nii                = nifti(path);
nii.dat.permission = permission;
%==========================================================================


%==========================================================================
function obj = read_mat(path, folder, json)

path = check_path(path, folder, json, '.mat');
if isempty(path)
    obj = [];
    return
end%==========================================================================

% Read nifti file
obj = load(path);
%==========================================================================


%==========================================================================
function metadata = check_metadata_model(metadata)
% FORMAT metadata = check_metadata_model(metadata)
%
% Check metadata in the model case.

if ~isfield(metadata,'type')
    error('[spm_json_manager] field ''type'' is mandatory.')        
end
if ~isfield(metadata,'modality')
    metadata.modality = '';
end
%==========================================================================


%==========================================================================
function metadata = check_metadata_dat(metadata)
% FORMAT metadata = check_metadata_dat(metadata)
%
% Check metadata in the dat case.

if ~isfield(metadata,'name')
    error('[spm_json_manager] Field ''name'' is mandatory.')        
end
if ~isfield(metadata,'population')
    metadata.population = '';       
end
if ~isfield(metadata,'modality')
    metadata.modality = '';
end
if ~isfield(metadata,'channel')
    metadata.channel = '';   
end
if ~isfield(metadata,'rater')
    metadata.rater = '';   
end
if ~isfield(metadata,'protocol')
    metadata.protocol = 'unknown';   
end
if ~isfield(metadata,'hospital')
    metadata.hospital = '';   
end
if ~isfield(metadata,'tissue')
    metadata.tissue = '';   
end
if ~isfield(metadata,'pth')
    metadata.pth = '';
else
    metadata.pth = strtrim(metadata.pth);
end
%==========================================================================

%==========================================================================
function ndat = get_population(dat,population)
ndat = {};
for s=1:numel(dat)
    if strcmpi(dat{s}.population,population)
        ndat{end + 1} = dat{s};
    end
end
%==========================================================================