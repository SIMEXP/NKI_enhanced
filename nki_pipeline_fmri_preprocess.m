function [] = nki_pipeleine_fmri_preprocess(opt);
%%% NKI_enhanced preprocessing pipeline
% Script to run a preprocessing pipeline analysis on the HCP database.
%
% SYNTAX:
% []= NKI_PIPELEINE_FMRI_PREPROCESS(OPT);
%
% _________________________________________________________________________
% INPUTS:
%
% OPT
%   (structure, optional) with the following fields :
%
%   TASK
%       (string, default 'checkerBoard') type of tasks that would be preprocessed. 
%        Possibles tasks are: 'checkerBoard','breathHold','rest'.
%
%   EXP
%       (string, default 'TR1400') TR used for the experiment.
%       Possibles experiment : 'TR1400', 'TR645'
%
% _________________________________________________________________________
% Copyright (c) Pierre Bellec, Yassine Benhajali
% Research Centre of the Montreal Geriatric Institute
% & Department of Computer Science and Operations Research
% University of Montreal, Québec, Canada, 2010-2012
% Maintainer : pierre.bellec@criugm.qc.ca
% See licensing information in the code.
% Keywords : fMRI, FIR, clustering, BASC
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in
% all copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
% THE SOFTWARE.

clear all
%%%%%%%%%%%%%%%%%%%%%
%% Parameters
%%%%%%%%%%%%%%%%%%%%%
list_fields   = { 'task'        , 'exp'    };
list_defaults = { 'checkerBoard', 'TR1400' };
if nargin > 1
    opt = psom_struct_defaults(opt,list_fields,list_defaults);
else
    opt = psom_struct_defaults(struct(),list_fields,list_defaults);
end

task  = opt.task;
exp   = opt.exp;

%% Setting input/output files 
[status,cmdout] = system ('uname -n');
server          = strtrim(cmdout);
if strfind(server,'lg-1r') % This is guillimin
    root_path = '/gs/scratch/yassinebha/NKI_enhanced/';
    fprintf ('server: %s (Guillimin) \n ',server)
    my_user_name = getenv('USER');
elseif strfind(server,'ip05') % this is mammouth
    root_path = '/mnt/parallel_scratch_ms2_wipe_on_april_2015/pbellec/benhajal/NKI_enhanced/';
    fprintf ('server: %s (Mammouth) \n',server)
    my_user_name = getenv('USER');
else
    switch server
        case 'peuplier' % this is peuplier
        root_path = '/media/database4/nki_enhanced/';
        fprintf ('server: %s\n',server)
        my_user_name = getenv('USER');
        
        case 'noisetier' % this is noisetier
        root_path = '/media/database1/';
        fprintf ('server: %s\n',server)
        my_user_name = getenv('USER');
    end
end

% path_raw = '/media/database3/nki_enhanced/raw_mnc/';
path_raw = [root_path 'raw_mnc/';
path_preproc = '/home/bellecp1/database/nki_enhanced/fmri_preprocess/';

%% Grab the raw data
list_subject = dir(path_raw);
list_subject = {list_subject.name};
list_subject = list_subject(~ismember(list_subject,{'.','..'}));

for num_s = 1:length(list_subject)
    subject = list_subject{num_s};
    id = ['X' subject];
    files_in.(id).anat = [path_raw subject filesep 'anat' filesep 'mprage.mnc.gz'];
    files_in.(id).fmri.sess1.std2500 = [path_raw subject filesep 'session_1' filesep 'RfMRI_std_2500' filesep 'rest.mnc.gz'];
    files_in.(id).fmri.sess1.mx1400 = [path_raw subject filesep 'session_1' filesep 'RfMRI_mx_1400' filesep 'rest.mnc.gz'];
    files_in.(id).fmri.sess1.mx645 = [path_raw subject filesep 'session_1' filesep 'RfMRI_mx_645' filesep 'rest.mnc.gz'];    
    
    files_c = psom_files2cell(files_in.(id));
    for num_f = 1:length(files_c)
        if ~psom_exist(files_c{num_f})
            warning ('The file %s does not exist, I suppressed subject %s',files_c{num_f},subject);
            files_in = rmfield(files_in,id);
            break
        end        
    end
end
    
%warning: The file /media/database3/nki_enhanced/raw_mnc/0103714/session_1/RfMRI_std_2500/rest.mnc.gz does not exist, I suppressed subject 0103714
%warning: The file /media/database3/nki_enhanced/raw_mnc/0120538/session_1/RfMRI_std_2500/rest.mnc.gz does not exist, I suppressed subject 0120538
%warning: The file /media/database3/nki_enhanced/raw_mnc/0120652/session_1/RfMRI_std_2500/rest.mnc.gz does not exist, I suppressed subject 0120652
%warning: The file /media/database3/nki_enhanced/raw_mnc/0121498/anat/mprage.mnc.gz does not exist, I suppressed subject 0121498
%warning: The file /media/database3/nki_enhanced/raw_mnc/0123657/session_1/RfMRI_mx_1400/rest.mnc.gz does not exist, I suppressed subject 0123657
%warning: The file /media/database3/nki_enhanced/raw_mnc/0136018/session_1/RfMRI_std_2500/rest.mnc.gz does not exist, I suppressed subject 0136018
%warning: The file /media/database3/nki_enhanced/raw_mnc/0141473/session_1/RfMRI_std_2500/rest.mnc.gz does not exist, I suppressed subject 0141473
%warning: The file /media/database3/nki_enhanced/raw_mnc/0144344/anat/mprage.mnc.gz does not exist, I suppressed subject 0144344
%warning: The file /media/database3/nki_enhanced/raw_mnc/0144495/session_1/RfMRI_std_2500/rest.mnc.gz does not exist, I suppressed subject 0144495
%warning: The file /media/database4/nki_enhanced/raw_mnc/0182604/session_1/RfMRI_mx_645/rest.mnc.gz does not exist, I suppressed subject 0182604

%% Pipeline options  %%
%% General
opt.folder_out  = [root_path 'fmri_preprocess_' upper(task) '_' exp];    % Where to store the results
opt.size_output = 'quality_control';                             % The amount of outputs that are generated by the pipeline. 'all' will keep intermediate outputs, 'quality_control' will only keep the quality control outputs.

%% General
opt.folder_out  = path_preproc;    % Where to store the results
opt.size_output = 'quality_control';                             % The amount of outputs that are generated by the pipeline. 'all' will keep intermediate outputs, 'quality_control' will only keep the quality control outputs. 

%% Pipeline manager 
opt.psom.qsub_options = '-q qwork@ms -l nodes=1:m32G,walltime=05:00:00';

%% Slice timing correction (niak_brick_slice_timing)
opt.slice_timing.type_acquisition = 'interleaved ascending'; % Slice timing order (available options : 'sequential ascending', 'sequential descending', 'interleaved ascending', 'interleaved descending')
opt.slice_timing.type_scanner     = 'Bruker';                % Scanner manufacturer. Only the value 'Siemens' will actually have an impact
opt.slice_timing.delay_in_tr      = 0;                       % The delay in TR ("blank" time between two volumes)
opt.slice_timing.suppress_vol     = 0;                       % Number of dummy scans to suppress.
opt.slice_timing.flag_nu_correct  = 1;                       % Apply a correction for non-uniformities on the EPI volumes (1: on, 0: of). This is particularly important for 32-channels coil.
opt.slice_timing.arg_nu_correct   = '-distance 200';         % The distance between control points for non-uniformity correction (in mm, lower values can capture faster varying slow spatial drifts).
opt.slice_timing.flag_center      = 0;                       % Set the origin of the volume at the center of mass of a brain mask. This is useful only if the voxel-to-world transformation from the DICOM header has somehow been damaged. This needs to be assessed on the raw images.
opt.slice_timing.flag_skip        = true;                    % Skip the slice timing (0: don't skip, 1 : skip). Note that only the slice timing corretion portion is skipped, not all other effects such as FLAG_CENTER or FLAG_NU_CORRECT
 
% resampling in stereotaxic space
opt.resample_vol.interpolation = 'trilinear'; % The resampling scheme. The fastest and most robust method is trilinear. 
opt.resample_vol.voxel_size    = [3 3 3];     % The voxel size to use in the stereotaxic space
opt.resample_vol.flag_skip     = 0;           % Skip resampling (data will stay in native functional space after slice timing/motion correction) (0: don't skip, 1 : skip)

% Linear and non-linear fit of the anatomical image in the stereotaxic
% space (niak_brick_t1_preprocess)
opt.t1_preprocess.nu_correct.arg = '-distance 75'; % Parameter for non-uniformity correction. 200 is a suggested value for 1.5T images, 75 for 3T images. If you find that this stage did not work well, this parameter is usually critical to improve the results.

% Temporal filtering (niak_brick_time_filter)
opt.time_filter.hp = 0.01; % Cut-off frequency for high-pass filtering, or removal of low frequencies (in Hz). A cut-off of -Inf will result in no high-pass filtering.
opt.time_filter.lp = Inf;  % Cut-off frequency for low-pass filtering, or removal of high frequencies (in Hz). A cut-off of Inf will result in no low-pass filtering.

% Regression of confounds and scrubbing (niak_brick_regress_confounds)
opt.regress_confounds.flag_wm = true;            % Turn on/off the regression of the average white matter signal (true: apply / false : don't apply)
opt.regress_confounds.flag_vent = true;          % Turn on/off the regression of the average of the ventricles (true: apply / false : don't apply)
opt.regress_confounds.flag_motion_params = true; % Turn on/off the regression of the motion parameters (true: apply / false : don't apply)
opt.regress_confounds.flag_gsc = false;          % Turn on/off the regression of the PCA-based estimation of the global signal (true: apply / false : don't apply)
opt.regress_confounds.flag_scrubbing = true;     % Turn on/off the scrubbing of time frames with excessive motion (true: apply / false : don't apply)
opt.regress_confounds.thre_fd = 0.5;             % The threshold on frame displacement that is used to determine frames with excessive motion in the scrubbing procedure

% Correction of physiological noise (niak_pipeline_corsica)
opt.corsica.sica.nb_comp             = 60;    % Number of components estimated during the ICA. 20 is a minimal number, 60 was used in the validation of CORSICA.
opt.corsica.threshold                = 0.15;  % This threshold has been calibrated on a validation database as providing good sensitivity with excellent specificity.
opt.corsica.flag_skip                = 1;     % Skip CORSICA (0: don't skip, 1 : skip). Even if it is skipped, ICA results will be generated for quality-control purposes. The method is not currently considered to be stable enough for production unless it is manually supervised.

% Spatial smoothing (niak_brick_smooth_vol)
opt.smooth_vol.fwhm      = 6;  % Full-width at maximum (FWHM) of the Gaussian blurring kernel, in mm.
opt.smooth_vol.flag_skip = 0;  % Skip spatial smoothing (0: don't skip, 1 : skip)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Tune the parameters for specific subjects %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opt.tune(1).subject = 'X0109819';
opt.tune(1).param.slice_timing = opt.slice_timing;
opt.tune(1).param.slice_timing.flag_center = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run the fmri_preprocess pipeline  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[pipeline,opt] = niak_pipeline_fmri_preprocess(files_in,opt);
