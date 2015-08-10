
% Script to run a STABILITY_FIR pipeline analysis on the NKI_enhanced database.
%
% Copyright (c) Pierre Bellec, Yassine Benhajali
% Research Centre of the Montreal Geriatric Institute
% & Department of Computer Science and Operations Research
% University of Montreal, Qubec, Canada, 2010-2014
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
% _________________________________________________________________________
%


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
        root_path = '/media/database8/NKI_enhanced/';
        fprintf ('server: %s\n',server)
        my_user_name = getenv('USER');
        
        case 'noisetier' % this is noisetier
        root_path = '/media/yassinebha/database2/nki_enhanced/';
        fprintf ('server: %s\n',server)
        my_user_name = getenv('USER');
    end
end

%for niak in docker 
%root_path = '/media/yassinebha/database2/nki_enhanced/';


%% create the csv model files
fmri_path = [root_path 'fmri_preprocess_ALL_task_noscrub/'];
mkdir(fmri_path,'onset');
path_folder = [ fmri_path 'onset/'];


opt_csv_read.separator= sprintf(',');
data.covariates_intrarun_names = {'times','duration'};
data.covariates_intrarun_cond  = {'baseline','breathhold',};
data.covariates_intrarun_values(1,1) = 0;
data.covariates_intrarun_values(1,2) = 10;
data.covariates_intrarun_values(2,1) = 0;
data.covariates_intrarun_values(2,2) = 250;

opt_csv_write.labels_y = data.covariates_intrarun_names;
opt_csv_write.labels_x = data.covariates_intrarun_cond;
opt_csv_write.precision = 2;
data.dir_output         = path_folder;
data.name_csv_intrarun  = 'nki_model_intrarun_breathhold_qc_timeseries';

niak_write_csv(strcat(data.dir_output,data.name_csv_intrarun,'.csv'),data.covariates_intrarun_values,opt_csv_write);



%%%%%%%%%%%%%%%%%%%%
%% Grabbing the results from the NIAK fMRI preprocessing pipeline
%%%%%%%%%%%%%%%%%%%%%
opt_g.min_nb_vol = 1;     % The minimum number of volumes for an fMRI dataset to be included. This option is useful when scrubbing is used, and the resulting time series may be too short.
opt_g.min_xcorr_func = 0.5; % The minimum xcorr score for an fMRI dataset to be included. This metric is a tool for quality control which assess the quality of non-linear coregistration of functional images in stereotaxic space. Manual inspection of the values during QC is necessary to properly set this threshold.
opt_g.min_xcorr_anat = 0.5; % The minimum xcorr score for an fMRI dataset to be included. This metric is a tool for quality control which assess the quality of non-linear coregistration of the anatomical image in stereotaxic space. Manual inspection of the values during QC is necessary to properly set this threshold.
opt_g.type_files = 'fir'; % Specify to the grabber to prepare the files for the STABILITY_FIR pipeline

%%Temporary grabber for debugging
liste_exclude = dir ([fmri_path 'anat']);
liste_exclude = liste_exclude(43:end -2);
liste_exclude = {liste_exclude.name};
opt_g.exclude_subject = liste_exclude;





opt_g.filter.run = {'breathHold1400'};
files_in = niak_grab_fmri_preprocess(fmri_path,opt_g); % replace the folder by the path where the results of the fmri preprocessing pipeline were stored. 
%% event times
data.covariates_group_subs = fieldnames(files_in.fmri);
for list = 1:length(data.covariates_group_subs)    
    files_in.timing.(data.covariates_group_subs{list}).sess1.('breathhold1400') = [fmri_path 'onset/nki_model_intrarun_breathhold_qc_timeseries.csv'];
end


%%%%%%%%%%%%%
%% Options %%
%%%%%%%%%%%%%

%% BASC
opt.folder_out = [ root_path '/stability_fir_perc_breathhold_1400_qc_timeserie' ]; % Where to store the results
opt.grid_scales = [5:5:50 60:10:200 220:20:400 500:100:900]; % Search in the range 2-900 clusters
opt.scales_maps = [15 15 15]; % Usually, this is initially left empty. After the pipeline ran a first time, the results of the MSTEPS procedure are used to select the final scales 
opt.stability_fir.nb_samps = 1;    % Number of bootstrap samples at the individual level. 100: the CI on indidividual stability is +/-0.1
opt.stability_fir.std_noise = 0;     % The standard deviation of the judo noise. The value 0 will not use judo noise. 
opt.stability_group.nb_samps = 500;  % Number of bootstrap samples at the group level. 500: the CI on group stability is +/-0.05
opt.nb_min_fir = 1;    % the minimum response windows number. By defaut is set to 1
opt.stability_group.min_subject = 2; % (integer, default 3) the minimal number of subjects to start the group-level stability analysis. An error message will be issued if this number is not reached.
%% FIR estimation 
opt.name_condition = 'breathhold'
opt.name_baseline = 'baseline';
opt.fir.type_norm     = 'fir';       % The type of normalization of the FIR.
opt.fir.time_window   = 189;        % The size (in sec) of the time window to evaluate the response
opt.fir.max_interpolation = 5;    % --> max 5 vols conscutifs manquants, sinon bloc rejet, mais a devrait tre irrelevant comme pas de scrubbing ici
opt.fir.time_sampling = 1.4;           % The time between two samples for the estimated response. Do not go below 1/2 TR unless there is a very large number of trials.
opt.fir.nb_min_baseline = 1 ;

%% FDR estimation
opt.nb_samps_fdr = 10000; % The number of samples to estimate the false-discovery rate

%% Multi-level options
opt.flag_ind = false;   % Generate maps/FIR at the individual level
opt.flag_mixed = false; % Generate maps/FIR at the mixed level (group-level networks mixed with individual stability matrices).
opt.flag_group = true;  % Generate maps/FIR at the group level

%%%%%%%%%%%%%%%%%%%%%%
%% Run the pipeline %%
%%%%%%%%%%%%%%%%%%%%%%
opt.flag_test = false; % Put this flag to true to just generate the pipeline without running it. Otherwise the pipeline will start.
%opt.psom.qsub_options = '-q sw -l nodes=1:ppn=4,walltime=05:00:00';
%opt.psom.max_queued = 16;
pipeline = niak_pipeline_stability_fir(files_in,opt);

%%extra
system(['cp ' mfilename('fullpath') '.m ' opt.folder_out '/.']); % make a copie of this script to output folder
