function [] = nki_pipeline_stability_fir(opt)
% Script to run a preprocessing pipeline analysis on the nki_enhanced database.
%
% SYNTAX:
% []= NKI_PIPELEINE_STABILITY_FIR(OPT);
%
% _________________________________________________________________________
% INPUTS:
%
% OPT
%   (structure, optional) with the following fields :
%
%   TASK
%       (string, default 'visualCheckerboard') type of tasks that would be extracted. Possibles tasks are: 'visualCheckerboard',
%       'BreathHold'.
%
%   TR
%       (string, default '1400') type of TR used. Possibles TR : '1400', '645'
%
%
% _________________________________________________________________________
%
% Script to run a STABILITY_FIR pipeline analysis on the NKI_enhanced database.
%
% Copyright (c) Pierre Bellec, Yassine Benhajali
% Research Centre of the Montreal Geriatric Institute
% & Department of Computer Science and Operations Research
% University of Montreal, Québec, Canada, 2010-2014
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

%%%%%%%%%%%%%%%%%%%%%
%% Parameters
%%%%%%%%%%%%%%%%%%%%%
%% set experimentent
list_fields   = { 'task' , 'tr' };
list_defaults = { 'visualCheckerboard', '1400' };
if ischar (opt.task ) &&  ischar(opt.exp)
   opt.task = upper(opt.task);
   if ismember(opt.task,{'visualCheckerboard','BreathHold'}) && ismember(opt.tr,{'1400','645'})
      opt = psom_struct_defaults(opt,list_fields,list_defaults);
   else
      error('wrong task or TR , see help nki_pipeline_stability_fir')
   end
else 
   error ( 'you must specify the task and the TR')
end

task  = opt.task;
tr   = opt.tr;
fprintf ('script to run nki_stability_fir pipeline \n Task: %s \n TR: %s\n ',task,tr)

%% Setting input/output files 
[status,cmdout] = system ('uname -n');
server          = strtrim(cmdout);
if strfind(server,'lg-1r') % This is guillimin
    root_path = '/gs/scratch/yassinebha//';
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
        root_path = '/media/database1/';
        fprintf ('server: %s\n',server)
        my_user_name = getenv('USER');
    end
end

%% create the csv model files
opt_model.task = task;
opt_model.tr  = tr;
eval([ 'nki_model_' lower(opt_model.task) '_csv(root_path,opt_model)']);

%%%%%%%%%%%%%%%%%%%%
%% Grabbing the results from the NIAK fMRI preprocessing pipeline
%%%%%%%%%%%%%%%%%%%%%
opt_g.min_nb_vol = 1;     % The minimum number of volumes for an fMRI dataset to be included. This option is useful when scrubbing is used, and the resulting time series may be too short.
opt_g.min_xcorr_func = 0.5; % The minimum xcorr score for an fMRI dataset to be included. This metric is a tool for quality control which assess the quality of non-linear coregistration of functional images in stereotaxic space. Manual inspection of the values during QC is necessary to properly set this threshold.
opt_g.min_xcorr_anat = 0.5; % The minimum xcorr score for an fMRI dataset to be included. This metric is a tool for quality control which assess the quality of non-linear coregistration of the anatomical image in stereotaxic space. Manual inspection of the values during QC is necessary to properly set this threshold.
opt_g.type_files = 'fir'; % Specify to the grabber to prepare the files for the STABILITY_FIR pipeline

%%Temporary grabber for debugging
%liste_exclude = dir ([root_path 'fmri_preprocess_' upper(task) '_' exp '/anat']);
%liste_exclude = liste_exclude(23:end -1);
%liste_exclude = {liste_exclude.name};
%opt_g.exclude_subject = liste_exclude;
opt_g.exclude_subject ={'HCP168139'}; % to be investigated later , it make the pipelne crash, strange artifact in functional images
files_in = niak_grab_fmri_preprocess([root_path 'fmri_preprocess_' upper(task) '_' exp],opt_g); % Replace the folder by the path where the results of the fMRI preprocessing pipeline were stored. 

%% Event times
data.covariates_group_subs = fieldnames(files_in.fmri);
for list = 1:length(data.covariates_group_subs)    
    files_in.timing.(data.covariates_group_subs{list}).session1.([lower(task)(1:2) 'RL']) = [root_path 'fmri_preprocess_' upper(task) '_' exp '/EVs/hcp_model_intrarunRL.csv'];
    files_in.timing.(data.covariates_group_subs{list}).session1.([lower(task)(1:2) 'LR']) = [root_path 'fmri_preprocess_' upper(task) '_' exp '/EVs/hcp_model_intrarunLR.csv'];
end


%%%%%%%%%%%%%
%% Options %%
%%%%%%%%%%%%%

%% BASC
opt.folder_out = [ root_path '/stability_fir_perc_' upper(task) trial '_' exp ]; % Where to store the results
opt.grid_scales = [5:5:50 60:10:200 220:20:400 500:100:900]; % Search in the range 2-900 clusters
opt.scales_maps = [ 10   7   7 ; 
                    80  80  83]; % Usually, this is initially left empty. After the pipeline ran a first time, the results of the MSTEPS procedure are used to select the final scales
opt.stability_fir.nb_samps = 100;    % Number of bootstrap samples at the individual level. 100: the CI on indidividual stability is +/-0.1
opt.stability_fir.std_noise = 0;     % The standard deviation of the judo noise. The value 0 will not use judo noise. 
opt.stability_group.nb_samps = 500;  % Number of bootstrap samples at the group level. 500: the CI on group stability is +/-0.05
opt.nb_min_fir = 1;    % the minimum response windows number. By defaut is set to 1
opt.stability_group.min_subject = 2; % (integer, default 3) the minimal number of subjects to start the group-level stability analysis. An error message will be issued if this number is not reached.
%% FIR estimation 
opt.name_condition = 'rh';
opt.name_baseline = 'baseline';
opt.fir.type_norm     = 'fir';       % The type of normalization of the FIR.
opt.fir.time_window   = 16.5;        % The size (in sec) of the time window to evaluate the response
opt.fir.max_interpolation = 7.2;    % --> max 10 vols consécutifs manquants (TR = 0.72s), sinon bloc rejeté, mais ça devrait être irrelevant comme pas de scrubbing ici
opt.fir.time_sampling = 0.72;           % The time between two samples for the estimated response. Do not go below 1/2 TR unless there is a very large number of trials.
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
%opt.psom.qsub_options = 'q lm -l nodes=1:ppn=12,walltime=05:00:00';
opt.psom.qsub_options = '-q sw -l nodes=1:ppn=4,walltime=05:00:00';
pipeline = niak_pipeline_stability_fir(files_in,opt);

%%extra
system(['cp ' mfilename('fullpath') '.m ' opt.folder_out '/.']); % make a copie of this script to output folder
