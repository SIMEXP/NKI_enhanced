function [] = nki_pipeline_glm_fir(opt)
% Script to run a glm analysis on the nki_enhanced fir responses database.
%
% SYNTAX:
% []= NKI_PIPELEINE_GLM_FIR(OPT);
%
% _________________________________________________________________________
% INPUTS:
%
% OPT
%   (structure, optional) with the following fields :
%
%   TASK
%       (string, default 'checkerboard') type of tasks that would be extracted. Possibles tasks are: 'checkerboard',
%       'breathhold'.
%
%   EXP
%       (string, default '1400') type of TR used. Possibles TR : '1400', '645'
%
%   TST
%       (string, default '') type of test used. Warning: put the prefix "_" before the test name (ex: "_noscrub")
%
%   GLM
%       (string, default '01') number of GLM experiment.
%
%   TYPE_NORM
%       (string, default 'fir') type of fir estimate normalisation. Possibles types are: 'fir',
%       'fir_shape'. see niak_normalize_fir for explanation
%
%   MODEL
%       (structure) see the OPT argument of  NKI_MODEL_<TASK-NAME>. 
%       The default parameters may work.
%
%
% _________________________________________________________________________
%
% Script to run a GLM_FIR pipeline analysis on the NKI_enhanced database.
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

%%%% Setting default input/output root path 
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
        root_path = '/media/database1/';
        fprintf ('server: %s\n',server)
        my_user_name = getenv('USER');
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%
%% Parameters
%%%%%%%%%%%%%%%%%%%%%
%% set experimentent
list_fields   = { 'task'         , 'exp'  , 'glm' ,'type_norm' ,'tst'      , 'model' };
list_defaults = { 'checkerboard' , '1400' , '01'  ,'fir'       ,'_noscrub' , struct() };

if ismember(opt.task,{'checkerboard','breathhold'}) && ismember(opt.exp,{'1400','645'})
   opt = psom_struct_defaults(opt,list_fields,list_defaults);
else
   error('wrong task or TR/EXP , see help nki_pipeline_glm_fir')
end

% Simplyfing and formatting input variables 
task  = opt.task;
exp   = opt.exp;
tst   = opt.tst;
glm   = opt.glm;
if ismember(opt.type_norm,{'fir'})
   type_norm = 'perc';
elseif ismember(opt.type_norm,{'fir_shape'})
   type_norm = 'shape';
else
   error('wrong normalisation type')
end
if isempty(tst); tst_tmp = 'on';else tst_tmp = 'off'; end
fprintf ('script to run nki_glm_fir pipeline \n Task: %s \n TR: %s\n normalisation: fir %s\n scrubbing: %s\n ',task,exp,type_norm,tst_tmp)

%%%%%%%%%%%%%%%%%%%%
%% Grabbing the results and models from fMRI preprocessing and fir pipelines
%%%%%%%%%%%%%%%%%%%%%
fmri_path = [root_path 'fmri_preprocess_ALL_task' tst '/'];
stability_fir_path = [root_path 'stability_fir_' type_norm  '_' lower(task) '_' exp tst];
opt_g.min_nb_vol = 1;     % The minimum number of volumes for an fMRI dataset to be included. This option is useful when scrubbing is used, and the resulting time series may be too short.
opt_g.min_xcorr_func = 0.5; % The minimum xcorr score for an fMRI dataset to be included. This metric is a tool for quality control which assess the quality of non-linear coregistration of functional images in stereotaxic space. Manual inspection of the values during QC is necessary to properly set this threshold.
opt_g.min_xcorr_anat = 0.5; % The minimum xcorr score for an fMRI dataset to be included. This metric is a tool for quality control which assess the quality of non-linear coregistration of the anatomical image in stereotaxic space. Manual inspection of the values during QC is necessary to properly set this threshold.
opt_g.type_files = 'fir'; % Specify to the grabber to prepare the files for the STABILITY_FIR pipeline

%%Temporary grabber for debugging
%liste_exclude = dir ([fmri_path 'anat']);
%liste_exclude = liste_exclude(43:end -1);
%liste_exclude = {liste_exclude.name};
%opt_g.exclude_subject = liste_exclude;

switch lower(task)
      case 'checkerboard'
      task_tmp = 'checBoard';
      opt_g.filter.run = {[task_tmp exp]};
      files_in = rmfield(niak_grab_fmri_preprocess(fmri_path,opt_g),{'mask','areas'}); % Replace the folder by the path where the results of the fMRI preprocessing pipeline were stored.
      case 'breathhold'
      task_tmp = 'breathHold';
      opt_g.filter.run = {[task_tmp exp]};
      files_in = rmfield(niak_grab_fmri_preprocess(fmri_path,opt_g),{'mask','areas'}); % Replace the folder by the path where the results of the fMRI preprocessing pipeline were stored.
end

%% Set the indiviual models
data.covariates_group_subs = fieldnames(files_in.fmri);
for list = 1:length(data.covariates_group_subs)    
    files_in.model.individual.(data.covariates_group_subs{list}).sess1.([task_tmp exp]) = [fmri_path 'onset/nki_model_intrarun_' lower(opt.task) '.csv'];
end

%% set the networks
files_in.networks = niak_grab_stability_fir(stability_fir_path).networks ;

%%%%%%%%%%%%%%%%%%%%
%% Create group model csv file
%%%%%%%%%%%%%%%%%%%%%

%%Load phenotypes and scrubbing data then combine them in one csv file
list_remove_pheno = { 'Download Group','frames_OK','frames_scrubbed'};
pheno_raw = niak_read_csv_cell([root_path 'nki-rs_lite_r1-2-3-4-5_phenotypic_v1.csv']);
master_cell = pheno_raw;
files_out  = niak_grab_all_preprocess([root_path 'fmri_preprocess_ALL_task' tst]);
slave_cell = niak_read_csv_cell(files_out.quality_control.group_motion.scrubbing);
ly = slave_cell(1,:);
slave_cell = slave_cell(2:end,:);
mask_slave_cell = strfind(slave_cell(:,1),[task_tmp exp]);%mask selected task_tmp and tr
mask_slave_cell = cellfun(@isempty,mask_slave_cell);
slave_cell(mask_slave_cell,:) = [];
slave_cell = [ly; slave_cell];
for cc = 1:length(slave_cell)-1;
    slave_cell{cc+1,1} = slave_cell{cc+1,1}(2:8);
end
pheno = combine_cell_tab(master_cell,slave_cell);

%%cleannig data
%remove unused pheno
mask_remove_pheno = ones(1,size(pheno,2));
for cc = 1: length(list_remove_pheno)
    mask_tmp = strfind(pheno(1,:),list_remove_pheno{cc});
    mask_tmp = cellfun(@isempty,mask_tmp);
    mask_remove_pheno = mask_remove_pheno & mask_tmp ;
end
pheno(:,~mask_remove_pheno)=[];
pheno(:,3) = strrep(pheno(:,3),'M','1'); %replace male 'M' by '1'
pheno(:,3) = strrep(pheno(:,3),'F','0'); %replace male 'M' by '0'
pheno(:,4) = strrep(pheno(:,4),'Right','1'); %replace 'Right' by '1'
pheno(:,4) = strrep(pheno(:,4),'Left','0'); %replace 'Left' by '0'
pheno(:,4) = strrep(pheno(:,4),'None','NaN'); %replace 'None' by 'NaN'
mask_pheno = cellfun(@(x) str2num(x)>100, pheno(2:end,2));%create mask for wrong age cells (higher than 99)
lx = pheno(2:end,1);
lx(mask_pheno,:) = [];%remove wrong age cells ID
for ii = 1:length(lx)
    lx{ii} = ['X' lx{ii}];%add X prefix in subject ID
end
ly = pheno(1,:);
pheno_tmp = pheno(2:end,2:end);
pheno_tmp(mask_pheno,:) = [];%remove wrong age cells data
pheno_clean = [lx  pheno_tmp];
pheno_clean = [ly ; pheno_clean];
name_save = [stability_fir_path filesep 'group_model_' type_norm  '_' lower(task) '_' exp tst '_' glm '.csv'];
niak_write_csv_cell(name_save,pheno_clean);

%goupe model 
files_in.model.group = name_save;

%%%%%%%%%%%%%
%% Options %%
%%%%%%%%%%%%%

opt.folder_out = [stability_fir_path filesep 'glm_fir_' type_norm  '_' lower(task) '_' exp tst '_' glm '/'];; % Where to store the results
opt.fdr        = 0.05; % The maximal false-discovery rate that is tolerated both for individual (single-seed) maps and whole-connectome discoveries, at each particular scale (multiple comparisons across scales are addressed via permutation testing)
opt.fwe        = 0.05; % The overall family-wise error, i.e. the probablity to have the observed number of discoveries, agregated across all scales, under the global null hypothesis of no association.
opt.nb_samps   = 1000; % The number of samples in the permutation test. This number has to be multiplied by OPT.NB_BATCH below to get the effective number of samples
opt.nb_batch   = 10; % The permutation tests are separated into NB_BATCH independent batches, which can run on parallel if sufficient computational resources are available

%% FIR 
opt.fir.name_condition = lower(task);
opt.fir.name_baseline = 'baseline';
opt.fir.type_norm     = opt.type_norm;       % The type of normalization of the FIR.
opt.fir.time_window   = opt.model.trial_duration;        % The size (in sec) of the time window to evaluate the response
opt.fir.max_interpolation = (str2num(exp)/1000)*5;    % --> max 5 vols consécutifs manquants, sinon bloc rejeté, mais ça devrait être irrelevant comme pas de scrubbing ici
opt.fir.time_sampling = str2num(exp)/1000;           % The time between two samples for the estimated response. Do not go below 1/2 TR unless there is a very large number of trials.
opt.fir.nb_min_baseline = 1 ;

%%%%%%%%%%%
%% TESTS %%
%%%%%%%%%%%  

% Regressing Age Sex and FD

%%Age
opt.test.Age.contrast.intercept = 1
opt.test.Age.contrast.Age       = 0;
opt.test.Age.contrast.Sex       = 0;
opt.test.Age.contrast.FD        = 0;


%%%%%%%%%%%%%%%%%%%%%%
%% Run the pipeline %%
%%%%%%%%%%%%%%%%%%%%%%
opt.flag_test = false; % Put this flag to true to just generate the pipeline without running it. Otherwise the pipeline will start.
opt.psom.qsub_options = '-q sw -l nodes=1:ppn=4,walltime=05:00:00';
[pipeline,opt_pipe] = niak_pipeline_glm_fir(files_in,opt);

%%extra
system(['cp ' mfilename('fullpath') '.m ' opt.folder_out '/.']); % make a copie of this script to output folder
