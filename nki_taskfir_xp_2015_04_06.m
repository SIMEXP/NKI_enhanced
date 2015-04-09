%  this script is to run on peuplier

%  NKI Breath hold subgoups fir
%  Experiments Summary:
%  
%  EXP1a: Scrubbing on - Normalisation 'perc' - scale 176
%  EXP1b: Scrubbing off - Normalisation 'perc' - scale 153
%  EXP1c: Scrubbing on - Normalisation 'shape' - scale 190
%  EXP1d: Scrubbing off - Normalisation 'shape' - scale 198


%  EXP1a:
%  In this EXP I'm using Breath-hold task individual fir from scale sci190_scg152_scf176, the fir normalisation is perc. Data are scrubbed.

clear
path_root = '/media/database10/nki_enhanced/';
list_scale = { 'sci190_scg152_scf176'};
task = 'breathhold' ;
tr = {'1400'};
fir = 'fir_perc';
scrub = '';

%%Load phenotypes
pheno = niak_read_csv_cell([path_root 'nki-rs_lite_r1-2-3-4-5_phenotypic_v1.csv']);
lx = pheno(2:end,1);
ly = pheno(1,2:end)';
pheno = pheno(2:end,2:end);

%% Load data
for tt = 1:length(tr)
    path_read  = [path_root 'stability_' fir '_' task '_' tr{tt} scrub '/stability_group/fir/'];
    path_fmri  = [path_root 'fmri_preprocess_ALL_task' scrub '/fmri/'];
    files_out  = niak_grab_all_preprocess([path_root 'fmri_preprocess_ALL_task' scrub]);
    list_files = dir([path_read 'fir_group_level_*']);
    list_files = {list_files.name};
    
    for ff = 1:length(list_files);
        subject = list_files{ff}(end-10:end-4);
        ind_s = find(ismember(lx,subject));
        if isempty(ind_s)
            warning('Could not find subject %s',subject)
            list_files{ff}= [];
        end
    end
    list_files(cellfun(@isempty,list_files)) = [];   %remove empty cells 
    pheno_r = cell(length(list_files),size(pheno,2));
    for ff = 1:length(list_files);
        subject = list_files{ff}(end-10:end-4);
        ind_s = find(ismember(lx,subject));
        pheno_r(ff,:) = pheno(ind_s,:);
        data = load([path_read list_files{ff}],list_scale{tt});
        fir_all{tt}(:,:,ff) = data.(list_scale{tt}).fir_mean;
    end
end


%  in this step we need to select the clusters based the maximum response to the task. First lets creat a max effect map then visualise it with mricron
%%matlab
path_scales =  [path_root 'stability_' fir '_' task '_' tr{1} scrub '/stability_group/' list_scale{1}];
niak_brick_mnc2nii(path_scales,[path_scales '_nii']);
%%matlab
%  system('source /opt/minc/minc-toolkit-config.sh')
cd([path_scales '_nii'])

max_effect_vol(['brain_partition_consensus_group_' list_scale{1} '.nii.gz'],['fdr_group_average_' list_scale{1} '.mat']);
system('mricron  ~/database/white_template.nii.gz -c -0 -o max_abs_eff.nii.gz -c "5redyell" -l 0.005 -h 0.5 -z  &');
system(['mricron ~/database/white_template.nii.gz -c -0 -o ' path_scales '_nii/brain_partition_consensus_group_' list_scale{1} '.nii.gz -c NIH -l 1 -h 177 -z &']);

%%matlab
%% More parameters
list_ind = [ 166 , 140, 51 ];
list_color = {'r','b','g','k','p'};

%% Hierarchical clustering
for tt = 1:length(tr)
    clf
    hold off
    for ii = 1:length(list_ind)
        % Clustering of subtypes
        figure(ii)
        fir_td = squeeze(fir_all{tt}(:,list_ind(ii),:));
        fir_td = fir_td./repmat(sqrt(sum(fir_td.^2,1)),[size(fir_td,1) 1]);
        fir_td(isnan(fir_td)) = 0;
        fir_td = fir_td - repmat(mean(fir_td,2),[1 size(fir_td,2)]);
        D = niak_build_distance (fir_td);
        hier = niak_hierarchical_clustering (-D);
        sil = niak_build_avg_silhouette(-D,hier);
        %[val,nb_clust(tt)] = max(sil(1:50));
        val = NaN;
        nb_clust(tt) = 5;
        fprintf('Task %s, max silhouette %1.2f at scale %i\n',tr{tt},val,nb_clust(tt))
        part = niak_threshold_hierarchy (hier,struct('thresh',nb_clust(tt)));
        order = niak_hier2order (hier);
        subplot(3,length(tr),tt)
        niak_visu_matrix(D(order,order));
        subplot(3,length(tr),tt+length(tr))
        niak_visu_part(part(order))
        subplot(3,length(tr),tt+2*length(tr))
        plot(sil)
        hold on
        plot(nb_clust(tt),val,'rx')
        
        % Show the subtypes
        figure(ii+length(list_ind))
        subplot(1,length(tr),tt)
        title(sprintf('Task %s',tr{tt}));
        for cc = 1:nb_clust(tt)        
            hold on 
            plot(mean(fir_all{tt}(:,list_ind(ii),part==cc),3),list_color{cc})
        end
        
        %% Build distance scores for all subtypes
        for cc = 1:nb_clust(tt)        
            avg_clust(:,cc) = mean(fir_td(:,part==cc),2);
            %weights(:,cc) = sum((fir_td-repmat(avg_clust(:,cc),[1 size(fir_td,2)])).^2);
            weights(:,cc) = corr(fir_td,avg_clust(:,cc));
        end
     end
end

%      %% GLM analysis 
%        for cc = 1:length(list_cov)
%            covar = tab2(:,list_cov(cc));
%            fd = tab2(:,21);
%            age = tab2(:,1);
%            sex = tab2(:,2);
%            mask_covar = ~isnan(covar);
%            model_covar.x = [ones(sum(mask_covar),1) niak_normalize_tseries([covar(mask_covar) fd(mask_covar)],'mean')];
%            model_covar.y = weights1(mask_covar,:);
%            model_covar.c = [0 ; 1 ; 0];
%            opt_glm.test = 'ttest';
%            opt_glm.flag_beta = true;
%            res_covar1 = niak_glm(model_covar,opt_glm);
%            model_covar.y = weights2(mask_covar,:);
%            res_covar2 = niak_glm(model_covar,opt_glm);
%            fprintf('%s\n',ly{list_cov(cc)});
%            pce1(cc,:,ind_net) = res_covar1.pce;
%            pce2(cc,:,ind_net) = res_covar2.pce;
%        end


%  EXP1b:
%  In this EXP I'm using Breath-hold task individual fir from scale sci150_scg150_scf153, the fir normalisation is perc. Data are not scrubbed.

clear
path_root =  '/media/database10/nki_enhanced/';
list_scale = { 'sci150_scg150_scf153'};
task = 'breathhold' ;
tr = {'1400'};
fir = 'fir_perc';
scrub = '_noscrub';

%%Load phenotypes
pheno = niak_read_csv_cell([path_root 'nki-rs_lite_r1-2-3-4-5_phenotypic_v1.csv']);
lx = pheno(2:end,1);
ly = pheno(1,2:end)';
pheno = pheno(2:end,2:end);

%% Load data
for tt = 1:length(tr)
    path_read  = [path_root 'stability_' fir '_' task '_' tr{tt} scrub '/stability_group/fir/'];
    path_fmri  = [path_root 'fmri_preprocess_ALL_task/fmri/'];
    list_files = dir([path_read 'fir_group_level_*']);
    list_files = {list_files.name};
    
    for ff = 1:length(list_files);
        subject = list_files{ff}(end-10:end-4);
        ind_s = find(ismember(lx,subject));
        if isempty(ind_s)
            warning('Could not find subject %s',subject)
            list_files{ff}= [];
        end
    end
    list_files(cellfun(@isempty,list_files)) = [];   %remove empty cells 
    pheno_r = cell(length(list_files),size(pheno,2));
    for ff = 1:length(list_files);
        niak_progress( ff , length(list_files))
        subject = list_files{ff}(end-10:end-4);
        ind_s = find(ismember(lx,subject));
        pheno_r(ff,:) = pheno(ind_s,:);
        data = load([path_read list_files{ff}],list_scale{tt});
        fir_all{tt}(:,:,ff) = data.(list_scale{tt}).fir_mean;
    end
end
%%matlab
pheno_num = zeros(size(pheno_r,1),1);
for xx = 1:size(pheno_r,1)
    %for yy = 1:size(pheno_r,2)
    for yy = 1
        if isempty(pheno_r{xx,yy})
            pheno_num(xx,yy) = NaN;
        else
            pheno_num(xx,yy) = str2num(pheno_r{xx,yy});
        end
    end
end
% visualise the partition
path_scales =  [path_root 'stability_' fir '_' task '_' tr{1} scrub '/stability_group/' list_scale{1} ];
opt.flag_zip = true;
niak_brick_mnc2nii(path_scales,[path_scales '_nii'],opt)
cd([path_scales '_nii'])
max_effect_vol(['brain_partition_consensus_group_' list_scale{1} '.nii.gz'],['fdr_group_average_' list_scale{1} '.mat']);
system('mricron  ~/database/white_template.nii.gz -c -0 -o max_abs_eff.nii.gz -c "5redyell" -l 0.005 -h 0.5 -z  &');
system(['mricron ~/database/white_template.nii.gz -c -0 -o ' path_scales '_nii/brain_partition_consensus_group_' list_scale{1} '.nii.gz -c NIH -l 1 -h 154 -z &']);

%% More parameters
list_ind = [ 84 , 128 , 59 , 149 , 110 , 83 ];
list_color = {'r','b','g','k','p'};

%% Hierarchical clustering
for tt = 1:length(tr)
    clf
    hold off
    for ii = 1:length(list_ind)
        % Clustering of subtypes
        figure(ii)
        fir_td = squeeze(fir_all{tt}(:,list_ind(ii),:));
        fir_td = fir_td./repmat(sqrt(sum(fir_td.^2,1)),[size(fir_td,1) 1]);
        fir_td(isnan(fir_td)) = 0;
        fir_td = fir_td - repmat(mean(fir_td,2),[1 size(fir_td,2)]);
        D = niak_build_distance (fir_td);
        hier = niak_hierarchical_clustering (-D);
        sil = niak_build_avg_silhouette(-D,hier);
        %[val,nb_clust(tt)] = max(sil(1:50));
        val = NaN;
        nb_clust(tt) = 5;
        fprintf('Task %s, max silhouette %1.2f at scale %i\n',tr{tt},val,nb_clust(tt))
        part = niak_threshold_hierarchy (hier,struct('thresh',nb_clust(tt)));
        order = niak_hier2order (hier);
        subplot(3,length(tr),tt)
        niak_visu_matrix(D(order,order));
        subplot(3,length(tr),tt+length(tr))
        niak_visu_part(part(order))
        subplot(3,length(tr),tt+2*length(tr))
        plot(sil)
        hold on
        plot(nb_clust(tt),val,'rx')
        
        % Show the subtypes
        figure(ii+length(list_ind))
        subplot(1,length(tr),tt)
        title(sprintf('Task %s',tr{tt}));
        for cc = 1:nb_clust(tt)        
            hold on 
            plot(mean(fir_all{tt}(:,list_ind(ii),part==cc),3),list_color{cc})
        end
        
        %% Build distance scores for all subtypes
        for cc = 1:nb_clust(tt)        
            avg_clust(:,cc) = mean(fir_td(:,part==cc),2);
            %weights(:,cc) = sum((fir_td-repmat(avg_clust(:,cc),[1 size(fir_td,2)])).^2);
            weights(:,cc) = corr(fir_td,avg_clust(:,cc));
        end

        %% GLM analysis 
        list_cov = { 'Age' };
        for cco = 1:length(list_cov)
            ind_cov = find(ismember(ly,list_cov{cco}));
            covar = pheno_num(:,ind_cov);
            mask_covar = ~isnan(covar);
            %%% Ugly fix
            mask_covar = mask_covar & (covar<999);
            model_covar.x = [ones(sum(mask_covar),1) niak_normalize_tseries(covar(mask_covar),'none')];
            model_covar.y = weights(mask_covar,:);
            model_covar.c = [0 ; 1 ];
            opt_glm.test = 'ttest';
            opt_glm.flag_beta = true;
            res_covar = niak_glm(model_covar,opt_glm);
            fprintf('%s\n',ly{ind_cov});
            pce(cco,:,ii) = res_covar.pce;
        end
        %plot(model_covar.x(:,2),model_covar.y(:,2),'.')
    end
end
[fdr,test] = niak_fdr(pce(:),'BH',0.05);