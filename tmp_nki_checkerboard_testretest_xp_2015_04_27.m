%________session1____________
clear
path_root =  '/media/database10/nki_enhanced/';
list_scale = { 'sci150_scg135_scf134'};
task = 'checkerboard' ;
tr = {'645'};
fir = 'fir_perc';
scrub = '_noscrub';
list_remove_pheno = { 'Download Group','frames_OK','frames_scrubbed'};

%%Load phenotypes and scrubbing data
%combine pheno and scrubbing
pheno_raw = niak_read_csv_cell([path_root 'nki-rs_lite_r1-2-3-4-5_phenotypic_v1.csv']);
master_cell = pheno_raw;
files_out  = niak_grab_all_preprocess([path_root 'fmri_preprocess_ALL_task' scrub]);
slave_cell = niak_read_csv_cell(files_out.quality_control.group_motion.scrubbing);
ly = slave_cell(1,:);
slave_cell = slave_cell(2:end,:);
mask_slave_cell = strfind(slave_cell(:,1),['checBoard' tr{1}]);%mask selected task and tr
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
mask_pheno = cellfun(@(x) str2num(x)>100, pheno(2:end,2));%create mask for wrong age cells less than 100
lx = pheno(2:end,1);
lx(mask_pheno,:) = [];%remove wrong age cells ID
ly = pheno(1,2:end)';
pheno = pheno(2:end,2:end);
pheno(mask_pheno,:) = [];%remove wrong age cells data

%% Load datapheno
for tt = 1:length(tr)
    path_read  = [path_root 'stability_' fir '_' lower(task) '_' tr{tt} scrub '/stability_group/fir/'];
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

%%convert all pheno from string to numeric
pheno_num = zeros(size(pheno_r));
for xx = 1:size(pheno_r,1)
    %for yy = 1:size(pheno_r,2)
    for yy = 1:size(pheno_r,2)
        if isempty(pheno_r{xx,yy})
            pheno_num(xx,yy) = NaN;
        else
            pheno_num(xx,yy) = str2num(pheno_r{xx,yy});
        end
    end
end

%% More parameters
list_ind = [  97 , 72 , 53 ];
%% Hierarchical clustering
for tt = 1:length(tr)
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
        
        %% Build distance scores for all subtypes
        for cc = 1:nb_clust(tt)        
            avg_clust(:,cc) = mean(fir_td(:,part==cc),2);
            weights(:,cc) = corr(fir_td,avg_clust(:,cc));
        end
        
        %% GLM analysis 
        list_cov = { 'Age','Sex','FD' };
        mask_covar = [];
        ind=[];
        covar = [];
        for cco = 1:length(list_cov)
            ind_cov = find(ismember(ly,list_cov{cco}));
            covar = [covar pheno_num(:,ind_cov)];
            mask_covar =[mask_covar ~isnan(covar)];   
        end
        [y,x]=find(mask_covar == 0);
        ind = ones(size(mask_covar),1);
        ind(unique(y)) = 0;
        model_tmp = [];
        % load model
        for ccx = 1 : length(list_cov)
            model_tmp = [model_tmp niak_normalize_tseries(covar(logical(ind),ccx),'none')];
        end    
        model_covar.x = [ones(sum(ind),1) model_tmp];
        model_covar.y = weights(logical(ind),:);
        for cco = 1:length(list_cov) 
            ind_cov = find(ismember(ly,list_cov{cco}));
            model_covar.c = zeros(1,size(model_covar.x,2))';           
            model_covar.c(cco+1) = 1;
            opt_glm.test = 'ttest';
            opt_glm.flag_beta = true;
            res_covar = niak_glm(model_covar,opt_glm);
            fprintf('%s\n',ly{ind_cov});
            pce1(cco,:,ii) = res_covar.pce;
        end
    end
end


%___________session2_____________________

path_root =  '/media/database10/nki_enhanced/';
list_scale = { 'sci150_scg135_scf134'};
task = 'checkerboard' ;
tr = {'645'};
fir = 'fir_perc';
scrub = '_noscrub';
list_cov = { 'Age','Sex','FD' };
list_remove_pheno = { 'Download Group','frames_OK','frames_scrubbed'};

%%Load phenotypes and scrubbing data
%combine pheno and scrubbing
pheno_raw = niak_read_csv_cell([path_root 'nki-rs_lite_r1-2-3-4-5_phenotypic_v1.csv']);
master_cell = pheno_raw;
files_out  = niak_grab_all_preprocess([path_root 'fmri_preprocess_ALL_task' scrub]);
slave_cell = niak_read_csv_cell(files_out.quality_control.group_motion.scrubbing);
ly = slave_cell(1,:);
slave_cell = slave_cell(2:end,:);
mask_slave_cell = strfind(slave_cell(:,1),['checBoard' tr{1}]);%mask selected task and tr
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
mask_pheno = cellfun(@(x) str2num(x)>100, pheno(2:end,2));%create mask for wrong age cells less than 100
lx = pheno(2:end,1);
lx(mask_pheno,:) = [];%remove wrong age cells ID
ly = pheno(1,2:end)';
pheno = pheno(2:end,2:end);
pheno(mask_pheno,:) = [];%remove wrong age cells data

%% Load datapheno
for tt = 1:length(tr)
    path_read  = [path_root 'stability_' fir '_' lower(task) '_' tr{tt} scrub '/stability_group/fir/'];
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

%%convert all pheno from string to numeric
pheno_num = zeros(size(pheno_r));
for xx = 1:size(pheno_r,1)
    %for yy = 1:size(pheno_r,2)
    for yy = 1:size(pheno_r,2)
        if isempty(pheno_r{xx,yy})
            pheno_num(xx,yy) = NaN;
        else
            pheno_num(xx,yy) = str2num(pheno_r{xx,yy});
        end
    end
end

%% More parameters
list_ind = [  97 , 72 , 53 ];
list_color = {'r','b','g','k','p'};

%% Hierarchical clustering
for tt = 1:length(tr)
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
        
        %% Build distance scores for all subtypes
        for cc = 1:nb_clust(tt)        
            avg_clust(:,cc) = mean(fir_td(:,part==cc),2);
            %weights(:,cc) = sum((fir_td-repmat(avg_clust(:,cc),[1 size(fir_td,2)])).^2);
            weights(:,cc) = corr(fir_td,avg_clust(:,cc));
        end
        
        %% GLM analysis 
        list_cov = { 'Age','Sex','FD' };
        mask_covar = [];
        ind=[];
        covar = [];
        for cco = 1:length(list_cov)
            ind_cov = find(ismember(ly,list_cov{cco}));
            covar = [covar pheno_num(:,ind_cov)];
            mask_covar =[mask_covar ~isnan(covar)];   
        end
        [y,x]=find(mask_covar == 0);
        ind = ones(size(mask_covar),1);
        ind(unique(y)) = 0;
        model_tmp = [];
        % load model
        for ccx = 1 : length(list_cov)
            model_tmp = [model_tmp niak_normalize_tseries(covar(logical(ind),ccx),'none')];
        end    
        model_covar.x = [ones(sum(ind),1) model_tmp];
        model_covar.y = weights(logical(ind),:);
        for cco = 1:length(list_cov) 
            ind_cov = find(ismember(ly,list_cov{cco}));
            model_covar.c = zeros(1,size(model_covar.x,2))';           
            model_covar.c(cco+1) = 1;
            opt_glm.test = 'ttest';
            opt_glm.flag_beta = true;
            res_covar = niak_glm(model_covar,opt_glm);
            fprintf('%s\n',ly{ind_cov});
            pce2(cco,:,ii) = res_covar.pce;
        end
    end
end


%  [fdr,test1] = niak_fdr(pce1(:),'BH',0.05);
%  [fdr,test2] = niak_fdr(pce2(:),'BH',0.05);


%% Test for associations session1
[fdr,test1] = niak_fdr(pce1(:),'BH',0.2);
test1 = reshape(test1,size(pce1))

%% Test for associations session2
[fdr,test2] = niak_fdr(pce2(:),'BH',0.2);
test2 = reshape(test2,size(pce2))

%% Test for associations combining runs 1 and 2 as if they were independent
chi2_trt = -2*(log(pce1)+log(pce2));
pce_comb = 1-chi2cdf(chi2_trt , 4);

[fdr,test_comb] = niak_fdr(pce_comb(:),'BH',0.2);
test_comb = reshape(test_comb,size(pce_comb))

