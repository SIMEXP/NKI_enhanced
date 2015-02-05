clear

%% Parameters
%stability_fir_perc_MOTORrh_hcp/stability_group/fir';
path_root =  '/media/database10/nki_enhanced/';
list_scale = { 'sci180_scg162_scf159' ; ...
             };
task = 'checkerboard' ;
tr = {'645'};

%% Load phenotypes
pheno = niak_read_csv_cell([path_root 'nki-rs_lite_r1-2-3-4-5_phenotypic_v1.csv']);
lx = pheno(2:end,1);
ly = pheno(1,2:end)';
pheno = pheno(2:end,2:end);

%% Load data
for tt = 1:length(tr)
    path_read  = [path_root 'stability_fir_perc_' task '_' tr{tt} '/stability_group/fir/'];
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
        subject = list_files{ff}(end-10:end-4);
        ind_s = find(ismember(lx,subject));
        pheno_r(ff,:) = pheno(ind_s,:);
        data = load([path_read list_files{ff}],list_scale{tt});
        fir_all{tt}(:,:,ff) = data.(list_scale{tt}).fir_mean;
    end
end

%% More parameters    0104892
%list_ind = [ 19 , 33 ];
list_ind = [ 113 ];
list_color = {'r','b','g','k','p'};
clf

%% Hierarchical clustering
hc = figure;
hp = figure;
for tt = 1:length(tr)
    for ii = 1:length(list_ind)
        % Clustering of subtypes
        figure(hc)
        fir_td = squeeze(fir_all{tt}(:,list_ind(ii),:));
        fir_td = fir_td./repmat(sqrt(sum(fir_td.^2,1)),[size(fir_td,1) 1]);
        fir_tf(isnan(fir_td)) = 0;
        D = niak_build_distance (fir_td);
        hier = niak_hierarchical_clustering (-D);
        %sil = niak_build_avg_silhouette(-D,hier);
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
        figure(hp)
        subplot(1,length(tr),tt)
        title(sprintf('Task %s',tr{tt}));
        for cc = 1:nb_clust(tt)        
            hold on 
            plot(mean(fir_all{tt}(:,list_ind(ii),part==cc),3),list_color{cc})
        end
     end
end

