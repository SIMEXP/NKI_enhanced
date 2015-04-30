clear

path_data = '/media/database10/nki_enhanced/';
scale = 'sci150_scg135_scf134';
task = 'breathhold' ;
tr = '1400';
fir = 'fir_perc';
scrub = '_noscrub';

path_run1 = [path_data 'stability_' fir '_' task '_' tr scrub];
path_run2 = [path_data 'stability_' fir '_' task '_' tr scrub '_resampl1400'];

%% Load data
for tt = 1:length(tr)
    path_read  = [path_root 'stability_' fir '_' lower(task) '_' tr{tt} scrub '/stability_group/fir/'];
    list_files = dir([path_read 'fir_group_level_*']);
    list_files = {list_files.name};
    for ff = 1:length(list_files);
        niak_progress( ff , length(list_files))
        data = load([path_read list_files{ff}],list_scale{tt});
        fir_all{tt}(:,:,ff) = data.(list_scale{tt}).fir_mean;
    end
end


%% More parameters
list_ind = [1:134];
%% Hierarchical clustering
for tt = 1:length(tr)
    for ii = 1:length(list_ind)
        % Clustering of subtypes
        fir_td1 = squeeze(fir_all{tt}(:,list_ind(ii),:));
        fir_td1 = fir_td1./repmat(sqrt(sum(fir_td1.^2,1)),[size(fir_td1,1) 1]);
        fir_td1(isnan(fir_td1)) = 0;
        fir_td1 = fir_td1 - repmat(mean(fir_td1,2),[1 size(fir_td1,2)]);
        
        fir_td2 = squeeze(fir_all{tt}(:,list_ind(ii),:));
        fir_td2 = fir_td2./repmat(sqrt(sum(fir_td2.^2,1)),[size(fir_td2,1) 1]);
        fir_td2(isnan(fir_td2)) = 0;
        fir_td2 = fir_td2 - repmat(mean(fir_td2,2),[1 size(fir_td2,2)]);
        
        D = niak_build_distance (fir_td1);
        hier = niak_hierarchical_clustering (-D);
        nb_clust(tt) = 5;
        part = niak_threshold_hierarchy (hier,struct('thresh',nb_clust(tt)));
        
        %% Build distance scores for all subtypes
        for cc = 1:nb_clust(tt)        
            avg_clust(ii,:,cc) = mean(fir_td1(:,part==cc),2);
            weights1(ii,:,cc) = corr(fir_td1,avg_clust(ii,:,cc));
            weights2(ii,:,cc) = corr(fir_td2,avg_clust(ii,:,cc));
        end
end

%% Reproc of weights
%tmp = corr([weights1,weights2]);
for cc = 1:nb_clust
    repro_weights(cc,ind_net) = IPN_icc([weights1(:,cc) weights2(:,cc)],2,'single'); 
end
