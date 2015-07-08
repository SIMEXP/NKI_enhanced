clear
%% parameters
path_data = '/peuplier/database10/nki_enhanced/';
scale = 'sci150_scg150_scf153';
task = 'breathhold' ;
tr = '1400';
fir = 'fir_perc';
scrub = '_noscrub';
path_run1 = [path_data 'stability_' fir '_' task '_' tr scrub '_trt1'];
path_run2 = [path_data 'stability_' fir '_' task '_' tr scrub '_trt2'];

%% Load data
path_read1  = [path_run1 '/stability_group/fir/'];
path_read2  = [path_run2 '/stability_group/fir/'];
list_files = dir([path_read1 'fir_group_level_*']);
list_files = {list_files.name};
for ff = 1:length(list_files);
    niak_progress( ff , length(list_files))
    data1 = load([path_read1 list_files{ff}],scale);
    data2 = load([path_read2 list_files{ff}],scale);
    fir_all1(:,:,ff) = data1.(scale).fir_mean;
    fir_all2(:,:,ff) = data2.(scale).fir_mean;
end

%% Hierarchical clustering
list_ind = [1:164];
for ii = 1:length(list_ind)
    % Clustering of subtypes for run1
    fir_td1 = squeeze(fir_all1(:,list_ind(ii),:));
    fir_td1 = fir_td1./repmat(sqrt(sum(fir_td1.^2,1)),[size(fir_td1,1) 1]);
    fir_td1(isnan(fir_td1)) = 0;
    fir_td1 = fir_td1 - repmat(mean(fir_td1,2),[1 size(fir_td1,2)]);
    % Clustering of subtypes for run2
    fir_td2 = squeeze(fir_all2(:,list_ind(ii),:));
    fir_td2 = fir_td2./repmat(sqrt(sum(fir_td2.^2,1)),[size(fir_td2,1) 1]);
    fir_td2(isnan(fir_td2)) = 0;
    fir_td2 = fir_td2 - repmat(mean(fir_td2,2),[1 size(fir_td2,2)]);
    % Build distance
    D = niak_build_distance (fir_td1);
    hier = niak_hierarchical_clustering (-D);
    nb_clust = 5;
    part = niak_threshold_hierarchy (hier,struct('thresh',nb_clust));
    % Build distance scores for all subtypes
    for cc = 1:nb_clust        
        avg_clust(ii,:,cc) = mean(fir_td1(:,part==cc),2);
        weights1(ii,:,cc) = corr(fir_td1,avg_clust(ii,:,cc));
        weights2(ii,:,cc) = corr(fir_td2,avg_clust(ii,:,cc));
        repro_weights(cc,ii) = IPN_icc([squeeze(weights1(ii,:,cc))' squeeze(weights2(ii,:,cc))'],2,'single'); 
    end
end
