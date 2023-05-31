tic
clear all
close all
addpath('/data/matlab_functions/')
terminate(pyenv)
pyenv("ExecutionMode","OutOfProcess")

savefig_flag = 0;
savefig_pdf = 0;

load PE_QC800_400_afterloading_30-Jan-2022.mat

cellidclsuter_tb = loadCellFile('stromal_step3_cellid_cluster_beforeclean_20-Nov-2022.txt');

tropho_clustername = loadCellFile('cluster_name_stromal_221120.txt');
tropho_clustername(:,2) = cellfun(@(x) ['STROMAL_',x], tropho_clustername(:,2),'UniformOutput',0);
tropho_clustername(:,2) = regexprep(tropho_clustername(:,2),'STROMAL_EXC','EXC');
tropho_clustername = tropho_clustername(2:end,:);
cellidclsuter_tb = [cellidclsuter_tb,cell(length(cellidclsuter_tb),1)];
c = cell2mat(cellidclsuter_tb(:,2));
for i=1:length(tropho_clustername)
    ind = find(c==tropho_clustername{i});
    cellidclsuter_tb(ind,3) = repmat(tropho_clustername(i,2),length(ind),1);
end


cellid_clusters = cell(length(cellid),3);
[~,loc]= ismember(cellidclsuter_tb(:,1),cellid);
cellid_clusters(loc(loc>0),:) = cellidclsuter_tb;
rmv = find(strcmpi(cellid_clusters(:,3),'EXC'));
cellid_clusters(rmv,:) = cell(length(rmv),3);

validcells = cell2mat(cellfun(@(x) ~isempty(x), cellid_clusters(:,1),'UniformOutput',0));
sum(validcells)

tot_mol = sum(data);
tot_mol(tot_mol>5e4) = 5e4;
tot_genes = sum(data>0);
loc = find(validcells);
% % % % % % 
data = data(:,loc);
cellid = cellid(loc);
sample = sample(loc);
tot_mol = tot_mol(loc);
v3_flag = v3_flag(loc);
ctrl_flag = ctrl_flag(loc);
peearly_flag = peearly_flag(loc);
pelate_flag = pelate_flag(loc);
female_flag = female_flag(loc);
iugr_flag = iugr_flag(loc);
ctrlearly_flag = ctrlearly_flag(loc);
csection_flag = csection_flag(loc);
vaginalbirth_flag = vaginalbirth_flag(loc);
induction_flag = induction_flag (loc);
noinduction_flag = noinduction_flag(loc);
magnesium_flag = magnesium_flag(loc);
spinal_flag = spinal_flag(loc);
epidural_flag = epidural_flag(loc);
generalanesthesia_flag = generalanesthesia_flag(loc);
week_flag = week_flag(loc);
weight_flag = weight_flag(loc);
weightprct_flag = weightprct_flag(loc);
donorage_flag = donorage_flag(loc);
donor_id_flag = donor_id_flag(loc);
cellid_clusters = cellid_clusters(loc,:);

% validcells = false(size(cellid));
% validcells(loc) = true;
sample_uni = {'48-1','53-1','61-1','63-1','63-2','67-1','67-2','68-1'...
    ,'68-2','77-1','78-1','86-1','87-1','88-1','92-1','99-1','101-1'...
    ,'103-1','104-1','105-1','106-1','107-1','107-2','108-1','113-2'...
    ,'114-2','116-1','117-1','119-1','120-1','120-2'};%
for i=1:length(sample_uni)
    fprintf(['valid cells in ',sample_uni{i},' = ',num2str(sum((strcmpi(sample,sample_uni{i})))),'\n']);
end


% % % % % % % % % % 
data = ceil(data./repmat(sum(data),length(data(:,1)),1)*10e3);

sex_genes = {'XIST','MTRNR2L8', 'EIF1AY', 'DDX3Y', 'RPS4Y1', 'KDM5D','MTRNR2L12','MTRNR2L10','MTRNR2L12'};
in = find(sum(data>0,2)>5 & sum(data>0,2)<length(data(1,:))*0.5 & ~ismember(geneid(:,1),sex_genes)); % & ~ismember(geneid(:,1),batch_genes)

corr_filt = cv_vs_m_selection(data(in,:),geneid(in),[],1,0);
% % % % % % % % % % % 


z = linkage(data(in(corr_filt),:),'ward','correlation');
tic
idx = cluster(z,'maxclust',50);
toc
d = corr_mat(data(in(corr_filt),:)');
dsum = sum(d>0.2,2);
leaforder = optimalleaforder(z,squareform(1-d,'tovector'));

[~,xi] = sort(idx);
figure;
set(gcf,'color','w','position',[20,20,900,800])
h1 = axes('position',[0.2,0.1,0.7,0.85]);
imagesc(d(leaforder,leaforder),[0,0.3]);
set(gca,'ytick',[1:length(leaforder)],'YTickLabel',geneid(in(corr_filt(leaforder))))
freezeColors(h1)
text([1:length(idx)],[1:length(idx)],num2str(idx(leaforder)))
h2 = axes('position',[0.91,0.1,0.04,0.85]);
imagesc(idx(leaforder));
colormap('prism');
h3 = axes('position',[0.95,0.1,0.05,0.85]);
plot(sum(d(leaforder,leaforder)>0.2,2),[1:length(leaforder)],'.');
set(h3,'ydir','reverse');
linkaxes([h1,h2,h3],'y')
table1 = [geneid(in(corr_filt(leaforder))),m2c(idx(leaforder))];
% saveCellFile(table1,['gene_correlation_step4_stromal_',date,'.txt'])
% % % % % % % % % % % 

data_orig_all = data;
geneid_all = geneid;
data = data(in(corr_filt),:);
geneid = geneid(in(corr_filt));

%% pca analysis
moldata = data;

datalog_tmp = cent_norm([(log2(moldata+1))]);
data_tsne = datalog_tmp;
% data_tsne = cent_norm(log2(datamarkers+1));
initial_dims = length(corr_filt);
[prj,m,D,V,Q] = pca_wis(data_tsne',initial_dims);
D = diag(D);
initial_dims = findknee(D);
figure;
subplot(1,2,1)
plot(cumsum(D)); hold on;
plot(initial_dims,sum(D(1:initial_dims)),'sk');
subplot(1,2,2)
plot(D); hold on;
plot(initial_dims,D(initial_dims),'sk');
title(['opt PC = ',num2str(initial_dims)]);
% initial_dims = 10;
prj = prj(:,1:initial_dims);
init = prj(:,1:2)/std(prj(:,1))*1e-4;
init = init-repmat(mean(init),length(init),1);
%% test harmony py
s = prj;
batchid = donor_id_flag;
usepylib = 1;
[sout]=harmonypy(s,batchid,usepylib);
prj = sout;
%% tsne 
perplexity = 200;
options = statset('MaxIter',1000);
mapped_xy = tsne(prj,'Algorithm','barneshut','Distance','correlation','NumDimensions',2,'NumPCAComponents',0,.....
    'Perplexity',perplexity,'Standardize',true,'InitialY',init,'LearnRate',length(cellid)/12,'Theta',0.5,'Verbose',1,'Options',options,'Exaggeration',20);
toc


% % % % % % % % % % % % % % % % % % % 
clusteruni = unique(cellid_clusters(:,3));
T_cells_tmp = zeros(length(cellid),1);
clustersize = zeros(length(clusteruni),1);
for i=1:length(clusteruni)
    ind = find(strcmpi(cellid_clusters(:,3),clusteruni{i}));
    T_cells_tmp(ind) = i;
    clustersize(i) = length(ind);
end
T_cells_tmp_uni = unique(T_cells_tmp);

%% color per cell type 
colors = loadCellFile_turbo('color_celltypes_ordered_HEX_221120.txt',1);
[~,l] = ismember(clusteruni,colors(:,1));
colors = colors(l,:);
% colors = distinguishable_colors(length(T_cells_tmp_uni)+1);
%% cellid tsne data
[~,loc]= ismember(cellid,cellid_clusters(:,1));
cellid_tsne = cell(length(cellid),6);
cellid_tsne(:,1) = cellid; 
% cellid_tsne(:,2) = cellid_clusters(loc,2);
cellid_tsne(:,3) = cellid_clusters(loc,3);
[~,loc]= ismember(cellid_tsne(:,3),colors(:,1));
cellid_tsne(:,2) = num2cell(loc);
cellid_tsne(:,4) = colors(loc,2);
cellid_tsne(:,[5,6]) = num2cell(mapped_xy);
if savefig_flag==1
    writecell(cellid_tsne,['Stromal_cellid-tsne_data_',date,'.csv']);
end
%% tsne by cluster
figure;
set(gcf,'color','w','position',[20,20,900,800])
for i=1:length(T_cells_tmp_uni)
    ii=find(T_cells_tmp==i); h=plot(mapped_xy(ii,1),mapped_xy(ii,2),'.','color',hex2rgb(colors{i,2}),'markersize',4); hold on;
    xx(i,:) = [clusteruni{i},m2c(length(ii))]
end
axis tight;
xl = get(gca,'xlim');
yl = get(gca,'ylim');
rr = 0.05*diff(xl);
for i=1:length(T_cells_tmp_uni)
    if i>=0
        in = T_cells_tmp==i;
        th = (rand*360)*pi/180;
        xo = rr*cos(th); yo = rr*sin(th);
        plot([median(mapped_xy(in,1)),median(mapped_xy(in,1))+xo] , ...
            [median(mapped_xy(in,2)),median(mapped_xy(in,2))+yo],'-','color',0.7*[1,1,1]);
        ht = text(median(mapped_xy(in,1))+xo,median(mapped_xy(in,2))+yo,regexprep(clusteruni{i},'_','-'));
        set(ht,'fontsize',8,'FontName','arial');%,'BackgroundColor',0.8*[1,1,1]
    end
end
axis tight;
axis equal
axis off

title(['#C=',num2str(max(T_cells_tmp)),', perplexity=',num2str(perplexity)],'fontsize',8);
if savefig_flag==1
    savefig(gcf,['tsne_final_step4_stromal_by_cluster_',date,'.fig'])    
end
if savefig_pdf==1
    eval(['export_fig tsne_final_step4_stromal_by_cluster_',date,'.pdf']);
end
%% flags sorting 
[T_cells_tmp,xi] = sort(T_cells_tmp);
xi(T_cells_tmp==0) = [];
T_cells_tmp(T_cells_tmp==0) = [];
data_sorted_all = data(:,xi);
data_orig_all_sorted = data_orig_all(:,xi);
cellid_sorted = cellid((xi));
sample_sorted = sample((xi));
mapped_xy = mapped_xy(xi,:);
tot_mol_sorted = tot_mol(xi);
v3_flag_sorted = v3_flag((xi));
ctrl_flag_sorted = ctrl_flag((xi));
peearly_flag_sorted = peearly_flag((xi));
pelate_flag_sorted = pelate_flag((xi));
female_flag_sorted = female_flag((xi));
iugr_flag_sorted = iugr_flag((xi));
ctrlearly_flag_sorted = ctrlearly_flag((xi));
csection_flag_sorted = csection_flag((xi));
vaginalbirth_flag_sorted = csection_flag((xi));
induction_flag_sorted = induction_flag((xi));
noinduction_flag_sorted = noinduction_flag((xi));
magnesium_flag_sorted = magnesium_flag((xi));
spinal_flag_sorted = spinal_flag((xi));
epidural_flag_sorted = epidural_flag((xi));
generalanesthesia_flag_sorted = generalanesthesia_flag((xi));
week_flag_sorted = week_flag(xi);
weight_flag_sorted = weight_flag(xi);
weightprct_flag_sorted = weightprct_flag(xi);
donorage_flag_sorted = donorage_flag(xi);
donor_id_flag_sorted = donor_id_flag(xi);
cellid_clusters = cellid_clusters(xi,:);

% % % % % % % % % % 

no_dims = 1;
initial_dims = 10;
perplexity = 5;
epsilon = 100;
dist_flag = 2;
theta = 0.5;
rand_seed = 13;
data_tsne = cent_norm(log2(data_sorted_all+1));
xi = [1:length(T_cells_tmp)];
for i=1:length(T_cells_tmp_uni)
    i
    ind = find(T_cells_tmp==i);
    if length(ind)>20
        tmp1d = fast_tsne((data_tsne(:,ind))', no_dims, initial_dims, perplexity,theta, rand_seed);
        [~,xitmp] = sort(tmp1d);
        xi(ind) = xi(ind((xitmp)));
    end
end

data_sorted_all = data_sorted_all(:,xi);
data_orig_all_sorted = data_orig_all_sorted(:,xi);
cellid_sorted = cellid_sorted((xi));
sample_sorted = sample_sorted((xi));
tot_mol_sorted = tot_mol_sorted(xi);
v3_flag_sorted = v3_flag_sorted((xi));
ctrl_flag_sorted = ctrl_flag_sorted((xi));
peearly_flag_sorted = peearly_flag_sorted((xi));
pelate_flag_sorted = pelate_flag_sorted((xi));
female_flag_sorted = female_flag_sorted((xi));
iugr_flag_sorted = iugr_flag_sorted((xi));
ctrlearly_flag_sorted = ctrlearly_flag_sorted((xi));
csection_flag_sorted = csection_flag_sorted((xi));
vaginalbirth_flag_sorted = vaginalbirth_flag_sorted((xi));
induction_flag_sorted = induction_flag_sorted((xi));
noinduction_flag_sorted = noinduction_flag_sorted((xi));
magnesium_flag_sorted = magnesium_flag_sorted((xi));
spinal_flag_sorted = spinal_flag_sorted((xi));
epidural_flag_sorted = epidural_flag_sorted((xi));
generalanesthesia_flag_sorted = generalanesthesia_flag_sorted((xi));
week_flag_sorted = week_flag_sorted(xi);
weight_flag_sorted = weight_flag_sorted(xi);
weightprct_flag_sorted = weightprct_flag_sorted(xi);
donorage_flag_sorted = donorage_flag_sorted(xi);
donor_id_flag_sorted = donor_id_flag_sorted(xi);
cellid_clusters = cellid_clusters(xi,:);

mapped_xy = mapped_xy(xi,:);
%% pca analysis
meangr_mat = zeros(length(moldata(:,1)),length(T_cells_tmp_uni));
clust_cent = zeros(length(T_cells_tmp_uni),2);
for jjj=1:length(T_cells_tmp_uni)
    jjj
    meangr_mat(:,jjj) = mean(log2(data_sorted_all(:,T_cells_tmp==T_cells_tmp_uni(jjj))+1),2);
    clust_cent(jjj,:) = [median(mapped_xy(T_cells_tmp==T_cells_tmp_uni(jjj),1)),median(mapped_xy(T_cells_tmp==T_cells_tmp_uni(jjj),2))];
end
meangr_mat1 = meangr_mat;

[prj,m,D,V,Q] = pca_wis(meangr_mat1',5);
Zpca = linkage(prj,'ward','correlation');
Dpca = pdist(prj,'correlation');
leaforder_pca = optimalleaforder(Zpca,Dpca);
figure;
set(gcf,'position',[100,100,1000,1000],'color','w')
% h1 = axes('position',[0.03,0.03,0.23,0.93]);
h1 = axes('position',[0.03,0.1,0.2,0.75]);
hden = dendrogram(Zpca,length(leaforder_pca),'Reorder',leaforder_pca,'Orientation','left');
axis off
set(gca,'ylim',[0.5,length(leaforder_pca)+0.5])
h2 = axes('position',[0.3,0.1,0.7,0.75]);
% h2 = axes('position',[0.38,0.03,0.98,0.93]);
x=squareform(1-Dpca)+eye(length(leaforder_pca)); imagesc(x(leaforder_pca,leaforder_pca));
colormap('summer')
set(gca,'ytick',[1:length(clusteruni)],'YTickLabel',regexprep(clusteruni(leaforder_pca),'_','-'),'xtick',[],'fontsize',8,'ydir','normal')
linkaxes([h1,h2],'y');
if savefig_flag==1
    savefig(gcf,['tree_step4_stromal_',date,'.fig']) 
    eval(['export_fig tree_step4_stromal_',date,'.pdf']);
end
leaforder = leaforder_pca;
clusteruni = clusteruni(leaforder_pca);

Zpca_post = linkage(prj(leaforder_pca,:),'ward','correlation');

T_cells_tmp_new = zeros(size(T_cells_tmp));
for i=1:length(leaforder)
    T_cells_tmp_new(T_cells_tmp==T_cells_tmp_uni(leaforder(i))) = i;
end
idxuni_new = unique(T_cells_tmp_new);
[~,xi] = sort(T_cells_tmp_new);
T_cells_tmp_new = T_cells_tmp_new(xi);
data_sorted_all = data_sorted_all(:,xi);
data_orig_all_sorted = data_orig_all_sorted(:,xi);
cellid_sorted = cellid_sorted(xi);
% % tissue = tissue(xi);
mapped_xy = mapped_xy(xi,:);
cells_bor_2 = find(diff(T_cells_tmp_new)>0)+1;
sample_sorted = sample_sorted(xi);
tot_mol_sorted = tot_mol_sorted(xi);
v3_flag_sorted = v3_flag_sorted((xi));
ctrl_flag_sorted = ctrl_flag_sorted((xi));
peearly_flag_sorted = peearly_flag_sorted((xi));
pelate_flag_sorted = pelate_flag_sorted((xi));
female_flag_sorted = female_flag_sorted((xi));
iugr_flag_sorted = iugr_flag_sorted((xi));
ctrlearly_flag_sorted = ctrlearly_flag_sorted((xi));
csection_flag_sorted = csection_flag_sorted((xi));
vaginalbirth_flag_sorted = vaginalbirth_flag_sorted((xi));
induction_flag_sorted = induction_flag_sorted((xi));
noinduction_flag_sorted = noinduction_flag_sorted((xi));
magnesium_flag_sorted = magnesium_flag_sorted((xi));
spinal_flag_sorted = spinal_flag_sorted((xi));
epidural_flag_sorted = epidural_flag_sorted((xi));
generalanesthesia_flag_sorted = generalanesthesia_flag_sorted((xi));
week_flag_sorted = week_flag_sorted(xi);
weight_flag_sorted = weight_flag_sorted(xi);
weightprct_flag_sorted = weightprct_flag_sorted(xi);
donorage_flag_sorted = donorage_flag_sorted(xi);
donor_id_flag_sorted = donor_id_flag_sorted(xi);
cellid_clusters = cellid_clusters(xi,:);

%% tsne by sample
T_cells_tmp = T_cells_tmp_new;
T_cells_tmp_uni = unique(T_cells_tmp);

colorsvec = distinguishable_colors(length(sample_uni)+1);
figure;
set(gcf,'color','w','position',[20,20,900,800]);
[ha, pos] = tight_subplot(1, 1, [0.02,0.02], [0.02,0.02], [0.02,0.02]);
axes(ha(1))
for i=1:length(sample_uni)
    s = scatter(mapped_xy(strcmpi(sample_sorted,sample_uni{i}),1),...
        mapped_xy(strcmpi(sample_sorted,sample_uni{i}),2),5,colorsvec(i,:),'filled'); hold on;
%     alpha(s,0.4);
end
axis tight
axis equal
axis off
legend(sample_uni)

if savefig_flag==1
    savefig(gcf,['stromal_step4_tsne_by_sample_perplexity_',num2str(perplexity),'_',date,'.fig'])
 
end
%% dendrogram tree 
[table1, table2] = dendrogram_split_markers(cellfun(@(x) num2str(x), m2c(T_cells_tmp_uni),'UniformOutput',0).....
    ,cellfun(@(x) num2str(x), m2c(T_cells_tmp(:,1)),'UniformOutput',0),Zpca_post,data_sorted_all,geneid);
saveCellFile(table1,['FC_final_step4_stromal_dendrogram_junction_split_markers_',date,'.txt']);
saveCellFile(table2,['FC_final_step4_stromal_dendrogram_junction_split_markers_by_average_',date,'.txt']);
ind = [[3:53:length(table1(:,1))];[4:53:length(table1(:,1))];[5:53:length(table1(:,1))];[6:53:length(table1(:,1))];[7:53:length(table1(:,1))]];
ind = ind(:);
treemark1 = table1(ind,[1,7]);
treemark1 = [reshape(treemark1(:,1),5,[]);reshape(treemark1(:,2),5,[])];
treemark1 = flipud(treemark1(:));
rmv = [];
for i=2:length(treemark1)
    if sum(strcmpi(treemark1{i},treemark1(1:i-1)))>0
        rmv = [rmv,i];
    end
end
treemark1(rmv) = [];
%% tsne by experimental group 
figure;
set(gcf,'color','w','position',[20,20,1200,1000]);
[ha, pos] = tight_subplot(2, 3, [0.02,0.02], [0.02,0.02], [0.02,0.02]);
axes(ha(1))
ms = 2;
f = female_flag_sorted;
plot(mapped_xy(:,1),mapped_xy(:,2),'.','color',0.6*[1,1,1],'markersize',ms); hold on;
plot(mapped_xy(f,1),mapped_xy(f,2),'.','color',[1,0,0],'markersize',ms); hold on;
axis tight
axis equal
axis off
title('Female')
axes(ha(2))
f = ~female_flag_sorted;
plot(mapped_xy(:,1),mapped_xy(:,2),'.','color',0.6*[1,1,1],'markersize',ms); hold on;
plot(mapped_xy(f,1),mapped_xy(f,2),'.','color',[1,0,0],'markersize',ms); hold on;
axis tight
axis equal
axis off
title('Male')
axes(ha(3))
f = iugr_flag_sorted;
plot(mapped_xy(:,1),mapped_xy(:,2),'.','color',0.6*[1,1,1],'markersize',ms); hold on;
plot(mapped_xy(f,1),mapped_xy(f,2),'.','color',[1,0,0],'markersize',ms); hold on;
axis tight
axis equal
axis off
title('IUGR')
axes(ha(4))
f = peearly_flag_sorted;
plot(mapped_xy(:,1),mapped_xy(:,2),'.','color',0.6*[1,1,1],'markersize',ms); hold on;
plot(mapped_xy(f,1),mapped_xy(f,2),'.','color',[1,0,0],'markersize',ms); hold on;
axis tight
axis equal
axis off
title('PE early')
axes(ha(5))
f = pelate_flag_sorted;
plot(mapped_xy(:,1),mapped_xy(:,2),'.','color',0.6*[1,1,1],'markersize',ms); hold on;
plot(mapped_xy(f,1),mapped_xy(f,2),'.','color',[1,0,0],'markersize',ms); hold on;
axis tight
axis equal
axis off
title('PE late')
axes(ha(6))
f = ctrl_flag_sorted;
plot(mapped_xy(:,1),mapped_xy(:,2),'.','color',0.6*[1,1,1],'markersize',ms); hold on;
plot(mapped_xy(f,1),mapped_xy(f,2),'.','color',[1,0,0],'markersize',ms); hold on;
axis tight
axis equal
axis off
title('Control')

%% median per condition

condition_sorted = zeros(6,length(peearly_flag_sorted));
condition_sorted(1,ctrlearly_flag_sorted) = 1;
condition_sorted(2,peearly_flag_sorted) = 1;
condition_sorted(3,ctrl_flag_sorted) = 1;
condition_sorted(4,pelate_flag_sorted) = 1;
condition_sorted(5,iugr_flag_sorted) = 1;
condition_sorted(6,(peearly_flag_sorted | pelate_flag_sorted) & ~iugr_flag_sorted) = 1;
donor_id_uni = unique(donor_id_flag_sorted);

% tsne by samples
figure;
set(gcf,'color','w','position',[20,20,1400,960])
[ha, pos] = tight_subplot(4, 8, [0.02,0.02], [0.02,0.02], [0.02,0.02]);
xy_range = [min(mapped_xy),max(mapped_xy)];
colorvec = distinguishable_colors(length((condition_sorted(:,1))));
for i=1:length(sample_uni)
    axes(ha(i));
    [~,c] = max(condition_sorted(1:4,strcmpi(sample_sorted,sample_uni{i})));
    c = unique(c);
    plot(mapped_xy(strcmpi(sample_sorted,sample_uni{i}),1),mapped_xy(strcmpi(sample_sorted,sample_uni{i}),2),'.','color',colorvec(c,:),'markersize',4);
    title(sample_uni{i})
    set(gca,'xlim',[xy_range(1),xy_range(3)],'ylim',[xy_range(2),xy_range(4)])
    axis off
    axis equal
    %     alpha(s,0.4);
end
% % % % % % % % % % % % % % % % % % 
% calculate the fraction of each cell type per sa,pe or donor
per_cluster_insample = zeros(length(sample_uni),length(T_cells_tmp_uni));
for i=1:length(sample_uni)
    in = find(strcmpi(sample_sorted,sample_uni{i}));
    per_cluster_insample(i,:) = histcounts(T_cells_tmp(in),[1:length(T_cells_tmp_uni)+1]);
end
per_cluster_insample = per_cluster_insample./repmat(sum(per_cluster_insample,2),1,length(T_cells_tmp_uni));
condition_sorted_uni = unique(condition_sorted);


per_cluster_indonor = zeros(length(donor_id_uni),length(T_cells_tmp_uni));
for i=1:length(donor_id_uni)
    in = find(strcmpi(donor_id_flag_sorted,donor_id_uni{i}));
    per_cluster_indonor(i,:) = histcounts(T_cells_tmp(in),[1:length(T_cells_tmp_uni)+1]);
end
per_cluster_indonor = per_cluster_indonor./repmat(sum(per_cluster_indonor,2),1,length(T_cells_tmp_uni));
% % % % % % % % % % % % 

%% cell types abundances

% [ha, pos] = tight_subplot(4, 4, [0.05,0.04], [0.02,0.02], [0.02,0.02]);
for j=1:length(T_cells_tmp_uni)
    hf2 = figure('color','w','position',[20,20,300,700],'name',clusteruni{j});
    ax1 = axes('position',[0.2,0.4,0.75,0.5]);
    t_cell = cell(6,1);
    for i=[1:4]
        s = unique(sample_sorted(condition_sorted(i,:)==1));
        [~,loc] = ismember(s,sample_uni);
        t = per_cluster_insample(loc,j);
        t_cell{i} = t;
        t_med(i) = mean(t);
        plot(i,100*t,'o','markerfacecolor',['#',colors{j,2}],'markeredgecolor','none','markersize',10); hold on;
        plot(i,100*median(t),'s','markerfacecolor','k','markersize',20,'markeredgecolor','none');
        plot(i,100*mean(t),'o','markerfacecolor','k','markersize',20,'markeredgecolor','none');
    end
    [~,p_ct_pl] = ttest2(t_cell{3},t_cell{4});
    [~,p_ce_pe] = ttest2(t_cell{1},t_cell{2});
    box off
    axis tight
    yl = get(gca,'ylim');
    plot([3,4],yl(2)*0.9*[1,1],'k');
    text(3,yl(2)*0.9,num2str(p_ct_pl,2),'fontsize',10,'VerticalAlignment','bottom');
    plot([1,2],yl(2)*0.9*[1,1],'k');
    text(1,yl(2)*0.9,num2str(p_ce_pe,2),'fontsize',10,'VerticalAlignment','bottom');
    set(ax1,'xlim',[0.5,4.5],'xtick',[1:4],'xticklabel',{'ce','pe','ct','pl'},'fontname','arial','fontsize',12)
    title(regexprep(clusteruni((j)),'_','-'),'fontsize',6)
    ax2 = axes('position',[0.2,0.07,0.75,0.2]);
    m_ratio = [(t_med(2)-t_med(1))/(min(t_med(2),t_med(1)));(t_med(4)-t_med(3))/(min(t_med(4),t_med(3)))];
    hb(1) = bar(1,100*m_ratio(1));hold on;
    hb(2) = bar(2,100*m_ratio(2));
    set(hb(1),'facecolor',211/255*[1,1,1],'barwidth',0.3)
    set(hb(2),'facecolor',[140,137,140]/255,'barwidth',0.3)
    axis tight
    ylabel('diff %');
    set(ax2,'ylim',[-200,200],'xlim',[0.5,2.5],'xtick',[1:2],'xticklabel',{'pe/ce','pl/ct'},'ytick',[-200,0,200],'fontname','arial','fontsize',12)
    eval(['export_fig stromal_abundance_celltype_221120/',clusteruni{j},'_sample_celltype_abundance_per_condition_',date,'.pdf']);
    close(hf2)
end
%donor
for j=1:length(T_cells_tmp_uni)
    hf2 = figure('color','w','position',[20,20,300,700],'name',clusteruni{j});
    ax1 = axes('position',[0.2,0.4,0.75,0.5]);
    t_cell = cell(6,1);
    for i=[1:4]
        s = unique(donor_id_flag_sorted(condition_sorted(i,:)==1));
        [~,loc] = ismember(s,donor_id_uni);
        t = per_cluster_indonor(loc,j);
        t_cell{i} = t;
        t_med(i) = mean(t);
        plot(i,100*t,'o','markerfacecolor',['#',colors{j,2}],'markeredgecolor','none','markersize',10); hold on;
        plot(i,100*median(t),'s','markerfacecolor','k','markersize',20,'markeredgecolor','none');
        plot(i,100*mean(t),'o','markerfacecolor','k','markersize',20,'markeredgecolor','none');
    end
    [~,p_ct_pl] = ttest2(t_cell{3},t_cell{4});
    [~,p_ce_pe] = ttest2(t_cell{1},t_cell{2});
    box off
    axis tight
    yl = get(gca,'ylim');
    plot([3,4],yl(2)*0.9*[1,1],'k');
    text(3,yl(2)*0.9,num2str(p_ct_pl,2),'fontsize',10,'VerticalAlignment','bottom');
    plot([1,2],yl(2)*0.9*[1,1],'k');
    text(1,yl(2)*0.9,num2str(p_ce_pe,2),'fontsize',10,'VerticalAlignment','bottom');
    set(ax1,'xlim',[0.5,4.5],'xtick',[1:4],'xticklabel',{'ce','pe','ct','pl'},'fontname','arial','fontsize',12)
    title(regexprep(clusteruni((j)),'_','-'),'fontsize',6)
    ax2 = axes('position',[0.2,0.07,0.75,0.2]);
    m_ratio = [(t_med(2)-t_med(1))/(min(t_med(2),t_med(1)));(t_med(4)-t_med(3))/(min(t_med(4),t_med(3)))];
    hb(1) = bar(1,100*m_ratio(1));hold on;
    hb(2) = bar(2,100*m_ratio(2));
    set(hb(1),'facecolor',211/255*[1,1,1],'barwidth',0.3)
    set(hb(2),'facecolor',[140,137,140]/255,'barwidth',0.3)
    axis tight
    text(1,100,num2str(m_ratio(1)*100));
    text(2,100,num2str(m_ratio(2)*100));
    ylabel('diff %');
    set(ax2,'ylim',[-200,200],'xlim',[0.5,2.5],'xtick',[1:2],'xticklabel',{'pe/ce','pl/ct'},'ytick',[-200,0,200],'fontname','arial','fontsize',12)
    eval(['export_fig stromal_abundance_celltype_221120/',clusteruni{j},'_donor_celltype_abundance_per_condition_',date,'.pdf']);
    close(hf2)
end

%% stackbar cell type per sample
figure('color','w','position',[20,20,1000,400]);
k = 0;
colorvec = distinguishable_colors(length(T_cells_tmp_uni));
[~,xi] = sort(median(per_cluster_insample));
sample_uni_order = [];
for i=[1,3,2,4]

    s = unique(sample_sorted(condition_sorted(i,:)==1));
    sample_uni_order = [sample_uni_order;s];
    [~,loc] = ismember(s,sample_uni);
    t = per_cluster_insample(loc,xi);
    %     axes(ha(k))
    hb = bar([k+1:k+length(loc)],t,'stacked'); hold on;
    plot([1,1]*k+length(loc)+0.5,[0,1],'-k','linewidth',3)
    k = k+length(loc);
    for jj=1:length(hb)
        set(hb(jj),'FaceColor',colorvec(jj,:));
    end
end
hleg = legend(clusteruni(xi),'interpreter','none');
set(hleg,'position',[0.78,0.3,0.2,0.6],'fontsize',6)
axis tight
set(gca,'position',[0.05,0.1,0.7,0.85],'xtick',[1:length(sample_uni_order)],'xticklabel',sample_uni_order)
if savefig_flag==1
%     savefig(gcf,['immune_step4_median_per_condition_',date,'.fig'])
    eval(['export_fig stromal_step4_stackedbar_celltype_per_sample_',date,'.pdf']);
end
%% markertable

[ind_gr_tmp_mark,cells_bor,gr_center] = markertablefeatures(T_cells_tmp,data_sorted_all,5);
% % % % % % % % % 
datamarkers = data_sorted_all(ind_gr_tmp_mark,:);
datamarkers_cn = cent_norm(log2(datamarkers+1));
% gr_tmp_mark = gr_tmp_mark(xi);
gr_tmp_mark = geneid(ind_gr_tmp_mark);

figure;
set(gcf,'position',[100,100,1400,770],'color','w')
ax1 = axes('position',[0.1,0.025,0.88,0.73]);
imagesc(datamarkers_cn,[prctile(datamarkers_cn(:),1),prctile(datamarkers_cn(:),99)]);
hold on;
linewid =0.5;
bor_color = 'grey11';%'green1';%
for jj=1:length(cells_bor)
    plot(cells_bor(jj)*[1,1]-0.5,[1,length(gr_tmp_mark)],'-','linewidth',linewid,'color',get_RGB(bor_color))
end
set(gca,'xtick',gr_center,'xticklabel',[1:length(gr_center)],'ytick',[1:length(gr_tmp_mark)],'yticklabel',gr_tmp_mark, 'fontsize', 10)
colormap('summer');
freezeColors(gca);
text([cells_bor;cells_bor(end)],linspace(5,length(gr_tmp_mark),length(gr_center)),regexprep(clusteruni(),'_','-'),'color','k')

xist = data_orig_all_sorted(strcmpi('xist',geneid_all),:);
% sample_uni = {'40-1','44-1','48-1','53-1'};
samples_num = false(length(sample_uni),length(sample_sorted));
for i=1:length(sample_uni)
    samples_num(i, strcmpi(sample_sorted,sample_uni{i})) = true;
end

ax2 = axes('position',[0.1,0.7580,0.88,0.01]);
imagesc(~ctrlearly_flag_sorted'); hold on;
grid on
colormap('gray');
freezeColors(gca);
set(gca,'xtick',[],'ytick',1,'yticklabel','Ctrlearly','fontsize',6);

ax3 = axes('position',[0.1,0.7680,0.88,0.01]);
imagesc(~peearly_flag_sorted'); hold on;
grid on
colormap('gray');
freezeColors(gca);
set(gca,'xtick',[],'ytick',1,'yticklabel','PEearly','fontsize',6);

ax4 = axes('position',[0.1,0.7780,0.88,0.01]);
imagesc(~pelate_flag_sorted'); hold on;
grid on
colormap('gray');
freezeColors(gca);
set(gca,'xtick',[],'ytick',1,'yticklabel','PElate','fontsize',6);

ax5 = axes('position',[0.1,0.7860,0.88,0.01]);
imagesc(~ctrl_flag_sorted'); hold on;
grid on
colormap('gray');
freezeColors(gca);
set(gca,'xtick',[],'ytick',1,'yticklabel','ctrl','fontsize',6);

ax6 = axes('position',[0.1,0.7960,0.88,0.01]);
imagesc(~iugr_flag_sorted'); hold on;
grid on
colormap('gray');
freezeColors(gca);
set(gca,'xtick',[],'ytick',1,'yticklabel','IUGR','fontsize',6);

ax7 = axes('position',[0.1,0.8060,0.88,0.01]);
imagesc(~female_flag_sorted'); hold on;
grid on
colormap('gray');
freezeColors(gca);
set(gca,'xtick',[],'ytick',1,'yticklabel','female','fontsize',6);

ax8 = axes('position',[0.1,0.8160,0.88,0.17]);
imagesc(~samples_num); hold on;
grid on
colormap('gray');
freezeColors(gca);
set(gca,'xtick',[],'ytick',[1:length(sample_uni)],'yticklabel',sample_uni,'fontsize',4.8);

ax9 = axes('position',[0.1,0.988,0.88,0.01]);
imagesc(~xist); hold on;
grid on
colormap('gray');
freezeColors(gca);
set(gca,'xtick',[],'ytick',1,'yticklabel','XIST','fontsize',6);

linkaxes([ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8,ax9],'x');

if savefig_flag==1
    savefig(gcf,['stromal_step4_markertable_',date,'.fig']) 
    eval(['export_fig stromal_step4_markertable_',date,'.pdf']);
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% table2 = [cellid_sorted,m2c(T_cells_tmp)];
% saveCellFile(table2,['trophoblasts_step4_cellid_cluster_beforeclean_',date,'.txt']);


%% volcano scatter plot
top_g = 50;
gr1name = 2;
gr2name = 3;
c1 = 1;
gr1 = find((T_cells_tmp == 2 ));
gr2 = find((T_cells_tmp ~= 2 ));
% in = find(mean(data_orig_all_sorted(:,gr1)>0,2)>0.03 | mean(data_orig_all_sorted(:,gr2)>0,2)<0.7);
in = find(mean(data_orig_all_sorted>0,2)>0.05 & mean(data_orig_all_sorted>0,2)<0.5 & ~ismember(geneid_all(:,1),sex_genes));
ptt = zeros(length(in),1);
for s=1:length(in)
    if mod(s,100)==0
        s
    end
    ptt(s) = ranksum(data_orig_all_sorted(in(s),gr1),data_orig_all_sorted(in(s),gr2),'tail','both');
end
ptt(isnan(ptt)) = 1;
ptt(ptt<1) = qval_from_pval(ptt(ptt<1));
ptt(ptt<1e-300) = 1e-300;
x1 = mean(log2(data_orig_all_sorted(:,gr1)+1),2);
x2 = mean(log2(data_orig_all_sorted(:,gr2)+1),2);
d = x1-x2 ;
figure('position',[200,200,1400,580],'color','w');
[ha, pos] = tight_subplot(1, 3, [0.05,0.05], [0.1,0.05], [0.05,0.05]);
axes(ha(1))
plot(d(in),-log10(ptt),'.'); hold on;
[~,xi] = sort(ptt);
plot(d(in(xi(1:100))),-log10(ptt(xi(1:100))),'or');
text(d(in(xi(1:100))),-log10(ptt(xi(1:100))),geneid_all(in(xi(1:100))))
xlabel([gr1name,'-',gr2name]);
ylabel(['-log10(q)']);
grid on;
title(['comparing ',gr1name,' to ',gr2name,' in cluster ',num2str(c1),', trophoblast'])
[~,xi] = sort(d);
axes(ha(2))
plot(x1, x2, '.');hold on;
xmax = max(x1);
plot([0,xmax],[0,xmax],'-k'); grid on
plot([0,xmax],[0,xmax]+1,'--k'); grid on
plot([1,xmax],[0,xmax-1],'--k'); grid on
plot(x1(xi(1:top_g)), x2(xi(1:top_g)),'.r'); hold on;
text(x1(xi(1:top_g)), x2(xi(1:top_g)),geneid_all(xi(1:top_g)),'fontsize',6);
xi = flipud(xi);
plot(x1(xi(1:top_g)), x2(xi(1:top_g)),'.r'); hold on;
text(x1(xi(1:top_g)), x2(xi(1:top_g)),geneid_all(xi(1:top_g)),'fontsize',6);
% title(['v3-nexgem (average)'])
xlabel(gr1name)
ylabel(gr2name)
axis tight

top_g = 50;
x1 = mean(data_orig_all_sorted(:,gr1)>0,2);
x2 = mean(data_orig_all_sorted(:,gr2)>0,2);
d = x1-x2 ;
[~,xi] = sort(d);
axes(ha(3))
plot(x1, x2, '.');hold on;
xmax = max(x1);
plot([0,xmax],[0,xmax],'-k'); grid on
plot([0,xmax],[0,xmax]+0.4,'--k'); grid on
plot([0.4,xmax],[0,xmax-0.4],'--k'); grid on
plot(x1(xi(1:top_g)), x2(xi(1:top_g)),'.r'); hold on;
text(x1(xi(1:top_g)), x2(xi(1:top_g)),geneid_all(xi(1:top_g)),'fontsize',6);
xi = flipud(xi);
plot(x1(xi(1:top_g)), x2(xi(1:top_g)),'.r'); hold on;
text(x1(xi(1:top_g)), x2(xi(1:top_g)),geneid_all(xi(1:top_g)),'fontsize',6);
% title(['v3-nexgem (precentage)'])
xlabel(gr1name)
ylabel(gr2name)
axis tight
% eval(['export_fig scatter_volcano_',gr1name,'_',gr2name,'_tropho_cluster',num2str(c1),'_',date,'.pdf']);
% % % % % % % 
 % % % % % % % % % % % % % % % % % % % 
 %%

mean_mat = zeros(length(geneid_all),length(T_cells_tmp_uni)*6);
median_mat = zeros(length(geneid_all),length(T_cells_tmp_uni)*6);
fracpos_mat = zeros(length(geneid_all),length(T_cells_tmp_uni)*6);
p75_mat = zeros(length(geneid_all),length(T_cells_tmp_uni)*6);
condgrsize = zeros(length(T_cells_tmp_uni)*6,1);
i= 0;
for k=1:length(T_cells_tmp_uni)
    for jj=1:6
        i = i+1;
        gr1 = find(T_cells_tmp==k & condition_sorted(jj,:)'==1);
        condgrsize(i) = length(gr1);
        mean_mat(:,i) = mean(data_orig_all_sorted(:,gr1),2);
        median_mat(:,i) = median(data_orig_all_sorted(:,gr1),2);
        fracpos_mat(:,i) = mean(data_orig_all_sorted(:,gr1)>0,2);
        p75_mat(:,i) = prctile(data_orig_all_sorted(:,gr1),75,2);
    end
end


%% gene expression on tsne
list = {'comp'};
figure;
set(gcf,'color','w','position',[20,20,1100,960])
[ha, pos] = tight_subplot(1, 1, [0.02,0.02], [0.02,0.02], [0.02,0.02]);
for i=1:length(list)
    genePlot = list{i};
    markergene = (data_orig_all_sorted(strcmpi(geneid_all,genePlot),:));
    inpos = markergene>0;
    tmpthlow = prctile(markergene(markergene>0),1);
    tmpthhigh = prctile(markergene(markergene>0),99);
    if tmpthlow==tmpthhigh
        tmpthlow = 0;
    end
    markergene(markergene>tmpthhigh) = tmpthhigh;
    markergene(markergene<tmpthlow) = tmpthlow;
    c_rgb = [1,0,0];rand([1,3]);
    markergene_color = [interp1([min(markergene),max(markergene)],[0.7,c_rgb(1)],markergene'),...
        interp1([min(markergene),max(markergene)],[0.7,c_rgb(2)],markergene')...
        ,interp1([min(markergene),max(markergene)],[0.7,c_rgb(3)],markergene')];
    axes(ha(i));
    scatter(mapped_xy(~inpos,1),mapped_xy(~inpos,2),20,markergene_color(~inpos,:),'.'); hold on;
    scatter(mapped_xy(inpos,1),mapped_xy(inpos,2),20,markergene_color(inpos,:),'.'); hold on;
    set(gca,'xlim',[-150,150],'ylim',[-150,150])
    title(genePlot);
    axis tight
    axis equal
    axis off
end

%% markers violins
typeplot = 1:11;
list = {'OGN','TOP2A'};
% list = {'C7','HGF','LTN1','PPM1N','DIO3OS','STEAP4','CXCL2','HSPA1A','COX4I2','MX2','GALNT15','MT1M'};
fnorm = normpdf(norminv([0.001:0.01:0.999]));
for kk = 1:length(list)
    hf = figure('color','w','position',[20,20,520*11/9,70*11/9],'name',list{kk});
    % [ha, pos] = tight_subplot(length(list{kk}), 1, [0.01,0.01], [0.01,0.01], [0.09,0.01]);
    ha = axes('position',[0,0,1,1]);
    logflag = 0;
    % for jjj=1:length(list)
    %     gn = list{jjj};
    gn = list{kk};
    g = find(strcmpi(gn,geneid_all));
    p99 = 0;%prctile(data_orig_all_sorted(g,:),100);
    %     axes(ha(jjj))
    p90 = max(ceil(prctile(data_orig_all_sorted(g,:),90)),1);
    t_ed = zeros(length(typeplot),1);
    t_av = zeros(length(typeplot),1);
    t_75 = zeros(length(typeplot),1);
    for k=1:length(typeplot)
        %         for k=1:length(T_cells_tmp_uni)
        c=k;
        gr2 = find(T_cells_tmp(:,1)==typeplot(c));
        y = (data_orig_all_sorted(g,gr2));
        if logflag ==1
            y = y+1;
            yscale = 'log';
        else
            yscale = 'linear';
        end
        if length(y)>5
%             [f,xi] = ksdensity(y);
%             fi = interp1(xi,f,y);
%             fi = fi/max(fi);

            xi = prctile(y,[0.1:1:99.9]);
            [xi,ia] = unique(xi);
            if length(xi)>1
            fi = interp1(xi,fnorm(ia),y);
            fi = fi/max(fi);
            else
                fi = zeros(size(y));
            end
            plot(k + fi'.*(0.5*rand(length(gr2),1)-0.25), 0.5*rand(length(gr2),1)-0.1+y'....
                ,'.','color',['#',colors{typeplot(k),2}],'markersize',7); hold on;
        else
            plot(k+0.7*rand(length(gr2),1)-0.1,0.5*rand(length(gr2),1)-0.1+data_orig_all_sorted(g,gr2)'....
                ,'.','color',['#',colors{typeplot(k),2}],'markersize',7); hold on;
        end
        t_ed(k) = median(data_orig_all_sorted(g,gr2));
        t_av(k) = mean(data_orig_all_sorted(g,gr2));
        t_75(k) = prctile(data_orig_all_sorted(g,gr2),75);
        t_95(k) = prctile(data_orig_all_sorted(g,gr2),95);
    end
    axis tight
    %     plot(t_ed,'sk');

    %         plot(t_av,'or');
    plot(t_av,'or','MarkerFaceColor','k');
    p99 = max(t_95);
    yl = get(gca,'ylim');
    yt = get(gca,'ytick');
    set(gca,'xlim',[0,1+length(T_cells_tmp_uni)],'xtick',[1:length(T_cells_tmp_uni)]...
        ,'xticklabel',cell(size(T_cells_tmp_uni)),'ylim',[-1,p99+1],'YScale',yscale,'ytick',yt,'yticklabel',cell(length(yt),1))
    axis off

    % end
    % set(gcf, 'units','normalized','outerposition',[0.1 0.1 0.3 0.11]);
    % eval(['export_fig markersViolin_Immune_Lymphoid_',char(list),'_',num2str(p99),'_',date,'.pdf']);
%     save2png(['TB_marker_violin/',list{kk},'_MarkersViolin_TB','_',num2str(round(p99)),'_',date],gcf,1000)
    save2png(['MarkersViolin_Stromal_',list{kk},'_',num2str(round(p99)),'_',date],gcf,1000)
    close(hf)
end



