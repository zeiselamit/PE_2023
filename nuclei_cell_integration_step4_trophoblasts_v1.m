tic
clear all
close all
addpath('/data/matlab_functions/')
terminate(pyenv)
pyenv("ExecutionMode","OutOfProcess")

savefig_flag = 0;
savefig_pdf = 0;

load PE-nuclei_QC800_400_afterloading_10-Apr-2023.mat
nucleiidclsuter_tb = loadCellFile('trophoblasts_step3_cellid_cluster_beforeclean_16-Apr-2023.txt');

% cellidclsuter_tb = loadCellFile('cellidCluster_step2_GABA_FC_07-Oct-2021.txt');
% cellidclsuter_tb = [cellidclsuter_tb, cellfun(@(x,y) ['gaba_',num2str(x),'_',num2str(y)],cellidclsuter_tb(:,2),cellidclsuter_tb(:,3),'UniformOutput',0)]; 

% typenum = cell2mat(cellid_type(:,2));
% clean_list = loadCellFile('cleaning_trophoblasts_final_211013.txt');
% clean_list = cell2mat(clean_list(3:end,:));
% for i=1:length(clean_list(:,1))
%     if clean_list(i,2)==0
%         typenum(typenum==clean_list(i,1)) = 0;
%     end
% end
% cellid_type(typenum==0,:) = [];% selecting only immune
% [~,loc] = ismember(cellid_type(:,1),cellid);
% sum(loc>0)

tropho_clustername = loadCellFile('nuceli_cluster_name_trophoblast_230423.txt');
tropho_clustername(:,2) = cellfun(@(x) ['TB_',x], tropho_clustername(:,2),'UniformOutput',0);
tropho_clustername(:,2) = regexprep(tropho_clustername(:,2),'TB_EXC','EXC');
tropho_clustername = tropho_clustername(2:end,:);
nucleiidclsuter_tb = [nucleiidclsuter_tb,cell(length(nucleiidclsuter_tb),1)];
c = cell2mat(nucleiidclsuter_tb(:,2));
for i=1:length(tropho_clustername)
    ind = find(c==tropho_clustername{i});
    nucleiidclsuter_tb(ind,3) = repmat(tropho_clustername(i,2),length(ind),1);
end


cellid_clusters = cell(length(cellid),3);
[~,loc]= ismember(nucleiidclsuter_tb(:,1),cellid);
cellid_clusters(loc(loc>0),:) = nucleiidclsuter_tb;

rmv = find(strcmpi(cellid_clusters(:,3),'EXC'));
cellid_clusters(rmv,:) = cell(length(rmv),3);


validcells = cell2mat(cellfun(@(x) ~isempty(x), cellid_clusters(:,1),'UniformOutput',0));
sum(validcells)

tot_mol = sum(data);
tot_mol(tot_mol>5e4) = 5e4;
tot_genes = sum(data>0);
loc = find(validcells);

 
n_data = data(:,loc);
n_cellid = cellid(loc);
n_sample = sample(loc);
n_tot_mol = tot_mol(loc);
n_v3_flag = v3_flag(loc);
n_ctrl_flag = ctrl_flag(loc);
n_peearly_flag = peearly_flag(loc);
n_pelate_flag = pelate_flag(loc);
n_female_flag = female_flag(loc);
n_iugr_flag = iugr_flag(loc);
n_ctrlearly_flag = ctrlearly_flag(loc);
n_csection_flag = csection_flag(loc);
n_vaginalbirth_flag = vaginalbirth_flag(loc);
n_induction_flag = induction_flag (loc);
n_noinduction_flag = noinduction_flag(loc);
n_magnesium_flag = magnesium_flag(loc);
n_spinal_flag = spinal_flag(loc);
n_epidural_flag = epidural_flag(loc);
n_generalanesthesia_flag = generalanesthesia_flag(loc);
n_week_flag = week_flag(loc);
n_weight_flag = weight_flag(loc);
n_weightprct_flag = weightprct_flag(loc);
n_donorage_flag = donorage_flag(loc);
n_donor_id_flag = donor_id_flag(loc);
n_donor_id_flag = append('n_',n_donor_id_flag);
n_cellid_clusters = cellid_clusters(loc,:);

%% single cell

load PE_QC800_400_afterloading_30-Jan-2022.mat
cellidclsuter_tb = loadCellFile('trophoblasts_step3_cellid_cluster_beforeclean_17-Oct-2021.txt');

tropho_clustername = loadCellFile('sc_cluster_name_trophoblast_221120.txt');
tropho_clustername(:,2) = cellfun(@(x) ['TB_',x], tropho_clustername(:,2),'UniformOutput',0);
tropho_clustername(:,2) = regexprep(tropho_clustername(:,2),'TB_EXC','EXC');
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


 
c_data = data(:,loc);
c_cellid = cellid(loc);
c_sample = sample(loc);
c_tot_mol = tot_mol(loc);
c_v3_flag = v3_flag(loc);
c_ctrl_flag = ctrl_flag(loc);
c_peearly_flag = peearly_flag(loc);
c_pelate_flag = pelate_flag(loc);
c_female_flag = female_flag(loc);
c_iugr_flag = iugr_flag(loc);
c_ctrlearly_flag = ctrlearly_flag(loc);
c_csection_flag = csection_flag(loc);
c_vaginalbirth_flag = vaginalbirth_flag(loc);
c_induction_flag = induction_flag (loc);
c_noinduction_flag = noinduction_flag(loc);
c_magnesium_flag = magnesium_flag(loc);
c_spinal_flag = spinal_flag(loc);
c_epidural_flag = epidural_flag(loc);
c_generalanesthesia_flag = generalanesthesia_flag(loc);
c_week_flag = week_flag(loc);
c_weight_flag = weight_flag(loc);
c_weightprct_flag = weightprct_flag(loc);
c_donorage_flag = donorage_flag(loc);
c_donor_id_flag = donor_id_flag(loc);
c_donor_id_flag = append('c_',c_donor_id_flag);
c_cellid_clusters = cellid_clusters(loc,:);

%% integration

data = [n_data,c_data];
cellid = [n_cellid;c_cellid];
sample = [n_sample;c_sample];
tot_mol = [n_tot_mol,c_tot_mol];
v3_flag = [n_v3_flag;c_v3_flag];
ctrl_flag = [n_ctrl_flag;c_ctrl_flag];
peearly_flag = [n_peearly_flag;c_peearly_flag];
pelate_flag = [n_pelate_flag;c_pelate_flag];
female_flag = [n_female_flag;c_female_flag];
iugr_flag = [n_iugr_flag;c_iugr_flag];
ctrlearly_flag = [n_ctrlearly_flag;c_ctrlearly_flag];
csection_flag = [n_csection_flag;c_csection_flag];
vaginalbirth_flag = [n_vaginalbirth_flag;c_vaginalbirth_flag];
induction_flag = [n_induction_flag;c_induction_flag];
noinduction_flag = [n_noinduction_flag;c_noinduction_flag];
magnesium_flag = [n_magnesium_flag;c_magnesium_flag];
spinal_flag = [n_spinal_flag;c_spinal_flag];
epidural_flag = [n_epidural_flag;c_epidural_flag];
generalanesthesia_flag = [n_generalanesthesia_flag;c_generalanesthesia_flag];
week_flag = [n_week_flag;c_week_flag];
weight_flag = [n_weight_flag;c_weight_flag];
weightprct_flag = [n_weightprct_flag;c_weightprct_flag];
donorage_flag = [n_donorage_flag;c_donorage_flag];
donor_id_flag = [n_donor_id_flag;c_donor_id_flag];
cellid_clusters = [n_cellid_clusters;c_cellid_clusters];
%%

% validcells = false(size(cellid));
% validcells(loc) = true;
sample_uni = {'132-1','134-2','134-3','135-1','135-2'...
    ,'48-1','53-1','61-1','63-1','63-2','67-1','67-2','68-1'...
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

data_orig_all = data;
geneid_all = geneid;
data = data(in(corr_filt),:);
geneid = geneid(in(corr_filt));
% % % % % % % % % % % % % % % % % 
colors = loadCellFile_turbo('cn_color_celltypes_ordered_HEX_230430.txt',1);
celltypes_ordered_hex = colors(contains(colors(:,1),'TB_'),:);
colors = celltypes_ordered_hex(:,1:2);
clusteruni = celltypes_ordered_hex(:,1);
T_cells_tmp = zeros(length(cellid),1);
clustersize = zeros(length(clusteruni),1);
for i=1:length(clusteruni)
    ind = find(strcmpi(cellid_clusters(:,3),clusteruni{i}));
    T_cells_tmp(ind) = i;
    clustersize(i) = length(ind);
end
T_cells_tmp_uni = unique(T_cells_tmp);

% % % % % % % % % % % % % % % % % % % 
%% pre-harmony
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
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
% initial_dims = 30;
ph_prj = prj(:,1:initial_dims);
init = ph_prj(:,1:2)/std(ph_prj(:,1))*1e-4;
init = init-repmat(mean(init),length(init),1);

perplexity = 200;
options = statset('MaxIter',1000);
ph_mapped_xy = tsne(ph_prj,'Algorithm','barneshut','Distance','correlation','NumDimensions',2,'NumPCAComponents',0,.....
    'Perplexity',perplexity,'Standardize',true,'InitialY',init,'LearnRate',length(cellid)/12,'Theta',0.5,'Verbose',1,'Options',options,'Exaggeration',20);
ph_prj = ph_prj';
% % % % % % % % % % % % 
figure;
set(gcf,'color','w','position',[20,20,900,800])
for i=1:length(T_cells_tmp_uni)
    ii=find(T_cells_tmp==i); h=plot(ph_mapped_xy(ii,1),ph_mapped_xy(ii,2),'.','color',hex2rgb(colors{i,2}),'markersize',4); hold on;
end
% plot(ph_mapped_xy([1:length(n_donor_id_flag)],1),ph_mapped_xy([1:length(n_donor_id_flag)],2),'^r'); hold on
% plot(ph_mapped_xy([length(n_donor_id_flag)+1:end],1),ph_mapped_xy([length(n_donor_id_flag)+1:end],2),'.b');
axis tight;
xl = get(gca,'xlim');
yl = get(gca,'ylim');
rr = 0.05*diff(xl);
for i=1:length(T_cells_tmp_uni)
    if i>=0
        in = T_cells_tmp==i;
        th = (rand*360)*pi/180;
        xo = rr*cos(th); yo = rr*sin(th);
        plot([median(ph_mapped_xy(in,1)),median(ph_mapped_xy(in,1))+xo] , ...
            [median(ph_mapped_xy(in,2)),median(ph_mapped_xy(in,2))+yo],'-','color',0.7*[1,1,1]);
        ht = text(median(ph_mapped_xy(in,1))+xo,median(ph_mapped_xy(in,2))+yo,regexprep(clusteruni{i},'_','-'));
        set(ht,'fontsize',8,'FontName','arial');%,'BackgroundColor',0.8*[1,1,1]
    end
end
axis tight;
axis equal
axis off
title(['#C=',num2str(max(T_cells_tmp)),', perplexity=',num2str(perplexity)],'fontsize',8);
if savefig_flag==1
    savefig(gcf,['cell_nuclei_preHarmony_tsne_final_TB_by_cluster_',date,'.fig'])    
    eval(['export_fig cell_nuclei_preHarmony_tsne_final_TB_by_cluster_',date,'.pdf']);
end
% % % % % % % % % % % % % % % % % % % % % test harmony py
%%
s = ph_prj;
batchid = donor_id_flag;
batchid([1:length(n_donor_id_flag)],2) = {'nuclei'}; %nuclei flag = 1
batchid([length(n_donor_id_flag)+1:end],2) = {'cell'}; % cell flag = 2
usepylib = 1;
[sout]=harmonypy(s,batchid(:,2),usepylib);
prj = sout;
init = prj(:,1:2)/std(prj(:,1))*1e-4;
init = init-repmat(mean(init),length(init),1);
h_prj = prj';
perplexity = 200;
options = statset('MaxIter',1000);
mapped_xy = tsne(prj,'Algorithm','barneshut','Distance','correlation','NumDimensions',2,'NumPCAComponents',0,.....
    'Perplexity',perplexity,'Standardize',true,'InitialY',init,'LearnRate',length(cellid)/12,'Theta',0.5,'Verbose',1,'Options',options,'Exaggeration',20);
toc



%% integrated tsne
figure;
set(gcf,'color','w','position',[20,20,900,800])
for i=1:length(T_cells_tmp_uni)
    ii=find(T_cells_tmp==i); h=plot(mapped_xy(ii,1),mapped_xy(ii,2),'.','color',hex2rgb(colors{i,2}),'markersize',4); hold on;
end
% plot(mapped_xy([1:length(n_donor_id_flag)],1),mapped_xy([1:length(n_donor_id_flag)],2),'^r'); hold on
% plot(mapped_xy([length(n_donor_id_flag)+1:end],1),mapped_xy([length(n_donor_id_flag)+1:end],2),'.b');
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
    savefig(gcf,['cell_nuclei_tsne_final_TB_by_cluster_',date,'.fig'])    
    eval(['export_fig cell_nuclei_tsne_final_TB_by_cluster_',date,'.pdf']);
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
ph_mapped_xy = ph_mapped_xy(xi,:);
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
h_prj = h_prj(:,xi);
ph_prj = ph_prj(:,xi);

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
h_prj = h_prj(:,xi);
ph_prj = ph_prj(:,xi);
mapped_xy = mapped_xy(xi,:);
ph_mapped_xy = ph_mapped_xy(xi,:);
% % % % % % % % % % % % % % % % 
meangr_mat = zeros(length(moldata(:,1)),length(T_cells_tmp_uni));
meangr_prj_harmony = zeros(length(h_prj(:,1)),length(T_cells_tmp_uni));
clust_cent = zeros(length(T_cells_tmp_uni),2);
for jjj=1:length(T_cells_tmp_uni)
    jjj
    meangr_mat(:,jjj) = mean(log2(data_sorted_all(:,T_cells_tmp==T_cells_tmp_uni(jjj))+1),2);
    meangr_prj_harmony(:,jjj) = mean(h_prj(:,T_cells_tmp==T_cells_tmp_uni(jjj)),2);
    clust_cent(jjj,:) = [median(mapped_xy(T_cells_tmp==T_cells_tmp_uni(jjj),1)),median(mapped_xy(T_cells_tmp==T_cells_tmp_uni(jjj),2))];
end
meangr_mat1 = meangr_mat;
% meangr_mat1(loc,:) = [];

% meangr_mat1 = cent_norm(meangr_mat1(:,leaforder));
initial_dims = 10;
[prj,m,D,V,Q] = pca_wis(meangr_mat1',initial_dims);
Zpca = linkage(prj,'ward','correlation');
Dpca = pdist(prj,'correlation');
leaforder_pca = optimalleaforder(Zpca,Dpca);
%%
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

%% linkage on integrated dimensions
Zpca_h = linkage(meangr_prj_harmony','ward','correlation');
Dpca_h = pdist(meangr_prj_harmony','correlation');

leaforder_pca_h = optimalleaforder(Zpca_h,Dpca_h);
figure;
set(gcf,'position',[100,100,1000,1000],'color','w')
% h1 = axes('position',[0.03,0.03,0.23,0.93]);
h1 = axes('position',[0.03,0.1,0.2,0.75]);
hden = dendrogram(Zpca_h,length(leaforder_pca_h),'Reorder',leaforder_pca_h,'Orientation','left');
axis off
set(gca,'ylim',[0.5,length(leaforder_pca_h)+0.5])
h2 = axes('position',[0.3,0.1,0.7,0.75]);
% h2 = axes('position',[0.38,0.03,0.98,0.93]);
x=squareform(1-Dpca_h)+eye(length(leaforder_pca_h)); imagesc(x(leaforder_pca_h,leaforder_pca_h));
colormap('summer')
set(gca,'ytick',[1:length(clusteruni)],'YTickLabel',regexprep(clusteruni(leaforder_pca_h),'_','-'),'xtick',[],'fontsize',8,'ydir','normal')
linkaxes([h1,h2],'y');
if savefig_flag==1
    savefig(gcf,['cell_nuclei_linkage_on_integrated_dimension_TB_by_cluster_',date,'.fig'])    
    eval(['export_fig cell_nuclei_linkage_on_integrated_dimension_TB_by_cluster_',date,'.pdf']);
end

figure;
set(gcf,'position',[100,100,1000,1000],'color','w')
imagesc(x([9:end],[1:8]));
colormap('summer')
set(gca,'xtick',[1:8],'XTickLabel',regexprep(clusteruni([1:8]),'_','-'),'ytick',[1:5],'YTickLabel',regexprep(clusteruni([9:13]),'_','-'),'fontsize',8,'ydir','normal')
linkaxes([h1,h2],'y');

if savefig_flag==1
    savefig(gcf,['cell_nuclei_correlation_on_integrated_dimension_TB_by_cluster_',date,'.fig'])    
    eval(['export_fig cell_nuclei_correlation_on_integrated_dimension_TB_by_cluster_',date,'.pdf']);
end
%%
leaforder = leaforder_pca_h;
clusteruni = clusteruni(leaforder_pca_h);

Zpca_post = linkage(prj(leaforder_pca_h,:),'ward','correlation');

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
mapped_xy = mapped_xy(xi,:);
ph_mapped_xy = ph_mapped_xy(xi,:);

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
    savefig(gcf,['cell_nuclei_trophoblasts_step4_tsne_by_sample_perplexity_',num2str(perplexity),'_',date,'.fig'])
end

% % % % % % % % % % % % % % 
[table1, table2] = dendrogram_split_markers(cellfun(@(x) num2str(x), m2c(T_cells_tmp_uni),'UniformOutput',0).....
    ,cellfun(@(x) num2str(x), m2c(T_cells_tmp(:,1)),'UniformOutput',0),Zpca_post,data_sorted_all,geneid);
saveCellFile(table1,['FC_final_cell_nuclei_step4_TB_dendrogram_junction_split_markers_',date,'.txt']);
saveCellFile(table2,['FC_final_cell_nuclei_step4_TB_dendrogram_junction_split_markers_by_average_',date,'.txt']);
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
%% plot by different flags
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
f = ctrlearly_flag_sorted;
plot(mapped_xy(:,1),mapped_xy(:,2),'.','color',0.6*[1,1,1],'markersize',ms); hold on;
plot(mapped_xy(f,1),mapped_xy(f,2),'.','color',[1,0,0],'markersize',ms); hold on;
axis tight
axis equal
axis off
title('Early control')


%%

condition_sorted = zeros(6,length(peearly_flag_sorted));
condition_sorted(1,ctrlearly_flag_sorted) = 1;
condition_sorted(2,peearly_flag_sorted) = 1;
condition_sorted(3,ctrl_flag_sorted) = 1;
condition_sorted(4,pelate_flag_sorted) = 1;
condition_sorted(5,iugr_flag_sorted) = 1;
condition_sorted(6,(peearly_flag_sorted | pelate_flag_sorted) & ~iugr_flag_sorted) = 1;
donor_id_uni = unique(donor_id_flag_sorted);
% calculate the fraction of each cell type per sample or donor
per_cluster_insample = zeros(length(sample_uni),length(T_cells_tmp_uni));
for i=1:length(sample_uni)
    in = find(strcmpi(sample_sorted,sample_uni{i}));
    per_cluster_insample(i,:) = histcounts(T_cells_tmp(in),[1:length(T_cells_tmp_uni)+1]);
end
per_cluster_insample = per_cluster_insample./repmat(sum(per_cluster_insample,2),1,length(T_cells_tmp_uni));
condition_sorted_uni = unique(condition_sorted);

sn_per_cluster_insample = zeros(5,length(T_cells_tmp_uni));
sc_per_cluster_insample = zeros(31,length(T_cells_tmp_uni));
for i=1:length(sample_uni)
    if i <= 5
        in = find(strcmpi(sample_sorted,sample_uni{i}));
        sn_per_cluster_insample(i,:) = histcounts(T_cells_tmp(in),[1:length(T_cells_tmp_uni)+1]);
    else
        in = find(strcmpi(sample_sorted,sample_uni{i}));
        sc_per_cluster_insample(i-5,:) = histcounts(T_cells_tmp(in),[1:length(T_cells_tmp_uni)+1]);
    end
end


per_cluster_indonor = zeros(length(donor_id_uni),length(T_cells_tmp_uni));
for i=1:length(donor_id_uni)
    in = find(strcmpi(donor_id_flag_sorted,donor_id_uni{i}));
    per_cluster_indonor(i,:) = histcounts(T_cells_tmp(in),[1:length(T_cells_tmp_uni)+1]);
end
per_cluster_indonor = per_cluster_indonor./repmat(sum(per_cluster_indonor,2),1,length(T_cells_tmp_uni));
% % % % % % % % % % % % 




[ind_gr_tmp_mark,cells_bor,gr_center] = markertablefeatures(T_cells_tmp,data_sorted_all,5);

%% markertable
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
    savefig(gcf,['cell_nuclei_trophoblats_step4_markertable_',date,'.fig']) 
    eval(['export_fig cell_nuclei_trophoblats_step4_markertable_',date,'.pdf']);
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% table2 = [cellid_sorted,m2c(T_cells_tmp)];
% saveCellFile(table2,['trophoblasts_step4_cellid_cluster_beforeclean_',date,'.txt']);

% % % % % % % % % % % % % % 

%% volcano scatter plot
top_g = 50;
gr1name = 'sc-pre-fused-SCT';
gr2name = 'sn--pre-fused-SCT';
% j = 4;
gr1 = find(T_cells_tmp == 7);
gr2 = find(T_cells_tmp == 6);
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
ptmp = ptt;
ptmp(d(in)>0) = 1;
[~,xi1] = sort(ptmp);
ptmp = ptt;
ptmp(d(in)<0) = 1;
[~,xi2] = sort(ptmp);
xi = [(xi1(1:20));(xi2(1:20))];
plot(d(in(xi)),-log10(ptt(xi)),'.r');
v1 = d(in(xi));
v2 = -log10(ptt(xi));
e = 1;
for i=1:length(v1)
    th = atan(v1(i)/v2(i));
    plot([v1(i),v1(i)+e*tan(th)],[v2(i),v2(i)+e*tan(th)],'-k');
    text(v1(i)+e*tan(th),v2(i)+e*tan(th),geneid_all(in(xi(i))));
end
% % text(d(in(xi(1:100))),-log10(ptt(xi(1:100))),geneid_all(in(xi(1:100))))
xlabel([gr1name,'-',gr2name]);
ylabel(['-log10(q)']);
grid on;
axis tight
title(['comparing ',gr1name,' to ',gr2name])
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
if savefig_flag==1
%     eval(['export_fig scatter_volcano_',gr1name,'_',gr2name,'_',clusteruni((j)),'_',date,'.pdf']);
    eval(['export_fig cell_nuclei_scatter_volcano_',gr1name,'_',gr2name,'_',date,'.pdf']);
end
%%


% % % % % % % % % % % % % % % %
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

%% tsne marker genes
list = {'NOTUM','DLK1','HSPB1','CD63','VIM','ACTB','TMSB10','CYP19A1','TNNI2','PAPPA','KISS1','PSG2','PSG9','ERVFRD-1','CGA','CGB','DLK1','SPARC','TMSB10','TFPI2','ADAM12'};
figure;
set(gcf,'color','w','position',[20,20,1100,960])
[ha, pos] = tight_subplot(5, 4, [0.02,0.02], [0.02,0.02], [0.02,0.02]);
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
%%
% list = {'XAGE2','XAGE3','ITGA5','EBI3','NOTUM','TPM1','PRG2','IL1R2',....
%     'HAPLN3','ACAN','THSD7A','CYP19A1','PSG3','PSG6'...
%     ,'ERVFRD-1','NUPR1','FAM118A','DDIT3','IFIT2','IFIT3','HMGB2','TOP2A'};
list = {'XAGE2','EGFR','TNNI2','HPGD','PRG2','FLT1','PAPPA','PSG3','KISS1',....
    'CYP19A1','ERVFRD-1','IFIT2','TOP2A','TPM1','FAR2','TENM3'};
% list = {'FLT1','KISS1','PAPPA','HPGD','TNNI2'};
typeplot = 1:13;
fnorm = normpdf(norminv([0.001:0.01:0.999]));
for kk = 1:length(list)
    hf = figure('color','w','position',[20,20,520*13/9,70],'name',list{kk});
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
    %     if jjj>length(list)
    %         set(gca,'xtick',[1:length(T_cells_tmp_uni)],'ylim',[-0.5,p99+1],'YScale',yscale,'ytick',yt,'yticklabel',cell(length(yt),1));
    %     else
    set(gca,'xlim',[0,1+length(T_cells_tmp_uni)],'xtick',[1:length(T_cells_tmp_uni)],'xticklabel',cell(size(T_cells_tmp_uni)),'ylim',[-1,p99+1],'YScale',yscale,'ytick',yt,'yticklabel',cell(length(yt),1))
    %     end
    %     text(0.5,p99/2,[gn],'horizontalalignment','right','verticalalignment','middle')
    %     box off
    axis off

    % end
    % set(gcf, 'units','normalized','outerposition',[0.1 0.1 0.3 0.11]);
    % eval(['export_fig markersViolin_Immune_Lymphoid_',char(list),'_',num2str(p99),'_',date,'.pdf']);
    %     save2png(['TB_marker_violin/',list{kk},'_MarkersViolin_TB','_',num2str(round(p99)),'_',date],gcf,1000)
    if savefig_flag==1
        save2png(['MarkersViolin_cell_nuclei_TB_',list{kk},'_',num2str(round(p99)),'_',date],gcf,1000)
    end
    close(hf)
end
%% 2 genes correlation scatter plot
typeplot = 1:13;
g1 = {'FLT1'};
g2 = {'PGF'};
% g1_exp = log2(1+data(strcmpi(geneid,g1),:));
% g2_exp = log2(1+data(strcmpi(geneid,g2),:));
g1_exp = data_orig_all_sorted(strcmpi(geneid_all,g1),:);
g2_exp = data_orig_all_sorted(strcmpi(geneid_all,g2),:);
for k = 1:length(typeplot)
    gr1 = find(T_cells_tmp==typeplot(k) & ctrlearly_flag_sorted == 1);
    gr2 = find(T_cells_tmp==typeplot(k) & peearly_flag_sorted == 1);
%     gr3 = find(T_cells_tmp==typeplot(k) & ctrl_flag == 1);
%     gr4 = find(T_cells_tmp==typeplot(k) & pelate_flag == 1);

   hf = figure('color','w','position',[100,20,675,400],'name',clusteruni{typeplot(k)});
%     [ha, pos] = tight_subplot(1, 4, [0.08,0.08], [0.04,0.04], [0.05,0.08]);
%     axes(ha(1))
    h1 = axes('position',[0.075,0.1,0.4,0.7]);
%     plot(g1_exp(gr1),g2_exp(gr1),'.')
% %     b = robustfit(g1_exp(gr1),g2_exp(gr1));
%     b = polyfit(g1_exp(gr1),g2_exp(gr1),1);
%     hold on
%     mn = min(g1_exp(gr1));
%     mx = max(g1_exp(gr1));
%     plot([mn,mx],b(2)+b(1)*[mn,mx],'-k')
%     [r,pv] = corr(g1_exp(gr1)',g2_exp(gr1)');
%     title([string(['CE, r=',num2str(r)]),['-log(pv) = ',num2str(-log10(pv))],string(['N = ',num2str(length(g1_exp(gr1)))])]);
%     ylabel(g1)
%     xlabel(g2)
    g1_gr1 = log2(1+g1_exp(gr1)+rand(size(g1_exp(gr1)))*0.5); %jitter
    g2_gr1 = log2(1+g2_exp(gr1)+rand(size(g2_exp(gr1)))*0.5);
    plot(g1_gr1,g2_gr1,'.')
    b = polyfit(g1_gr1,g2_gr1,1);
    hold on
    mn = min(g1_gr1);
    mx = max(g1_gr1);
    plot([mn,mx],b(2)+b(1)*[mn,mx],'-k')
    [r,pv] = corr(g1_gr1',g2_gr1');
    title([string(['ECT, r=',num2str(r)]),['-log(pv) = ',num2str(-log10(pv))],string(['N = ',num2str(length(g1_exp(gr1)))])]);
    ylabel(upper(g2))
    xlabel(upper(g1))
    box off
%     axes(ha(2))
    h2 = axes('position',[0.55,0.1,0.4,0.7]);
%     plot(g1_exp(gr2),g2_exp(gr2),'.')
% %     b = robustfit(g1_exp(gr2),g2_exp(gr2));
%     b = polyfit(g1_exp(gr2),g2_exp(gr2),1);
%     hold on
%     mn = min(g1_exp(gr2));
%     mx = max(g1_exp(gr2));
%     plot([mn,mx],b(2)+b(1)*[mn,mx],'-k')
%     [r,pv] = corr(g1_exp(gr2)',g2_exp(gr2)');
%     title([string(['PE, r=',num2str(r)]),['-log(pv) = ',num2str(-log10(pv))],string(['N = ',num2str(length(g1_exp(gr2)))])]);
    g1_gr2 = log2(1+g1_exp(gr2)+rand(size(g1_exp(gr2)))*0.5);
    g2_gr2 = log2(1+g2_exp(gr2)+rand(size(g2_exp(gr2)))*0.5);
    plot(g1_gr2,g2_gr2,'.')
    b = polyfit(g1_gr2,g2_gr2,1);
    hold on
    mn = min(g1_gr2);
    mx = max(g1_gr2);
    plot([mn,mx],b(2)+b(1)*[mn,mx],'-k')
    [r,pv] = corr(g1_gr2',g2_gr2');
    title([string(['EPE, r=',num2str(r)]),['-log(pv) = ',num2str(-log10(pv))],string(['N = ',num2str(length(g1_exp(gr2)))])]);
    xlabel(upper(g1))
    box off
    %     axes(ha(3))

    linkaxes([h1,h2],'xy')
    if savefig_flag==1
        eval(['export_fig Two_genes_correlation_cell_nuclei_TB/',clusteruni{typeplot(k)},'_',char(g1),'-',char(g2),'_',date,'.pdf']);
    end
    close(hf)
end
%% violin 
gn = upper('notch1');
g = find(strcmpi(gn,geneid_all));
fnorm = normpdf(norminv([0.001:0.01:0.999]));

figure('color','w','position',[20,20,1500,500],'Name',gn);
ha = axes('Position',[0.03,0.3,0.95,0.55]);
[~,cond] = max(condition_sorted(1:2,:));
cond_uni = unique(cond);
logflag = 0;
p99 = 0;% prctile(data(g,:),100);
markervec = 'opds<>';
xt = [];
% donorid_uni = unique(donor_id_flag_sorted);
notext_flag = 1;
for k=1:length(typeplot)
%     axes(ha(k))
    c=k;
    t_ed = zeros(2,1);%choose the desired column
    t_av = zeros(2,1);
    t_75 = zeros(2,1);
    for i=[1:2]
        gr2 = find(T_cells_tmp==typeplot(c) & condition_sorted(i,:)'==1);%find(fc_time_sorted==fc_time_uni(i));%
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

            y = 0.5*rand(length(gr2),1)-0.1+y';
            x = fi'.*(0.2*rand(length(gr2),1)-0.1);
           
            plot(k+i/3 + x-0.5, y,'.','color','none','marker',markervec(i)...
                ,'markerfacecolor',['#',colors{typeplot(k),2}],'markersize',3); hold on;
        else
            plot(-0.5+k+i/3+0.1*rand(length(gr2),1)-0.1,0.2*rand(length(gr2),1)-0.1+data_orig_all_sorted(g,gr2)'....
                ,'.','color','none','marker',markervec(i),'markerfacecolor',['#',colors{typeplot(k),2}],'markersize',3); hold on;
        end
        xt = [xt,k+i/3-0.5];
        %         plot(i+0.5*rand(length(gr2),1)-0.1,0.5*rand(length(gr2),1)-0.1+data(g,gr2)','.'); hold on;
%         t_ed(i) = median(y);
        t_av(i) = mean(data_orig_all_sorted(g,gr2));
%         t_75(i) = prctile(y,75);
        t_98(i) = prctile(y,95);
    end
    p99 = max([p99,t_98]);
    plot(k+[1:2]/3-0.5,t_av,'sr','MarkerFaceColor','k'); hold on;
%     plot(t_av,'or'); hold on;
%     plot(t_75,'dg'); hold on;
    axis tight
%     yl = get(gca,'ylim');
%     set(gca,'xtick',[k+[1:4]/5],'XTickLabel',[{'ce','pe','ct','pl',}],'ylim',[0,p99])
    set(gca,'ylim',[0,p99])
    %     axis tight
%     title([gn,', ',regexprep(clusteruni{typeplot(c)},'_','-')])
    comp = [[1,2]];
    pcomp = ones(length(comp(:,1)),1);
%     for jj=1:length(comp)
        gr1 = find(T_cells_tmp==typeplot(c) & condition_sorted(comp(1,1),:)'==1);
        gr2 = find(T_cells_tmp==typeplot(c) & condition_sorted(comp(1,2),:)'==1);
        y1 = (data_orig_all_sorted(g,gr1));
        y2 = (data_orig_all_sorted(g,gr2));
        pcomp = ranksum(y1,y2);
        if pcomp<1e-4
            sigstar='***';
        elseif pcomp<1e-3
            sigstar='**';
        elseif pcomp<1e-2
            sigstar='*';
        else
            sigstar='';
        end        
        plot(k + comp(1,:)/3-0.5,p99*(0.93+1*0.01)*[1,1],'-','color',0.5*[1,1,1]);
        text(mean(k + comp(1,:)/3-0.5),p99*(0.93+1*0.02),sigstar,'fontsize',7,'VerticalAlignment','bottom','HorizontalAlignment','center');
%     end        
end
set(gca,'ylim',[0,p99]);
yl = get(gca,'ylim');
text(xt,-(yl(2)-yl(1))*0.05*ones(size(xt)),repmat({'ce','pe'},1,length(typeplot)),'fontsize',5)
set(gca,'xtick',[1:length(typeplot)],'XTickLabel',clusteruni(typeplot),'TickLabelInterpreter'....
    ,'none','XTickLabelRotation',45,'ylim',[-(yl(2)-yl(1))*0.1,yl(2)]);
box off
% eval(['export_fig violin_',gn,'_',date,'.pdf']);

%%
figure('color','w','position',[20,20,700,400]);
k = 0;
colorvec = distinguishable_colors(length(T_cells_tmp_uni));
sn_per_cluster = sum(sn_per_cluster_insample);
[~,xi_sn] = sort(sn_per_cluster);
sc_per_cluster = sum(sc_per_cluster_insample);
[~,xi_sc] = sort(sc_per_cluster);
hb = bar(0.5,sn_per_cluster(xi_sn)./(sum(sn_per_cluster(xi_sn))),0.5,'stacked'); hold on
% for jj=1:length(hb)
%     set(hb(jj),'FaceColor',colorvec(jj,:));
% end
hb = bar(1.5,sc_per_cluster(xi_sc)./(sum(sc_per_cluster(xi_sc))),0.5,'stacked');
for jj=1:length(hb)
    set(hb(jj),'FaceColor',colorvec(jj,:));
end
% hleg = legend(clusteruni,'interpreter','none');
hleg = legend(clusteruni([xi_sn(end-4:end),xi_sc(end-7:end)]),'interpreter','none');
set(hleg,'position',[0.78,0.3,0.2,0.6],'fontsize',6)
% axis tight
set(gca,'position',[0.05,0.1,0.7,0.85],'xlim',[0,2],'xtick',[0.5,1.5],'xticklabel',{'sn','sc'})
%% abundaces comparison

figure('color','w','position',[20,20,700,400]);
k = 0;
colorvec = distinguishable_colors(length(T_cells_tmp_uni));
sn_per_cluster = sum(sn_per_cluster_insample);
prc_sn = sn_per_cluster./(sum(sn_per_cluster));
sc_per_cluster = sum(sc_per_cluster_insample);
prc_sc = sc_per_cluster./(sum(sc_per_cluster));
hb = bar(0.5,prc_sn(find(prc_sn>0)),0.5,'stacked'); hold on
hb = bar(1.5,prc_sc(find(prc_sc>0)),0.5,'stacked');
for jj=1:length(hb)
    set(hb(jj),'FaceColor',colorvec(jj,:));
end
% hleg = legend(clusteruni,'interpreter','none');
hleg = legend([clusteruni(find(prc_sn>0));clusteruni(find(prc_sc>0))],'interpreter','none');
set(hleg,'position',[0.78,0.3,0.2,0.6],'fontsize',6)
% axis tight
set(gca,'position',[0.05,0.1,0.7,0.85],'xlim',[0,2],'xtick',[0.5,1.5],'xticklabel',{'sn','sc'})

if savefig_flag==1
    eval(['export_fig cell_nuclei_abundances_comparison_',date,'.pdf']);
end
%% 
datalog = log2(data_orig_all_sorted+1);
i= 0;
for k=1:length(T_cells_tmp_uni)
    for jj=1:4
        i = i+1;
        gr1 = find(T_cells_tmp==k & condition_sorted(jj,:)'==1);
        condgrsize(i) = length(gr1);
        mean_mat(:,i) = mean(datalog(:,gr1),2);
        mean_lin_mat(:,i) = mean(data_orig_all_sorted(:,gr1),2);
        median_mat(:,i) = median(datalog(:,gr1),2);
        fracpos_mat(:,i) = mean(data_orig_all_sorted(:,gr1)>0,2);
        p75_mat(:,i) = prctile(datalog(:,gr1),75,2);
    end
end
%% bin anaysis

figure;
set(gcf,'color','w','position',[20,20,600,600]);

%     [ha, pos] = tight_subplot(4, ceil(1/4), [0.08,0.08], [0.04,0.04], [0.05,0.08]);

sig_gene_cluster = cell(1,1);
th = 0.7;
list = [];
geneper_bin = zeros(1,8);

gr1 = find( (T_cells_tmp == 3) & peearly_flag_sorted==1);%sc-EVT1-epe
gr2 = find( (T_cells_tmp == 3) & ctrlearly_flag_sorted==1); %sc-EVT1-ect
gr3 = find( (T_cells_tmp == 2) & peearly_flag_sorted==1);  %sn-EVT-epe   
gr4 = find( (T_cells_tmp == 2) & ctrlearly_flag_sorted==1); %sn-EVT-ect

x1 = mean(log2(data_orig_all_sorted(:,gr1)+1),2);
x2 = mean(log2(data_orig_all_sorted(:,gr2)+1),2);
x3 = mean(log2(data_orig_all_sorted(:,gr3)+1),2);
x4 = mean(log2(data_orig_all_sorted(:,gr4)+1),2);
dx = x1-x2 ;
dy = x3-x4 ;

inexc = find(mean(data_orig_all_sorted>0,2)>0.9 | sum(data_orig_all_sorted>0,2)<100);
dx(inexc) = 0;
dy(inexc) = 0;
dis = sqrt(dx.^2 + dy.^2);

indsig = find(dis>th);
psig_x = ones(size(indsig));
psig_y = ones(size(indsig));
for jj=1:length(indsig)
    psig_x(jj) = ranksum(data_orig_all_sorted(indsig(jj),gr1),data_orig_all_sorted(indsig(jj),gr3));
    psig_y(jj) = ranksum(data_orig_all_sorted(indsig(jj),gr2),data_orig_all_sorted(indsig(jj),gr4));
end
q = qval_from_pval([psig_x;psig_y])';
q = min(reshape(q,length(indsig),2),[],2);
indsig = indsig(q<0.01);
sig_gene_cluster = geneid_all(indsig);%[dx(dis>th),dy(dis>th)];
list = [list;sig_gene_cluster];
[~,xi] = sort(dis,'descend');
theta = atan2(dy(indsig),dx(indsig))*180/pi;
theta(theta<0) = 360+theta(theta<0);
anglebins = [[0:45:(360)]-22.5,360];
[tmp,~,bins] = histcounts(theta,anglebins);
tmp = [sum(tmp([1,9])),tmp(2:8)];
geneper_bin(1,:) = [tmp];
bins(bins==9) = 1;
tmp2 = cell(max(tmp),8);
for jj=1:8
    tmp2(1:tmp(jj),jj) = sig_gene_cluster(bins==jj);
end
%     polarhistogram(atan2(dy(dis>th),dx(dis>th)),100);
plot(dx(indsig),dy(indsig),'.'); hold on;
axis tight
xl = get(gcf,'xlim');
plot([-2.5,2.5],[-2.5,2.5],'-k');
plot([-2.5,2.5],tan(22.5/180*pi)*[-2.5,2.5],'-r');
plot([-2.5,2.5],tan(67.5/180*pi)*[-2.5,2.5],'-r');
plot([-2.5,2.5],tan(-22.5/180*pi)*[-2.5,2.5],'-r');
plot([-2.5,2.5],tan(-67.5/180*pi)*[-2.5,2.5],'-r');
plot([-2.5,2.5],0*[-2.5,2.5],'-k');
plot(0*[-2.5,2.5],[-2.5,2.5],'-k');
% axis equal
% set(gcf,'xlim',[-1.5,1.5],'ylim',[-3,3])
xlim([-2,2]); xlabel('sc-EVT epe-ect')
ylim([-3.5,3.5]); ylabel('sn-EVT epe-ect')
%     text(dx(xi(1:ntop)),dy(xi(1:ntop)),geneid_all(xi(1:ntop)),'fontsize',8);
text(dx(indsig),dy(indsig),geneid_all(indsig),'fontsize',8);
tmp = m2c(geneper_bin(1,:));
tmp = cellfun(@num2str,tmp,'UniformOutput',0);
text(th*cos((anglebins(1:8)+22.5)*pi/180),th*sin((anglebins(1:8)+22.5)*pi/180),tmp,'fontsize',8)

eval(['export_fig sc_EVT1_sn_EVT_pe_ct_bin_analysis_',date,'.pdf']);
