tic
clear all
close all
addpath('/data/matlab_functions/')
terminate(pyenv)
pyenv("ExecutionMode","OutOfProcess")

savefig_flag = 1;
savefig_pdf = 1;

load PE-nuclei_QC800_400_afterloading_30-May-2023.mat

cellid_type = loadCellFile('trophoblasts_step2_cellid_cluster_beforeclean_30-May-2023.txt');
typenum = cell2mat(cellid_type(:,2));
clean_list = loadCellFile('cleaning_trophoblasts_final_230413.txt');
clean_list = cell2mat(clean_list(3:end,:));
for i=1:length(clean_list(:,1))
    if clean_list(i,2)==0
        typenum(typenum==clean_list(i,1)) = 0;
    end
end
cellid_type(typenum==0,:) = [];% selecting only TB
[~,loc] = ismember(cellid_type(:,1),cellid);
sum(loc>0)

tot_mol = sum(data);
tot_mol(tot_mol>5e4) = 5e4;
tot_genes = sum(data>0);
% validcells = (tot_mol>1300 & tot_mol<5e4 & tot_genes>800);
% sum(validcells)
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


% validcells = false(size(cellid));
% validcells(loc) = true;

sample_uni = {'132-1','134-2','134-3','135-1','135-2'};%
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
saveCellFile(table1,['gene_correlation_trophoblasts_',date,'.txt'])
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
% initial_dims = 30;
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

D = squareform(pdist(prj,'correlation'),'tomatrix');
[Dsort,XI] = sort(D,'ascend');
per_range = 500;
x = 1:per_range;%length(D);
optk = zeros(length(D),1);
for i=1:length(D)
    y = Dsort(1:per_range,i);
    x = x(:);
    y = y(:);
    a = atan((y(end)-y(1))/(x(end)-x(1)));
    xn = x*cos(a) + y*sin(a);
    yn = -x*sin(a) + y*cos(a);
    [~,imax] = max(yn);
    optk(i) = round((x(imax)));
end

perplexity = median(optk);

options = statset('MaxIter',1000);
mapped_xy = tsne(prj,'Algorithm','barneshut','Distance','correlation','NumDimensions',2,'NumPCAComponents',0,.....
    'Perplexity',perplexity,'Standardize',true,'InitialY',init,'LearnRate',length(cellid)/12,'Theta',0.5,'Verbose',1,'Options',options,'Exaggeration',20);
toc

% this is just the initial tsne, can be commented later
figure;
set(gcf,'color','w','position',[20,20,900,800])
plot(mapped_xy(:,1),mapped_xy(:,2),'.'); axis tight; axis off
% plot by sample
% sample_uni = {'40-1','44-1','48-1','53-1'};
colors = distinguishable_colors(length(sample_uni)+1);
figure;
set(gcf,'color','w','position',[20,20,900,800]);
[ha, pos] = tight_subplot(1, 1, [0.02,0.02], [0.02,0.02], [0.02,0.02]);
axes(ha(1))
for i=1:length(sample_uni)
    s = scatter(mapped_xy(strcmpi(sample,sample_uni{i}),1),mapped_xy(strcmpi(sample,sample_uni{i}),2),10,colors(i,:),'filled'); hold on;
%     alpha(s,0.4);
end
axis tight
axis equal
axis off
legend(sample_uni)

%% clustering with dbscan
MinPts = 60;
eps_prc = 50;
[idx, isnoise] = dbscan_epsprc_mipts(mapped_xy,eps_prc,MinPts);

colors = distinguishable_colors(length(unique(idx))+1);
figure;
set(gcf,'color','w','position',[20,20,900,800])
for i=unique(idx)'
    if i>=0
        ii=find(idx==i); h=plot(mapped_xy(ii,1),mapped_xy(ii,2),'o','color',colors(i+1,:),'markersize',3); hold on;
    elseif i>0
        ii=find(idx==i); h=plot(mapped_xy(ii,1),mapped_xy(ii,2),'.','color',colors(i+1,:),'markersize',4); hold on;
    end
end

axis tight;
axis equal
axis off

title(['perplexity=',num2str(perplexity),', PC=',num2str(initial_dims),',#C=',num2str(max(idx)),',#out=',num2str(sum(idx==0))],'fontsize',8);

%% flags sorting 
[idx,xi] = sort(idx);
xi(idx==0) = [];
idx(idx==0) = [];
idxuni = unique(idx);
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


no_dims = 1;
initial_dims = 8;
perplexity = 5;
epsilon = 100;
dist_flag = 2;
theta = 0.5;
rand_seed = 13;
data_tsne = cent_norm(log2(data_sorted_all+1));
xi = [1:length(idx)];
for i=1:length(idxuni)
    i
    ind = find(idx==i);
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

mapped_xy = mapped_xy(xi,:);
%% dendrogram tree

meangr_mat = zeros(length(moldata(:,1)),length(idxuni));
clust_cent = zeros(length(idxuni),2);
for jjj=1:length(idxuni)
    jjj
    meangr_mat(:,jjj) = mean(log2(data_sorted_all(:,idx==idxuni(jjj))+1),2);
    clust_cent(jjj,:) = [median(mapped_xy(idx==idxuni(jjj),1)),median(mapped_xy(idx==idxuni(jjj),2))];
end
meangr_mat1 = meangr_mat;
% meangr_mat1(loc,:) = [];

% meangr_mat1 = cent_norm(meangr_mat1(:,leaforder));
[prj,m,D,V,Q] = pca_wis(meangr_mat1',initial_dims);
Zpca = linkage(prj,'ward','correlation');
Dpca = pdist(prj,'correlation');
leaforder_pca = optimalleaforder(Zpca,Dpca);
figure;
set(gcf,'position',[100,100,1000,1000],'color','w')
ax1 = axes('position',[0.03,0.03,0.3,0.93]);
% hden = dendrogram(Zpca,length(leaforder_pca),'Orientation','left');
hden = dendrogram(Zpca,length(leaforder_pca),'Reorder',leaforder_pca,'Orientation','left');
axis off
set(gca,'ylim',[0.5,length(leaforder_pca)+0.5])
ax2 = axes('position',[0.35,0.03,0.63,0.93]);
x=squareform(Dpca); imagesc(x(leaforder_pca,leaforder_pca));
colormap('summer')
set(gca,'ytick',[1:length(leaforder_pca)],'xtick',[],'fontsize',8,'ydir','normal')
linkaxes([ax1,ax2],'y')
if savefig_flag==1
    savefig(gcf,['tree_step3_trophoblasts_',date,'.fig'])    
end
leaforder = leaforder_pca;


T_cells_tmp_new = zeros(size(idx));
for i=1:length(leaforder)
    T_cells_tmp_new(idx==idxuni(leaforder(i))) = i;
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


%% tsne by cluster
T_cells_tmp = T_cells_tmp_new;
T_cells_tmp_uni = unique(T_cells_tmp);
idx = T_cells_tmp_new;
idxuni = idxuni_new;


colors = distinguishable_colors(length(unique(idx))+1);
figure;
set(gcf,'color','w','position',[20,20,900,800])
for i=unique(idx)'
    if i==-1
        ii=find(idx==i); h=plot(mapped_xy(ii,1),mapped_xy(ii,2),'.','color',[0,0,1],'markersize',3); hold on;
    else
        ii=find(idx==i); h=plot(mapped_xy(ii,1),mapped_xy(ii,2),'.','color',colors(i+1,:),'markersize',5); hold on;
    end
end
for i=idxuni'
    in = idx==i;
    ht = text(median(mapped_xy(in,1))+1,median(mapped_xy(in,2))+1,num2str(i));    
%     set(ht,'BackgroundColor',0.8*[1,1,1],'fontsize',8);
    set(ht,'BackgroundColor','none','fontsize',8)
    plot([median(mapped_xy(in,1)),median(mapped_xy(in,1))+1],[median(mapped_xy(in,2)),median(mapped_xy(in,2))+1]....
        ,'linewidth',0.5,'color',0.8*[1,1,1]);
end
axis tight;
axis equal
axis off
title(['MinPts=',num2str(MinPts),', epsprc=',num2str(eps_prc),',#C=',num2str(max(idx)),',#out=',num2str(sum(idx==0))],'fontsize',8);
if savefig_flag==1
    savefig(gcf,['tsne_by_cluster_step3_trophoblasts_perplexity_',num2str(perplexity),'_',date,'.fig'])
    eval(['export_fig tsne_by_cluster_step3_trophoblasts','_',date,'.pdf']);
end
%% tsne by sample

colors = distinguishable_colors(length(sample_uni)+1);
figure;
set(gcf,'color','w','position',[20,20,900,800]);
[ha, pos] = tight_subplot(1, 1, [0.02,0.02], [0.02,0.02], [0.02,0.02]);
axes(ha(1))
for i=1:length(sample_uni)
    s = scatter(mapped_xy(strcmpi(sample_sorted,sample_uni{i}),1),mapped_xy(strcmpi(sample_sorted,sample_uni{i}),2),5,colors(i,:),'filled'); hold on;
%     alpha(s,0.4);
end
axis tight
axis equal
axis off
legend(sample_uni)

if savefig_flag==1
    savefig(gcf,['trophoblasts_step3_tsne_by_sample_perplexity_',num2str(perplexity),'_',date,'.fig'])
 
end

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
f = ctrlearly_flag_sorted;
plot(mapped_xy(:,1),mapped_xy(:,2),'.','color',0.6*[1,1,1],'markersize',ms); hold on;
plot(mapped_xy(f,1),mapped_xy(f,2),'.','color',[1,0,0],'markersize',ms); hold on;
axis tight
axis equal
axis off
title('Early control')



%% median per condition 
condition_sorted = zeros(6,length(peearly_flag_sorted));
condition_sorted(1,ctrlearly_flag_sorted) = 1;
condition_sorted(3,ctrl_flag_sorted) = 1;
condition_sorted(2,peearly_flag_sorted) = 1;
condition_sorted(4,pelate_flag_sorted) = 1;
condition_sorted(5,iugr_flag_sorted) = 1;
condition_sorted(6,(peearly_flag_sorted | pelate_flag_sorted) & ~iugr_flag_sorted) = 1;

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
%%
per_cluster_insample = zeros(length(sample_uni),length(T_cells_tmp_uni));
for i=1:length(sample_uni)
    in = find(strcmpi(sample_sorted,sample_uni{i}));
    per_cluster_insample(i,:) = histcounts(T_cells_tmp(in),[1:length(T_cells_tmp_uni)+1]);
end
per_cluster_insample = per_cluster_insample./repmat(sum(per_cluster_insample,2),1,length(T_cells_tmp_uni));
condition_sorted_uni = unique(condition_sorted);
colorvec = distinguishable_colors(2);

figure;
set(gcf,'color','w','position',[20,20,1200,1000]);
[ha, pos] = tight_subplot(8, ceil(length(T_cells_tmp_uni)/8), [0.02,0.02], [0.02,0.02], [0.02,0.02]);
for j=1:length(T_cells_tmp_uni)
    axes(ha(j))
    t_cell = cell(2,1);
    for i=1:2
        s = unique(sample_sorted(condition_sorted(i,:)==1));
        [~,loc] = ismember(s,sample_uni);
        t = per_cluster_insample(loc,j);
        t_cell{i} = t;
        plot(i,100*t,'.','color',colorvec(i,:)); hold on;
        plot(i,100*median(t),'sk','markersize',12)
    end
%     [~,p_ct_pe] = ttest2(t_cell{2},t_cell{3});
%     [~,p_ct_pl] = ttest2(t_cell{2},t_cell{4});
%     [~,p_pe_pl] = ttest2(t_cell{3},t_cell{4});
%     [~,p_up_um] = ttest2(t_cell{5},t_cell{6});
    [~,p_ce_pe] = ttest2(t_cell{1},t_cell{2});
    box off
    axis tight
    yl = get(gca,'ylim');
%     plot([2,3],yl(2)*0.7*[1,1],'k');
%     text(2.5,yl(2)*0.7,num2str(p_ct_pe,2),'fontsize',5,'VerticalAlignment','bottom');
%     plot([3,4],yl(2)*0.8*[1,1],'k');
%     text(3.5,yl(2)*0.8,num2str(p_pe_pl,2),'fontsize',5,'VerticalAlignment','bottom');
%     plot([2,4],yl(2)*0.9*[1,1],'k');
%     text(3,yl(2)*0.9,num2str(p_ct_pl,2),'fontsize',5,'VerticalAlignment','bottom');
%     plot([5,6],yl(2)*0.7*[1,1],'k');
%     text(5.5,yl(2)*0.7,num2str(p_up_um,2),'fontsize',5,'VerticalAlignment','bottom');
    plot([1,2],yl(2)*0.9*[1,1],'k');
    text(2,yl(2)*0.9,num2str(p_ce_pe,2),'fontsize',5,'VerticalAlignment','bottom');
    set(ha(j),'xlim',[0.5,2.5],'xtick',[1:2],'xticklabel',{'ce','pe'})
    title(num2str(j))
end

if savefig_flag==1
    savefig(gcf,['trophoblasts_step3_median_per_condition_',date,'.fig'])
    eval(['export_fig trophoblasts_step3_median_per_condition_',date,'.pdf']);
end
% % % % % % % % % % % % % % % % % % 
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
    savefig(gcf,['trophoblats_step3_markertable_',date,'.fig']) 
    eval(['export_fig trophoblats_step3_markertable_',date,'.pdf']);
end


% % % % % % % % % % % % % % % % % % % % % % % % % % 
table2 = [cellid_sorted,m2c(T_cells_tmp)];
saveCellFile(table2,['trophoblasts_step3_cellid_cluster_beforeclean_',date,'.txt']);

%%  gene expression on tsne
cat_markers = {'PTPRC','TYROBP','FCER1G','CLDN5','CD34','SPARCL1','HBA1','ACTA2','APOD','TCF21','CXCL14','XAGE2','EGFR','XAGE3','CSH2','HLA-G'};
% exclude_markers = {'C1qc','C1qa','C1qb','Gja1','Cx3cr1','Acta2','Ly6c1','Mfge8','Plp1'....
%     ,'Aqp4','Vtn','Cldn5','Pdgfrb','Flt1','Slc1a3','Pdgfra','Foxj1','Olig1','Olig2','Sox10','Hbb-bs','Hbb-bt','Hba-a2'};
[~,loc] = ismember(cat_markers,geneid_all);

nontype = zeros(5,length(T_cells_tmp_uni));
nontype_sc = zeros(size(T_cells_tmp));
for jjj=1:length(T_cells_tmp_uni)
    jjj
    tmp = mean(data_orig_all_sorted(loc,T_cells_tmp==jjj)>0,2);
    tmp = tmp>0.5;
    
    nontype(:,jjj) = [max(tmp(1:3)),max(tmp(4:6)),max(tmp(7:8)),max(tmp(9:11)),max(tmp(12:16))]';
    if sum(nontype(:,jjj))>1
        nontype_sc(T_cells_tmp==jjj) = -1;
    elseif sum(nontype(:,jjj))==0
        [~,imax] = max(mean(data_orig_all_sorted(loc,T_cells_tmp==jjj)>0,2));
        c = [1,1,1,2,2,2,3,3,4,4,4,5,5,5,5,5];
        nontype_sc(T_cells_tmp==jjj) = c(imax);
    else
        nontype_sc(T_cells_tmp==jjj) = find(nontype(:,jjj));
    end
end

figure;
set(gcf,'color','w','position',[20,20,800,800])
colors = distinguishable_colors(length(unique(nontype_sc)));
k = 0;
for idx=unique(nontype_sc)'
    k = k+1;
    ii=find(nontype_sc==idx); 
    h=plot(mapped_xy(ii,1),mapped_xy(ii,2),'.','color',colors(k,:),'markersize',3); hold on;
end
axis tight
axis equal
axis off
title('main categories')
legend('doublets','immune','vascular','eryth','stromal','trophblast')

% if savefig_flag==1
%     savefig(gcf,['tsne_by_big_clusters_v8',num2str(perplexity),'_',date,'.fig'])
%     eval(['export_fig tsne_by_big_clusters_v8','_',date,'.pdf']);
% end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%%
list = {'PAGE4','PEG10','PHLDA2','XAGE3','DSP','ITGA6','FBN2','NOTUM','PAPPA2','PSG3','PSG2','PSG4','PSG5','ERVFRD-1','PLAC8','LVRN'};
figure;
set(gcf,'color','w','position',[20,20,1100,960])
[ha, pos] = tight_subplot(4, 4, [0.02,0.02], [0.02,0.02], [0.02,0.02]);
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




list = {'PTPRC','NOTUM','DIO2','HLA-G','PRG2','LAIR2','TAC3','PSG1','PAPPA','PAPPA2','PRG2','ENG','CYP19A1','FLT1','DHRS9','FAR2'};
figure;
set(gcf,'color','w','position',[20,20,1100,960])
[ha, pos] = tight_subplot(4, 4, [0.02,0.02], [0.02,0.02], [0.02,0.02]);
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
% % % % % % % % % % % % % % % % % % % % % % 
list = {'EBI3','IGF2','ORMDL3','GATA2','KIR2DL4','NOTUM','PSG9','HLA-G','DIO2','TCF21','CXCL14','XAGE2','EGFR','XAGE3','CSH2','IFIT2'};;
figure;
set(gcf,'color','w','position',[20,20,1100,960])
[ha, pos] = tight_subplot(4, 4, [0.02,0.02], [0.02,0.02], [0.02,0.02]);
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
%% top genes - scatter comparison
top_g = 50;
gr1name = '3';
gr2name = '4';
gr1 = find( (T_cells_tmp == 3 ));
gr2 = find( (T_cells_tmp == 4 ));
in = find(mean(data_orig_all_sorted(:,gr1)>0,2)>0.03 | mean(data_orig_all_sorted(:,gr2)>0,2)>0.03);
ptt = zeros(length(in),1);
for s=1:length(in)
    ptt(s) = ranksum(data_orig_all_sorted(in(s),gr1),data_orig_all_sorted(in(s),gr2),'tail','both');
end
ptt(isnan(ptt)) = 1;
x1 = mean(log2(data_orig_all_sorted(:,gr1)+1),2);
x2 = mean(log2(data_orig_all_sorted(:,gr2)+1),2);
d = x1-x2 ;
figure('position',[200,200,1000,580],'color','w');
plot(d(in),-log10(ptt),'.'); hold on;
[~,xi] = sort(ptt);
plot(d(in(xi(1:100))),-log10(ptt(xi(1:100))),'or');
text(d(in(xi(1:100))),-log10(ptt(xi(1:100))),geneid_all(in(xi(1:100))))
[~,xi] = sort(d);
figure('position',[200,200,1000,580],'color','w');
[ha, pos] = tight_subplot(1, 2, [0.05,0.05], [0.1,0.05], [0.05,0.05]);
axes(ha(1))
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
axes(ha(2))
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
% % % % % % % % % % % % % % % % % % % % % % % % % 

flag_mat = [female_flag_sorted,~female_flag_sorted,ctrl_flag_sorted,ctrlearly_flag_sorted,.....
    peearly_flag_sorted,pelate_flag_sorted,iugr_flag_sorted]';
flag_title = {'female','male','ctrl','ctrl early','PE early','PE late','IUGR'};
t1 = [];
t2 = [];
for ii=1:length(flag_mat(:,1))
    marker = flag_mat(ii,:);
    marker_percent = zeros(1, length(T_cells_tmp_uni));
    p_marker = zeros(1, length(T_cells_tmp_uni));
    M = length(marker);
    K = sum(marker>0);
    for j=1:length(T_cells_tmp_uni)
        c1 = sum( T_cells_tmp==T_cells_tmp_uni(j));
        c2 = sum( T_cells_tmp==T_cells_tmp_uni(j) & marker'>0);
        marker_percent(j) = 100*c2/c1;
        p_marker(j) = hygecdf(c2,M,K,c1,'upper');
    end
    marker_percent(isnan(marker_percent)) = 0;
%     eval([markers{ii},'_percent = marker_percent;'])
%     eval([markers{ii},'_p = p_marker;']);
    t1 = [t1;marker_percent];
    t2 = [t2;p_marker];
end
figure;
set(gcf,'color','w','position',[20,20,900,1800])
[ha, pos] = tight_subplot(1, 7, [0.02,0.00], [0.02,0.02], [0.05,0.02]);
for i=1:7
    axes(ha(i));
    barh(t1(i,:));
    set(gca,'ytick',[1:length(idxuni)],'YTickLabel',cell(length(idxuni),1),'XLim',[0,100])
    title([flag_title{i},' %>0'])
end
axes(ha(1));
set(gca,'ytick',[1:length(idxuni)],'YTickLabel',[1:length(idxuni)]);
linkaxes([ha(1),ha(2),ha(3),ha(4),ha(5),ha(6),ha(7)],'y');

figure;
set(gcf,'color','w','position',[20,20,900,1800])
[ha, pos] = tight_subplot(1, 7, [0.02,0.00], [0.02,0.02], [0.05,0.02]);
for i=1:7
    axes(ha(i));
    barh(-log10(t2(i,:)));
    set(gca,'ytick',[1:length(idxuni)],'YTickLabel',cell(length(idxuni),1),'XLim',[0,10])
    title([flag_title{i},' -log10(p)'])
end
axes(ha(1));
set(gca,'ytick',[1:length(idxuni)],'YTickLabel',[1:length(idxuni)]);
linkaxes([ha(1),ha(2),ha(3),ha(4),ha(5),ha(6),ha(7)],'y');
%% gene expression - violin 
gn = 'dio2';
g = find(strcmpi(gn,geneid_all));

figure;
set(gcf,'color','w','position',[20,20,1200,1000]);
b = 1;
[ha, pos] = tight_subplot(4, ceil(length(unique(T_cells_tmp(:,b)))/4), [0.04,0.04], [0.04,0.04], [0.04,0.04]);
[~,cond] = max(condition_sorted(1:4,:));
cond_uni = unique(cond);
for k=1:length(unique(T_cells_tmp(:,b)))
    axes(ha(k))
    c=k;
    t_ed = zeros(2,1);
    t_av = zeros(2,1);
    t_75 = zeros(2,1);
    for i=[1:2]
        gr2 = find(T_cells_tmp(:,b)==c & cond'==cond_uni(i));%find(fc_time_sorted==fc_time_uni(i));%
        plot(i+0.5*rand(length(gr2),1)-0.1,0.5*rand(length(gr2),1)-0.1+data_orig_all_sorted(g,gr2)','.'); hold on;
        t_ed(i) = median(data_orig_all_sorted(g,gr2));
        t_av(i) = mean(data_orig_all_sorted(g,gr2));
        t_75(i) = prctile(data_orig_all_sorted(g,gr2),75);
    end
    plot(t_ed,'-sk'); hold on;
    plot(t_av,'-or'); hold on;
    plot(t_75,'-dg'); hold on; 
    axis tight
    yl = get(gca,'ylim');
    set(gca,'xtick',[1:4],'XTickLabel',[{'ce','ct','pe','pl'}],'yscale','linear','ylim',[-0.5,ceil(yl(2))+1])
    axis tight
    title([gn,',c=',num2str(c)])
end

