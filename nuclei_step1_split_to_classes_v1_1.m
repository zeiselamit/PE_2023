tic
clear all
close all
addpath('/data/matlab_functions/')
terminate(pyenv)
pyenv("ExecutionMode","OutOfProcess")

savefig_flag = 0;
savefig_pdf = 0;

load PE-nuclei_QC800_400_afterloading_30-May-2023.mat

tot_mol = sum(data);
tot_genes = sum(data>0);


figure('position',[100,100,800,800],'color','w');
subplot(2,2,1);
[f,xi] = ksdensity(tot_mol);
plot(xi,f);
axis tight;
set(gca,'xlim',[0,prctile(tot_mol,98)]);
xlabel('total mol')
subplot(2,2,3);
[f,xi] = ecdf(tot_mol);
plot(xi,f);
axis tight;
set(gca,'xlim',[0,prctile(tot_mol,98)]);
xlabel('total mol')
subplot(2,2,2);
[f,xi] = ksdensity(tot_genes);
plot(xi,f);
axis tight;
set(gca,'xlim',[0,prctile(tot_genes,98)]);
xlabel('total genes')
subplot(2,2,4);
[f,xi] = ecdf(tot_mol);
plot(xi,f);
axis tight;
set(gca,'xlim',[0,prctile(tot_genes,100)]);
xlabel('total genes')

% % % % % % % % % % 

data = ceil(data./repmat(sum(data),length(data(:,1)),1)*10e3);

sex_genes = {'XIST','MTRNR2L8', 'EIF1AY', 'DDX3Y', 'RPS4Y1', 'KDM5D','MTRNR2L12','MTRNR2L10','MTRNR2L12'};
in = find(sum(data>0,2)>100 & sum(data>0,2)<length(data(1,:))*0.5 & ~ismember(geneid(:,1),sex_genes)); % & ~ismember(geneid(:,1),batch_genes)


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
saveCellFile(table1,['gene_correlation_general_',date,'.txt'])
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
batchid = sample;%donor_id_flag;
usepylib = 1;
[sout]=harmonypy(s,batchid,usepylib);
prj = sout;
%% tsne 
rng(13);
if length(sample)>10000
    subsam = randperm (length(sample));
    subsam = subsam(1:1e4);
    D = squareform(pdist(prj(subsam,:),'correlation'),'tomatrix');
else
    D = squareform(pdist(prj,'correlation'),'tomatrix');
end
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

% plot by sample
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
MinPts = 30;
eps_prc = 85;
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
vaginalbirth_flag_sorted = vaginalbirth_flag((xi));
induction_flag_sorted = induction_flag((xi));
noinduction_flag_sorted = noinduction_flag((xi));
magnesium_flag_sorted = magnesium_flag((xi));
spinal_flag_sorted = spinal_flag((xi));
epidural_flag_sorted = epidural_flag((xi));
generalanesthesia_flag_sorted = generalanesthesia_flag((xi));
donor_id_flag_sorted = donor_id_flag(xi);
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
[prj,m,D,V,Q] = pca_wis(meangr_mat1',10);
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
    savefig(gcf,['tree_QC2000_',date,'.fig'])    
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
    savefig(gcf,['general_tsne_by_cluster_v7_perplexity',num2str(perplexity),'_',date,'.fig'])
    eval(['export_fig general_tsne_by_cluster_v7','_',date,'.pdf']);
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
    savefig(gcf,['general_tsne_by_sample_v7_perplexity',num2str(perplexity),'_',date,'.fig'])
    eval(['export_fig general_tsne_by_sample_v7_perplexity_','_',date,'.pdf']);
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
f = ~ctrlearly_flag_sorted;
plot(mapped_xy(:,1),mapped_xy(:,2),'.','color',0.6*[1,1,1],'markersize',ms); hold on;
plot(mapped_xy(f,1),mapped_xy(f,2),'.','color',[1,0,0],'markersize',ms); hold on;
axis tight
axis equal
axis off
title('Ctrl early')
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


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
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
    savefig(gcf,['markertable_general_v7_',date,'.fig'])    
end

%% gene expression on tsne 
% list = {'NOTUM','KISS1','XAGE2','ECM1','ICAM1','HGF','CDH5','CNN1','COX4I2','CD34','MYH11','ACTA2','ESAM','CXCL14','PECAM1','APOD',};
list = {'PTPRC','TYROBP','FCER1G','APLNR','CD34','SPARCL1','NOTCH3','ACTA2','APOD','TCF21','CXCL14','XAGE2','EGFR','XAGE3','CSH2'};
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

% cat_markers = {'PTPRC','TYROBP','FCER1G','APLNR','CD34','SPARCL1','HBA1','HBA2','APOD','TCF21','CXCL14','XAGE2','EGFR','XAGE3','CSH2'};
cat_markers = {'PTPRC','TYROBP','FCER1G','APLNR','CD34','SPARCL1','NOTCH3','ACTA2','APOD','TCF21','CXCL14','XAGE2','EGFR','NOTUM','CSH2'};
% exclude_markers = {'C1qc','C1qa','C1qb','Gja1','Cx3cr1','Acta2','Ly6c1','Mfge8','Plp1'....
%     ,'Aqp4','Vtn','Cldn5','Pdgfrb','Flt1','Slc1a3','Pdgfra','Foxj1','Olig1','Olig2','Sox10','Hbb-bs','Hbb-bt','Hba-a2'};
[~,loc] = ismember(cat_markers,geneid_all);

nontype = zeros(5,length(T_cells_tmp_uni));
nontype_sc = zeros(size(T_cells_tmp));
for jjj=1:length(T_cells_tmp_uni)
    jjj
    tmp = mean(data_orig_all_sorted(loc,T_cells_tmp==jjj)>0,2);
    tmp = tmp>0.7;
    
    nontype(:,jjj) = [max(tmp(1:3)),max(tmp(4:6)),max(tmp(7:8)),max(tmp(9:11)),max(tmp(12:15))]';
    if sum(nontype(:,jjj))>1
        nontype_sc(T_cells_tmp==jjj) = -1;
    elseif sum(nontype(:,jjj))==0
        [~,imax] = max(mean(data_orig_all_sorted(loc,T_cells_tmp==jjj)>0,2));
        c = [1,1,1,2,2,2,3,3,4,4,4,5,5,5,5];
        nontype_sc(T_cells_tmp==jjj) = c(imax);
    else
        nontype_sc(T_cells_tmp==jjj) = find(nontype(:,jjj));
    end
end
% nontype_sc(T_cells_tmp==14) = 5;%manually set the PSG cluster to trophoblasts

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
legend('doublets','immune','vascular','mural','stromal','trophblast')

if savefig_flag==1
    savefig(gcf,['tsne_by_big_clusters_v7',num2str(perplexity),'_',date,'.fig'])
    eval(['export_fig tsne_by_big_clusters_v7','_',date,'.pdf']);
end

table1 = [cellid_sorted,m2c([nontype_sc,T_cells_tmp])];
saveCellFile(table1,['cellid_type_PE_QC800_400_',date,'.txt']);

% PTPRC,TYROBP,FCER1G  - immune cells, ESAM,CAV1 - vascular, HBA1 -
% erythrocytes, APOD,TCF21 - stromal, XAGE2,EGFR - trophoblasts
% % 
%% gene expression 
g = 'spp1';
gx = data_orig_all_sorted(strcmpi(geneid_all,g),:);
figure('Position',[100,100,1400,300],'Color','w');
bar(gx); axis tight; hold on;
yl = max(get(gca,'ylim'));
linewid =0.5;
for jj=1:length(cells_bor)
    plot(cells_bor(jj)*[1,1]-0.5,[0,yl],'--','linewidth',linewid,'color',[0.6,0.6,0.6])
end
set(gca,'xtick',gr_center,'xticklabel',[1:length(gr_center)], 'fontsize', 8)
ylabel(g);
% % % % % % % % % % 





