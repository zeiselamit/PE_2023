tic
clear all
close all
addpath('/data/matlab_functions/')
terminate(pyenv)
pyenv("ExecutionMode","OutOfProcess")

savefig_flag = 0;
savefig_pdf = 0;

load PE_QC800_400_afterloading_19-Oct-2021.mat

%%
cellidclsuter_imm = loadCellFile('immune_step3_cellid_cluster_beforeclean_04-Aug-2022.txt');
cellidclsuter_str = loadCellFile('stromal_step3_cellid_cluster_beforeclean_20-Nov-2022.txt');
cellidclsuter_tb = loadCellFile('trophoblasts_step3_cellid_cluster_beforeclean_17-Oct-2021.txt');
cellidclsuter_vas = loadCellFile('vascular_step3_cellid_cluster_beforeclean_17-Oct-2021.txt');

% geneselection_TB = loadCellFile('gene_correlation_step4_TB_15-Aug-2022.txt');
% geneselection_immune = loadCellFile('gene_correlation_step4_immune_15-Aug-2022.txt');
% geneselection_stromal = loadCellFile('gene_correlation_step4_stromal_15-Aug-2022.txt');
% geneselection_vascular = loadCellFile('gene_correlation_step4_vascular_15-Aug-2022.txt');


imm_clustername = loadCellFile('cluster_name_immune_221120.txt');
imm_clustername(:,2) = cellfun(@(x) ['IMMUNE_',x], imm_clustername(:,2),'UniformOutput',0);
imm_clustername(:,2) = regexprep(imm_clustername(:,2),'IMMUNE_EXC','EXC');
imm_clustername = imm_clustername(2:end,:);
str_clustername = loadCellFile('cluster_name_stromal_221120.txt');
str_clustername(:,2) = cellfun(@(x) ['STROMAL_',x], str_clustername(:,2),'UniformOutput',0);
str_clustername(:,2) = regexprep(str_clustername(:,2),'STROMAL_EXC','EXC');
str_clustername = str_clustername(2:end,:);
tb_clustername = loadCellFile('cluster_name_trophoblast_221120.txt');
tb_clustername(:,2) = cellfun(@(x) ['TB_',x], tb_clustername(:,2),'UniformOutput',0);
tb_clustername(:,2) = regexprep(tb_clustername(:,2),'TB_EXC','EXC');
tb_clustername = tb_clustername(2:end,:);
vas_clustername = loadCellFile('cluster_name_vascular_221120.txt');
vas_clustername(:,2) = cellfun(@(x) ['VASCULAR_',x], vas_clustername(:,2),'UniformOutput',0);
vas_clustername(:,2) = regexprep(vas_clustername(:,2),'VASCULAR_EXC','EXC');
vas_clustername = vas_clustername(2:end,:);


cellidclsuter_imm = [cellidclsuter_imm,cell(length(cellidclsuter_imm),1)];
c = cell2mat(cellidclsuter_imm(:,2));
for i=1:length(imm_clustername)
    ind = find(c==imm_clustername{i});
    cellidclsuter_imm(ind,3) = repmat(imm_clustername(i,2),length(ind),1);
end

cellidclsuter_str = [cellidclsuter_str,cell(length(cellidclsuter_str),1)];
c = cell2mat(cellidclsuter_str(:,2));
for i=1:length(str_clustername)
    ind = find(c==str_clustername{i});
    cellidclsuter_str(ind,3) = repmat(str_clustername(i,2),length(ind),1);
end

cellidclsuter_tb = [cellidclsuter_tb,cell(length(cellidclsuter_tb),1)];
c = cell2mat(cellidclsuter_tb(:,2));
for i=1:length(tb_clustername)
    ind = find(c==tb_clustername{i});
    cellidclsuter_tb(ind,3) = repmat(tb_clustername(i,2),length(ind),1);
end

cellidclsuter_vas = [cellidclsuter_vas,cell(length(cellidclsuter_vas),1)];
c = cell2mat(cellidclsuter_vas(:,2));
for i=1:length(vas_clustername)
    ind = find(c==vas_clustername{i});
    cellidclsuter_vas(ind,3) = repmat(vas_clustername(i,2),length(ind),1);
end


cellid_clusters = cell(length(cellid),3);
[~,loc]= ismember(cellidclsuter_imm(:,1),cellid);
cellid_clusters(loc(loc>0),:) = cellidclsuter_imm;
[~,loc]= ismember(cellidclsuter_str(:,1),cellid);
cellid_clusters(loc(loc>0),:) = cellidclsuter_str;
[~,loc]= ismember(cellidclsuter_tb(:,1),cellid);
cellid_clusters(loc(loc>0),:) = cellidclsuter_tb;
[~,loc]= ismember(cellidclsuter_vas(:,1),cellid);
cellid_clusters(loc(loc>0),:) = cellidclsuter_vas;
rmv = find(strcmpi(cellid_clusters(:,3),'EXC'));
% ind_heavy_metal_exc = find(strcmpi(sample,'105-1') & strcmpi(cellid_clusters(:,3),'STROMAL_HEAVY_METAL'));
% rmv = [rmv;ind_heavy_metal_exc];% excluding cells from heavy metal cluster sample 105-1
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



data = ceil(data./repmat(sum(data),length(data(:,1)),1)*10e3); %normalization 

donorid_uni = {'control_early_1','control_early_2','control_early_3',...
    'PE_early_1','PE_early_2','PE_early_3','PE_early_4','PE_early_5','PE_early_6','PE_early_7','PE_early_8','PE_early_9','PE_early_10',....
    'control_3','control_4','control_5','control_6','control_7','control_8',....
    'PE_late_1','PE_late_2','PE_late_3','PE_late_4','PE_late_5','PE_late_6','PE_late_7'};
donorid_uni_group = [1,1,1,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4,4];

%% color per cell type
% color_celltypes_ordered = loadCellFile_turbo('color_celltypes_ordered_HEX_221120.txt',1);
color_celltypes_ordered = loadCellFile_turbo('treeOrder_color_celltypes_ordered_HEX_20-Nov-2022.txt',1);
clusteruni = color_celltypes_ordered(1:46,1);%unique(cellid_clusters(:,3));
T_cells_tmp = zeros(length(cellid),1);
clustersize = zeros(length(clusteruni),1);
for i=1:length(clusteruni)
    
    ind = find(strcmpi(cellid_clusters(:,3),clusteruni{i}));
    T_cells_tmp(ind) = i;
    clustersize(i) = length(ind);
end
T_cells_tmp_uni = unique(T_cells_tmp);
[clusteruni,m2c(clustersize)]
immune_clusters = find(contains(clusteruni,'IMMUNE_') & ~contains(clusteruni,'STROMAL_'));
stromal_clusters = find(contains(clusteruni,'STROMAL_'));
tb_clusters = find(contains(clusteruni,'TB_'));
vas_clusters = find(contains(clusteruni,'VASCULAR_') & ~contains(clusteruni,'STROMAL_'));
%%
mean_mat = zeros(length(geneid),length(T_cells_tmp_uni)*4);
mean_lin_mat = zeros(length(geneid),length(T_cells_tmp_uni)*4);
median_mat = zeros(length(geneid),length(T_cells_tmp_uni)*4);
fracpos_mat = zeros(length(geneid),length(T_cells_tmp_uni)*4);
p75_mat = zeros(length(geneid),length(T_cells_tmp_uni)*4);
condgrsize = zeros(length(T_cells_tmp_uni)*4,1);
condition_sorted = zeros(4,length(peearly_flag));
condition_sorted(1,ctrlearly_flag) = 1;
condition_sorted(2,peearly_flag) = 1;
condition_sorted(3,ctrl_flag) = 1;
condition_sorted(4,pelate_flag) = 1;
% condition_sorted(5,iugr_flag) = 1;
% condition_sorted(4,(peearly_flag | pelate_flag) & ~iugr_flag) = 1;
datalog = log2(data+1);
i= 0;
for k=1:length(T_cells_tmp_uni)
    for jj=1:4
        i = i+1;
        gr1 = find(T_cells_tmp==k & condition_sorted(jj,:)'==1);
        condgrsize(i) = length(gr1);
        mean_mat(:,i) = mean(datalog(:,gr1),2);
        mean_lin_mat(:,i) = mean(data(:,gr1),2);
        median_mat(:,i) = median(datalog(:,gr1),2);
        fracpos_mat(:,i) = mean(data(:,gr1)>0,2);
        p75_mat(:,i) = prctile(datalog(:,gr1),75,2);
    end
end
%% bin anaysis
class_type = 3; %1=immune, 2=stromal, 3=tb, 4= vascular
if class_type==1
    typeplot = immune_clusters;
elseif class_type==2
    typeplot = stromal_clusters;
elseif class_type==3
    typeplot = tb_clusters;
elseif class_type==4
    typeplot = vas_clusters;
elseif class_type==5
    typeplot = T_cells_tmp_uni;
end
figure;
set(gcf,'color','w','position',[20,20,1200,1000]);
if class_type==5
    [ha, pos] = tight_subplot(6, ceil(length(typeplot)/6), [0.08,0.08], [0.04,0.04], [0.05,0.08]);
else
    [ha, pos] = tight_subplot(4, ceil(length(typeplot)/4), [0.08,0.08], [0.04,0.04], [0.05,0.08]);
end
sig_gene_cluster = cell(length(typeplot),1);
th = 0.5;
list = [];
geneper_bin = zeros(length(typeplot),8);
for j=1:length(typeplot)
    
        gr1 = find( (T_cells_tmp == typeplot(j)) & peearly_flag==1);
        gr2 = find( (T_cells_tmp == typeplot(j) ) & pelate_flag==1);
        gr3 = find( (T_cells_tmp == typeplot(j)) & ctrlearly_flag==1);
        gr4 = find( (T_cells_tmp == typeplot(j) ) & ctrl_flag==1);
        x1 = mean(log2(data(:,gr1)+1),2);
        x2 = mean(log2(data(:,gr2)+1),2);
        x3 = mean(log2(data(:,gr3)+1),2);
        x4 = mean(log2(data(:,gr4)+1),2);
        dx = x1-x3 ;
        dy = x2-x4 ;
        dis = sqrt(dx.^2 + dy.^2);
        
        indsig = find(dis>th);
        psig_x = ones(size(indsig));
        psig_y = ones(size(indsig));
        for jj=1:length(indsig)
            psig_x(jj) = ranksum(data(indsig(jj),gr1),data(indsig(jj),gr3));
            psig_y(jj) = ranksum(data(indsig(jj),gr2),data(indsig(jj),gr4));
        end
        q = qval_from_pval([psig_x;psig_y])';
        q = min(reshape(q,length(indsig),2),[],2);
        indsig = indsig(q<0.1);
        sig_gene_cluster{j} = geneid(indsig);%[dx(dis>th),dy(dis>th)];
        list = [list;sig_gene_cluster{j}];
        [~,xi] = sort(dis,'descend');
        theta = atan2(dy(indsig),dx(indsig))*180/pi;
        theta(theta<0) = 360+theta(theta<0);
        anglebins = [[0:45:(360)]-22.5,360];
        [tmp,~,bins] = histcounts(theta,anglebins);
        tmp = [sum(tmp([1,9])),tmp(2:8)];
        geneper_bin(j,:) = [tmp];
        bins(bins==9) = 1;
        tmp2 = cell(max(tmp),8);
        for jj=1:8
            tmp2(1:tmp(jj),jj) = sig_gene_cluster{j}(bins==jj);
        end
        axes(ha(j))
        %     polarhistogram(atan2(dy(dis>th),dx(dis>th)),100);
        plot(dx(indsig),dy(indsig),'.'); hold on;
        axis tight
        xl = get(ha(j),'xlim');
        plot([-2.5,2.5],[-2.5,2.5],'-k');
        plot([-2.5,2.5],tan(22.5/180*pi)*[-2.5,2.5],'-r');
        plot([-2.5,2.5],tan(67.5/180*pi)*[-2.5,2.5],'-r');
        plot([-2.5,2.5],tan(-22.5/180*pi)*[-2.5,2.5],'-r');
        plot([-2.5,2.5],tan(-67.5/180*pi)*[-2.5,2.5],'-r');
        plot([-2.5,2.5],0*[-2.5,2.5],'-k');
        plot(0*[-2.5,2.5],[-2.5,2.5],'-k');
        axis equal
        set(ha(j),'xlim',[-1.5,1.5],'ylim',[-1.5,1.5])
        %     text(dx(xi(1:ntop)),dy(xi(1:ntop)),geneid_all(xi(1:ntop)),'fontsize',8);
        text(dx(indsig),dy(indsig),geneid(indsig),'fontsize',8);
        tmp = m2c(geneper_bin(j,:));
        tmp = cellfun(@num2str,tmp,'UniformOutput',0);
        text(th*cos((anglebins(1:8)+22.5)*pi/180),th*sin((anglebins(1:8)+22.5)*pi/180),tmp,'fontsize',8)
        [r,pr] = corr(dx,dy)
        c = strfind(clusteruni{typeplot(j)},'_');
        tl = [regexprep(clusteruni{typeplot(j)}(c(1)+1:end),'_','-'),',N=',num2str(clustersize(typeplot(j)))]; 
        title(tl,'fontsize',8,'fontname','arial')
       
        %     xlabel('pe-ce');
        %     ylabel('pl-ct');
    
end
figure('color','w','position',[20,20,800,1000]);
imagesc(cent_norm(geneper_bin));
set(gca,'xtick',[1:8],'ytick',[1:length(typeplot)],'YTickLabel',clusteruni(typeplot),'TickLabelInterpreter','none')

list = unique(list);
[~,loc] = ismember(list,geneid);
tmpmat = p75_mat(loc,[(typeplot(1)-1)*4+1:4*typeplot(end)]);
for i=1:length(typeplot)
    ind = (i-1)*4+1:i*4;
    tmpmat(:,ind) = tmpmat(:,ind)./repmat(mean(tmpmat(:,ind),2),1,4);
end
tmpmat(isnan(tmpmat)) = 0;

z = linkage(tmpmat,'average','euclidean');
% d = corr_mat(tmpmat');
leaforder = optimalleaforder(z,pdist(tmpmat));
% leaforder = optimalleaforder(z,squareform(1-d,'tovector'));

figure('color','w','position',[20,20,800,1000]);
imagesc(tmpmat(leaforder,:));
set(gca,'xtick',[1:4*length(typeplot)],'ytick',[1:length(list)],'YTickLabel',list(leaforder),'TickLabelInterpreter','none')

%% gene expression bars
class_type = 3; %1=immune, 2=stromal, 3=tb, 4= vascular
if class_type==1
    typeplot = immune_clusters;
elseif class_type==2
    typeplot = stromal_clusters;
elseif class_type==3
    typeplot = tb_clusters;
elseif class_type==4
    typeplot = vas_clusters;
elseif class_type==5
    typeplot = T_cells_tmp_uni;
end
% typeplot = tb_clusters;
% gn_list = {'CXCL9','IDO1','SEPP1','RPS26','RPS10'};
gn_list = {'flt1','pgf'};
% gn_list = {'SPP1','CTSD','CD163','C1QA','C1QB','C1QC','CCL3L3','FTL','STAB1','CD74'};%{'IGFBP5','IGFBP3','COL4A1','COL4A2','COL6A1','COL6A2'};%
figure('color','w','position',[20,20,800,length(gn_list)*200]);
log_lin_flag = 1;
[ha, pos] = tight_subplot(length(gn_list),1, [0.01,0.01], [0.1,0.04], [0.1,0.02]);
for jj=1:length(gn_list)
    gn = gn_list{jj};
    g = find(strcmpi(geneid,gn));
    if log_lin_flag==1
    y = mean_mat(g,[(typeplot(1)-1)*4+1:4*typeplot(end)]);
    t = [y([2:4:end])-y([1:4:end]);y([4:4:end])-y([3:4:end])];
    elseif log_lin_flag==2
        y = mean_lin_mat(g,[(typeplot(1)-1)*4+1:4*typeplot(end)]);
    t = log2([y([2:4:end])./y([1:4:end]);y([4:4:end])./y([3:4:end])]);
    end

    axes(ha(jj))
    for k=1:length(typeplot)
%     for k=1:9
        bar(k,t(:,k),'FaceColor',['#',color_celltypes_ordered{typeplot(k),2}]);hold on;
    end
    axis tight
    p99 = get(gca,'ylim');
    set(gca,'ylim',[p99(1),p99(2)*1.1])
    p99 = p99(2);
    comp = [[1,2];[3,4]];
    pcomp = ones(length(comp(:,1)),1);
    for c=1:length(typeplot)
        for jjj=1:length(comp)
            gr1 = find(T_cells_tmp==typeplot(c) & condition_sorted(comp(jjj,1),:)'==1);
            gr2 = find(T_cells_tmp==typeplot(c) & condition_sorted(comp(jjj,2),:)'==1);
            y1 = (data(g,gr1));
            y2 = (data(g,gr2));
            pcomp(jjj) = ranksum(y1,y2);
            if pcomp(jjj)<1e-4
                sigstar='***';
            elseif pcomp(jjj)<1e-3
                sigstar='**';
            elseif pcomp(jjj)<1e-2
                sigstar='*';
            else
                sigstar='';
            end
            %         plot(k ,p99*(0.93+jjj*0.01)*[1,1],'-','color',0.5*[1,1,1]);
            if jjj==1
                text(c-0.1 ,p99*1,sigstar,'fontsize',7,'VerticalAlignment','bottom','HorizontalAlignment','center','color','g');
            elseif jjj==2
                text(c+0.1 ,p99*1,sigstar,'fontsize',7,'VerticalAlignment','bottom','HorizontalAlignment','center','color','r');
            end
        end
    end


    ylabel([gn,', diff log'],'fontsize',18,'FontWeight','bold')
    set(gca,'xtick',[1:length(typeplot)],'XTickLabel',cell(length(typeplot),1),'fontsize',6)
end
set(gca,'xtick',[1:length(typeplot)],'XTickLabel',clusteruni(typeplot),'TickLabelInterpreter','none','XTickLabelRotation',45,'fontsize',6)
% legend({'PE-CE','PL-CT'})
% ylabel('FC')
% linkaxes(ha,'y')
eval(['export_fig ST_2_FC_bars_',date,'.pdf']);   

%%

mean_mat_donor = repmat({zeros(length(geneid),length(donorid_uni))},length(T_cells_tmp_uni),1);
grsize_donor = zeros(length(donorid_uni),length(T_cells_tmp_uni));%repmat({zeros(length(donorid_uni),1)},length(T_cells_tmp_uni),1);
    
for k=1:length(T_cells_tmp_uni)
    k
    tmp = zeros(length(geneid),length(donorid_uni));
    condgrsize = zeros(length(donorid_uni),1);
    for kk=1:length(donorid_uni)
        kk
        in = strcmpi(donor_id_flag,donorid_uni{kk});
        gr1 = find(T_cells_tmp==k & in);
        if ~isempty(gr1)
            condgrsize(kk) = length(gr1);
            tmp(:,kk) = mean(log2(data(:,gr1)+1),2);
        end
    end
    mean_mat_donor{k} = tmp;
    grsize_donor(:,k) = condgrsize;
end


%%
class_type = 5; %1=immune, 2=stromal, 3=tb, 4= vascular
if class_type==1
    typeplot = immune_clusters;
elseif class_type==2
    typeplot = stromal_clusters;
elseif class_type==3
    typeplot = tb_clusters;
elseif class_type==4
    typeplot = vas_clusters;
elseif class_type==5
    typeplot = T_cells_tmp_uni;
end
% typeplot = tb_clusters;
gn = 'FLT1';
g = find(strcmpi(geneid,gn));
figure('color','w','position',[20,20,800,500]);
xt = [];
donorid_cell = m2c([1:26]);
donorid_cell = cellfun(@num2str,donorid_cell,'UniformOutput',0);
for k=1:length(typeplot)
    y = mean_mat_donor{typeplot(k)}(g,:);
    for kk=1:4
        yy(kk) = median(y(donorid_uni_group==kk));
    end
    bar(k+[1:4]/5-0.5,yy,'FaceColor',['#',color_celltypes_ordered{typeplot(k),2}]);hold on;
    xt = [xt,k+[1:4]/5-0.5];
    plot(k+donorid_uni_group/5-0.5,y,'+','color','k','markersize',2);hold on;
    text(k+donorid_uni_group/5-0.5,y,donorid_cell,'fontsize',7);hold on;
    %     plot(k+donorid_uni_group/5-0.5,y,'+','color',['#',color_celltypes_ordered{typeplot(k),2}]);hold on;
end
axis tight
yl = get(gca,'ylim');
set(gca,'xtick',[1:length(typeplot)],'XTickLabel',clusteruni(typeplot),'TickLabelInterpreter'....
    ,'none','XTickLabelRotation',45,'ylim',[-(yl(2)-yl(1))*0.1,yl(2)]);
text(xt,-(yl(2)-yl(1))*0.05*ones(size(xt)),repmat({'ce','pe','ct','pl'},1,length(typeplot)),'fontsize',5)
% legend({'PE-CE','PL-CT'})
ylabel('Expression (log2+1)')
title(gn)

%% gene expression 
class_type = 3; %1=immune, 2=stromal, 3=tb, 4= vascular
if class_type==1
    typeplot = immune_clusters;
elseif class_type==2
    typeplot = stromal_clusters;
elseif class_type==3
    typeplot = tb_clusters;
elseif class_type==4
    typeplot = vas_clusters;
elseif class_type==5
    typeplot = T_cells_tmp_uni;
end
% typeplot = tb_clusters;
gn = 'FLT1';
g = find(strcmpi(geneid,gn));
figure('color','w','position',[20,20,800,500]);
xt = [];
for k=1:length(typeplot)
    y = mean_mat_donor{typeplot(k)}(g,:);
    for kk=1:4
        yy(kk) = median(y(donorid_uni_group==kk));
    end
    bar(k+[1:2]/4-0.5,[yy(2)-yy(1),yy(4)-yy(3)],'FaceColor',['#',color_celltypes_ordered{typeplot(k),2}]);hold on;
%     bar(k+[1:4]/5-0.5,yy,'FaceColor',['#',color_celltypes_ordered{typeplot(k),2}]);hold on;
%     xt = [xt,k+[1:4]/5-0.5];
%     plot(k+donorid_uni_group/5-0.5,y,'+','color','k');hold on;
    %     plot(k+donorid_uni_group/5-0.5,y,'+','color',['#',color_celltypes_ordered{typeplot(k),2}]);hold on;
end
axis tight
yl = get(gca,'ylim');
set(gca,'xtick',[1:length(typeplot)],'XTickLabel',clusteruni(typeplot),'TickLabelInterpreter'....
    ,'none','XTickLabelRotation',45,'ylim',[-(yl(2)-yl(1))*0.1,yl(2)]);
% text(xt,-(yl(2)-yl(1))*0.05*ones(size(xt)),repmat({'ce','pe','ct','pl'},1,length(typeplot)),'fontsize',5)
% legend({'PE-CE','PL-CT'})
ylabel('FC')
title(gn)


%% bin analysis
class_type = 2; %1=immune, 2=stromal, 3=tb, 4= vascular
if class_type==1
    typeplot = immune_clusters;
elseif class_type==2
    typeplot = stromal_clusters;
elseif class_type==3
    typeplot = tb_clusters;
elseif class_type==4
    typeplot = vas_clusters;
elseif class_type==5
    typeplot = T_cells_tmp_uni;
end
figure;
set(gcf,'color','w','position',[20,20,1200,1000]);
if class_type==5
    [ha, pos] = tight_subplot(6, ceil(length(typeplot)/6), [0.08,0.08], [0.04,0.04], [0.05,0.08]);
else
    [ha, pos] = tight_subplot(4, ceil(length(typeplot)/4), [0.08,0.08], [0.04,0.04], [0.05,0.08]);
end
sig_gene_cluster = cell(length(typeplot),1);
th = 0.7;
list = [];
geneper_bin = zeros(length(typeplot),8);
pmat_early_up = ones(length(func_genes(:,1)),length(typeplot));
pmat_early_d = pmat_early_up;
pmat_late_up = pmat_early_up;
pmat_late_d = pmat_early_up;
for j=1:length(typeplot)

    gr1 = find( (T_cells_tmp == typeplot(j)) & peearly_flag==1);
    gr2 = find( (T_cells_tmp == typeplot(j) ) & pelate_flag==1);
    gr3 = find( (T_cells_tmp == typeplot(j)) & ctrlearly_flag==1);
    gr4 = find( (T_cells_tmp == typeplot(j) ) & ctrl_flag==1);
    x1 = mean(log2(data(:,gr1)+1),2);
    x2 = mean(log2(data(:,gr2)+1),2);
    x3 = mean(log2(data(:,gr3)+1),2);
    x4 = mean(log2(data(:,gr4)+1),2);
    %         x11 = mean((data(:,gr1)),2);
    %         x22 = mean((data(:,gr2)),2);
    %         x33 = mean((data(:,gr3)),2);
    %         x44 = mean((data(:,gr4)),2);
    dx = x1-x3 ;
    dy = x2-x4 ;
    dis = sqrt(dx.^2 + dy.^2);
    m = max([x1,x2,x3,x4],[],2);

    indsig = find(dis>th & m>3);
    psig_x = ones(size(indsig));
    psig_y = ones(size(indsig));
    for jj=1:length(indsig)
        psig_x(jj) = ranksum(data(indsig(jj),gr1),data(indsig(jj),gr3));
        psig_y(jj) = ranksum(data(indsig(jj),gr2),data(indsig(jj),gr4));
    end
    q = qval_from_pval([psig_x;psig_y])';
    q = min(reshape(q,length(indsig),2),[],2);
    indsig = indsig(q<0.1);
    sig_gene_cluster{j} = geneid(indsig);%[dx(dis>th),dy(dis>th)];
    list = [list;sig_gene_cluster{j}];
    [~,xi] = sort(dis,'descend');
    theta = atan2(dy(indsig),dx(indsig))*180/pi;
    theta(theta<0) = 360+theta(theta<0);
    anglebins = [[0:45:(360)]-22.5,360];
    [tmp,~,bins] = histcounts(theta,anglebins);
    tmp = [sum(tmp([1,9])),tmp(2:8)];
    geneper_bin(j,:) = [tmp];
    bins(bins==9) = 1;
    tmp2 = cell(max(tmp),8);
    for jj=1:8
        tmp2(1:tmp(jj),jj) = sig_gene_cluster{j}(bins==jj);
    end
    axes(ha(j))
    %     polarhistogram(atan2(dy(dis>th),dx(dis>th)),100);
    plot(dx(indsig),dy(indsig),'.'); hold on;
    axis tight
    xl = get(ha(j),'xlim');
    plot([-2.5,2.5],[-2.5,2.5],'-k');
    plot([-2.5,2.5],tan(22.5/180*pi)*[-2.5,2.5],'-r');
    plot([-2.5,2.5],tan(67.5/180*pi)*[-2.5,2.5],'-r');
    plot([-2.5,2.5],tan(-22.5/180*pi)*[-2.5,2.5],'-r');
    plot([-2.5,2.5],tan(-67.5/180*pi)*[-2.5,2.5],'-r');
    plot([-2.5,2.5],0*[-2.5,2.5],'-k');
    plot(0*[-2.5,2.5],[-2.5,2.5],'-k');
    axis equal
    set(ha(j),'xlim',[-1.5,1.5],'ylim',[-1.5,1.5])
    %     text(dx(xi(1:ntop)),dy(xi(1:ntop)),geneid_all(xi(1:ntop)),'fontsize',8);
    text(dx(indsig),dy(indsig),geneid(indsig),'fontsize',8);
    tmp = m2c(geneper_bin(j,:));
    tmp = cellfun(@num2str,tmp,'UniformOutput',0);
    text(th*cos((anglebins(1:8)+22.5)*pi/180),th*sin((anglebins(1:8)+22.5)*pi/180),tmp,'fontsize',8)
    [r,pr] = corr(dx,dy)
    c = strfind(clusteruni{typeplot(j)},'_');
    tl = [regexprep(clusteruni{typeplot(j)}(c(1)+1:end),'_','-'),',N=',num2str(clustersize(typeplot(j)))];
    title(tl,'fontsize',8,'fontname','arial')

    %     xlabel('pe-ce');
    %     ylabel('pl-ct');
    ind1 = find(x1>1 | x3>1); %early
    ind2 = find(x2>1 | x4>1); %late
    [rs1,xi1] = sort(dx(ind1),'descend');
    rk1 = geneid(ind1(xi1));
    rk1_up = rk1(rs1>0);
    rs1_up = rs1(rs1>0);
    rk1_d = flipud(rk1(rs1<0));
    rs1_d = flipud(-rs1(rs1<0));
    [rs2,xi2] = sort(dy(ind2),'descend');
    rk2 = geneid(ind2(xi2));
    rk2_up = rk2(rs2>0);
    rs2_up = rs2(rs2>0);
    rk2_d = flipud(rk2(rs2<0));
    rs2_d = flipud(-rs2(rs2<0));
    func_pv = [];
    np = 10000;
    w = 1;
    for ii = 1:size(func_genes,1)
        ii
        tmpset = func_genes(ii,2:end);
        tmpset = tmpset(cellfun(@(x) ~isempty(x),tmpset));
        %             early
        [~,gs1_pos_up] = ismember(tmpset, rk1_up);
        gs1_pos_up(gs1_pos_up==0) = [];
        gs1_neg_up = setdiff(rk1_up,rk1_up(gs1_pos_up));
        [~,gs1_pos_d] = ismember(tmpset, rk1_d);
        gs1_pos_d(gs1_pos_d==0) = [];
        gs1_neg_d = setdiff(rk1_d,rk1_d(gs1_pos_d));
        if length(gs1_pos_up)<5
            pv1_up = 1;
        else
            [es1_up,nes1_up,pv1_up,ledge1_up] = gsea2(rk1_up,rs1_up,rk1_up(gs1_pos_up),gs1_neg_up,np,w);
        end
        if length(gs1_pos_d)<5
            pv1_d = 1;
        else
            [es1_d,nes1_d,pv1_d,ledge1_d] = gsea2(rk1_d,rs1_d,rk1_d(gs1_pos_d),gs1_neg_d,np,w);
        end
        %             late
        [~,gs2_pos_up] = ismember(tmpset, rk2_up);
        gs2_pos_up(gs2_pos_up==0) = [];
        gs2_neg_up = setdiff(rk2_up,rk2_up(gs2_pos_up));
        [~,gs2_pos_d] = ismember(tmpset, rk2_d);
        gs2_pos_d(gs2_pos_d==0) = [];
        gs2_neg_d = setdiff(rk2_d,rk2_d(gs2_pos_d));
        if length(gs2_pos_up)<5
            pv2_up = 1;
        else
            [es2_up,nes2_up,pv2_up,ledge2_up] = gsea2(rk2_up,rs2_up,rk2_up(gs2_pos_up),gs2_neg_up,np,w);
        end
        if length(gs2_pos_d)<5
            pv2_d = 1;
        else
            [es2_d,nes2_d,pv2_d,ledge2_d] = gsea2(rk2_d,rs2_d,rk2_d(gs2_pos_d),gs2_neg_d,np,w);
        end


        %             [~,gs11] = ismember(tmpset, rk1);
        %             gs11(gs11==0) = [];
        %             gs21 = setdiff(rk1,rk1(gs11));
        %             [~,gs12] = ismember(tmpset, rk2);
        %             gs12(gs12==0) = [];
        %             gs22 = setdiff(rk2,rk2(gs12));
        %
        %             [es1,nes1,pv1,ledge1] = gsea2(rk1,rs1,rk1(gs11),gs21,np,w); %rk- list of 1530.. genes, rs -
        %             [es2,nes2,pv2,ledge2] = gsea2(rk2,rs2,rk2(gs12),gs22,np,w);


        %             func_pv = [func_pv;func_genes(ii,1),pv1,pv2,nes1,nes2];
        pmat_early_up(ii,j) = pv1_up;
        pmat_early_d(ii,j) = pv1_d;
        pmat_late_up(ii,j) = pv2_up;
        pmat_late_d(ii,j) = pv2_d;
    end
    func_pv_cell{j} = func_pv;
end
table_proc = [[{''};func_genes(:,1)],[clusteruni(typeplot)';m2c(pmat_early_up)]];
figure('position',[100,100,800,1000],'color','w');
h1 = axes('position',[0.2,0.2,0.35,0.75]);
imagesc(round(-log10(pmat_early_up+1e-5)),[0,5]);
colormap('jet')
colorbar
set(h1,'YTick',[1:length(func_genes(:,1))],'YTickLabel',func_genes(:,1)....
    ,'xTick',[1:length(typeplot)],'xTickLabel',clusteruni(typeplot),'TickLabelInterpreter','none','fontsize',8)
title('early-up')
% xline([17.5,30.5,38.5],'w','LineWidth',1); % for all only!
h2 = axes('position',[0.6,0.2,0.35,0.75]);
imagesc(round(-log10(pmat_late_up+1e-5)),[0,5]);
colormap('jet')
colorbar
set(h2,'YTick',[1:length(func_genes(:,1))],'YTickLabel',cell(size(func_genes(:,1)))....
    ,'xTick',[1:length(typeplot)],'xTickLabel',clusteruni(typeplot),'TickLabelInterpreter','none','fontsize',8)
title('late-up')
% xline([17.5,30.5,38.5],'w','LineWidth',1); % for all only!
linkaxes([h1,h2],'xy')

figure('position',[100,100,800,1000],'color','w');
h1 = axes('position',[0.2,0.2,0.35,0.75]);
imagesc(round(-log10(pmat_early_d+1e-5)),[0,5]);
colormap('jet')
colorbar
set(h1,'YTick',[1:length(func_genes(:,1))],'YTickLabel',func_genes(:,1)....
    ,'xTick',[1:length(typeplot)],'xTickLabel',clusteruni(typeplot),'TickLabelInterpreter','none','fontsize',8)
title('early-down')
% xline([17.5,30.5,38.5],'w','LineWidth',1); % for all only!
h2 = axes('position',[0.6,0.2,0.35,0.75]);
imagesc(round(-log10(pmat_late_d+1e-5)),[0,5]);
colormap('jet')
colorbar
set(h2,'YTick',[1:length(func_genes(:,1))],'YTickLabel',cell(size(func_genes(:,1)))....
    ,'xTick',[1:length(typeplot)],'xTickLabel',clusteruni(typeplot),'TickLabelInterpreter','none','fontsize',8)
title('late-down')
% xline([17.5,30.5,38.5],'w','LineWidth',1); % for all only!
linkaxes([h1,h2],'xy')

figure('position',[100,100,800,1000],'color','w');
h1 = axes('position',[0.2,0.2,0.35,0.75]);
imagesc(-log10(pmat_early_up+1e-3)+log10(pmat_late_up+1e-3));
colormap('summer')
colorbar
set(h1,'YTick',[1:length(func_genes(:,1))],'YTickLabel',func_genes(:,1)....
    ,'xTick',[1:length(typeplot)],'xTickLabel',clusteruni(typeplot),'TickLabelInterpreter','none','fontsize',8)
title('early-late (up)')
h2 = axes('position',[0.6,0.2,0.35,0.75]);
imagesc(-log10(pmat_early_d+1e-3)+log10(pmat_late_d+1e-3));
colormap('summer')
colorbar
set(h2,'YTick',[1:length(func_genes(:,1))],'YTickLabel',cell(size(func_genes(:,1)))....
    ,'xTick',[1:length(typeplot)],'xTickLabel',clusteruni(typeplot),'TickLabelInterpreter','none','fontsize',8)
title('early-late (down)')
linkaxes([h1,h2],'xy')


figure('position',[100,100,800,1000],'color','w');
h1 = axes('position',[0.2,0.2,0.35,0.75]);
imagesc(-log10(pmat_early_up),[0,5]);
colormap('summer')
colorbar
set(h1,'YTick',[1:length(func_genes(:,1))],'YTickLabel',func_genes(:,1)....
    ,'xTick',[1:length(typeplot)],'xTickLabel',clusteruni(typeplot),'TickLabelInterpreter','none','fontsize',8)
title('early')
h2 = axes('position',[0.6,0.2,0.35,0.75]);
imagesc(-log10(pmat_late_up),[0,5]);
colormap('summer')
colorbar
set(h2,'YTick',[1:length(func_genes(:,1))],'YTickLabel',cell(size(func_genes(:,1)))....
    ,'xTick',[1:length(typeplot)],'xTickLabel',clusteruni(typeplot),'TickLabelInterpreter','none','fontsize',8)
title('late')
linkaxes([h1,h2],'xy')


%%
figure('color','w','position',[20,20,600,1000]);
h1 = axes('position',[0.05,0.05,0.65,0.8]);
imagesc(cent_norm(geneper_bin));
colormap('summer')
set(gca,'xtick',[1:8],'ytick',[1:length(typeplot)],'YTickLabel',clusteruni(typeplot),'TickLabelInterpreter','none','YDir','normal')
h2 = axes('position',[0.72,0.05,0.01,0.8]);
scatter(ones(size(typeplot)),[1:length(typeplot)],clustersize*0.05,[0.6,0.6,0.6],'filled');
axis off
h3 = axes('position',[0.75,0.05,0.2,0.8]);
barh(sum(geneper_bin,2),'BarWidth',0.5);
set(gca,'ytick',[1:length(typeplot)],'YTickLabel',cell(size(typeplot)),'TickLabelInterpreter','none')
box off
h4 = axes('position',[0.05,0.86,0.65,0.13]);
bar(sum(geneper_bin),'BarWidth',0.5);
set(gca,'xtick',[1:8],'XTickLabel',cell(1,8),'TickLabelInterpreter','none')
box off
linkaxes([h1,h2,h3],'y')
linkaxes([h1,h4],'x')
axes(h1)
axis tight

%% figure 1h
class_type = 5; %1=immune, 2=stromal, 3=tb, 4= vascular
if class_type==1
    typeplot = immune_clusters;
elseif class_type==2
    typeplot = stromal_clusters;
elseif class_type==3
    typeplot = tb_clusters;
elseif class_type==4
    typeplot = vas_clusters;
elseif class_type==5
    typeplot = T_cells_tmp_uni;
end

% tree_order = loadCellFile('celltype_order_tree_fig1.txt');
% tree_order = regexprep(tree_order,'-','_');
% clusteruni1 = regexprep(clusteruni,'-','_');
% [~,loc] = ismember(tree_order,clusteruni1);
% typeplot = T_cells_tmp_uni(loc);
% typeplot = typeplot(end:-1:1);

figure;
set(gcf,'color','w','position',[20,20,1200,1000]);
if class_type==5
    [ha, pos] = tight_subplot(6, ceil(length(typeplot)/6), [0.08,0.08], [0.04,0.04], [0.05,0.08]);
else
    [ha, pos] = tight_subplot(4, ceil(length(typeplot)/4), [0.08,0.08], [0.04,0.04], [0.05,0.08]);
end
sig_gene_cluster = cell(length(typeplot),1);
th = 0.5;
list = [];
geneper_bin = zeros(length(typeplot),8);
% pmat_early_up = ones(length(func_genes(:,1)),length(typeplot));
% pmat_early_d = pmat_early_up;
% pmat_late_up = pmat_early_up;
% pmat_late_d = pmat_early_up;
for j=1:length(typeplot)
    
        gr1 = find( (T_cells_tmp == typeplot(j)) & peearly_flag==1);
        gr2 = find( (T_cells_tmp == typeplot(j) ) & pelate_flag==1);
        gr3 = find( (T_cells_tmp == typeplot(j)) & ctrlearly_flag==1);
        gr4 = find( (T_cells_tmp == typeplot(j) ) & ctrl_flag==1);
        x1 = mean(log2(data(:,gr1)+1),2);
        x2 = mean(log2(data(:,gr2)+1),2);
        x3 = mean(log2(data(:,gr3)+1),2);
        x4 = mean(log2(data(:,gr4)+1),2);
%         x11 = mean((data(:,gr1)),2);
%         x22 = mean((data(:,gr2)),2);
%         x33 = mean((data(:,gr3)),2);
%         x44 = mean((data(:,gr4)),2);
        dx = x1-x3 ;
        dy = x2-x4 ;
        dis = sqrt(dx.^2 + dy.^2);
        m = max([x1,x2,x3,x4],[],2);
        
        indsig = find(dis>th & m>3);
        psig_x = ones(size(indsig));
        psig_y = ones(size(indsig));
        for jj=1:length(indsig)
            psig_x(jj) = ranksum(data(indsig(jj),gr1),data(indsig(jj),gr3));
            psig_y(jj) = ranksum(data(indsig(jj),gr2),data(indsig(jj),gr4));
        end
        q = qval_from_pval([psig_x;psig_y])';
        q = min(reshape(q,length(indsig),2),[],2);
        indsig = indsig(q<0.1);
        sig_gene_cluster{j} = geneid(indsig);%[dx(dis>th),dy(dis>th)];
        list = [list;sig_gene_cluster{j}];
        [~,xi] = sort(dis,'descend');
        theta = atan2(dy(indsig),dx(indsig))*180/pi;
        theta(theta<0) = 360+theta(theta<0);
        anglebins = [[0:45:(360)]-22.5,360];
        [tmp,~,bins] = histcounts(theta,anglebins);
        tmp = [sum(tmp([1,9])),tmp(2:8)];
        geneper_bin(j,:) = [tmp];
        bins(bins==9) = 1;
        tmp2 = cell(max(tmp),8);
        for jj=1:8
            tmp2(1:tmp(jj),jj) = sig_gene_cluster{j}(bins==jj);
        end
        axes(ha(j))
        %     polarhistogram(atan2(dy(dis>th),dx(dis>th)),100);
        plot(dx(indsig),dy(indsig),'.'); hold on;
        axis tight
        xl = get(ha(j),'xlim');
        plot([-2.5,2.5],[-2.5,2.5],'-k');
        plot([-2.5,2.5],tan(22.5/180*pi)*[-2.5,2.5],'-r');
        plot([-2.5,2.5],tan(67.5/180*pi)*[-2.5,2.5],'-r');
        plot([-2.5,2.5],tan(-22.5/180*pi)*[-2.5,2.5],'-r');
        plot([-2.5,2.5],tan(-67.5/180*pi)*[-2.5,2.5],'-r');
        plot([-2.5,2.5],0*[-2.5,2.5],'-k');
        plot(0*[-2.5,2.5],[-2.5,2.5],'-k');
        axis equal
        set(ha(j),'xlim',[-1.5,1.5],'ylim',[-1.5,1.5])
        %     text(dx(xi(1:ntop)),dy(xi(1:ntop)),geneid_all(xi(1:ntop)),'fontsize',8);
        text(dx(indsig),dy(indsig),geneid(indsig),'fontsize',8);
        tmp = m2c(geneper_bin(j,:));
        tmp = cellfun(@num2str,tmp,'UniformOutput',0);
        text(th*cos((anglebins(1:8)+22.5)*pi/180),th*sin((anglebins(1:8)+22.5)*pi/180),tmp,'fontsize',8)
        [r,pr] = corr(dx,dy)
        c = strfind(clusteruni{typeplot(j)},'_');
        tl = [regexprep(clusteruni{typeplot(j)}(c(1)+1:end),'_','-'),',N=',num2str(clustersize(typeplot(j)))]; 
        title(tl,'fontsize',8,'fontname','arial')
       
end

figure('color','w','position',[20,20,500,960]);
h1 = axes('position',[0.2,0.05,0.55,0.8]);
imagesc(cent_norm(geneper_bin));
colormap('summer')
set(gca,'xtick',[1:8],'ytick',[1:length(typeplot)],'YTickLabel',clusteruni(typeplot)....
    ,'TickLabelInterpreter','none','YDir','normal','fontsize',7)
h2 = axes('position',[0.78,0.05,0.01,0.8]);
scatter(ones(size(typeplot)),[1:length(typeplot)],clustersize*0.05,[0.6,0.6,0.6],'filled');
axis off
h3 = axes('position',[0.82,0.05,0.15,0.8]);
barh(sum(geneper_bin,2),'BarWidth',0.5); hold on;
text(sum(geneper_bin,2),[1:length(typeplot)],cellfun(@num2str, m2c(sum(geneper_bin,2)),'uniformoutput',0))
set(gca,'ytick',[1:length(typeplot)],'YTickLabel',cell(size(typeplot)),'TickLabelInterpreter','none')
box off
h4 = axes('position',[0.2,0.86,0.55,0.13]);
bar(sum(geneper_bin),'BarWidth',0.5);hold on;
text([1:8],sum(geneper_bin),cellfun(@num2str, m2c(sum(geneper_bin)),'uniformoutput',0))
set(gca,'xtick',[1:8],'XTickLabel',cell(1,8),'TickLabelInterpreter','none')
box off
linkaxes([h1,h2,h3],'y')
linkaxes([h1,h4],'x')
axes(h1)
axis tight
eval(['export_fig PE_ALL_num_genes_perbin_th05_pece_plc_',date,'.pdf']);
%% anti/inflammatory cytokines analysis.
class_type = 1; %1=immune, 2=stromal, 3=tb, 4= vascular
if class_type==1
    typeplot = immune_clusters;
elseif class_type==2
    typeplot = stromal_clusters;
elseif class_type==3
    typeplot = tb_clusters;
elseif class_type==4
    typeplot = vas_clusters;
elseif class_type==5
    typeplot = T_cells_tmp_uni;
end
typeplot = typeplot([1:9]); %for myeloid
gn_list = loadCellFile('pro_inflammatory_niv_2.csv'); %pro
% gn_list = loadCellFile('anti_inflammatory_niv.csv'); %anti
% gn_list = {'IL1B','CCL2','CXCL8','TNF','SPP1','REL','NFKBIA','EREG'};
% gn_list = {'IL10','CD74','KLF2','IDO1'}
gn_list = upper(gn_list);
[~,loc] = ismember(gn_list,geneid);
gn_list = gn_list(loc>0);
loc = loc(loc>0);
th = 0.3;
y = mean_mat(loc,[(typeplot(1)-1)*4+1:4*typeplot(end)]);
y_early = y(:,[2:4:end])-y(:,[1:4:end]);
y_early(y_early<th & y_early>-th) = 0;
loc_early =  find(sum(y_early,2)~=0);
y_late = y(:,[4:4:end])-y(:,[3:4:end]);
y_late(y_late<th & y_late>-th) = 0;
loc_late =  find(sum(y_late,2)~=0);
[~,loc_el] = ismember(loc_late,loc_early);
loc_el = loc_early(loc_el(loc_el>0));
% tmp = [y_early(loc_el),y_late(loc_el)];
% tmp = cent_norm(tmp);
Y = pdist(y_early(loc_el,:),'correlation');
z = linkage(Y,'ward');
sorted_ind = optimalleaforder(z,Y);
[ii,jj] = find(sign(y_early(loc_el(sorted_ind),:)) == sign(y_late(loc_el(sorted_ind),:)) & sign(y_early(loc_el(sorted_ind),:))~=0);
figure('color','w','position',[20,20,800,600]);
h1 = axes('position',[0.12,0.25,0.38,0.6]);
% imagesc(tmp(:,[1:length(typeplot)]));
imagesc(y_early(loc_el(sorted_ind),:),[-1,1]); hold on
plot(jj,ii,'.k');
% colormap('summer')
cm = ([[linspace(0,1,500)',linspace(0,1,500)',linspace(1,1,500)'];[linspace(1,1,500)',linspace(1,0,500)',linspace(1,0,500)']]);
colormap((cm))
set(gca,'ytick',[1:length(gn_list(loc_el))],'YTickLabel',gn_list(loc_el(sorted_ind)),'TickLabelInterpreter','none')
set(gca,'xtick',[1:length(typeplot)],'XTickLabel',clusteruni(typeplot),'TickLabelInterpreter','none')
title('PE-CE')
h2 = axes('position',[0.52,0.25,0.38,0.6]);
imagesc(y_late(loc_el(sorted_ind),:),[-1,1]); hold on
% imagesc(tmp(:,[length(typeplot)+1:end]));
% set(gca,'ytick',[1:length(gn_list(loc_el))],'YTickLabel',gn_list(loc_el),'TickLabelInterpreter','none')
set(gca,'ytick',[1:length(gn_list(loc_el))],'YTickLabel',[],'TickLabelInterpreter','none')
set(gca,'xtick',[1:length(typeplot)],'XTickLabel',clusteruni(typeplot),'TickLabelInterpreter','none')
title('PL-CT')
c = colorbar;
w = get(c,'Position');
w([1,3]) = [0.91,0.02];
set(c,'Position',w)
linkaxes([h1,h2],'y')
plot(jj,ii,'.k');
sgtitle('cytoskeleton','fontweight','bold')
% eval(['export_fig lymphoid_pro-inflammatory_cytokines_',date,'.pdf']);
%% 2 genes correlation scatter plot
class_type = 1; %1=immune, 2=stromal, 3=tb, 4= vascular
if class_type==1
    typeplot = immune_clusters;
elseif class_type==2
    typeplot = stromal_clusters;
elseif class_type==3
    typeplot = tb_clusters;
elseif class_type==4
    typeplot = vas_clusters;
elseif class_type==5
    typeplot = T_cells_tmp_uni;
end
typeplot = 33;

g1 = {'FLT1'};
g2 = {'PGF'};
% g1_exp = log2(1+data(strcmpi(geneid,g1),:));
% g2_exp = log2(1+data(strcmpi(geneid,g2),:));
g1_exp = data(strcmpi(geneid,g1),:);
g2_exp = data(strcmpi(geneid,g2),:);
for k = 1:length(typeplot)
    gr1 = find(T_cells_tmp==typeplot(k) & ctrlearly_flag == 1);
    gr2 = find(T_cells_tmp==typeplot(k) & peearly_flag == 1);
    gr3 = find(T_cells_tmp==typeplot(k) & ctrl_flag == 1);
    gr4 = find(T_cells_tmp==typeplot(k) & pelate_flag == 1);

   hf = figure('color','w','position',[100,20,1350,400],'name',clusteruni{typeplot(k)});
%     [ha, pos] = tight_subplot(1, 4, [0.08,0.08], [0.04,0.04], [0.05,0.08]);
%     axes(ha(1))
    h1 = axes('position',[0.075,0.1,0.2,0.7]);
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
    title([string(['CT, r=',num2str(r)]),['-log(pv) = ',num2str(-log10(pv))],string(['N = ',num2str(length(g1_exp(gr1)))])]);
    ylabel(upper(g2))
    xlabel(upper(g1))
    box off
%     axes(ha(2))
    h2 = axes('position',[0.3,0.1,0.2,0.7]);
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
    title([string(['CT, r=',num2str(r)]),['-log(pv) = ',num2str(-log10(pv))],string(['N = ',num2str(length(g1_exp(gr2)))])]);
    xlabel(upper(g1))
    box off
%     axes(ha(3))
    h3 = axes('position',[0.525,0.1,0.2,0.7]);
%     plot(g1_exp(gr3),g2_exp(gr3),'.')
%     b = polyfit(g1_exp(gr3),g2_exp(gr3),1);
%     hold on
%     mn = min(g1_exp(gr3));
%     mx = max(g1_exp(gr3));
%     plot([mn,mx],b(2)+b(1)*[mn,mx],'-k')
%     [r,pv] = corr(g1_exp(gr3)',g2_exp(gr3)');
%     title([string(['CT, r=',num2str(r)]),['-log(pv) = ',num2str(-log10(pv))],string(['N = ',num2str(length(g1_exp(gr3)))])]);
    g1_gr3 = log2(1+g1_exp(gr3)+rand(size(g1_exp(gr3)))*0.5);
    g2_gr3 = log2(1+g2_exp(gr3)+rand(size(g2_exp(gr3)))*0.5);
    plot(g1_gr3,g2_gr3,'.')
    b = polyfit(g1_gr3,g2_gr3,1);
    hold on
    mn = min(g1_gr3);
    mx = max(g1_gr3);
    plot([mn,mx],b(2)+b(1)*[mn,mx],'-k')
    [r,pv] = corr(g1_gr3',g2_gr3');
    title([string(['CT, r=',num2str(r)]),['-log(pv) = ',num2str(-log10(pv))],string(['N = ',num2str(length(g1_exp(gr3)))])]);
    xlabel(upper(g1))
    box off
%     axes(ha(4))
    h4 = axes('position',[0.75,0.1,0.2,0.7]);
%     plot(g1_exp(gr4),g2_exp(gr4),'.')
%     b = polyfit(g1_exp(gr4),g2_exp(gr4),1);
%     hold on
%     mn = min(g1_exp(gr4));
%     mx = max(g1_exp(gr4));
%     plot([mn,mx],b(2)+b(1)*[mn,mx],'-k')
%     [r,pv] = corr(g1_exp(gr4)',g2_exp(gr4)');
%     title([string(['PL, r=',num2str(r)]),['-log(pv) = ',num2str(-log10(pv))],string(['N = ',num2str(length(g1_exp(gr4)))])]);
%     sgtitle([regexprep(clusteruni{typeplot(k)},'_','-')])
    g1_gr4 = log2(1+g1_exp(gr4)+rand(size(g1_exp(gr4)))*0.5);
    g2_gr4 = log2(1+g2_exp(gr4)+rand(size(g2_exp(gr4)))*0.5);
    plot(g1_gr4,g2_gr4,'.')
    b = polyfit(g1_gr4,g2_gr4,1);
    hold on
    mn = min(g1_gr4);
    mx = max(g1_gr4);
    plot([mn,mx],b(2)+b(1)*[mn,mx],'-k')
    [r,pv] = corr(g1_gr4',g2_gr4');
    title([string(['CT, r=',num2str(r)]),['-log(pv) = ',num2str(-log10(pv))],string(['N = ',num2str(length(g1_exp(gr4)))])]);
    xlabel(upper(g1))
    box off
    linkaxes([h1,h2,h3,h4],'xy')
    eval(['export_fig angiogenesis_correlation_genes/',clusteruni{typeplot(k)},'_',char(g1),'-',char(g2),'_',date,'.pdf']);
    close(hf)
end
%% violin 
class_type = 1; %1=immune, 2=stromal, 3=tb, 4= vascular
if class_type==1
    typeplot = immune_clusters;
elseif class_type==2
    typeplot = stromal_clusters;
elseif class_type==3
    typeplot = tb_clusters;
elseif class_type==4
    typeplot = vas_clusters;
elseif class_type==5
    typeplot = T_cells_tmp_uni;
end
% typeplot = [2,8,9];
gn = upper('spp1');
g = find(strcmpi(gn,geneid));
fnorm = normpdf(norminv([0.001:0.01:0.999]));

figure('color','w','position',[20,20,1500,500],'Name',gn);
ha = axes('Position',[0.03,0.3,0.95,0.55]);
[~,cond] = max(condition_sorted(1:4,:));
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
    t_ed = zeros(4,1);%choose the desired column
    t_av = zeros(4,1);
    t_75 = zeros(4,1);
    for i=[1:4]
        gr2 = find(T_cells_tmp==typeplot(c) & condition_sorted(i,:)'==1);%find(fc_time_sorted==fc_time_uni(i));%
        y = (data(g,gr2));
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
           
            plot(k+i/5 + x-0.5, y,'.','color','none','marker',markervec(i)...
                ,'markerfacecolor',['#',color_celltypes_ordered{typeplot(k),2}],'markersize',3); hold on;
        else
            plot(-0.5+k+i/5+0.1*rand(length(gr2),1)-0.1,0.2*rand(length(gr2),1)-0.1+data(g,gr2)'....
                ,'.','color','none','marker',markervec(i),'markerfacecolor',['#',color_celltypes_ordered{typeplot(k),2}],'markersize',3); hold on;
        end
        xt = [xt,k+i/5-0.5];
        %         plot(i+0.5*rand(length(gr2),1)-0.1,0.5*rand(length(gr2),1)-0.1+data_orig_all_sorted(g,gr2)','.'); hold on;
%         t_ed(i) = median(y);
        t_av(i) = mean(data(g,gr2));
%         t_75(i) = prctile(y,75);
        t_98(i) = prctile(y,95);
    end
    p99 = max([p99,t_98]);
    plot(k+[1:4]/5-0.5,t_av,'sr','MarkerFaceColor','k'); hold on;
%     plot(t_av,'or'); hold on;
%     plot(t_75,'dg'); hold on;
    axis tight
%     yl = get(gca,'ylim');
%     set(gca,'xtick',[k+[1:4]/5],'XTickLabel',[{'ce','pe','ct','pl',}],'ylim',[0,p99])
    set(gca,'ylim',[0,p99])
    %     axis tight
%     title([gn,', ',regexprep(clusteruni{typeplot(c)},'_','-')])
    comp = [[1,2];[3,4]];
    pcomp = ones(length(comp(:,1)),1);
    for jj=1:length(comp)
        gr1 = find(T_cells_tmp==typeplot(c) & condition_sorted(comp(jj,1),:)'==1);
        gr2 = find(T_cells_tmp==typeplot(c) & condition_sorted(comp(jj,2),:)'==1);
        y1 = (data(g,gr1));
        y2 = (data(g,gr2));
        pcomp(jj) = ranksum(y1,y2);
        if pcomp(jj)<1e-4
            sigstar='***';
        elseif pcomp(jj)<1e-3
            sigstar='**';
        elseif pcomp(jj)<1e-2
            sigstar='*';
        else
            sigstar='';
        end        
        plot(k + comp(jj,:)/5-0.5,p99*(0.93+jj*0.01)*[1,1],'-','color',0.5*[1,1,1]);
        text(mean(k + comp(jj,:)/5-0.5),p99*(0.93+jj*0.02),sigstar,'fontsize',7,'VerticalAlignment','bottom','HorizontalAlignment','center');
    end        
end
set(gca,'ylim',[0,p99]);
yl = get(gca,'ylim');
text(xt,-(yl(2)-yl(1))*0.05*ones(size(xt)),repmat({'ce','pe','ct','pl'},1,length(typeplot)),'fontsize',5)
set(gca,'xtick',[1:length(typeplot)],'XTickLabel',clusteruni(typeplot),'TickLabelInterpreter'....
    ,'none','XTickLabelRotation',45,'ylim',[-(yl(2)-yl(1))*0.1,yl(2)]);
box off
% eval(['export_fig violin_',gn,'_',date,'.pdf']);

%% volcano-scatter plot
sex_genes = {'XIST','MTRNR2L8', 'EIF1AY', 'DDX3Y', 'RPS4Y1', 'KDM5D','MTRNR2L12','MTRNR2L10','MTRNR2L12'};
top_g = 50;
gr1name = 2;
gr2name = 3;
c1 = 1;
% gr1 = find(sum((T_cells_tmp ==[19, 21:29,46]),2));

gr1 = find(T_cells_tmp>=36 &  pelate_flag==1);%find(T_cells_tmp==24 );
gr2 = find(T_cells_tmp>=36  & ctrl_flag==1);%find(T_cells_tmp>=19 & T_cells_tmp<=29 & T_cells_tmp~=24);
% in = find(mean(data_orig_all_sorted(:,gr1)>0,2)>0.03 | mean(data_orig_all_sorted(:,gr2)>0,2)<0.7);
in = find(mean(data>0,2)>0.05 & mean(data>0,2)<0.5 & ~ismember(geneid(:,1),sex_genes));
% data=data';
ptt = zeros(length(in),1);
for s=1:length(in)
    if mod(s,100)==0
        s
    end
    ptt(s) = ranksum(data(in(s),gr1),data(in(s),gr2),'tail','both');
end
ptt(isnan(ptt)) = 1;
ptt(ptt<1) = qval_from_pval(ptt(ptt<1));
ptt(ptt<1e-300) = 1e-300;
x1 = mean(log2(data(:,gr1)+1),2);
x2 = mean(log2(data(:,gr2)+1),2);
d = x1-x2 ;
figure('position',[200,200,1400,580],'color','w');
[ha, pos] = tight_subplot(1, 3, [0.05,0.05], [0.1,0.05], [0.05,0.05]);
axes(ha(1))
plot(d(in),-log10(ptt),'.'); hold on;
[~,xi] = sort(ptt);
plot(d(in(xi(1:100))),-log10(ptt(xi(1:100))),'or');
text(d(in(xi(1:100))),-log10(ptt(xi(1:100))),geneid(in(xi(1:100))))
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
text(x1(xi(1:top_g)), x2(xi(1:top_g)),geneid(xi(1:top_g)),'fontsize',6);
xi = flipud(xi);
plot(x1(xi(1:top_g)), x2(xi(1:top_g)),'.r'); hold on;
text(x1(xi(1:top_g)), x2(xi(1:top_g)),geneid(xi(1:top_g)),'fontsize',6);
% title(['v3-nexgem (average)'])
xlabel(gr1name)
ylabel(gr2name)
axis tight

top_g = 50;
x1 = mean(data(:,gr1)>0,2);
x2 = mean(data(:,gr2)>0,2);
d = x1-x2 ;
[~,xi] = sort(d);
axes(ha(3))
plot(x1, x2, '.');hold on;
xmax = max(x1);
plot([0,xmax],[0,xmax],'-k'); grid on
plot([0,xmax],[0,xmax]+0.4,'--k'); grid on
plot([0.4,xmax],[0,xmax-0.4],'--k'); grid on
plot(x1(xi(1:top_g)), x2(xi(1:top_g)),'.r'); hold on;
text(x1(xi(1:top_g)), x2(xi(1:top_g)),geneid(xi(1:top_g)),'fontsize',6);
xi = flipud(xi);
plot(x1(xi(1:top_g)), x2(xi(1:top_g)),'.r'); hold on;
text(x1(xi(1:top_g)), x2(xi(1:top_g)),geneid(xi(1:top_g)),'fontsize',6);
% title(['v3-nexgem (precentage)'])
xlabel(gr1name)
ylabel(gr2name)
axis tight
% eval(['export_fig scatter_volcano_',gr1name,'_',gr2name,'_tropho_cluster',num2str(c1),'_',date,'.pdf']);
eval(['export_fig scatter_volcano_stromal_LPE_vs_LCT_',date,'.pdf']);


gr1 = find(T_cells_tmp>=36 &  peearly_flag==1);%find(T_cells_tmp==24 );
gr2 = find(T_cells_tmp>=36  & ctrlearly_flag==1);%find(T_cells_tmp>=19 & T_cells_tmp<=29 & T_cells_tmp~=24);
% in = find(mean(data_orig_all_sorted(:,gr1)>0,2)>0.03 | mean(data_orig_all_sorted(:,gr2)>0,2)<0.7);
in = find(mean(data>0,2)>0.05 & mean(data>0,2)<0.5 & ~ismember(geneid(:,1),sex_genes));
% data=data';
ptt = zeros(length(in),1);
for s=1:length(in)
    if mod(s,100)==0
        s
    end
    ptt(s) = ranksum(data(in(s),gr1),data(in(s),gr2),'tail','both');
end
ptt(isnan(ptt)) = 1;
ptt(ptt<1) = qval_from_pval(ptt(ptt<1));
ptt(ptt<1e-300) = 1e-300;
x1 = mean(log2(data(:,gr1)+1),2);
x2 = mean(log2(data(:,gr2)+1),2);
d = x1-x2 ;
figure('position',[200,200,1400,580],'color','w');
[ha, pos] = tight_subplot(1, 3, [0.05,0.05], [0.1,0.05], [0.05,0.05]);
axes(ha(1))
plot(d(in),-log10(ptt),'.'); hold on;
[~,xi] = sort(ptt);
plot(d(in(xi(1:100))),-log10(ptt(xi(1:100))),'or');
text(d(in(xi(1:100))),-log10(ptt(xi(1:100))),geneid(in(xi(1:100))))
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
text(x1(xi(1:top_g)), x2(xi(1:top_g)),geneid(xi(1:top_g)),'fontsize',6);
xi = flipud(xi);
plot(x1(xi(1:top_g)), x2(xi(1:top_g)),'.r'); hold on;
text(x1(xi(1:top_g)), x2(xi(1:top_g)),geneid(xi(1:top_g)),'fontsize',6);
% title(['v3-nexgem (average)'])
xlabel(gr1name)
ylabel(gr2name)
axis tight

top_g = 50;
x1 = mean(data(:,gr1)>0,2);
x2 = mean(data(:,gr2)>0,2);
d = x1-x2 ;
[~,xi] = sort(d);
axes(ha(3))
plot(x1, x2, '.');hold on;
xmax = max(x1);
plot([0,xmax],[0,xmax],'-k'); grid on
plot([0,xmax],[0,xmax]+0.4,'--k'); grid on
plot([0.4,xmax],[0,xmax-0.4],'--k'); grid on
plot(x1(xi(1:top_g)), x2(xi(1:top_g)),'.r'); hold on;
text(x1(xi(1:top_g)), x2(xi(1:top_g)),geneid(xi(1:top_g)),'fontsize',6);
xi = flipud(xi);
plot(x1(xi(1:top_g)), x2(xi(1:top_g)),'.r'); hold on;
text(x1(xi(1:top_g)), x2(xi(1:top_g)),geneid(xi(1:top_g)),'fontsize',6);
% title(['v3-nexgem (precentage)'])
xlabel(gr1name)
ylabel(gr2name)
axis tight
% eval(['export_fig scatter_volcano_',gr1name,'_',gr2name,'_tropho_cluster',num2str(c1),'_',date,'.pdf']);
eval(['export_fig scatter_volcano_stromal_EPE_vs_ECT_',date,'.pdf']);
%% stromal fraction enrichment

gr1 = find(T_cells_tmp>=36 );%stromal
gr2 = find(T_cells_tmp<36);%all rest

in = find(mean(data(:,[gr1(:);gr2(:)])>0,2)>0.05 & mean(data>0,2)<0.5 & ~ismember(geneid(:,1),sex_genes));
x1 = mean(log2(data(in,gr1)+1),2);
x2 = mean(log2(data(in,gr2)+1),2);
d1 = x1-x2;
[~,xi] = sort(d1,'descend');
table1 = [geneid(in(xi(1:200))),m2c([d1(xi(1:200)),x1(xi(1:200)),x2(xi(1:200))])];
saveCellFile(table1,['stromal_markers_mean_top200_',date,'.txt'])

x1 = mean(data(in,gr1)>0,2);
x2 = mean(data(in,gr2)>0,2);
d2 = x1./x2;
d2(isnan(d2)) = 0;
d2(x1<0.3) = 0;
[~,xi] = sort(d2,'descend');
table1 = [geneid(in(xi(1:200))),m2c([d2(xi(1:200)),x1(xi(1:200)),x2(xi(1:200))])];
saveCellFile(table1,['stromal_markers_frac_',date,'.txt'])

figure('position',[200,200,600,1000],'color','w');
[ha, pos] = tight_subplot(1, 2, [0.1,0.1], [0.1,0.05], [0.05,0.05]);
axes(ha(1))
barh([x1(xi(1:50)),x2(xi(1:50))])
xlabel('frac. cells >0')
set(ha(1),'ytick',[1:50],'yticklabel',geneid(in(xi(1:50))),'xdir','reverse','ydir','reverse','ylim',[0,51])
axes(ha(2))
barh([d2(xi(1:50))])
xlabel('fold enrichment (ST/rest)')
set(ha(2),'ytick',[1:50],'yticklabel',geneid(in(xi(1:50))),'ydir','reverse','ylim',[0,51])
eval(['export_fig stromal_vs_rest_frac_enrichment_',date,'.pdf']);

%% stromal bin plots 

gr1 = find(T_cells_tmp>=36 & peearly_flag==1);%EPE stromal
gr2 = find(T_cells_tmp>=36 & ctrlearly_flag==1);%ECT stromal
gr3 = find(T_cells_tmp>=36 & pelate_flag==1);%LPE stromal
gr4 = find(T_cells_tmp>=36 & ctrl_flag==1);%LCT stromal

% x1 = mean(data(:,gr1)>0,2);
% x2 = mean(data(:,gr2)>0,2);
% x3 = mean(data(:,gr3)>0,2);
% x4 = mean(data(:,gr4)>0,2);
% d1 = x1-x2;
% d2 = x3-x4;
% t = find( sqrt(d1.^2+d2.^2)>0.15);
% figure('position',[100,100,800,800],'color','w');
% plot(d1,d2,'.'); hold on;
% text(d1(t),d2(t),geneid(t))
% plot(0.5*[1,-1],0*[1,-1],'-k');
% plot(0*[1,-1],0.5*[1,-1],'-k');
% axis equal
% 
x1 = mean(log2(data(:,gr1)+1),2);
x2 = mean(log2(data(:,gr2)+1),2);
x3 = mean(log2(data(:,gr3)+1),2);
x4 = mean(log2(data(:,gr4)+1),2);
in = find(max([x1,x2,x3,x4],[],2)>0.3);
dx = x1-x2;
dy = x3-x4;
r = sqrt(dx(in).^2+dy(in).^2);
indsig = in(r>0.5);
% theta = atan2(d2(in),d1(in))*180/pi;
theta = atan2(dy(indsig),dx(indsig))*180/pi;

theta(theta<0) = 360+theta(theta<0);
anglebins = [[0:45:(360)]-22.5,360];
[tmp,~,bins] = histcounts(theta,anglebins);
tmp = [sum(tmp([1,9])),tmp(2:8)];


bin1 = indsig((theta<22.5 | theta>337.5));
bin5 = indsig((theta>157.5 & theta<202.5));
binrest = setdiff(indsig,[bin1;bin5]);

% t = find( sqrt(d1(in).^2+d2(in).^2)>0.5);
figure('position',[100,100,800,800],'color','w');
% plot(d1(in),d2(in),'.'); hold on;
plot(dx(bin1),dy(bin1),'o','markerfacecolor','r','markeredgecolor','none'); hold on;
plot(dx(bin5),dy(bin5),'o','markerfacecolor','b','markeredgecolor','none'); hold on;
plot(dx(binrest),dy(binrest),'.k'); hold on;
text(dx(indsig),dy(indsig),geneid(indsig));

plot(1.5*[1,-1],0*[1,-1],'-k');
plot(0*[1,-1],1.5*[1,-1],'-k');
plot([-2.5,2.5],[-2.5,2.5],'-k');
plot([-2.5,2.5],tan(22.5/180*pi)*[-2.5,2.5],'-r');
plot([-2.5,2.5],tan(67.5/180*pi)*[-2.5,2.5],'-r');
plot([-2.5,2.5],tan(-22.5/180*pi)*[-2.5,2.5],'-r');
plot([-2.5,2.5],tan(-67.5/180*pi)*[-2.5,2.5],'-r');
plot([-2.5,2.5],0*[-2.5,2.5],'-k');
plot(0*[-2.5,2.5],[-2.5,2.5],'-k');
for i=1:8
    alpha = 45*(i-1)/180*pi;
    text(0.3*cos(alpha),0.3*sin(alpha),num2str(tmp(i)));
end
axis equal
set(gca,'xlim',[-1.5,1.5],'ylim',[-1,1])
xlabel('FC early')
ylabel('FC late')
eval(['export_fig stromal_FCearly_FC_late_',date,'.pdf']);


% r = sqrt(d1(in).^2+d2(in).^2);
% theta = atan2(d2(in),d1(in))*180/pi;
% theta(theta<0) = 360+theta(theta<0);
indsig = in(r>0.3);
theta = atan2(dy(indsig),dx(indsig))*180/pi;
theta(theta<0) = 360+theta(theta<0);
bin1 = indsig((theta<22.5 | theta>337.5));
bin5 = indsig((theta>157.5 & theta<202.5));

pehigh = geneid(bin1);
pelow = geneid(bin5);
saveCellFile(pehigh,['stromal_high_PEearly_',date,'.txt']);
saveCellFile(pelow,['stromal_low_PEearly_',date,'.txt']);

%%
gr1 = find(T_cells_tmp>=38 & T_cells_tmp<=43 & peearly_flag==1);%find(T_cells_tmp==24 );
gr2 = find(T_cells_tmp>=38 & T_cells_tmp<=43 & ctrlearly_flag==1);%find(T_cells_tmp>=19 & T_cells_tmp<=29 & T_cells_tmp~=24);
gr3 = find(T_cells_tmp>=38 & T_cells_tmp<=43 & pelate_flag==1);%find(T_cells_tmp==24 );
gr4 = find(T_cells_tmp>=38 & T_cells_tmp<=43 & ctrl_flag==1);%find(T_cells_tmp>=19 & T_cells_tmp<=29 & T_cells_tmp~=24);

x1 = mean(data(:,gr1)>0,2);
x2 = mean(data(:,gr2)>0,2);
x3 = mean(data(:,gr3)>0,2);
x4 = mean(data(:,gr4)>0,2);
d1 = x1-x2;
d2 = x3-x4;
t = find( sqrt(d1.^2+d2.^2)>0.15);
figure('position',[100,100,800,800],'color','w');
plot(d1,d2,'.'); hold on;
text(d1(t),d2(t),geneid(t))
plot(0.5*[1,-1],0*[1,-1],'-k');
plot(0*[1,-1],0.5*[1,-1],'-k');
axis equal
% 
x1 = mean(log2(data(:,gr1)+1),2);
x2 = mean(log2(data(:,gr2)+1),2);
x3 = mean(log2(data(:,gr3)+1),2);
x4 = mean(log2(data(:,gr4)+1),2);
in = find(max([x1,x2,x3,x4],[],2)>0.5);
d1 = x1-x2;
d2 = x3-x4;
t = find( sqrt(d1(in).^2+d2(in).^2)>0.5);
figure('position',[100,100,800,800],'color','w');
plot(d1(in),d2(in),'.'); hold on;
plot(d1(in(t)),d2(in(t)),'or')
text(d1(in(t)),d2(in(t)),geneid(in(t)))
plot(1.5*[1,-1],0*[1,-1],'-k');
plot(0*[1,-1],1.5*[1,-1],'-k');
plot([-2.5,2.5],[-2.5,2.5],'-k');
plot([-2.5,2.5],tan(22.5/180*pi)*[-2.5,2.5],'-r');
plot([-2.5,2.5],tan(67.5/180*pi)*[-2.5,2.5],'-r');
plot([-2.5,2.5],tan(-22.5/180*pi)*[-2.5,2.5],'-r');
plot([-2.5,2.5],tan(-67.5/180*pi)*[-2.5,2.5],'-r');
plot([-2.5,2.5],0*[-2.5,2.5],'-k');
plot(0*[-2.5,2.5],[-2.5,2.5],'-k');
axis equal
set(gca,'xlim',[-1,1.5],'ylim',[-1,1])
xlabel('FC early')
ylabel('FC late')
eval(['export_fig endot_FCearly_FC_late_',date,'.pdf']);


r = sqrt(d1(in).^2+d2(in).^2);
theta = atan2(d2(in),d1(in))*180/pi;
theta(theta<0) = 360+theta(theta<0);

pehigh = geneid(in( (theta<22.5 | theta>(360-22.5)) & r>0.3));
pelow = geneid(in( (theta>157.5 & theta<(180+22.5)) & r>0.3));
saveCellFile(pehigh,['endot_high_PEearly_',date,'.txt']);
saveCellFile(pelow,['endot_low_PEearly_',date,'.txt']);
%
gr1 = find(T_cells_tmp>=44 & T_cells_tmp<=46 & peearly_flag==1);%find(T_cells_tmp==24 );
gr2 = find(T_cells_tmp>=44 & T_cells_tmp<=46 & ctrlearly_flag==1);%find(T_cells_tmp>=19 & T_cells_tmp<=29 & T_cells_tmp~=24);
gr3 = find(T_cells_tmp>=44 & T_cells_tmp<=46 & pelate_flag==1);%find(T_cells_tmp==24 );
gr4 = find(T_cells_tmp>=44 & T_cells_tmp<=46 & ctrl_flag==1);%find(T_cells_tmp>=19 & T_cells_tmp<=29 & T_cells_tmp~=24);

x1 = mean(data(:,gr1)>0,2);
x2 = mean(data(:,gr2)>0,2);
x3 = mean(data(:,gr3)>0,2);
x4 = mean(data(:,gr4)>0,2);
d1 = x1-x2;
d2 = x3-x4;
t = find( sqrt(d1.^2+d2.^2)>0.15);
figure('position',[100,100,800,800],'color','w');
plot(d1,d2,'.'); hold on;
text(d1(t),d2(t),geneid(t))
plot(0.5*[1,-1],0*[1,-1],'-k');
plot(0*[1,-1],0.5*[1,-1],'-k');
axis equal
% 
x1 = mean(log2(data(:,gr1)+1),2);
x2 = mean(log2(data(:,gr2)+1),2);
x3 = mean(log2(data(:,gr3)+1),2);
x4 = mean(log2(data(:,gr4)+1),2);
in = find(max([x1,x2,x3,x4],[],2)>0.5);
d1 = x1-x2;
d2 = x3-x4;
t = find( sqrt(d1(in).^2+d2(in).^2)>0.5);
figure('position',[100,100,800,800],'color','w');
plot(d1(in),d2(in),'.'); hold on;
plot(d1(in(t)),d2(in(t)),'or')
text(d1(in(t)),d2(in(t)),geneid(in(t)))
plot(1.5*[1,-1],0*[1,-1],'-k');
plot(0*[1,-1],1.5*[1,-1],'-k');
plot([-2.5,2.5],[-2.5,2.5],'-k');
plot([-2.5,2.5],tan(22.5/180*pi)*[-2.5,2.5],'-r');
plot([-2.5,2.5],tan(67.5/180*pi)*[-2.5,2.5],'-r');
plot([-2.5,2.5],tan(-22.5/180*pi)*[-2.5,2.5],'-r');
plot([-2.5,2.5],tan(-67.5/180*pi)*[-2.5,2.5],'-r');
plot([-2.5,2.5],0*[-2.5,2.5],'-k');
plot(0*[-2.5,2.5],[-2.5,2.5],'-k');
axis equal
set(gca,'xlim',[-1,1.5],'ylim',[-1,1])
xlabel('FC early')
ylabel('FC late')
eval(['export_fig VSM_FCearly_FC_late_',date,'.pdf']);


r = sqrt(d1(in).^2+d2(in).^2);
theta = atan2(d2(in),d1(in))*180/pi;
theta(theta<0) = 360+theta(theta<0);

pehigh = geneid(in( (theta<22.5 | theta>(360-22.5)) & r>0.3));
pelow = geneid(in( (theta>157.5 & theta<(180+22.5)) & r>0.3));
saveCellFile(pehigh,['VSM_high_PEearly_',date,'.txt']);
saveCellFile(pelow,['VSM_low_PEearly_',date,'.txt']);
% % % % % % % 
 % % % % % % % % % % % % % % % % % % % 

 %%


class_type = 2; %1=immune, 2=stromal, 3=tb, 4= vascular
if class_type==1
    typeplot = immune_clusters;
elseif class_type==2
    typeplot = stromal_clusters;
elseif class_type==3
    typeplot = tb_clusters;
elseif class_type==4
    typeplot = vas_clusters;
elseif class_type==5
    typeplot = T_cells_tmp_uni;
end
typeplot = typeplot([2,4]);
figure;
set(gcf,'color','w','position',[20,20,1400,1000]);
if class_type==5
    [ha, pos] = tight_subplot(6, ceil(length(typeplot)/6), [0.08,0.08], [0.04,0.04], [0.05,0.08]);
else
    [ha, pos] = tight_subplot(1, ceil(length(typeplot)/1), [0.08,0.08], [0.04,0.04], [0.05,0.08]);
end
sig_gene_cluster = cell(length(typeplot),1);
th = 0.7;
list = [];
geneper_bin = zeros(length(typeplot),8);
for j=1:length(typeplot)
    
        gr1 = find( (T_cells_tmp == typeplot(j)) & peearly_flag==1);
        gr2 = find( (T_cells_tmp == typeplot(j) ) & pelate_flag==1);
        gr3 = find( (T_cells_tmp == typeplot(j)) & ctrlearly_flag==1);
        gr4 = find( (T_cells_tmp == typeplot(j) ) & ctrl_flag==1);
        x1 = mean(log2(data(:,gr1)+1),2);
        x2 = mean(log2(data(:,gr2)+1),2);
        x3 = mean(log2(data(:,gr3)+1),2);
        x4 = mean(log2(data(:,gr4)+1),2);
        dx = x1-x3 ;
        dy = x2-x4 ;
        dis = sqrt(dx.^2 + dy.^2);
        
        indsig = find(dis>th);
        psig_x = ones(size(indsig));
        psig_y = ones(size(indsig));
        for jj=1:length(indsig)
            psig_x(jj) = ranksum(data(indsig(jj),gr1),data(indsig(jj),gr3));
            psig_y(jj) = ranksum(data(indsig(jj),gr2),data(indsig(jj),gr4));
        end
        q = qval_from_pval([psig_x;psig_y])';
        q = min(reshape(q,length(indsig),2),[],2);
        indsig = indsig(q<0.1);
        sig_gene_cluster{j} = geneid(indsig);%[dx(dis>th),dy(dis>th)];
        list = [list;sig_gene_cluster{j}];
        [~,xi] = sort(dis,'descend');
        theta = atan2(dy(indsig),dx(indsig))*180/pi;
        theta(theta<0) = 360+theta(theta<0);
        anglebins = [[0:45:(360)]-22.5,360];
        [tmp,~,bins] = histcounts(theta,anglebins);
        tmp = [sum(tmp([1,9])),tmp(2:8)];
        geneper_bin(j,:) = [tmp];
        bins(bins==9) = 1;
        tmp2 = cell(max(tmp),8);
        for jj=1:8
            tmp2(1:tmp(jj),jj) = sig_gene_cluster{j}(bins==jj);
        end
        axes(ha(j))
        %     polarhistogram(atan2(dy(dis>th),dx(dis>th)),100);
        plot(dx(indsig),dy(indsig),'.'); hold on;
        axis tight
        xl = get(ha(j),'xlim');
        plot([-2.5,2.5],[-2.5,2.5],'-k');
        plot([-2.5,2.5],tan(22.5/180*pi)*[-2.5,2.5],'-r');
        plot([-2.5,2.5],tan(67.5/180*pi)*[-2.5,2.5],'-r');
        plot([-2.5,2.5],tan(-22.5/180*pi)*[-2.5,2.5],'-r');
        plot([-2.5,2.5],tan(-67.5/180*pi)*[-2.5,2.5],'-r');
        plot([-2.5,2.5],0*[-2.5,2.5],'-k');
        plot(0*[-2.5,2.5],[-2.5,2.5],'-k');
        axis equal
        set(ha(j),'xlim',[-1.5,2],'ylim',[-1.5,2])
        %     text(dx(xi(1:ntop)),dy(xi(1:ntop)),geneid_all(xi(1:ntop)),'fontsize',8);
        text(dx(indsig),dy(indsig),geneid(indsig),'fontsize',8);
        tmp = m2c(geneper_bin(j,:));
        tmp = cellfun(@num2str,tmp,'UniformOutput',0);
        text(th*cos((anglebins(1:8)+22.5)*pi/180),th*sin((anglebins(1:8)+22.5)*pi/180),tmp,'fontsize',8)
        [r,pr] = corr(dx,dy)
        c = strfind(clusteruni{typeplot(j)},'_');
        tl = [regexprep(clusteruni{typeplot(j)}(c(1)+1:end),'_','-'),',N=',num2str(clustersize(typeplot(j)))]; 
        title(tl,'fontsize',8,'fontname','arial')
       
        %     xlabel('pe-ce');
        %     ylabel('pl-ct');
    
end


% eval(['export_fig EVT2_SCT_FC_plane_',date,'.pdf']);

%% stromal panel e, heatmap of selected genes
class_type = 2; %1=immune, 2=stromal, 3=tb, 4= vascular
if class_type==1
    typeplot = immune_clusters;
elseif class_type==2
    typeplot = stromal_clusters;
elseif class_type==3
    typeplot = tb_clusters;
elseif class_type==4
    typeplot = vas_clusters;
elseif class_type==5
    typeplot = T_cells_tmp_uni;
end
% typeplot = typeplot([1:9]); %for myeloid
gn_list = loadCellFile('stromal_genes_modules4.txt'); %pro
% gn_list = loadCellFile('anti_inflammatory_niv.csv'); %anti
% gn_list = {'IL1B','CCL2','CXCL8','TNF','SPP1','REL','NFKBIA','EREG'};
% gn_list = {'IL10','CD74','KLF2','IDO1'}
% gn_list = upper(gn_list);
[~,loc] = ismember(gn_list,geneid);
gn_list = gn_list(loc>0);
loc = loc(loc>0);
th = 0.3;
y = mean_mat(loc,[(typeplot(1)-1)*4+1:4*typeplot(end)]);
y_early = y(:,[2:4:end])-y(:,[1:4:end]);
y_early(y_early<th & y_early>-th) = 0;
loc_early =  find(sum(y_early,2)~=0);
y_late = y(:,[4:4:end])-y(:,[3:4:end]);
y_late(y_late<th & y_late>-th) = 0;
loc_late =  find(sum(y_late,2)~=0);
% [~,loc_el] = ismember(loc_late,loc_early);
% loc_el = loc_early(loc_el(loc_el>0));
loc_el = union(loc_early,loc_late);
% tmp = [y_early(loc_el),y_late(loc_el)];
% tmp = cent_norm(tmp);
% Y = pdist(y_early(loc_el,:),'correlation');
% Y = pdist(y_early,'correlation');
% z = linkage(Y,'ward');
% sorted_ind = optimalleaforder(z,Y);
[ii,jj] = find(sign(y_early(loc_el,:)) == sign(y_late(loc_el,:)) & sign(y_early(loc_el,:))~=0);
% [ii,jj] = find(sign(y_early) == sign(y_late) & sign(y_early)~=0);
figure('color','w','position',[20,20,800,600]);
h1 = axes('position',[0.12,0.25,0.38,0.6]);
% imagesc(tmp(:,[1:length(typeplot)]));
imagesc(y_early(loc_el,:),[-1,1]); hold on
% imagesc(y_early,[-1,1]); hold on
plot(jj,ii,'.k');
% colormap('summer')
cm = ([[linspace(0,1,500)',linspace(0,1,500)',linspace(1,1,500)'];[linspace(1,1,500)',linspace(1,0,500)',linspace(1,0,500)']]);
colormap((cm))
% set(gca,'ytick',[1:length(gn_list(loc_el))],'YTickLabel',gn_list(loc_el),'TickLabelInterpreter','none')
set(gca,'ytick',[1:length(gn_list)],'YTickLabel',gn_list,'TickLabelInterpreter','none')
set(gca,'xtick',[1:length(typeplot)],'XTickLabel',clusteruni(typeplot),'TickLabelInterpreter','none')
title('PE-CE')
h2 = axes('position',[0.52,0.25,0.38,0.6]);
imagesc(y_late(loc_el,:),[-1,1]); hold on
% imagesc(y_late,[-1,1]); hold on
% imagesc(tmp(:,[length(typeplot)+1:end]));
% set(gca,'ytick',[1:length(gn_list(loc_el))],'YTickLabel',gn_list(loc_el),'TickLabelInterpreter','none')
% set(gca,'ytick',[1:length(gn_list(loc_el))],'YTickLabel',[],'TickLabelInterpreter','none')
set(gca,'ytick',[1:length(gn_list)],'YTickLabel',[],'TickLabelInterpreter','none')
set(gca,'xtick',[1:length(typeplot)],'XTickLabel',clusteruni(typeplot),'TickLabelInterpreter','none')
title('PL-CT')
c = colorbar;
w = get(c,'Position');
w([1,3]) = [0.91,0.02];
set(c,'Position',w)
linkaxes([h1,h2],'y')
plot(jj,ii,'.k');
sgtitle('','fontweight','bold')
eval(['export_fig stromal_panel_e_gene_groups_noleaf_',date,'.pdf']);

 
%% plot stress genes in TB
list = loadCellFile('stress_cytokines_TB.txt');
gr1 = find(T_cells_tmp>=25 & T_cells_tmp<=26 & peearly_flag==1);%EPE EVT
gr2 = find(T_cells_tmp>=25 & T_cells_tmp<=26 & ctrlearly_flag==1);%ECT EVT
gr3 = find(T_cells_tmp>=25 & T_cells_tmp<=26 & pelate_flag==1);%LPE EVT
gr4 = find(T_cells_tmp>=25 & T_cells_tmp<=26 & ctrl_flag==1);%LCT EVT

x1 = mean(log2(data(:,gr1)+1),2);
x2 = mean(log2(data(:,gr2)+1),2);
x3 = mean(log2(data(:,gr3)+1),2);
x4 = mean(log2(data(:,gr4)+1),2);
in = find(max([x1,x2,x3,x4],[],2)>0.3);
dx = x1-x2;% FC EPE
dy = x3-x4;% FC LPE

[~,loc] = ismember(list,geneid);
loc = intersect(loc,in);
list = geneid(loc);
t = [dx(loc),dy(loc)];
figure('position',[100,100,800,200],'color','w');
bar(t,'barwidth',0.5);
set(gca,'xtick',[1:length(list)],'xticklabel',list)
legend({'EPE','LPE'})
eval(['export_fig EVT_stressgenes_bar_EPE_LPE_',date,'.pdf']);

% % % % % syn TB and co.

gr1 = find(T_cells_tmp>=20 & T_cells_tmp<=24 & peearly_flag==1);%EPE EVT
gr2 = find(T_cells_tmp>=20 & T_cells_tmp<=24 & ctrlearly_flag==1);%ECT EVT
gr3 = find(T_cells_tmp>=20 & T_cells_tmp<=24 & pelate_flag==1);%LPE EVT
gr4 = find(T_cells_tmp>=20 & T_cells_tmp<=24 & ctrl_flag==1);%LCT EVT

x1 = mean(log2(data(:,gr1)+1),2);
x2 = mean(log2(data(:,gr2)+1),2);
x3 = mean(log2(data(:,gr3)+1),2);
x4 = mean(log2(data(:,gr4)+1),2);
in = find(max([x1,x2,x3,x4],[],2)>0.3);
dx = x1-x2;% FC EPE
dy = x3-x4;% FC LPE


[~,loc] = ismember(list,geneid);
loc = intersect(loc,in);
list = geneid(loc);
t = [dx(loc),dy(loc)];
figure('position',[100,100,800,200],'color','w');
bar(t,'barwidth',0.5);
set(gca,'xtick',[1:length(list)],'xticklabel',list)
legend({'EPE','LPE'})
eval(['export_fig SCTVCT_stressgenes_bar_EPE_LPE_',date,'.pdf']);


%%
