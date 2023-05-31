tic
clear all
close all
addpath('/data/matlab_functions/');
savefig_flag = 1;

sam_cell = {'134_3','134_2','132_1','135_1','135_2'};
sam_folder = {'out10X134_3_230220','out10X134_2_230220','out10X132_1_230219','out10X135_1_230404','out10X135_2_230404'};
for i=1:length(sam_cell)
    fprintf(['loading sample 10X',sam_cell{i},'\n']);
    mtx_file = ['/bigdata/runs_from_Aug_2022/samples/10X',sam_cell{i},'/',sam_folder{i},'/outs/filtered_feature_bc_matrix/matrix.mtx'];
    bc_file = ['/bigdata/runs_from_Aug_2022/samples/10X',sam_cell{i},'/',sam_folder{i},'/outs/filtered_feature_bc_matrix/barcodes.tsv'];
    gene_file = ['/bigdata/runs_from_Aug_2022/samples/10X',sam_cell{i},'/',sam_folder{i},'/outs/filtered_feature_bc_matrix/features.tsv'];
    eval(['[data_',sam_cell{i},', geneid_',sam_cell{i},', barcodes_',sam_cell{i},'] = load10xMtxFile(mtx_file,bc_file,gene_file,800,1e5);']);
%     eval(['data_',sam_cell{i},' = data_',sam_cell{i},' - repmat(g99_',sam_cell{i},',1,length(barcodes_',sam_cell{i},'));'])
end
str = [];
for i=1:length(sam_cell)
    str = [str,'repmat({''',sam_cell{i},'''},length(barcodes_',sam_cell{i},'),1);'];
end
eval(['sample = [',str,'];'])
sample = regexprep(sample,'_','-');
str = [];
for i=1:length(sam_cell)
    str = [str,'data_',sam_cell{i},','];
end
eval(['data = [',str,'];'])
data(data<0) = 0;
str = [];
for i=1:length(sam_cell)
    str = [str,'barcodes_',sam_cell{i},';'];
end
eval(['cellid = [',str,'];'])
geneid = geneid_132_1(:,1);

[num,txt,raw] = xlsread('PE-nuclei_metadata_230530.xlsx');
num(isnan(num)) = 0;
flagnames = txt(1,2:end);
samname = regexprep(txt(2:end,1),'_','-');
v3_flag = false(size(cellid));
ctrl_flag = v3_flag;
peearly_flag = v3_flag;
pelate_flag = v3_flag;
female_flag = v3_flag;
iugr_flag = v3_flag;
ctrlearly_flag = v3_flag;
csection_flag = v3_flag;
vaginalbirth_flag = v3_flag;
induction_flag = v3_flag;
noinduction_flag = v3_flag;
magnesium_flag = v3_flag;
spinal_flag = v3_flag;
epidural_flag = v3_flag;
generalanesthesia_flag = v3_flag;
donor_id_flag = cell(size(cellid));
week_flag = zeros(size(cellid));
weight_flag = zeros(size(cellid));
weightprct_flag = zeros(size(cellid));
donorage_flag = zeros(size(cellid));
deliverydate_flag = cell(size(cellid));
week_diagnosis_flag = zeros(size(cellid));
weightprct_hadloack_flag = zeros(size(cellid));
for i=1:length(samname)
    ind = find(strcmpi(samname{i},sample));
    v3_flag(ind) = num(i,find(strcmpi('v3_flag',flagnames))-1);
    ctrl_flag(ind) = num(i,find(strcmpi('ctrl_flag',flagnames))-1);
    peearly_flag(ind) = num(i,find(strcmpi('peearly_flag',flagnames))-1);
    pelate_flag(ind) = num(i,find(strcmpi('pelate_flag',flagnames))-1);
    female_flag(ind) = num(i,find(strcmpi('female_flag',flagnames))-1);
    iugr_flag(ind) = num(i,find(strcmpi('iugr_flag',flagnames))-1);
    ctrlearly_flag(ind) = num(i,find(strcmpi('ctrlearly_flag',flagnames))-1);
    csection_flag(ind) = num(i,find(strcmpi('csection_flag',flagnames))-1);
    vaginalbirth_flag(ind) = num(i,find(strcmpi('vaginalbirth_flag',flagnames))-1);
    induction_flag(ind) = num(i,find(strcmpi('induction_flag',flagnames))-1);
    noinduction_flag(ind) = num(i,find(strcmpi('noinduction_flag',flagnames))-1);
    magnesium_flag(ind) = num(i,find(strcmpi('magnesium_flag',flagnames))-1);
    spinal_flag(ind) = num(i,find(strcmpi('spinal_flag',flagnames))-1);
    epidural_flag(ind) = num(i,find(strcmpi('epidural_flag',flagnames))-1);
    generalanesthesia_flag(ind) = num(i,find(strcmpi('generalanesthesia_flag',flagnames))-1);
    donor_id_flag(ind) = txt(i+1,1+find(strcmpi('Donor_ID',flagnames)));
    deliverydate_flag(ind) = txt(i+1,1+find(strcmpi('Data_of_Delivery',flagnames)));
    week_flag(ind) = num(i,find(strcmpi('gestational_week',flagnames))-1);
    week_diagnosis_flag(ind) = num(i,find(strcmpi('Gestational_age_diagnosis',flagnames))-1);
    weight_flag(ind) = num(i,find(strcmpi('weight_gr',flagnames))-1);
    weightprct_flag(ind) = num(i,find(strcmpi('Weight_precentage_dolberg',flagnames))-1);
    weightprct_hadloack_flag(ind) = num(i,find(strcmpi('Weight_percentile_hadlock',flagnames))-1);
    donorage_flag(ind) = num(i,find(strcmpi('Mother_age',flagnames))-1);
end


clear geneid_*str = [];
clear data_*
clear barcodes_*

%  load PE_afterloading_04-May-2021.mat

cellid = cellfun(@(x,y) [x,'_',y], cellid, sample,'UniformOutput',0);

tot_mol = sum(data);
tot_mol(tot_mol>5e4) = 5e4;
tot_genes = sum(data>0);
validcells = (tot_mol>800 & tot_mol<5e4 & tot_genes>400);
sum(validcells)
data = data(:,validcells);
cellid = cellid(validcells);
sample = sample(validcells);
tot_mol = tot_mol(validcells);
v3_flag = v3_flag(validcells);
ctrl_flag = ctrl_flag(validcells);
peearly_flag = peearly_flag(validcells);
pelate_flag = pelate_flag(validcells);
female_flag = female_flag(validcells);
iugr_flag = iugr_flag(validcells);
ctrlearly_flag = ctrlearly_flag(validcells);
csection_flag = csection_flag(validcells);
vaginalbirth_flag = vaginalbirth_flag(validcells);
induction_flag = induction_flag(validcells);
noinduction_flag = noinduction_flag(validcells);
magnesium_flag = magnesium_flag(validcells);
spinal_flag = spinal_flag(validcells);
epidural_flag = epidural_flag(validcells);
generalanesthesia_flag = generalanesthesia_flag(validcells);
week_flag = week_flag(validcells);
weight_flag = weight_flag(validcells);
weightprct_flag = weightprct_flag(validcells);
donorage_flag = donorage_flag(validcells);
donor_id_flag = donor_id_flag(validcells);
deliverydate_flag = deliverydate_flag(validcells);
weightprct_hadloack_flag = weightprct_hadloack_flag(validcells);
week_diagnosis_flag = week_diagnosis_flag(validcells);


sample_uni = {'134-3','134-2','132-1','135-1','135-2'};%
for i=1:length(sample_uni)
    fprintf(['valid cells in ',sample_uni{i},' = ',num2str(sum((strcmpi(sample,sample_uni{i})))),'\n']);
end


save(['PE-nuclei_QC800_400_afterloading_',date],'data','cellid','sample','sample_uni','geneid'...
    ,'ctrl_flag','v3_flag'...
    ,'peearly_flag','pelate_flag'...
    ,'female_flag','iugr_flag','ctrlearly_flag'...
    ,'csection_flag','vaginalbirth_flag'...
    ,'induction_flag','noinduction_flag'...
    ,'magnesium_flag','spinal_flag'...
    ,'epidural_flag','generalanesthesia_flag'...
    ,'week_flag','weight_flag'...
    ,'weightprct_flag','donorage_flag'...
    ,'donor_id_flag','deliverydate_flag','week_diagnosis_flag','weightprct_hadloack_flag','-v7.3');


toc
