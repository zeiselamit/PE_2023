tic
clear all
close all
addpath('/data/matlab_functions/')
% cd('crop_for_paper/TB/CYP19A1-PGF/')
% foldername = 'crop_for_paper/TB/CYP19A1-PGF/';
cd('crop_for_paper/TB/ERVFRD1-FLT1/')
foldername = 'crop_for_paper/TB/ERVFRD1-FLT1/';
% foldername = [];
f = dir([foldername,'*tif*']);
f = dir('*tif');
fsize = [f.bytes];
f = {f.name};
nbit =16;
sizethresh = 1e7;
highprc_th = 99.99; %right LUT
lowprc_th = 96; %left LUT
saveflag = 0;
sumIpercell_tot = cell(length(f),2);
for i=[1:length(f)]
    % for i = 1:4 %ePE
    %     for i = 5:8 %ect
    %     i = 5;
    fname = f{i}
    %     if fsize(i)<sizethresh
    if ~contains(fname,'DAPIonly')
        tmp = strsplit(fname,'_');
        tmp = tmp(1);
        tmp = strsplit(tmp{1},'-');
        c1name = tmp{1};
        c2name = tmp{2};
        c3name = tmp{3};

        geneC1 = imread([foldername,fname],'index',1);
        low_in = double(prctile(geneC1(:),30))/(2^nbit-1);
        high_in = double(prctile(geneC1(:),98))/(2^nbit-1);
        geneC1_adj = imadjust(geneC1,[low_in,high_in],[0,1]);

        bw = im2bw(geneC1,double(prctile(geneC1(geneC1(:)>0),50))/(2^nbit-1));%im2bw(dapi,graythresh(dapi));
        bw1 = bwareaopen(bw,100);
        %         bw1 = imfill(bw1,'holes');
        bw1 = bwareaopen(bw1,100);
        [bw2] = segment_dapi_by_gradient(geneC1,bw1,100,500,0);
        bw80 = bwareaopen(bw2,100);
        bw80 = bw80 - bwareaopen(bw80,2000);
        bw = im2bw(geneC1,double(prctile(geneC1(geneC1(:)>0),90))/(2^nbit-1));%im2bw(dapi,graythresh(dapi));
        bw1 = bwareaopen(bw,100);
        bw1 = imfill(bw1,'holes');
        bw1 = bwareaopen(bw1,100);
        [bw2] = segment_dapi_by_gradient(geneC1,bw1,100,500,0);
        bw95 = bwareaopen(bw2,100);
        B80 = bwboundaries(bw80);
        B95 = bwboundaries(bw95);
        %             figure;
        %             imshow(geneC1_adj); hold on;
        %             for jjj=1:length(B80)
        %                 plot(B80{jjj}(:,2),B80{jjj}(:,1),'g');
        %             end
        %             for jjj=1:length(B95)
        %                 plot(B95{jjj}(:,2),B95{jjj}(:,1),'m');
        %             end
        bw = logical(bw80 - bw95);
        bw = bw | bw95;
        bw = bwareaopen(bw,500);
        x = regionprops(bw,"EulerNumber",'PixelIdxList');
        for zz=1:length(x)
            if x(zz).EulerNumber<1
                bw(x(zz).PixelIdxList) = false;
            end
        end
        %             B95 = bwboundaries(bw);
        %             for jjj=1:length(B95)
        %                 plot(B95{jjj}(:,2),B95{jjj}(:,1),'y');
        %             end


        dapicells = regionprops(bw,'Centroid','Area','PixelIdxList');

        geneC2 = imread([foldername,fname],'index',2);
        low_in = double(prctile(geneC2(:),lowprc_th))/(2^nbit-1);
        high_in = double(prctile(geneC2(:),highprc_th))/(2^nbit-1);
        geneC2_adj = imadjust(geneC2,[low_in,high_in],[0,1]);
        C2_rgb = uint16(zeros([size(geneC1),3]));
        C2_rgb(:,:,1) =   + geneC1_adj/2;
        C2_rgb(:,:,2) = geneC2_adj  + geneC1_adj/2;
        C2_rgb(:,:,3) = geneC2_adj  + geneC1_adj/2;
        [bw2,pos2] = binary_dot_per_channel_x20_v2(geneC2_adj,99.999,10000);
        bw2_rgb = label2rgb(bw2 + 1,[0,0,0;0,1,1]);0
        geneC3 = imread([foldername,fname],'index',3);
        low_in = double(prctile(geneC3(:),lowprc_th))/(2^nbit-1);
        high_in = double(prctile(geneC3(:),highprc_th))/(2^nbit-1);
        geneC3_adj = imadjust(geneC3,[low_in,high_in],[0,1]);
        C3_rgb = uint16(zeros([size(geneC1),3]));
        C3_rgb(:,:,1) = geneC3_adj  + geneC1_adj/2;
        C3_rgb(:,:,2) =  + geneC1_adj/2;
        C3_rgb(:,:,3) = geneC3_adj + geneC1_adj/2;
        [bw3,pos3] = binary_dot_per_channel_x20_v2(geneC3_adj,99.999,5000);
        bw3_rgb = label2rgb(bw3 + 1,[0,0,0;1,0,1]);

        pos_tot = zeros(max([length(pos2),length(pos3)]),4);
        pos_tot(1:length(pos2),1) = pos2(:,1); pos_tot(1:length(pos2),2) = pos2(:,2); ...
            pos_tot(1:length(pos3),3) = pos3(:,1); pos_tot(1:length(pos3),4) = pos3(:,2); ...

        sumIpercell = zeros(length(dapicells),2);
        C2mean = zeros(size(geneC1));
        C3mean = zeros(size(geneC1));

        for zzz=1:length(dapicells)
            ind = dapicells(zzz).PixelIdxList;
            sumIpercell(zzz,:) = [mean(geneC2(ind)),mean(geneC3(ind))];
            C2mean(ind) = sumIpercell(zzz,1);
            C3mean(ind) = sumIpercell(zzz,2);
        end
        sumIpercell(:,3) = [dapicells.Area]';
        sumIpercell_tot{i,1} = fname;
        sumIpercell_tot{i,2} = (sumIpercell);
        B95 = bwboundaries(bw);
        figure('Position',[100,100,1200,600],'color','w');
        [ha, pos] = tight_subplot(1,2, [0.01,0.01], [0.05,0.15], [0.01,0.01]);
        axes(ha(1))
        imagesc(C2mean,[min(sumIpercell(:,1)),max(sumIpercell(:,1))]); hold on;
        for jjj=1:length(B95)
            plot(B95{jjj}(:,2),B95{jjj}(:,1),'w','LineWidth',0.1);
        end
        axis equal, axis off
        title(c2name)
        axes(ha(2))
        imagesc(C3mean,[min(sumIpercell(:,2)),max(sumIpercell(:,2))]);hold on;
        for jjj=1:length(B95)
            plot(B95{jjj}(:,2),B95{jjj}(:,1),'w','LineWidth',0.1);
        end
        axis equal, axis off
        title(c3name)



        I_rgb = uint16(zeros([size(geneC1),3]));
        I_rgb(:,:,1) = geneC3_adj + geneC1_adj/2;
        I_rgb(:,:,2) = geneC2_adj + geneC1_adj/2;
        I_rgb(:,:,3) = geneC2_adj + geneC3_adj + geneC1_adj/2;
        figure('color','w','position',[20,20,1200,1000],'name',fname);
        [ha, pos] = tight_subplot(2,3, [0.01,0.01], [0.01,0.01], [0.01,0.01]);
        %             ax = zeros(1,8);
        for jj = 1:2
            for j=1:3
                if jj == 1 & j>1
                    axes(ha(j));
                    %                             ax(j) = subplot(2,4,j);
                    eval(['imshow(C',num2str(j),'_rgb)'])
                    eval(['title(c',num2str(j),'name);'])
                elseif jj == 1 & j==1
                    axes(ha(j));
                    %                             ax(1) = subplot(2,4,1);
                    eval(['imshow(geneC1_adj)'])
                    eval(['title(c',num2str(j),'name);'])
                elseif jj == 2 & j>1
                    axes(ha(j+jj+1));
                    %                             ax(j+2*jj) = subplot(2,4,j+2*jj);
                    eval(['imshow(bw',num2str(j),'_rgb)'])
                    %                             linkaxes([ax(j),ax(j+2*jj)],'xy')
                elseif jj == 2 & j==1
                    %                            ax(5) = subplot(2,4,5);
                    axes(ha(4));
                    eval(['imshow(geneC1_adj)'])
                end
            end
        end
        linkaxes(ha,'xy');
        % set(ax(1),'position',[0,0.6,0.25,0.25]);set(ax(2),'position',[0.25,0.6,0.25,0.25]);set(ax(3),'position',[0.50,0.6,0.25,0.25]);set(ax(4),'position',[0.75,0.6,0.25,0.25]);
        % set(ax(5),'position',[0,0.15,0.25,0.25]);set(ax(6),'position',[0.25,0.15,0.25,0.25]);set(ax(7),'position',[0.50,0.15,0.25,0.25]);set(ax(8),'position',[0.75,0.15,0.25,0.25]);
        f2 = figure('color','k','position',[20,20,800,700],'name',fname);
        h1 = axes('Position',[0,0,1,1]);
        imshow(I_rgb);
        axes(h1)
        text(-1,50,c3name,'horizontalalignment','right','color',[1,0,1]); %red
        text(-1,150,c2name,'horizontalalignment','right','color',[0,1,1]); %green

        if saveflag==1
            imwrite(uint8(C2_rgb/(2^8)),[c2name,'_',fname(1:end-4),'.jpg'],'jpg','BitDepth',8,'Quality',100);
            imwrite(uint8(C3_rgb/(2^8)),[c3name,'_',fname(1:end-4),'.jpg'],'jpg','BitDepth',8,'Quality',100);
            imwrite(uint8(I_rgb/(2^8)),['all_',fname(1:end-4),'.jpg'],'jpg');
        end
    else
        c1name = 'DAPI';
        geneC1 = imread([fname],'index',1);
        low_in = double(prctile(geneC1(:),30))/(2^nbit-1);
        high_in = double(prctile(geneC1(:),98))/(2^nbit-1);
        geneC1_adj = imadjust(geneC1,[low_in,high_in],[0,1]);
        imwrite(uint8(geneC1_adj/(2^8)),[c1name,'_',fname(1:end-4),'.jpg'],'jpg','BitDepth',8,'Quality',100);
    end
    %         pause
    close all

    figure; imshow(I_rgb)
    % Get the current axes handle.
    button = 1;
    while sum(button) <=3   % read ginputs until a mouse right-button occurs
        [x,y,button] = ginput;
        %    plot(x,y,'r*'); hold on
    end
    sumIpercell_tot{i,3} = ([x,y]);
    pause
    % end
end

f3 = figure('color','k','position',[20,20,800,700],'name',fname);
h1 = axes('Position',[0,0,1,1]);
imshow(geneC1_adj/2, 'InitialMAg', 'fit');
%             imshow(bw4_rgb, 'InitialMAg', 'fit');
yel = cat(3,ones(size(geneC1_adj)),ones(size(geneC1_adj)),zeros(size(geneC1_adj)));
cyn = cat(3,zeros(size(geneC1_adj)),ones(size(geneC1_adj)),ones(size(geneC1_adj)));
mgn = cat(3,ones(size(geneC1_adj)),zeros(size(geneC1_adj)),ones(size(geneC1_adj)));
hold on
h = imshow(cyn);
hold off
dc2 = set(h, 'AlphaData', bw2);
hold on
h = imshow(mgn);
hold off
dc2c3 = set(h, 'AlphaData', bw3);
hold on
h = imshow(yel);
hold off

axes(h1)
text(-1,50,c3name,'horizontalalignment','right','color',[1,0,1]); %red
text(-1,150,c2name,'horizontalalignment','right','color',[0,1,1]); %green
if saveflag==1
    imwrite(uint8(C2_rgb/(2^8)),[c2name,'_',fname(1:end-4),'.jpg'],'jpg','BitDepth',8,'Quality',100);
    imwrite(uint8(C3_rgb/(2^8)),[c3name,'_',fname(1:end-4),'.jpg'],'jpg','BitDepth',8,'Quality',100);
    imwrite(uint8(I_rgb/(2^8)),['all_',fname(1:end-4),'.jpg'],'jpg');
end

f4 = figure('color','k','position',[20,20,800,700],'name',fname);
imshow(geneC1_adj/2, 'InitialMAg', 'fit'); hold on
plot(pos2(:,1),pos2(:,2),'.c'); hold on
plot(pos3(:,1),pos3(:,2),'.m'); hold on
%          text(-1,50,c3name,'horizontalalignment','right','color',[1,0,1]); %red
%          text(-1,150,c2name,'horizontalalignment','right','color',[0,1,1]); %green
%          text(-1,250,c4name,'horizontalalignment','right','color',[1,1,0]); %yellow
if saveflag==1
    eval(['export_fig ','dots_all_', fname(1:end-4),'.pdf']);
end
f5 = figure('color','k','position',[20,20,800,700],'name',fname);
imshow(geneC1_adj/2, 'InitialMAg', 'fit'); hold on
plot(pos2(:,1),pos2(:,2),'.c');
if saveflag==1
    eval(['export_fig ','dots_',c2name,'_',fname(1:end-4),'.pdf']);
end

f6 = figure('color','k','position',[20,20,800,700],'name',fname);
imshow(geneC1_adj/2, 'InitialMAg', 'fit'); hold on
plot(pos3(:,1),pos3(:,2),'.m');
if saveflag==1
    eval(['export_fig ','dots_',c3name,'_',fname(1:end-4),'.pdf']);
end



% for i=1:2:5
%     f5 = figure('color','k','position',[20,20,800,700],'name',fname);
%     imshow(geneC1_adj/2, 'InitialMAg', 'fit'); hold on
%     plot(pos_tot(:,i),pos_tot(:,i+1),'.c');
% end
%% intensity analysis hist
loc_ctrl = contains(sumIpercell_tot(:,1),'ctrl');
sumIpercell_cntl = sumIpercell_tot(loc_ctrl,:);
sumIpercell_PE = sumIpercell_tot(~loc_ctrl,:);
f8 = figure('color','w','position',[20,20,800,400],'name','Intensity');
h1 = axes('Position',[0.1,0.1,0.35,0.7]);
h2 = axes('Position',[0.47,0.1,0.05,0.7]);
h3 = axes('Position',[0.1,0.82,0.35,0.05]);
xmin = 650;
ymin = 900;
c1t = [];
c1f = [];
c2f = [];
for jj = 1:length(sumIpercell_cntl)
    cintensity = cell2mat(sumIpercell_cntl(jj,2));
    c1_int = cintensity(:,1);
    c1t = [c1t;c1_int];
    c2_int = cintensity(:,2);
    %         c3_int = cintensity(:,3);
    c1_pass = find(c1_int >= xmin);
    c2_pass = find(c2_int >= ymin);
    c12 = intersect(c1_pass,c2_pass);
    c1f = [c1f;c1_int(c12)];
    c2f = [c2f;c2_int(c12)];
    axes(h1);
    plot(c1_int(c12),c2_int(c12),'.b'); hold on
    c1_int(c12) = []; c2_int(c12) = [];
    plot(c1_int,c2_int,'.', 'Color', [0.5,0.5,0.5]); hold on
    if jj == length(sumIpercell_cntl)
        xline(xmin, '--'); yline(ymin, '--');
        xlabel(c2name); ylabel(c3name);
        [ksd1,x1] = ksdensity(c1f,'support',[min(c1f)-0.1,4500],'function','cdf');
        [ksd2,x2] = ksdensity(c2f,'support',[min(c2f)-0.1,4500],'function','cdf');
        c1fsrt = sort(c1f);
        c2fsrt = sort(c2f);
        %         [ksd3,x3] = ksdensity(c3_int);
        medc1 = median(c1f); c175 = prctile(c1f,75);
        medc2 = median(c2f); c275 = prctile(c2f,75);
        xline(medc1,'-',num2str(medc1));  xline(c175,'--r',num2str(c175));
        yline(medc2,'-',num2str(medc2));  yline(c275,'--r',num2str(c275));
        axes(h2);
        %         plot(ksd2,x2); hold on
        plot([1:length(c2fsrt)]/length(c2fsrt),c2fsrt); hold on
        xline(0,'-',{'0'}); xline(1,'-',{'1'}); xline(0.5);
        %     set(gca,'xtick',[],'ytick',[]); box off;
        box off; axis tight; axis off;
        axes(h3);
        %         plot(x1,ksd1); hold on
        plot(c1fsrt,[1:length(c1fsrt)]/length(c1fsrt));
        yline(0,'-',{'0'}); yline(1,'-',{'1'});  yline(0.5);
        box off; axis tight; axis off;
        %    set(gca,'xtick',[],'ytick',[]); box off;
    end

end
% linkaxes([h1,h2],'y'); linkaxes([h1,h3],'x');
% set(h1, 'xlim', [500,1000], 'ylim', [500,4500], 'xtick', [500:100:1000], 'ytick', [500:500:4500])
set(h1, 'xlim', [500,1000], 'ylim', [500,4500], 'xtick', [500:100:1000], 'ytick', [500:500:4500])
set(h2, 'ylim', [500,4500],'xlim',[0,1])
set(h3, 'xlim', [500,1000],'ylim',[0,1])

title(['CNTL'],['n =', num2str(length(c1f)),'/', num2str(length(c1t))]);

h4 = axes('Position',[0.57,0.1,0.35,0.7]);
h5 = axes('Position',[0.94,0.1,0.05,0.7]);
h6 = axes('Position',[0.57,0.82,0.35,0.05]);
c1t = [];
c1f = [];
c2f = [];
for jj = 1:length(sumIpercell_PE)
    cintensity = cell2mat(sumIpercell_PE(jj,2));
    c1_int = cintensity(:,1);
    c1t = [c1t;c1_int];
    c2_int = cintensity(:,2);
    %         c3_int = cintensity(:,3);
    c1_pass = find(c1_int >= xmin);
    c2_pass = find(c2_int >= ymin);
    c12 = intersect(c1_pass,c2_pass);
    c1f = [c1f;c1_int(c12)];
    c2f = [c2f;c2_int(c12)];
    axes(h4);
    plot(c1_int(c12),c2_int(c12),'.b'); hold on
    c1_int(c12) = []; c2_int(c12) = [];
    plot(c1_int,c2_int,'.', 'Color', [0.5,0.5,0.5]); hold on
    if jj == length(sumIpercell_PE)
        xline(xmin, '--'); yline(ymin, '--');
        [ksd1,x1] = ksdensity(c1f,'support',[min(c1f)-0.1,4500],'function','cdf');
        [ksd2,x2] = ksdensity(c2f,'support',[min(c2f)-0.1,4500],'function','cdf');
        c1fsrt = sort(c1f);
        c2fsrt = sort(c2f);
        medc1 = median(c1f); c175 = prctile(c1f,75);
        medc2 = median(c2f); c275 = prctile(c2f,75);
        %         [ksd3,x3] = ksdensity(c3_int);
        xline(medc1,'-',num2str(medc1));  xline(c175,'--r',num2str(c175));
        yline(medc2,'-',num2str(medc2));  yline(c275,'--r',num2str(c275));
        axes(h5);
        %         plot(ksd2,x2); hold on
        plot([1:length(c2fsrt)]/length(c2fsrt),c2fsrt); hold on
        xline(0,'-',{'0'}); xline(1,'-',{'1'}); xline(0.5);
        %     set(gca,'xtick',[],'ytick',[]); box off;
        box off; axis tight; axis off;
        axes(h6);
        %         plot(x1,ksd2); hold on
        plot(c1fsrt,[1:length(c1fsrt)]/length(c1fsrt));
        yline(0,'-',{'0'}); yline(1,'-',{'1'}); yline(0.5);
        box off; axis tight; axis off;
        %    set(gca,'xtick',[],'ytick',[]); box off;
    end
end
% linkaxes([h4,h5],'y'); linkaxes([h4,h6],'x'); linkaxes([h1,h4],'xy'); linkaxes([h1,h2,h4,h5],'y');
set(h4, 'xlim', [500,1000], 'ylim', [500,4500], 'xtick', [500:100:1000], 'ytick', [500:500:4500])
set(h5, 'ylim', [500,4500],'xlim',[0,1])
set(h6, 'xlim', [500,1000],'ylim',[0,1])
title(['PE'],['n =', num2str(length(c1f)),'/', num2str(length(c1t))]);

eval(['export_fig ',c2name,'_',c3name,'_','intensity_analysis_CT-PE_',date,'.pdf']);
%% intensity analysis by area no hist
loc_ctrl = contains(sumIpercell_tot(:,1),'ctrl');
sumIpercell_cntl = sumIpercell_tot(loc_ctrl,:);
sumIpercell_PE = sumIpercell_tot(~loc_ctrl,:);
f8 = figure('color','w','position',[20,20,800,400],'name','Intensity');
h1 = axes('Position',[0.1,0.18,0.35,0.7]);
% h2 = axes('Position',[0.47,0.1,0.05,0.7]);
% h3 = axes('Position',[0.1,0.82,0.35,0.05]);
xmin = 4000;
ymin = 0;
areamin = 700;
c1t = [];
c1f = [];
c2f = [];
for jj = 1:length(sumIpercell_cntl)
    cintensity = cell2mat(sumIpercell_cntl(jj,2));
    c1_int = cintensity(:,1);
    c1t = [c1t;c1_int];
    c2_int = cintensity(:,2);
    area = cintensity(:,3);
    %         c3_int = cintensity(:,3);
    c1_pass = find(c1_int >= xmin);
    c2_pass = find(c2_int >= ymin);
    area_pass = find(area >= areamin);
    c12 = intersect(c1_pass,c2_pass);
    c12area = intersect(c12,area_pass)
    c1f = [c1f;c1_int(c12area)];
    c2f = [c2f;c2_int(c12area)];
    axes(h1);
    plot(c1_int(c12area),c2_int(c12area),'.b'); hold on
    c1_int(c12area) = []; c2_int(c12area) = [];
    plot(c1_int,c2_int,'.', 'Color', [0.5,0.5,0.5]); hold on
    if jj == length(sumIpercell_cntl)
        xline(xmin, '--'); yline(ymin, '--');
        xlabel(c2name); ylabel(c3name);
        %         [ksd1,x1] = ksdensity(c1f,'support',[min(c1f)-0.1,11000],'function','cdf');
        %         [ksd2,x2] = ksdensity(c2f,'support',[min(c2f)-0.1,6500],'function','cdf');
        %         [ksd3,x3] = ksdensity(c3_int);
        medc1 = median(c1f); c175 = prctile(c1f,75);
        medc2 = median(c2f); c275 = prctile(c2f,75)
        xline(medc1,'-',num2str(medc1));  xline(c175,'--r',num2str(c175));
        yline(medc2,'-',num2str(medc2));  yline(c275,'--r',num2str(c275));

    end

end
% linkaxes([h1,h2],'y'); linkaxes([h1,h3],'x');
% set(h1, 'xlim', [500,1000], 'ylim', [500,4500], 'xtick', [500:100:1000], 'ytick', [500:500:4500])
set(h1, 'xlim', [500,13000], 'ylim', [500,10000], 'xtick', [500:1500:13000], 'ytick', [500:1500:10000])
% set(h2, 'ylim', [500,4500],'xlim',[0,1])
% set(h3, 'xlim', [500,1000],'ylim',[0,1])

title(['ECT'],['n =', num2str(length(c1f)),'/', num2str(length(c1t))]);

h2 = axes('Position',[0.57,0.18,0.35,0.7]);
% h5 = axes('Position',[0.94,0.1,0.05,0.7]);
% h6 = axes('Position',[0.57,0.82,0.35,0.05]);
c1t = [];
c1f = [];
c2f = [];
for jj = 1:length(sumIpercell_PE)
    cintensity = cell2mat(sumIpercell_PE(jj,2));
    c1_int = cintensity(:,1);
    c1t = [c1t;c1_int];
    c2_int = cintensity(:,2);
    %         c3_int = cintensity(:,3);
    c1_pass = find(c1_int >= xmin);
    c2_pass = find(c2_int >= ymin);
    area_pass = find(area >= areamin);
    c12 = intersect(c1_pass,c2_pass);
    c12area = intersect(c12,area_pass)
    c1f = [c1f;c1_int(c12area)];
    c2f = [c2f;c2_int(c12area)];
    axes(h2);
    plot(c1_int(c12area),c2_int(c12area),'.b'); hold on
    c1_int(c12area) = []; c2_int(c12area) = [];
    plot(c1_int,c2_int,'.', 'Color', [0.5,0.5,0.5]); hold on
    if jj == length(sumIpercell_PE)
        xline(xmin, '--'); yline(ymin, '--');
        xlabel(c2name);
        %         [ksd1,x1] = ksdensity(c1f,'support',[min(c1f)-0.1,4500],'function','cdf');
        %         [ksd2,x2] = ksdensity(c2f,'support',[min(c2f)-0.1,4500],'function','cdf');
        %         [ksd3,x3] = ksdensity(c3_int);
        medc1 = median(c1f); c175 = prctile(c1f,75);
        medc2 = median(c2f); c275 = prctile(c2f,75)
        xline(medc1,'-',num2str(medc1));  xline(c175,'--r',num2str(c175));
        yline(medc2,'-',num2str(medc2));  yline(c275,'--r',num2str(c275));
    end
end
% linkaxes([h4,h5],'y'); linkaxes([h4,h6],'x'); linkaxes([h1,h4],'xy'); linkaxes([h1,h2,h4,h5],'y');
set(h2, 'xlim', [500,13000], 'ylim', [500,10000], 'xtick', [500:1500:13000], 'ytick', [500:1500:10000])
% set(h5, 'ylim', [500,4500],'xlim',[0,1])
% set(h6, 'xlim', [500,1000],'ylim',[0,1])
title(['EPE'],['n =', num2str(length(c1f)),'/', num2str(length(c1t))]);

% eval(['export_fig ',c2name,'_',c3name,'_','intensity_analysis_CT-PE_no_hist',date,'.pdf']);
%% Area vs gene1

loc_ctrl = contains(sumIpercell_tot(:,1),'ctrl');
sumIpercell_cntl = sumIpercell_tot(loc_ctrl,:);
sumIpercell_PE = sumIpercell_tot(~loc_ctrl,:);
f9 = figure('color','w','position',[20,20,800,400],'name','Intensity');
h1 = axes('Position',[0.1,0.18,0.35,0.7]);
% h2 = axes('Position',[0.47,0.1,0.05,0.7]);
% h3 = axes('Position',[0.1,0.82,0.35,0.05]);
xmin = 0;
c1t = [];
c1f = [];
c3f = [];
for jj = 1:length(sumIpercell_cntl)
    cintensity = cell2mat(sumIpercell_cntl(jj,2));
    c1_int = cintensity(:,1);
    c1t = [c1t;c1_int];
    c3_int = cintensity(:,3);
    %         c3_int = cintensity(:,3);
    c1_pass = find(c1_int >= xmin);
    c1f = [c1f;c1_int(c1_pass)];
    c3f = [c3f;c3_int(c1_pass)];
    axes(h1);
    plot(c1_int(c1_pass),c3_int(c1_pass),'.b'); hold on
    c1_int(c1_pass) = []; c3_int(c1_pass) = [];
    plot(c1_int,c3_int,'.', 'Color', [0.5,0.5,0.5]); hold on
    if jj == length(sumIpercell_cntl)
        xline(xmin, '--'); yline(ymin, '--');
        xlabel(c2name); ylabel('Area');
        %         [ksd1,x1] = ksdensity(c1f,'support',[min(c1f)-0.1,11000],'function','cdf');
        %         [ksd2,x2] = ksdensity(c2f,'support',[min(c2f)-0.1,6500],'function','cdf');
        %         [ksd3,x3] = ksdensity(c3_int);
        medc1 = median(c1f); c175 = prctile(c1f,75);
        medc3 = median(c3f); c375 = prctile(c3f,75)
        %         xline(medc1,'-',num2str(medc1));  xline(c175,'--r',num2str(c175));
        %         yline(medc3,'-',num2str(medc3));  yline(c375,'--r',num2str(c375));

    end

end
% linkaxes([h1,h2],'y'); linkaxes([h1,h3],'x');
% set(h1, 'xlim', [500,1000], 'ylim', [500,4500], 'xtick', [500:100:1000], 'ytick', [500:500:4500])
set(h1, 'xlim', [500,13000], 'ylim', [500,10000], 'xtick', [500:1500:13000], 'ytick', [500:1500:10000])
% set(h2, 'ylim', [500,4500],'xlim',[0,1])
% set(h3, 'xlim', [500,1000],'ylim',[0,1])

title(['ECT'],['n =', num2str(length(c1f)),'/', num2str(length(c1t))]);

h2 = axes('Position',[0.57,0.18,0.35,0.7]);
% h5 = axes('Position',[0.94,0.1,0.05,0.7]);
% h6 = axes('Position',[0.57,0.82,0.35,0.05]);
c1t = [];
c1f = [];
c3f = [];
for jj = 1:length(sumIpercell_PE)
    cintensity = cell2mat(sumIpercell_PE(jj,2));
    c1_int = cintensity(:,1);
    c1t = [c1t;c1_int];
    c3_int = cintensity(:,3);
    %         c3_int = cintensity(:,3);
    c1_pass = find(c1_int >= xmin);
    c1f = [c1f;c1_int(c1_pass)];
    c3f = [c3f;c3_int(c1_pass)];
    axes(h2);
    plot(c1_int(c1_pass),c3_int(c1_pass),'.b'); hold on
    c1_int(c1_pass) = []; c3_int(c1_pass) = [];
    plot(c1_int,c3_int,'.', 'Color', [0.5,0.5,0.5]); hold on
    if jj == length(sumIpercell_cntl)
        xline(xmin, '--'); yline(ymin, '--');
        xlabel(c2name); ylabel('Area');
        %         [ksd1,x1] = ksdensity(c1f,'support',[min(c1f)-0.1,11000],'function','cdf');
        %         [ksd2,x2] = ksdensity(c2f,'support',[min(c2f)-0.1,6500],'function','cdf');
        %         [ksd3,x3] = ksdensity(c3_int);
        medc1 = median(c1f); c175 = prctile(c1f,75);
        medc3 = median(c3f); c375 = prctile(c3f,75)
        %         xline(medc1,'-',num2str(medc1));  xline(c175,'--r',num2str(c175));
        %         yline(medc3,'-',num2str(medc3));  yline(c375,'--r',num2str(c375));

    end
end
% linkaxes([h4,h5],'y'); linkaxes([h4,h6],'x'); linkaxes([h1,h4],'xy'); linkaxes([h1,h2,h4,h5],'y');
set(h2, 'xlim', [500,13000], 'ylim', [500,10000], 'xtick', [500:1500:13000], 'ytick', [500:1500:10000])
% set(h5, 'ylim', [500,4500],'xlim',[0,1])
% set(h6, 'xlim', [500,1000],'ylim',[0,1])
title(['EPE'],['n =', num2str(length(c1f)),'/', num2str(length(c1t))]);

% eval(['export_fig ',c2name,'_',c3name,'_','intensity_analysis_CT-PE_no_hist',date,'.pdf']);
%% polygon circumference
for i = 1:length(f)
    dx = diff(sumIpercell_tot{i,3}(:,1));
    dy = diff(sumIpercell_tot{i,3}(:,2));
    for j = 1:length(dx)
        pc(j) = sqrt(dx(j).^2+dy(j).^2);
    end
    sumIpercell_tot{i,4} = sum(pc);
    pc = [];
end
%% #pre-fused TB per villi circumference - boxplot
loc_ctrl = contains(sumIpercell_tot(:,1),'ctrl');
loc_210609 = contains(sumIpercell_tot(:,1),'210609'); loc_230111 = contains(sumIpercell_tot(:,1),'230111'); loc_210901 = contains(sumIpercell_tot(:,1),'210901'); %ce
loc_210511 = contains(sumIpercell_tot(:,1),'210511'); loc_210803 = contains(sumIpercell_tot(:,1),'210803'); loc_211222 = contains(sumIpercell_tot(:,1),'211222'); %pe
loc_ce_per_sample = [loc_210609,loc_230111,loc_210901];
loc_pe_per_sample = [loc_210511,loc_210803,loc_211222];
syn2_per_circ = cell2mat(sumIpercell_tot(:,5))./cell2mat(sumIpercell_tot(:,4));
syn2_per_circ_cntl = syn2_per_circ(loc_ctrl); syn2_per_circ_PE = syn2_per_circ(~loc_ctrl);

t_95 = prctile(syn2_per_circ,95);
p99 = max(t_95);

hf = figure('color','w','position',[100,20,400,400]);
boxplot([syn2_per_circ_cntl,syn2_per_circ_PE]); hold on

for j = 1:size(loc_ce_per_sample,2)
    plot(ones(7,1),syn2_per_circ(loc_ce_per_sample(:,j)),'.'); hold on
end
for j = 1:size(loc_pe_per_sample,2)
    plot(2*ones(7,1),syn2_per_circ(loc_pe_per_sample(:,j)),'.'); hold on
end
pcomp = ranksum(syn2_per_circ_cntl,syn2_per_circ_PE);
if pcomp<1e-4
    sigstar='***';
elseif pcomp<1e-3
    sigstar='**';
elseif pcomp<1e-2
    sigstar='*';
else
    sigstar='';
end
plot(2 + pcomp/3-0.5,p99*(0.93+1*0.01)*[1,1],'-','color',0.5*[1,1,1]);
text(mean(2 + pcomp/3-0.5),p99*(0.93+1*0.02),sigstar,'fontsize',7,'VerticalAlignment','bottom','HorizontalAlignment','center');
legend('ect210609','ect230111','ect210901','pe210511','pe210803','pe211222')
xticklabels(['ce';'pe']); ylabel('#pre-fused TB/ villi circumference')
eval(['export_fig ','prefusedTB_per_villi_circumference_ce_pe_box_',date,'.pdf']);
%% #pre-fused TB per villi circumference - bar
loc_ctrl = contains(sumIpercell_tot(:,1),'ctrl');
loc_210609 = contains(sumIpercell_tot(:,1),'210609'); loc_230111 = contains(sumIpercell_tot(:,1),'230111'); loc_210901 = contains(sumIpercell_tot(:,1),'210901'); %ce
loc_210511 = contains(sumIpercell_tot(:,1),'210511'); loc_210803 = contains(sumIpercell_tot(:,1),'210803'); loc_211222 = contains(sumIpercell_tot(:,1),'211222'); %pe
loc_ce_per_sample = [loc_210609,loc_230111,loc_210901];
loc_pe_per_sample = [loc_210511,loc_210803,loc_211222];
syn2_per_circ = cell2mat(sumIpercell_tot(:,5))./cell2mat(sumIpercell_tot(:,4));
syn2_per_circ_cntl = syn2_per_circ(loc_ctrl); syn2_per_circ_PE = syn2_per_circ(~loc_ctrl);
err_ect = std(syn2_per_circ_cntl)/sqrt(length(syn2_per_circ_cntl)); err_epe = std(syn2_per_circ_PE)/sqrt(length(syn2_per_circ_PE));
t_95 = prctile(syn2_per_circ,95);
p99 = max(t_95);

hf = figure('color','w','position',[100,20,400,400]);
bar([median(syn2_per_circ_cntl),median(syn2_per_circ_PE)]); hold on
er = errorbar([median(syn2_per_circ_cntl);median(syn2_per_circ_PE)],[err_ect;err_epe]); hold on
er.Color = [0 0 0];
er.LineStyle = 'none';
for j = 1:size(loc_ce_per_sample,2)
    plot(ones(7,1),syn2_per_circ(loc_ce_per_sample(:,j)),'.'); hold on
end
for j = 1:size(loc_pe_per_sample,2)
    plot(2*ones(7,1),syn2_per_circ(loc_pe_per_sample(:,j)),'.'); hold on
end
pcomp = ranksum(syn2_per_circ_cntl,syn2_per_circ_PE);
if pcomp<1e-4
    sigstar='***';
elseif pcomp<1e-3
    sigstar='**';
elseif pcomp<1e-2
    sigstar='*';
else
    sigstar='';
end
plot(2 + pcomp/3-0.5,p99*(0.93+1*0.01)*[1,1],'-','color',0.5*[1,1,1]);
text(mean(2 + pcomp/3-0.5),p99*(0.93+1*0.02),sigstar,'fontsize',7,'VerticalAlignment','bottom','HorizontalAlignment','center');
legend('','','ect210609','ect230111','ect210901','pe210511','pe210803','pe211222')
xticklabels(['ce';'pe']); ylabel('#pre-fused TB/ villi circumference')
% eval(['export_fig ','prefusedTB_per_villi_circumference_ce_pe_err_',date,'.pdf']);
