tic
clear all
close all
addpath('/data/matlab_functions/')
cd('crop_for_paper/TB/CYP19A1-FLT1-PGF/')
foldername = 'crop_for_paper/TB/CYP19A1-FLT1-PGF/';
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
% tmp2 = cell(1,length(f));
sumIpercell_tot = cell(length(f),2);
for i=[1:length(f)]
% i = 14;
    fname = f{i}



%     if fsize(i)<sizethresh
        if ~contains(fname,'DAPIonly')
            tmp = strsplit(fname,'_');
            tmp = tmp(1);
            tmp = strsplit(tmp{1},'-');
            c1name = tmp{1};
            c2name = tmp{2};
            c3name = tmp{3};
            c4name = tmp{4};
            geneC1 = imread([foldername,fname],'index',1);
            low_in = double(prctile(geneC1(:),30))/(2^nbit-1);
            high_in = double(prctile(geneC1(:),98))/(2^nbit-1);    
            geneC1_adj = imadjust(geneC1,[low_in,high_in],[0,1]);

            bw = im2bw(geneC1,double(prctile(geneC1(geneC1(:)>0),50))/(2^nbit-1));%im2bw(dapi,graythresh(dapi));
            bw1 = bwareaopen(bw,100);
            bw1 = imfill(bw1,'holes');
            bw1 = bwareaopen(bw1,100);
            [bw2] = segment_dapi_by_gradient(geneC1,bw1,100,500,0);
            bw80 = bwareaopen(bw2,100);
            bw80 = bw80 - bwareaopen(bw80,2000);
            bw = im2bw(geneC1,double(prctile(geneC1(geneC1(:)>0),85))/(2^nbit-1));%im2bw(dapi,graythresh(dapi));
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
            x_epe = regionprops(bw,"EulerNumber",'PixelIdxList');
            for zz=1:length(x_epe)
                if x_epe(zz).EulerNumber<1
                    bw(x_epe(zz).PixelIdxList) = false;
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
            [bw3,pos3] = binary_dot_per_channel_x20_v2(geneC3_adj,99.999,15000);
            bw3_rgb = label2rgb(bw3 + 1,[0,0,0;1,0,1]);
            geneC4 = imread([foldername,fname],'index',4);
%             low_in = double(prctile(geneC4(:),lowprc_th))/(2^nbit-1);
%             high_in = double(prctile(geneC4(:),highprc_th))/(2^nbit-1);
            low_in = double(prctile(geneC4(:),lowprc_th))/(2^nbit-1);
            high_in = double(prctile(geneC4(:),99.95))/(2^nbit-1);
            geneC4_adj = imadjust(geneC4,[low_in,high_in],[0,1]);
            geneC4_adj = medfilt2(geneC4_adj,[3,3]); %background reduce
            C4_rgb = uint16(zeros([size(geneC1),3]));
            C4_rgb(:,:,1) =  + geneC4_adj + geneC1_adj/2;
            C4_rgb(:,:,2) =  + geneC4_adj + geneC1_adj/2;
            C4_rgb(:,:,3) = + geneC1_adj/2;
%             [bw4,pos4] = binary_dot_per_channel_x20(medfilt2(geneC4,[3,3]),99.999,8000); %for high noise 
            [bw4,pos4] = binary_dot_per_channel_x20_v2(geneC4_adj,99.999,5000);
            bw4_rgb = label2rgb(bw4 + 1,[0,0,0;1,1,0]);
            pos_tot = zeros(max([length(pos2),length(pos3),length(pos4)]),6);
            pos_tot(1:length(pos2),1) = pos2(:,1); pos_tot(1:length(pos2),2) = pos2(:,2); ...
                pos_tot(1:length(pos3),3) = pos3(:,1); pos_tot(1:length(pos3),4) = pos3(:,2); ...
                pos_tot(1:length(pos4),5) = pos4(:,1); pos_tot(1:length(pos4),6) = pos4(:,2);
            sumIpercell = zeros(length(dapicells),3);
            C2mean = zeros(size(geneC1));
            C3mean = zeros(size(geneC1));
            C4mean = zeros(size(geneC1));
            for zzz=1:length(dapicells)
                ind = dapicells(zzz).PixelIdxList;
                sumIpercell(zzz,:) = [mean(geneC2(ind)),mean(geneC3(ind)),mean(geneC4(ind))];
                C2mean(ind) = sumIpercell(zzz,1);
                C3mean(ind) = sumIpercell(zzz,2);
                C4mean(ind) = sumIpercell(zzz,3);
            end
            
%             tmp2{1,i} = {fname, sumIpercell};
            sumIpercell(:,4) = [dapicells.Area]';
            sumIpercell_tot{i,1} = fname;
            sumIpercell_tot{i,2} = (sumIpercell);
            B95 = bwboundaries(bw);

            c1min = 1000;
            figure('Position',[100,100,1200,600],'color','w');
            [ha, pos] = tight_subplot(1,3, [0.01,0.01], [0.05,0.15], [0.01,0.01]);
            axes(ha(1))
            imagesc(C2mean,[min(sumIpercell(:,1)),max(sumIpercell(:,1))]); hold on;
            for jjj=1:length(B95)
                if sumIpercell_tot{i,2}(jjj,1)>c1min
                    plot(B95{jjj}(:,2),B95{jjj}(:,1),'g','LineWidth',0.1);
                else
                plot(B95{jjj}(:,2),B95{jjj}(:,1),'w','LineWidth',0.1);
                end
            end
            axis equal, axis off
            title(c2name)
            axes(ha(2))
            imagesc(C3mean,[min(sumIpercell(:,2)),max(sumIpercell(:,2))]);hold on;
            for jjj=1:length(B95)
                if sumIpercell_tot{i,2}(jjj,1)>c1min
                    plot(B95{jjj}(:,2),B95{jjj}(:,1),'g','LineWidth',0.1);
                else
                    plot(B95{jjj}(:,2),B95{jjj}(:,1),'w','LineWidth',0.1);
                end
            end
            axis equal, axis off
            title(c3name)
            axes(ha(3))
            imagesc(C4mean,[min(sumIpercell(:,3)),max(sumIpercell(:,3))]);hold on;
            for jjj=1:length(B95)
                if sumIpercell_tot{i,2}(jjj,1)>c1min
                    plot(B95{jjj}(:,2),B95{jjj}(:,1),'g','LineWidth',0.1);
                else
                    plot(B95{jjj}(:,2),B95{jjj}(:,1),'w','LineWidth',0.1);
                end
            end
            axis equal, axis off
            colormap('copper')
            title(c4name)
            if saveflag == 1
%             imwrite(uint8(geneC1_adj/(2^8)),['contour_by_CYP19A1_',fname(1:end-4),'.pdf'],'pdf','BitDepth',8,'Quality',100);
%             saveas(gcf,['contour_by_CYP19A1_',fname(1:end-4),'.pdf']);
                   eval(['export_fig ','contour_by_CYP19A1_',num2str(c1min),'_',fname(1:end-4),'.pdf']);
            end
            
            
            I_rgb = uint16(zeros([size(geneC1),3]));
            I_rgb(:,:,1) = geneC3_adj + geneC4_adj + geneC1_adj/2;
            I_rgb(:,:,2) = geneC2_adj + geneC4_adj + geneC1_adj/2;
            I_rgb(:,:,3) = geneC2_adj + geneC3_adj + geneC1_adj/2;
            figure('color','w','position',[20,20,1200,1000],'name',fname);
            [ha, pos] = tight_subplot(2,4, [0.01,0.01], [0.01,0.01], [0.01,0.01]);
%             ax = zeros(1,8);
            for jj = 1:2
                for j=1:4
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
                            axes(ha(j+2*jj));
%                             ax(j+2*jj) = subplot(2,4,j+2*jj);
                            eval(['imshow(bw',num2str(j),'_rgb)'])
%                             linkaxes([ax(j),ax(j+2*jj)],'xy')
                    elseif jj == 2 & j==1
%                            ax(5) = subplot(2,4,5);
                           axes(ha(5));
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
            text(-1,250,c4name,'horizontalalignment','right','color',[1,1,0]); %yellow
            if saveflag==1
                imwrite(uint8(C2_rgb/(2^8)),[c2name,'_',fname(1:end-4),'.jpg'],'jpg','BitDepth',8,'Quality',100);
                imwrite(uint8(C3_rgb/(2^8)),[c3name,'_',fname(1:end-4),'.jpg'],'jpg','BitDepth',8,'Quality',100);
                imwrite(uint8(C4_rgb/(2^8)),[c4name,'_',fname(1:end-4),'.jpg'],'jpg','BitDepth',8,'Quality',100);
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
        pause
close all
%     end
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
            dc2c3cd = set(h, 'AlphaData', bw4);
            axes(h1)
            text(-1,50,c3name,'horizontalalignment','right','color',[1,0,1]); %red
            text(-1,150,c2name,'horizontalalignment','right','color',[0,1,1]); %green
            text(-1,250,c4name,'horizontalalignment','right','color',[1,1,0]); %yellow
            if saveflag==1
                imwrite(uint8(C2_rgb/(2^8)),[c2name,'_',fname(1:end-4),'.jpg'],'jpg','BitDepth',8,'Quality',100);
                imwrite(uint8(C3_rgb/(2^8)),[c3name,'_',fname(1:end-4),'.jpg'],'jpg','BitDepth',8,'Quality',100);
                imwrite(uint8(C4_rgb/(2^8)),[c4name,'_',fname(1:end-4),'.jpg'],'jpg','BitDepth',8,'Quality',100);
                imwrite(uint8(I_rgb/(2^8)),['all_',fname(1:end-4),'.jpg'],'jpg');
            end
            
  f4 = figure('color','k','position',[20,20,800,700],'name',fname);
         imshow(geneC1_adj/2, 'InitialMAg', 'fit'); hold on
         plot(pos2(:,1),pos2(:,2),'.c'); hold on
         plot(pos3(:,1),pos3(:,2),'.m'); hold on
         plot(pos4(:,1),pos4(:,2),'.y'); hold off
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
    
    f7 = figure('color','k','position',[20,20,800,700],'name',fname);
    imshow(geneC1_adj/2, 'InitialMAg', 'fit'); hold on
    plot(pos4(:,1),pos4(:,2),'.y');
    if saveflag==1
       eval(['export_fig ','dots_',c4name,'_', fname(1:end-4),'.pdf']);
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
c1min = 650;
c2min = 900;
c3min = 500;
c1t = [];
c1f_12 = []; c1f_13 = []; c1f_23 = [];
c2f_12 = []; c2f_13 = []; c2f_23 = [];
for jj = 1:length(sumIpercell_cntl)
    cintensity = cell2mat(sumIpercell_cntl(jj,2));
    c1_int = cintensity(:,1);
    c1t = [c1t;c1_int];
    c2_int = cintensity(:,2);
    c3_int = cintensity(:,3);
    c1_pass = find(c1_int >= c1min);
    c2_pass = find(c2_int >= c2min);
    c3_pass = find(c3_int >= c3min);
    c12 = intersect(c1_pass,c2_pass);
    c13 = intersect(c1_pass,c3_pass);
    c23 = intersect(c2_pass,c3_pass);
    c1f_12 = [c1f_12;c1_int(c12)];
    c1f_13 = [c1f_13;c1_int(c13)];
    c1f_23 = [c1f_23;c1_int(c23)];
    c2f_12 = [c2f_12;c2_int(c12)];
    c2f_13 = [c2f_13;c2_int(c13)];
    c2f_23 = [c2f_23;c2_int(c23)];
    axes(h1);
    plot(c1_int(c12),c2_int(c12),'.b'); hold on
    c1_int(c12) = []; c2_int(c12) = [];
    plot(c1_int,c2_int,'.', 'Color', [0.5,0.5,0.5]); hold on
    if jj == length(sumIpercell_cntl)
        xline(c1min, '--'); yline(c2min, '--');
        xlabel(c2name); ylabel(c3name);
        [ksd1,x1] = ksdensity(c1f_12,'support',[min(c1f_12)-0.1,4500],'function','cdf');
        [ksd2,x2] = ksdensity(c2f_12,'support',[min(c2f_12)-0.1,4500],'function','cdf');
        c1fsrt = sort(c1f_12);
        c2fsrt = sort(c2f_12);
        %         [ksd3,x3] = ksdensity(c3_int);
        medc1 = median(c1f_12); c175 = prctile(c1f_12,75);
        medc2 = median(c2f_12); c275 = prctile(c2f_12,75);
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

title(['CNTL'],['n =', num2str(length(c1f_12)),'/', num2str(length(c1t))]);

h4 = axes('Position',[0.57,0.1,0.35,0.7]);
h5 = axes('Position',[0.94,0.1,0.05,0.7]);
h6 = axes('Position',[0.57,0.82,0.35,0.05]);
c1t = [];
c1f_12 = [];
c2f_12 = [];
for jj = 1:length(sumIpercell_PE)
    cintensity = cell2mat(sumIpercell_PE(jj,2));
    c1_int = cintensity(:,1);
    c1t = [c1t;c1_int];
    c2_int = cintensity(:,2);
    %         c3_int = cintensity(:,3);
    c1_pass = find(c1_int >= c1min);
    c2_pass = find(c2_int >= c2min);
    c12 = intersect(c1_pass,c2_pass);
    c1f_12 = [c1f_12;c1_int(c12)];
    c2f_12 = [c2f_12;c2_int(c12)];
    axes(h4);
    plot(c1_int(c12),c2_int(c12),'.b'); hold on
    c1_int(c12) = []; c2_int(c12) = [];
    plot(c1_int,c2_int,'.', 'Color', [0.5,0.5,0.5]); hold on
    if jj == length(sumIpercell_PE)
        xline(c1min, '--'); yline(c2min, '--');
        [ksd1,x1] = ksdensity(c1f_12,'support',[min(c1f_12)-0.1,4500],'function','cdf');
        [ksd2,x2] = ksdensity(c2f_12,'support',[min(c2f_12)-0.1,4500],'function','cdf');
        c1fsrt = sort(c1f_12);
        c2fsrt = sort(c2f_12); 
        medc1 = median(c1f_12); c175 = prctile(c1f_12,75);
        medc2 = median(c2f_12); c275 = prctile(c2f_12,75);
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
title(['PE'],['n =', num2str(length(c1f_12)),'/', num2str(length(c1t))]);

% eval(['export_fig ',c2name,'_',c3name,'_','intensity_analysis_CT-PE_',date,'.pdf']);
%% intensity analysis by area no hist gene 1 vs 2
loc_ctrl = contains(sumIpercell_tot(:,1),'ctrl');
sumIpercell_cntl = sumIpercell_tot(loc_ctrl,:);
sumIpercell_PE = sumIpercell_tot(~loc_ctrl,:);
f8 = figure('color','w','position',[20,20,800,400],'name','Intensity');
h1 = axes('Position',[0.1,0.18,0.35,0.7]);
% h2 = axes('Position',[0.47,0.1,0.05,0.7]);
% h3 = axes('Position',[0.1,0.82,0.35,0.05]);
c1min = 0;
c2min = 0;
areamin = 700;
c1t = [];
c1f_12 = [];
c2f_12 = [];
for jj = 1:length(sumIpercell_cntl)
    cintensity = cell2mat(sumIpercell_cntl(jj,2));
    c1_int = cintensity(:,1);
    c1t = [c1t;c1_int];
    c2_int = cintensity(:,2);
    area = cintensity(:,4);
    %         c3_int = cintensity(:,3);
    c1_pass = find(c1_int >= c1min);
    c2_pass = find(c2_int >= c2min);
    area_pass = find(area >= areamin);
    c12 = intersect(c1_pass,c2_pass);
    c12area = intersect(c12,area_pass);
    c1f_12 = [c1f_12;c1_int(c12area)];
    c2f_12 = [c2f_12;c2_int(c12area)];
    axes(h1);
    plot(c1_int(c12area),c2_int(c12area),'.b'); hold on
    c1_int(c12area) = []; c2_int(c12area) = [];
    plot(c1_int,c2_int,'.', 'Color', [0.5,0.5,0.5]); hold on
    if jj == length(sumIpercell_cntl)
        xline(c1min, '--'); yline(c2min, '--');
        xlabel(c2name); ylabel(c3name);
%         [ksd1,x1] = ksdensity(c1f,'support',[min(c1f)-0.1,11000],'function','cdf');
%         [ksd2,x2] = ksdensity(c2f,'support',[min(c2f)-0.1,6500],'function','cdf');
        %         [ksd3,x3] = ksdensity(c3_int);
        medc1 = median(c1f_12); c175 = prctile(c1f_12,75);
        medc2 = median(c2f_12); c275 = prctile(c2f_12,75);
        xline(medc1,'-',num2str(medc1));  xline(c175,'--r',num2str(c175));
        yline(medc2,'-',num2str(medc2));  yline(c275,'--r',num2str(c275));
       
    end

end
% linkaxes([h1,h2],'y'); linkaxes([h1,h3],'x');
% set(h1, 'xlim', [500,1000], 'ylim', [500,4500], 'xtick', [500:100:1000], 'ytick', [500:500:4500])
set(h1, 'xlim', [500,4500], 'ylim', [500,50000], 'xtick', [500:1500:4500], 'ytick', [500:10000:50000])
% set(h2, 'ylim', [500,4500],'xlim',[0,1])
% set(h3, 'xlim', [500,1000],'ylim',[0,1])

title(['ECT'],['n =', num2str(length(c1f_12)),'/', num2str(length(c1t))]);

h2 = axes('Position',[0.57,0.18,0.35,0.7]);
% h5 = axes('Position',[0.94,0.1,0.05,0.7]);
% h6 = axes('Position',[0.57,0.82,0.35,0.05]);
c1t = [];
c1f_12 = [];
c2f_12 = [];
for jj = 1:length(sumIpercell_PE)
    cintensity = cell2mat(sumIpercell_PE(jj,2));
    c1_int = cintensity(:,1);
    c1t = [c1t;c1_int];
    c2_int = cintensity(:,2);
    %         c3_int = cintensity(:,3);
    c1_pass = find(c1_int >= c1min);
    c2_pass = find(c2_int >= c2min);
    area_pass = find(area >= areamin);
    c12 = intersect(c1_pass,c2_pass);
    c12area = intersect(c12,area_pass);
    c1f_12 = [c1f_12;c1_int(c12area)];
    c2f_12 = [c2f_12;c2_int(c12area)];
    axes(h2);
    plot(c1_int(c12area),c2_int(c12area),'.b'); hold on
    c1_int(c12area) = []; c2_int(c12area) = [];
    plot(c1_int,c2_int,'.', 'Color', [0.5,0.5,0.5]); hold on
    if jj == length(sumIpercell_PE)
        xline(c1min, '--'); yline(c2min, '--');
        xlabel(c2name);
%         [ksd1,x1] = ksdensity(c1f,'support',[min(c1f)-0.1,4500],'function','cdf');
%         [ksd2,x2] = ksdensity(c2f,'support',[min(c2f)-0.1,4500],'function','cdf');
    %         [ksd3,x3] = ksdensity(c3_int);
        medc1 = median(c1f_12); c175 = prctile(c1f_12,75);
        medc2 = median(c2f_12); c275 = prctile(c2f_12,75);
        xline(medc1,'-',num2str(medc1));  xline(c175,'--r',num2str(c175));
        yline(medc2,'-',num2str(medc2));  yline(c275,'--r',num2str(c275));
    end
end
% linkaxes([h4,h5],'y'); linkaxes([h4,h6],'x'); linkaxes([h1,h4],'xy'); linkaxes([h1,h2,h4,h5],'y');
set(h2, 'xlim', [500,4500], 'ylim', [500,50000], 'xtick', [500:1500:4500], 'ytick', [500:10000:50000])
% set(h5, 'ylim', [500,4500],'xlim',[0,1])
% set(h6, 'xlim', [500,1000],'ylim',[0,1])
title(['EPE'],['n =', num2str(length(c1f_12)),'/', num2str(length(c1t))]);

% eval(['export_fig ',c2name,'_',c3name,'_','intensity_analysis_CT-PE_no_hist',date,'.pdf']);
%% intensity analysis by area no hist gene 1 vs 3
loc_ctrl = contains(sumIpercell_tot(:,1),'ctrl');
sumIpercell_cntl = sumIpercell_tot(loc_ctrl,:);
sumIpercell_PE = sumIpercell_tot(~loc_ctrl,:);
f8 = figure('color','w','position',[20,20,800,400],'name','Intensity');
h1 = axes('Position',[0.1,0.18,0.35,0.7]);
% h2 = axes('Position',[0.47,0.1,0.05,0.7]);
% h3 = axes('Position',[0.1,0.82,0.35,0.05]);
c1min = 0;
c3min = 0;
areamin = 700;
c1t = [];
c1f_13 = [];
c3f_13 = [];
for jj = 1:length(sumIpercell_cntl)
    cintensity = cell2mat(sumIpercell_cntl(jj,2));
    c1_int = cintensity(:,1);
    c1t = [c1t;c1_int];
    c3_int = cintensity(:,3);
    area = cintensity(:,4);
    %         c3_int = cintensity(:,3);
    c1_pass = find(c1_int >= c1min);
    c3_pass = find(c3_int >= c3min);
    area_pass = find(area >= areamin);
    c13 = intersect(c1_pass,c3_pass);
    c13area = intersect(c13,area_pass);
    c1f_13 = [c1f_13;c1_int(c13area)];
    c3f_13 = [c3f_13;c3_int(c13area)];
    axes(h1);
    plot(c1_int(c13area),c3_int(c13area),'.b'); hold on
    c1_int(c13area) = []; c3_int(c13area) = [];
    plot(c1_int,c3_int,'.', 'Color', [0.5,0.5,0.5]); hold on
    if jj == length(sumIpercell_cntl)
        xline(c1min, '--'); yline(c3min, '--');
        xlabel(c2name); ylabel(c4name);
%         [ksd1,x1] = ksdensity(c1f,'support',[min(c1f)-0.1,11000],'function','cdf');
%         [ksd2,x2] = ksdensity(c2f,'support',[min(c2f)-0.1,6500],'function','cdf');
        %         [ksd3,x3] = ksdensity(c3_int);
        medc1 = median(c1f_13); c175 = prctile(c1f_13,75);
        medc3 = median(c3f_13); c375 = prctile(c3f_13,75);
        xline(medc1,'-',num2str(medc1));  xline(c175,'--r',num2str(c175));
        yline(medc3,'-',num2str(medc3));  yline(c375,'--r',num2str(c375));
       
    end

end
% linkaxes([h1,h2],'y'); linkaxes([h1,h3],'x');
% set(h1, 'xlim', [500,1000], 'ylim', [500,4500], 'xtick', [500:100:1000], 'ytick', [500:500:4500])
set(h1, 'xlim', [500,4500], 'ylim', [300,6500], 'xtick', [500:1500:4500], 'ytick', [300:1500:6500])
% set(h2, 'ylim', [500,4500],'xlim',[0,1])
% set(h3, 'xlim', [500,1000],'ylim',[0,1])

title(['ECT'],['n =', num2str(length(c1f_13)),'/', num2str(length(c1t))]);

h2 = axes('Position',[0.57,0.18,0.35,0.7]);
% h5 = axes('Position',[0.94,0.1,0.05,0.7]);
% h6 = axes('Position',[0.57,0.82,0.35,0.05]);
c1t = [];
c1f_13 = [];
c3f_13 = [];
for jj = 1:length(sumIpercell_PE)
    cintensity = cell2mat(sumIpercell_PE(jj,2));
    c1_int = cintensity(:,1);
    c1t = [c1t;c1_int];
    c3_int = cintensity(:,3);
    %         c3_int = cintensity(:,3);
    c1_pass = find(c1_int >= c1min);
    c3_pass = find(c3_int >= c3min);
    area_pass = find(area >= areamin);
    c13 = intersect(c1_pass,c3_pass);
    c13area = intersect(c13,area_pass);
    c1f_13 = [c1f_13;c1_int(c13area)];
    c3f_13 = [c3f_13;c3_int(c13area)];
    axes(h2);
    plot(c1_int(c13area),c3_int(c13area),'.b'); hold on
    c1_int(c13area) = []; c3_int(c13area) = [];
    plot(c1_int,c3_int,'.', 'Color', [0.5,0.5,0.5]); hold on
    if jj == length(sumIpercell_PE)
        xline(c1min, '--'); yline(c3min, '--');
        xlabel(c2name);
%         [ksd1,x1] = ksdensity(c1f,'support',[min(c1f)-0.1,4500],'function','cdf');
%         [ksd2,x2] = ksdensity(c2f,'support',[min(c2f)-0.1,4500],'function','cdf');
    %         [ksd3,x3] = ksdensity(c3_int);
        medc1 = median(c1f_13); c175 = prctile(c1f_13,75);
        medc3 = median(c3f_13); c375 = prctile(c3f_13,75);
        xline(medc1,'-',num2str(medc1));  xline(c175,'--r',num2str(c175));
        yline(medc3,'-',num2str(medc3));  yline(c375,'--r',num2str(c375));
    end
end
% linkaxes([h4,h5],'y'); linkaxes([h4,h6],'x'); linkaxes([h1,h4],'xy'); linkaxes([h1,h2,h4,h5],'y');
set(h2, 'xlim', [500,4500], 'ylim', [300,6500], 'xtick', [500:1500:4500], 'ytick', [300:1500:6500])
% set(h5, 'ylim', [500,4500],'xlim',[0,1])
% set(h6, 'xlim', [500,1000],'ylim',[0,1])
title(['EPE'],['n =', num2str(length(c1f_13)),'/', num2str(length(c1t))]);

% eval(['export_fig ',c2name,'_',c3name,'_','intensity_analysis_CT-PE_no_hist',date,'.pdf']);
%% intensity analysis by area no hist gene 2 vs 3
loc_ctrl = contains(sumIpercell_tot(:,1),'ctrl');
sumIpercell_cntl = sumIpercell_tot(loc_ctrl,:);
sumIpercell_PE = sumIpercell_tot(~loc_ctrl,:);
f8 = figure('color','w','position',[20,20,800,400],'name','Intensity');
h1 = axes('Position',[0.1,0.1,0.35,0.7]);
h2 = axes('Position',[0.47,0.1,0.05,0.7]);
h3 = axes('Position',[0.1,0.82,0.35,0.05]);

c1min = 1000; % thershold on gene 1
c2min = 0;
c3min = 0;
areamin = 0;
c2t = [];
c2f_123 = [];
c3f_123 = [];
c2_pass_cntl = [];
c3_pass_cntl = [];
for jj = 1:length(sumIpercell_cntl)
    cintensity = cell2mat(sumIpercell_cntl(jj,2));
    c1_int = cintensity(:,1);
    c2_int = cintensity(:,2);
    c2t = [c2t;c2_int];
    c3_int = cintensity(:,3);
    area = cintensity(:,4);
    %         c3_int = cintensity(:,3);
    c1_pass = find(c1_int >= c1min);
    c2_pass = find(c2_int >= c2min);
    c3_pass = find(c3_int >= c3min);
    area_pass = find(area >= areamin);
    c23 = intersect(c2_pass,c3_pass);
    c123 = intersect(c1_pass,c23);
    c123area = intersect(c123,area_pass);
    c2f_123 = [c2f_123;c2_int(c123area)];
    c3f_123 = [c3f_123;c3_int(c123area)];
    axes(h1);
    plot(c2_int(c123area),c3_int(c123area),'.b'); hold on
    c2_pass_cntl = [c2_pass_cntl;c2_int(c123area)];
    c3_pass_cntl = [c3_pass_cntl;c3_int(c123area)];
    c2_int(c123area) = []; c3_int(c123area) = [];
    plot(c2_int,c3_int,'.', 'Color', [0.5,0.5,0.5]); hold on
    if jj == length(sumIpercell_cntl)
        xline(c2min, '--'); yline(c3min, '--');
        xlabel(c3name); ylabel(c4name);

        medc2 = median(c2f_123); c275 = prctile(c2f_123,75);
        medc3 = median(c3f_123); c375 = prctile(c3f_123,75);
        xline(medc2,'-',num2str(medc2));  xline(c275,'--r',num2str(c275));
        yline(medc3,'-',num2str(medc3));  yline(c375,'--r',num2str(c375));
       
    end

end

c2fsrt = sort(c2_pass_cntl);
c3fsrt = sort(c3_pass_cntl);
axes(h2);
plot([1:length(c3fsrt)]/length(c3fsrt),c3fsrt); hold on
xline(0,'-',{'0'}); xline(1,'-',{'1'}); xline(0.5);
linkaxes([h1,h2],'y')
set(h2,'yscale','log')
box off; axis off;
axes(h3);
plot(c2fsrt,[1:length(c2fsrt)]/length(c2fsrt));
yline(0,'-',{'0'}); yline(1,'-',{'1'});  yline(0.5);
box off; axis off;
linkaxes([h1,h3],'x')
set(h3,'xscale','log')
set(h1, 'xlim', [500,50000], 'ylim', [300,6500], 'xtick', [1e2,1e3,1e4], 'ytick', [1e2,1e3],'xscale','log','yscale','log')


title(['ECT'],['n =', num2str(length(c2f_123)),'/', num2str(length(c2t))]);




h4 = axes('Position',[0.57,0.1,0.35,0.7]);
h5 = axes('Position',[0.94,0.1,0.05,0.7]);
h6 = axes('Position',[0.57,0.82,0.35,0.05]);
c2t = [];
c2f_123 = [];
c3f_123 = [];
c2_pass_PE = [];
c3_pass_PE = [];
for jj = 1:length(sumIpercell_PE)
    cintensity = cell2mat(sumIpercell_PE(jj,2));
    c1_int = cintensity(:,1);
    c2_int = cintensity(:,2);
    c2t = [c2t;c2_int];
    c3_int = cintensity(:,3);
    %         c3_int = cintensity(:,3);
    c1_pass = find(c1_int >= c1min);
    c2_pass = find(c2_int >= c2min);
    c3_pass = find(c3_int >= c3min);
    area_pass = find(area >= areamin);
    c23 = intersect(c2_pass,c3_pass);
    c123 = intersect(c1_pass,c23);
    c123area = intersect(c123,area_pass);
    c2f_123 = [c2f_123;c2_int(c123area)];
    c3f_123 = [c3f_123;c3_int(c123area)];
    axes(h4);
    plot(c2_int(c123area),c3_int(c123area),'.b'); hold on
    c2_pass_PE = [c2_pass_PE;c2_int(c123area)];
    c3_pass_PE = [c3_pass_PE;c3_int(c123area)];
    c2_int(c123area) = []; c3_int(c123area) = [];
    plot(c2_int,c3_int,'.', 'Color', [0.5,0.5,0.5]); hold on
    if jj == length(sumIpercell_PE)
        xline(c2min, '--'); yline(c3min, '--');
        xlabel(c3name);

        medc2 = median(c2f_123); c275 = prctile(c2f_123,75);
        medc3 = median(c3f_123); c375 = prctile(c3f_123,75);
        xline(medc2,'-',num2str(medc2));  xline(c275,'--r',num2str(c275));
        yline(medc3,'-',num2str(medc3));  yline(c375,'--r',num2str(c375));
    end
end

c2fsrt = sort(c2_pass_PE);
c3fsrt = sort(c3_pass_PE);
axes(h5);
plot([1:length(c3fsrt)]/length(c3fsrt),c3fsrt); hold on
xline(0,'-',{'0'}); xline(1,'-',{'1'}); xline(0.5);
linkaxes([h4,h5],'y')
set(h5,'yscale','log')
box off; axis off;
axes(h6);
plot(c2fsrt,[1:length(c2fsrt)]/length(c2fsrt));
yline(0,'-',{'0'}); yline(1,'-',{'1'});  yline(0.5);
box off; axis off;
linkaxes([h4,h6],'x')
set(h6,'xscale','log')

set(h4, 'xlim', [500,50000], 'ylim', [300,6500], 'xtick', [1e2,1e3,1e4], 'ytick', [1e2,1e3],'xscale','log','yscale','log')

title(['EPE'],['n =', num2str(length(c2f_123)),'/', num2str(length(c2t))]);

% eval(['export_fig ',c2name,'_',c3name,'_','intensity_analysis_CT-PE_no_hist',date,'.pdf']);
% eval(['export_fig ',c3name,'_',c4name,'_','intensity_analysis_ECT-EPE_CYP19A1th1000_',date,'.pdf']);

%% Area vs gene1

loc_ctrl = contains(sumIpercell_tot(:,1),'ctrl');
sumIpercell_cntl = sumIpercell_tot(loc_ctrl,:);
sumIpercell_PE = sumIpercell_tot(~loc_ctrl,:);
f9 = figure('color','w','position',[20,20,800,400],'name','Intensity');
h1 = axes('Position',[0.1,0.18,0.35,0.7]);
% h2 = axes('Position',[0.47,0.1,0.05,0.7]);
% h3 = axes('Position',[0.1,0.82,0.35,0.05]);
c1min = 0;
c1t = [];
c1f_12 = [];
c3f = [];
for jj = 1:length(sumIpercell_cntl)
    cintensity = cell2mat(sumIpercell_cntl(jj,2));
    c1_int = cintensity(:,1);
    c1t = [c1t;c1_int];
    c3_int = cintensity(:,3);
    %         c3_int = cintensity(:,3);
    c1_pass = find(c1_int >= c1min);
    c1f_12 = [c1f_12;c1_int(c1_pass)];
    c3f = [c3f;c3_int(c1_pass)];
    axes(h1);
    plot(c1_int(c1_pass),c3_int(c1_pass),'.b'); hold on
    c1_int(c1_pass) = []; c3_int(c1_pass) = [];
    plot(c1_int,c3_int,'.', 'Color', [0.5,0.5,0.5]); hold on
    if jj == length(sumIpercell_cntl)
        xline(c1min, '--'); yline(c2min, '--');
        xlabel(c2name); ylabel('Area');
%         [ksd1,x1] = ksdensity(c1f,'support',[min(c1f)-0.1,11000],'function','cdf');
%         [ksd2,x2] = ksdensity(c2f,'support',[min(c2f)-0.1,6500],'function','cdf');
        %         [ksd3,x3] = ksdensity(c3_int);
        medc1 = median(c1f_12); c175 = prctile(c1f_12,75);
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

title(['ECT'],['n =', num2str(length(c1f_12)),'/', num2str(length(c1t))]);

h2 = axes('Position',[0.57,0.18,0.35,0.7]);
% h5 = axes('Position',[0.94,0.1,0.05,0.7]);
% h6 = axes('Position',[0.57,0.82,0.35,0.05]);
c1t = [];
c1f_12 = [];
c3f = [];
for jj = 1:length(sumIpercell_PE)
   cintensity = cell2mat(sumIpercell_PE(jj,2));
    c1_int = cintensity(:,1);
    c1t = [c1t;c1_int];
    c3_int = cintensity(:,3);
    %         c3_int = cintensity(:,3);
    c1_pass = find(c1_int >= c1min);
    c1f_12 = [c1f_12;c1_int(c1_pass)];
    c3f = [c3f;c3_int(c1_pass)];
    axes(h2);
    plot(c1_int(c1_pass),c3_int(c1_pass),'.b'); hold on
    c1_int(c1_pass) = []; c3_int(c1_pass) = [];
    plot(c1_int,c3_int,'.', 'Color', [0.5,0.5,0.5]); hold on
    if jj == length(sumIpercell_cntl)
        xline(c1min, '--'); yline(c2min, '--');
        xlabel(c2name); ylabel('Area');
%         [ksd1,x1] = ksdensity(c1f,'support',[min(c1f)-0.1,11000],'function','cdf');
%         [ksd2,x2] = ksdensity(c2f,'support',[min(c2f)-0.1,6500],'function','cdf');
        %         [ksd3,x3] = ksdensity(c3_int);
        medc1 = median(c1f_12); c175 = prctile(c1f_12,75);
        medc3 = median(c3f); c375 = prctile(c3f,75)
%         xline(medc1,'-',num2str(medc1));  xline(c175,'--r',num2str(c175));
%         yline(medc3,'-',num2str(medc3));  yline(c375,'--r',num2str(c375));
       
    end
end
% linkaxes([h4,h5],'y'); linkaxes([h4,h6],'x'); linkaxes([h1,h4],'xy'); linkaxes([h1,h2,h4,h5],'y');
set(h2, 'xlim', [500,13000], 'ylim', [500,10000], 'xtick', [500:1500:13000], 'ytick', [500:1500:10000])
% set(h5, 'ylim', [500,4500],'xlim',[0,1])
% set(h6, 'xlim', [500,1000],'ylim',[0,1])
title(['EPE'],['n =', num2str(length(c1f_12)),'/', num2str(length(c1t))]);

% eval(['export_fig ',c2name,'_',c3name,'_','intensity_analysis_CT-PE_no_hist',date,'.pdf']);
%% violin
g1 = c2name;
g2 = c3name;
g3 = c4name;
fnorm = normpdf(norminv([0.001:0.01:0.999]));

figure('color','w','position',[20,20,150,500],'Name',g3); %change g# according to the name
ha = axes('Position',[0.1,0.3,0.85,0.55]);
logflag = 0;
p99 = 0;% prctile(data(g,:),100);
markervec = 'opds<>';
xt = [];
t_av = [];
% donorid_uni = unique(donor_id_flag_sorted);
notext_flag = 1;
for i=[1:2]
    if i == 1 %ect
        %     axes(ha(k))
        c=1;
        k = c;
        %             gr2 = find(T_cells_tmp==typeplot(c) & condition_sorted(i,:)'==1);%find(fc_time_sorted==fc_time_uni(i));%
        yc_ect = c3_pass_cntl; %change channels
        y_ect = yc_ect; 
        if logflag ==1
            y_ect = y_ect+1;
            yscale = 'log';
        else
            yscale = 'linear';
        end
        if length(y_ect)>5
            %             [f,xi] = ksdensity(y);
            %             fi = interp1(xi,f,y);
            %             fi = fi/max(fi);

            xi_ect = prctile(y_ect,[0.1:1:99.9]);
            [xi_ect,ia_ect] = unique(xi_ect);
            if length(xi_ect)>1
                fi_ect = interp1(xi_ect,fnorm(ia_ect),y_ect);
                fi_ect = fi_ect/max(fi_ect);
            else
                fi_ect = zeros(size(y_ect));
            end

            y_ect = 0.5*rand(length(yc_ect),1)-0.1+y_ect;
            x_ect = fi_ect.*(0.2*rand(length(yc_ect),1)-0.1);

            plot(k+i/3 + x_ect-0.5, y_ect,'.','marker',markervec(i),'markersize',3); hold on;
        else
            plot(-0.5+k+i/3+0.1*rand(length(yc_ect),1)-0.1,0.2*rand(length(yc_ect),1)-0.1+yc_ect....
                ,'.','marker',markervec(i),'markersize',3); hold on;
        end
%                 xt = [xt,k+i/3-0.5];
        %         plot(i+0.5*rand(length(gr2),1)-0.1,0.5*rand(length(gr2),1)-0.1+data(g,gr2)','.'); hold on;
        %         t_ed(i) = median(y);
        t_av(i) = mean(yc_ect);
        %         t_75(i) = prctile(y,75);
        t_98(i) = prctile(y_ect,99);
        xt = [xt,k+i/3-0.5];

    end
    if i == 2 %epe
        %     axes(ha(k))
       
        yc_epe = c3_pass_PE; %change channels
        y_epe = yc_epe;
        if logflag ==1
            y_epe = y_epe+1;
            yscale = 'log';
        else
            yscale = 'linear';
        end
        if length(y_epe)>5
            %             [f,xi] = ksdensity(y);
            %             fi = interp1(xi,f,y);
            %             fi = fi/max(fi);

            xi_epe = prctile(y_epe,[0.1:1:99.9]);
            [xi_epe,ia_epe] = unique(xi_epe);
            if length(xi_epe)>1
                fi_epe = interp1(xi_epe,fnorm(ia_epe),y_epe);
                fi_epe = fi_epe/max(fi_epe);
            else
                fi_epe = zeros(size(y_epe));
            end

            y_epe = 0.5*rand(length(yc_epe),1)-0.1+y_epe;
            x_epe = fi_epe.*(0.2*rand(length(yc_epe),1)-0.1);

            plot(k+i/3 + x_epe-0.5, y_epe,'.','marker',markervec(i),'markersize',3); hold on;
        else
            plot(-0.5+k+i/3+0.1*rand(length(yc_epe),1)-0.1,0.2*rand(length(yc_epe),1)-0.1+yc_epe....
                ,'.','marker',markervec(i),'markersize',3); hold on;
        end

    xt = [xt,k+i/3-0.5];
    %         plot(i+0.5*rand(length(gr2),1)-0.1,0.5*rand(length(gr2),1)-0.1+data(g,gr2)','.'); hold on;
    %         t_ed(i) = median(y);
    t_av(i) = mean(yc_epe);
    %         t_75(i) = prctile(y,75);
    t_98(i) = prctile(y_epe,99);
    end
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
pcomp = ranksum(y_ect,y_epe);
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

% set(gca,'ylim',[0,p99]);
yl = get(gca,'ylim');
text(xt,-(yl(2)-yl(1))*0.05*ones(size(xt)),{'ce','pe'},'fontsize',5)
set(gca,'xtick',[],'TickLabelInterpreter'....
    ,'none','XTickLabelRotation',45);
box off
axis tight
eval(['export_fig HCR_violin_',g2,'_',date,'.pdf']);

% text(xt,-(yl(2)-yl(1))*0.05*ones(size(xt)),repmat({'ce','pe','ct','pl'},1,length(typeplot)),'fontsize',5)
% set(gca,'xtick',[1:length(typeplot)],'XTickLabel',clusteruni(typeplot),'TickLabelInterpreter'....
%     ,'none','XTickLabelRotation',45,'ylim',[-(yl(2)-yl(1))*0.1,yl(2)]);
%% 2 genes correlation scatter plot
typeplot = 1:5;
g1 = {'FLT1'};
g2 = {'PGF'};
% g1_exp = log2(1+data(strcmpi(geneid,g1),:));
% g2_exp = log2(1+data(strcmpi(geneid,g2),:));

%     gr3 = find(T_cells_tmp==typeplot(k) & ctrl_flag == 1);
%     gr4 = find(T_cells_tmp==typeplot(k) & pelate_flag == 1);

   hf = figure('color','w','position',[100,20,675,400]);
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
    g1_gr1 = log2(1+c2_pass_cntl+rand(size(c2_pass_cntl))*0.5); %jitter
    g2_gr1 = log2(1+c3_pass_cntl+rand(size(c3_pass_cntl))*0.5);
    plot(g1_gr1,g2_gr1,'.')
    b = polyfit(g1_gr1,g2_gr1,1);
    hold on
    mn = min(g1_gr1);
    mx = max(g1_gr1);
    plot([mn,mx],b(2)+b(1)*[mn,mx],'-k')
    [r,pv] = corr(g1_gr1,g2_gr1);
    title([string(['ECT, r=',num2str(r)]),['-log(pv) = ',num2str(-log10(pv))],string(['N = ',num2str(length(c2_pass_cntl))])]);
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
    g1_gr2 = log2(1+c2_pass_PE+rand(size(c2_pass_PE))*0.5);
    g2_gr2 = log2(1+c3_pass_PE+rand(size(c3_pass_PE))*0.5);
    plot(g1_gr2,g2_gr2,'.')
    b = polyfit(g1_gr2,g2_gr2,1);
    hold on
    mn = min(g1_gr2);
    mx = max(g1_gr2);
    plot([mn,mx],b(2)+b(1)*[mn,mx],'-k')
    [r,pv] = corr(g1_gr2,g2_gr2);
    title([string(['EPE, r=',num2str(r)]),['-log(pv) = ',num2str(-log10(pv))],string(['N = ',num2str(length(c2_pass_PE))])]);
    xlabel(upper(g1))
    box off
%     axes(ha(3))
  
    linkaxes([h1,h2],'xy')
%     eval(['export_fig Two_genes_correlation_TB/',clusteruni{typeplot(k)},'_',char(g1),'-',char(g2),'_',date,'.pdf']);
    