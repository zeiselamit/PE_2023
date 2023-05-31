function [bw,pos] = binary_dot_per_channel_x20_v2(chanlname,enhan_prct,p_maxtranf)
%enhan_prct ~99
%p_maxtranf ~10000
nbit = 16;
disth = round(10/(6.5/20));
% T = medfilt2(chanlname,[5,5]);
% chanlname = chanlname-T;
% chanlname = imadjust(chanlname,[0,double(prctile(chanlname(:),enhan_prct))/(2^nbit-1)],[0,1]);

bw = imextendedmax(chanlname,p_maxtranf, 4);
bw = imfill(bw,'holes');
bw = bw - bwareaopen(bw,30);
bw = bwareaopen(bw,2);
pos = regionprops(bw>0,'Centroid','Area','PixelIdxList');
xy = reshape([pos.Centroid]',2,length(pos))';
d = squareform(pdist(xy,'euclidean'),'tomatrix');
sumd = sum(d<disth,2);
rm = find(sumd<1);
for i=1:length(rm)
    bw(pos(rm(i)).PixelIdxList) = 0;
end
pos = regionprops(bw>0,'Centroid','Area');
% pos = reshape([pos.Centroid]',2,length(pos))';
pos = cat(1,pos.Centroid);
