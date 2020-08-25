function [] = plotMaps...
    (T2img, mapFile, IndSlice, ROImask, ROI_img, colorRange, textColor, ...
    label_str, verticalFLAG)

if nargin < 9
    verticalFLAG = 0;
end

nii = load_untouch_nii(mapFile);
Map = double(nii.img(:,:,IndSlice));

% crop images
T2img = T2img(17:239,17:239,:);
Map = Map(17:239,17:239,:);
ROImask = ROImask(17:239,17:239,:);
ROI_img = ROI_img(17:239,17:239,:);

if verticalFLAG
    T2img = [T2img(:,:,1); T2img(:,:,2)];
    Map = [Map(:,:,1); Map(:,:,2)];
    ROImask = [ROImask(:,:,1); ROImask(:,:,2)];
    ROI_img = [ROI_img(:,:,1); ROI_img(:,:,2)];
end

if verticalFLAG
    figure1 = figure('rend','painters','pos',[10 10 450 700]);
else
    figure1 = figure('rend','painters','pos',[10 10 800 400]);
end
ax1 = axes('Parent',figure1);
ax2 = axes('Parent',figure1);
imshow(T2img(:,:),[0 prctile(T2img(:),99)],'Colormap',gray,'Parent',ax1);
Foreground = imshow(Map(:,:),colorRange,...
    'Colormap',cool,'Parent',ax2);
if verticalFLAG
    c = colorbar('eastoutside','AxisLocation','out','FontSize',20);
else
    c = colorbar('southoutside','AxisLocation','out','FontSize',20);
end
set(get(c,'label'),'string',label_str);
axPos = ax2.Position;
ax1.Position = axPos;
set(Foreground, 'AlphaData', ROImask(:,:));
axis tight; axis off; axis equal;

[text_str, text_pos] = textOnImg(ROI_img(:,:));
for ii = 1:length(text_str)
    text(text_pos(ii,1)-8,text_pos(ii,2),text_str{ii},...
        'FontSize',20,'Color',textColor);
end

end