function [] = createROItemplate(template_folder)

nii = load_untouch_nii([template_folder '/t1_template.nii.gz']);
ref_img = int16(nii.img);

nii = load_untouch_nii([template_folder '/ROI.nii.gz']);
ROI_img = int16(nii.img);

ROIedit_img = ROI_img*0;
for i = 1:18
    
    disp(['ROI #' num2str(i)]);
    ROI = ROI_img == i;
    ROI_activeContour = int16(activecontour(ref_img,ROI,12));
    
    se1 = strel('sphere',2);
    ROIedit = imerode(ROI_activeContour,se1);
    se2 = strel('disk',8);
    ROIedit = imopen(ROIedit,se2);
    ROIedit_img = ROIedit_img + i*ROIedit;
    
end

% figure; volshow(ROI_img);
% figure; volshow(ROIedit_img);
% figure; imagesc(ROIedit_img(:,:,98));

nii.img = ROIedit_img;
save_untouch_nii(nii,[template_folder  '/ROI_template.nii.gz']);

end