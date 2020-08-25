function [text_str, text_pos] = textOnImg(ROI_img)

Nsample = max(ROI_img(:));

text_str = cell(Nsample,1);
text_pos = zeros(Nsample,2);
for ii = 1:Nsample
    text_str{ii} = [num2str(ii,'%02d')];
    ROImask = ROI_img(:,:) == ii;
    s = regionprops(ROImask,'centroid');
    text_pos(ii,:)  = s.Centroid;
end

end