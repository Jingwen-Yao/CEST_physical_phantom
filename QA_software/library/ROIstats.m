function [ROIstats] = ROIstats(ROIfile,mapfile)

% ROIstats size: Nx4
% mean, std, median, mad

nii = load_untouch_nii(ROIfile);
ROI = double(nii.img);

nii = load_untouch_nii(mapfile);
Map = double(nii.img);

Nroi = max(ROI(:));

ROIstats = zeros(Nroi,4);
for i = 1:Nroi
    ROImask = double(ROI == i);
    ROIdata = nonzeros(ROImask.*Map);
    ROIdata(isnan(ROIdata)) = [];
    if strcmp(mapfile,'T2map.nii.gz') || strcmp(mapfile,'R2map.nii.gz') ...
            || strcmp(mapfile,'T1map.nii.gz') || strcmp(mapfile,'R1map.nii.gz')
        ROIdata(ROIdata <= 0) = [];
    end
    ROIdata(ROIdata > prctile(ROIdata,99) | ROIdata < prctile(ROIdata,1)) = [];
    ROIstats(i,1) = mean(ROIdata);
    ROIstats(i,2) = std(ROIdata);
    ROIstats(i,3) = median(ROIdata);
    ROIstats(i,4) = mad(ROIdata,1);
end

end