function [] = sort_me_cest()

scanfolders = dir('DICOM');
scanfolders(contains({scanfolders.name},'.')) = [];
scanfolders(~contains({scanfolders.name},'me_cest')) = [];
idx = strfind(scanfolders(1).name,'_');
scanfolders(1).idx = str2double(scanfolders(1).name((idx(end)+1):end));
scanfolders(2).idx = str2double(scanfolders(2).name((idx(end)+1):end));
[~,idx] = sort([scanfolders.idx]);
scanfolders = scanfolders(idx);

movefile(['DICOM/' scanfolders(1).name],['DICOM/me_cest_image']);
movefile(['DICOM/' scanfolders(2).name],['DICOM/me_cest_phase']);

end