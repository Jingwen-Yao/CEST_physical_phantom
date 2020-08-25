function [TEs] = extractTE(DicomFolder, DicomNameStr)

sequencefolders = dir(DicomFolder);
sequencefolders(~cellfun('isempty',strfind({sequencefolders.name},'.'))) = [];
sequencefolders(cellfun('isempty',strfind({sequencefolders.name},DicomNameStr))) = [];
seqFolderPath = [DicomFolder '/' sequencefolders(1).name];

dicomfile = dir(seqFolderPath);
dicomfile(cellfun('isempty',strfind({dicomfile.name},'.dcm'))) = [];

TEs = zeros(1,numel(dicomfile));
for i = 1:numel(dicomfile)
    try
        metadata = dicominfo([seqFolderPath '/' dicomfile(i).name], 'UseDictionaryVR', true);
        TEs(i) = metadata.EchoTime;
    catch
    end
end

TEs = sort(unique(TEs));
TEs(TEs == 0) = [];

end