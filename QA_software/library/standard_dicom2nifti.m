function [] = standard_dicom2nifti(pathToData)

dicom_folder = [pathToData '/DICOM'];
nifti_folder = [pathToData '/NIFTI'];

cd(pathToData);
if ~isfolder(dicom_folder)
    mkdir(dicom_folder);
end

dicomfolders = dir(pathToData);
dicomfolders(contains({dicomfolders.name},'.') | ...
    contains({dicomfolders.name},'NIFTI') | ...
    contains({dicomfolders.name},'Report') | ...
    contains({dicomfolders.name},'DICOM')) = [];
for d = 1:numel(dicomfolders)
    movefile(dicomfolders(d).name,dicom_folder);
end

if ~isfolder(nifti_folder)
    mkdir(nifti_folder);
end

if ~isfile([nifti_folder '/t2.nii.gz'])
    dicom2nifti(dicom_folder , 't2_tse' , '../NIFTI/t2.nii.gz', 'dcm2nii');
end
if ~isfile([nifti_folder '/t1.nii.gz'])
    dicom2nifti(dicom_folder , 't1_mprage' , '../NIFTI/t1.nii.gz', 'dcm2nii');
end
if ~isfile([nifti_folder '/MSE.nii.gz'])
    dicom2nifti(dicom_folder , 'ME_SE' , '../NIFTI/MSE.nii.gz', 'dcm2nii');
end
if ~isfile([nifti_folder '/MGE.nii.gz'])
    dicom2nifti(dicom_folder , 'ME_GRE' , '../NIFTI/MGE.nii.gz', 'dcm2nii');
end
if ~isfile([nifti_folder '/gre_fieldmap.nii.gz'])
    dicom2nifti(dicom_folder , 'gre_field' , '../NIFTI/gre_fieldmap.nii.gz', 'dcm2nii');
end

if ~isfile([nifti_folder '/me_cest.nii.gz'])
    if ~isfolder('DICOM/me_cest_image')
        sort_me_cest();
    end
    dicom2nifti(dicom_folder , 'me_cest_image' , '../NIFTI/me_cest.nii.gz', 'dcm2nii');
end

if ~isfolder([nifti_folder '/T1mapping'])
    scanfolders = dir('DICOM');
    scanfolders(contains({scanfolders.name},'.')) = [];
    scanfolders(~contains({scanfolders.name},'TI')) = [];
    idx1 = strfind(scanfolders(1).name,'TI');
    
    mkdir('DICOM/T1mapping');
    mkdir('NIFTI/T1mapping');
    for i = 1:length(scanfolders)
        movefile([dicom_folder '/' scanfolders(i).name],...
            [dicom_folder '/T1mapping/' scanfolders(i).name]);
        idx2 = strfind(scanfolders(i).name,'_');
        TI = scanfolders(i).name((idx1(end)+2):(idx2(end)-1));
        dicom2nifti([dicom_folder '/T1mapping'] , scanfolders(i).name, ...
            ['../../NIFTI/T1mapping/TI' TI '.nii.gz'], 'dcm2nii');
    end
end

end