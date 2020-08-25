function [ output ] = dicom2nifti(input_path , input_name , output_name, alg)
%dicom2nifti uses dcm2nii
%   Inputs:
%       input_path - path to folder with sequences
%       input_name - name of the sequence to convert
%       output_name - name of output NiFTI file name
%       alg - 'dcm2nii' or 'mri_convert'
%
%   Output:
%       output - 1 for sucess 0 for error;
%
% author: CRaymond@mednet.ucla.edu
% date: April 17th 2017

disp(['     Creating ' output_name]);
sequencesfolders = dir(input_path);
sequencesfolders(~cellfun('isempty',strfind({sequencesfolders.name},'.'))) = []; % keep just the folders on the list (removes '.' and '..')
try
    fullname = sequencesfolders(~cellfun('isempty',strfind({sequencesfolders.name},input_name))).name;
catch
    disp(['File ' input_name ' not found']);
    output = 0;
    return;
end

folder_temp = [input_path '/' fullname];
folder = strrep(folder_temp,' ' , '\ '); % replace all blank spaces in the names with \space
folder = strrep(folder,'&' , '\&'); % replace all & in the names with \&
folder = strrep(folder,'<' , '\<'); % replace all & in the names with \&
folder = strrep(folder,'(' , '\('); % replace all ( in the names with \(
folder = strrep(folder,')' , '\)'); % replace all ) in the names with \)

dicoms = dir([folder_temp '/']); % get files inside the folder
dicoms(cellfun('isempty',strfind({dicoms.name},'.dcm'))) = []; % remove all files that are not dicoms from list
try
    dicom_temp = [folder_temp '/' dicoms(1).name]; % create path to the first dicom in the series
    dicom = strrep(dicom_temp,' ' , '\ '); % replace all blank spaces in the names with \space
    
    outputpath_temp = [input_path '/' output_name];
    outputpath = strrep(outputpath_temp,' ' , '\ '); % replace all blank spaces in the names with \space
    outputpath = strrep(outputpath,'(' , '\('); % replace all & in the names with \&
    outputpath = strrep(outputpath,')' , '\)'); % replace all & in the names with \&
    
    if strcmp(alg, 'mri_convert')
        cmd = ['mri_convert ' dicom ' ' outputpath]; %create bash command
        status = system(cmd); %execute the command
    end
    
    if strcmp(alg, 'dcm2nii')
        [status, ~] = dcm2nii_all(folder, input_path, outputpath);
    end
    
    if ~status
        disp(['     ' output_name ' created']);
    else % try to decompress
        warning('on'); warning( ['     ERROR: unable to create ' output_name ] ); warning('off');
        disp(['     INFO: Trying to decompress with jpeg ' output_name ]);
        cmd = ['cd ' folder '; find . -name "*.dcm" -type f -exec dcmdjpeg "{}" "{}" \;']; %create bash command
        [~, ~] = system(cmd); %execute the command
        [status, ~] = dcm2nii_all(folder, input_path, outputpath);
        if ~status
            disp(['     ' output_name ' created']);
        else
            disp(['     INFO: Trying to decompress with rle ' output_name ]);
            cmd = ['cd ' folder '; find . -name "*.dcm" -type f -exec dcmdrle "{}" "{}" \;']; %create bash command
            [~, ~] = system(cmd); %execute the command
            [status, ~] = dcm2nii_all(folder, input_path, outputpath);
            if ~status
                disp(['     ' output_name ' created']);
            else
               warning('on'); warning( ['     ERROR: unable to decompress and create ' output_name ]); warning('off');
            end
        end
    end
    output = 1; % return sucess
catch
    warning('on'); warning(['     ERROR: unable to create path to dicom series ' output_name ]); warning('off');
end
end

function [status, cmdout] = dcm2nii_all(folder, input_path, outputpath)

input_path_nospace = strrep(input_path,' ' , '\ '); % replace all blank spaces in the names with \space
input_path_nospace = strrep(input_path_nospace,'&' , '\&'); % replace all & in the names with \&
input_path_nospace = strrep(input_path_nospace,'<' , '\<'); % replace all & in the names with \&
input_path_nospace = strrep(input_path_nospace,'(' , '\('); % replace all ( in the names with \(
input_path_nospace = strrep(input_path_nospace,')' , '\)'); % replace all ( in the names with \)
cmd = ['cd ' folder '; dcm2nii -d N -e N -p Y -o ' input_path_nospace ' *']; %create bash command
[status, cmdout] = system(cmd); %execute the command

if status
    warning('on'); warning(['     ERROR: unable to create ' outputpath ' files had problems']);warning('off');
    warning('on'); warning( cmdout ); warning('off');
else % try to decompress
    
    % get files created
    a = strsplit(cmdout,'GZip...');
    b = strsplit(char(a(2)),'.nii.gz');
    name_on_file = char(b(1));
    
    % rename nii
    cmd = ['mv ' input_path_nospace '/' name_on_file '.nii.gz ' outputpath ];
    system(cmd);
    
    % delete remaining niftis
    niftis = dir(input_path);
    niftis(cellfun('isempty',strfind({niftis.name},name_on_file))) = []; % keep just the nii we created
    for i = 1:1:numel(niftis)
        delete([input_path '/' niftis(i).name]);
    end
end
end
