function[TE] = T2_mapping(MSEfile, contrastFLAG)

% Load image
nii = load_untouch_nii(MSEfile);
SE = double(nii.img);
matSize = size(SE);

% echo time (change according to protocol)
if strcmp(contrastFLAG,'T2')
    try
        TE = extractTE('../DICOM', 'ME_SE');
    catch
        TE = 49:49:392; % ms
    end
elseif strcmp(contrastFLAG,'T2star')
    try
        TE = extractTE('../DICOM', 'ME_GRE');
    catch
        TE = [10 20 40 60 80]; % ms
    end
else
    msg = 'Contrast flag is neither T2 nor T2star';
    error(msg);
end

% save TEs in txt file
fileID = fopen([contrastFLAG '.txt'],'wt');
fprintf(fileID,'%.1f\n',TE);
fclose(fileID);

% Calculate T2 map
% Initialize matrices
R2map = zeros(matSize(1:3));
T2map = zeros(matSize(1:3));
GOFmap = zeros(matSize(1:3));

mask = SE(:,:,:,1) > max(SE(:))/10;

if strcmp(contrastFLAG,'T2')
    disp('T2 mapping calculation...')
else
    disp('T2star mapping calculation...')
end
for sli = 1:matSize(3)
    disp([' -- Slice ' num2str(sli)]);
    for row = 1:matSize(1)
        for col = 1:matSize(2)
            if mask(row,col,sli) < 1
                R2map(row,col,sli) = nan;
                T2map(row,col,sli) = nan;
                GOFmap(row,col,sli) = nan;
            else
                % linear fitting of ln(S) and TE
                mdl = fitlm(TE(1:end)',log(squeeze(SE(row,col,sli,1:end))));
                R2map(row,col,sli) = -mdl.Coefficients.Estimate(2);
                T2map(row,col,sli) = 1/R2map(row,col,sli);
                GOFmap(row,col,sli) = mdl.Rsquared.Ordinary;
            end
        end
    end
end

R2map = R2map*1000; % ms^-1 -> s^-1

% Save T2* map and R2* map to NIFTI
% change data type to double
system(['3dcalc -a ' MSEfile '[0] -expr a -datum float -prefix dummy.nii.gz']);

nii = load_untouch_nii('dummy.nii.gz');

if strcmp(contrastFLAG,'T2')
    nii.img = R2map;
    save_untouch_nii(nii,'R2map.nii.gz');
    nii.img = T2map;
    save_untouch_nii(nii,'T2map.nii.gz');
    nii.img = GOFmap;
    save_untouch_nii(nii,'GOFmap.nii.gz');
else
    nii.img = R2map;
    save_untouch_nii(nii,'R2starmap.nii.gz');
    nii.img = T2map;
    save_untouch_nii(nii,'T2starmap.nii.gz');
    nii.img = GOFmap;
    save_untouch_nii(nii,'GOFstarmap.nii.gz');
end

system('rm dummy.nii.gz');

end
