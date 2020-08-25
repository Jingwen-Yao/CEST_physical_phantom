function[TI] = T1_mapping(TIfolder)

% Load image
TIfile = dir(TIfolder);
TIfile(cellfun('isempty',strfind({TIfile.name},'.nii.gz'))) = [];
TIfile(cellfun('isempty',strfind({TIfile.name},'TI'))) = [];

nii = load_untouch_nii([TIfolder '/' TIfile(1).name]);
SE = zeros([size(nii.img) numel(TIfile)]);
TIvalue = zeros(1,numel(TIfile));
for i = 1:numel(TIfile)
    nii = load_untouch_nii([TIfolder '/' TIfile(i).name]);
    SE(:,:,:,i) = double(nii.img);
    TIvalue(i) = str2double(TIfile(i).name(3:end-7));
end
matSize = size(SE);
[TI,idx] = sort(TIvalue);
SE = SE(:,:,:,idx);

% save TIs in txt file
fileID = fopen(['T1.txt'],'wt');
fprintf(fileID,'%.1f\n',TI);
fclose(fileID);

% Calculate T1 map
% Initialize matrices
R1map = zeros(matSize(1:3));
T1map = zeros(matSize(1:3));
GOFmap = zeros(matSize(1:3));

mask = SE(:,:,:,end) > max(SE(:))/10;

disp('T1 mapping calculation...')
ft = fittype('abs(a*(1-2*exp(-b*x)))');
for sli = 1:matSize(3)
    disp([' -- Slice ' num2str(sli)]);
    for row = 1:matSize(1)
%         disp([' -- Row ' num2str(row)]);
        for col = 1:matSize(2)
            if mask(row,col,sli) < 1
                R1map(row,col,sli) = nan;
                T1map(row,col,sli) = nan;
                GOFmap(row,col,sli) = nan;
            else
                % correct values beyond inversion point
%                 if ndims(SE) == 4
%                     [~, i] = min(SE(row,col,sli,:));
%                     for it = 1:i-1
%                         SE(row,col,sli,it) = - SE(row,col,sli,it);
%                     end
%                 end
                % fitting of S and TI
                [F,gof] = fit(TI',squeeze(SE(row,col,sli,:)),ft,...
                    'StartPoint',[2000 0.001],...
                    'Upper',[4096 100],...
                    'Lower',[0 0]);
                R1map(row,col,sli) = F.b;
                T1map(row,col,sli) = 1/F.b;
                GOFmap(row,col,sli) = gof.rsquare;
            end
        end
    end
end

R1map = R1map*1000; % ms^-1 -> s^-1

% Save T1 map and R1 map to NIFTI
% change data type to double

system(['3dcalc -a ' TIfolder '/' TIfile(1).name ' -expr a -datum float -prefix dummy.nii.gz']);

nii = load_untouch_nii('dummy.nii.gz');

nii.img = R1map;
save_untouch_nii(nii,'R1map.nii.gz');
nii.img = T1map;
save_untouch_nii(nii,'T1map.nii.gz');
nii.img = GOFmap;
save_untouch_nii(nii,'GOFT1map.nii.gz');

% save multi TI 4D data
TIfile = TIfile(idx);
system(['3dcalc -a ' TIfolder '/' TIfile(1).name ' -expr a -prefix MTI.nii.gz']);
pause(10);
for i = 2:length(TI)
    system(['3dTcat MTI.nii.gz ' TIfolder '/' TIfile(i).name ' -prefix temp.nii.gz']);
    system('rm MTI.nii.gz');
    system('mv temp.nii.gz MTI.nii.gz');
end

system('rm dummy.nii.gz');

end