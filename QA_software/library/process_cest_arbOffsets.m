function [] = process_cest_arbOffsets(cest_flag, offsets)

if nargin < 2
    offsets = [-3.5:0.1:-2.5 -0.3:0.1:0.3 2.5:0.1:3.5];
end

if strcmp(cest_flag, 'me_cest')
    if ~isfile('me_cest_s0_TE1.nii.gz')
        nii = load_untouch_nii('me_cest.nii.gz');
        Nspec = length(offsets);
        Ns0 = size(nii.img,4)/4-Nspec;
        
        % motion correction, create mask
        [st, cou] =  system(['mcflirt -in me_cest.nii.gz -refvol 2 -out me_cest_all_s0_moco.nii.gz;'...
            '3dcalc -a me_cest_all_s0_moco.nii.gz[' num2str(Ns0) '..' num2str(Ns0+Nspec-1) '] -expr a -prefix me_cest_TE1_moco.nii.gz;'...
            '3dcalc -a me_cest_all_s0_moco.nii.gz[' num2str(Ns0*2+Nspec) '..' num2str(Ns0*2+Nspec*2-1) '] -expr a -prefix me_cest_TE2_moco.nii.gz;'...
            '3dTstat -mean -prefix me_cest_s0_TE1.nii.gz me_cest_all_s0_moco.nii.gz[0..' num2str(Ns0-1) '];'...
            '3dTstat -mean -prefix me_cest_s0_TE2.nii.gz me_cest_all_s0_moco.nii.gz[' num2str(Ns0+Nspec) '..' num2str(Ns0*2+Nspec-1) '];'...
            '3dTstat -mean -prefix me_cest_s0_TE3.nii.gz me_cest_all_s0_moco.nii.gz[' num2str(Ns0*2+Nspec*2) '..' num2str(Ns0*3+Nspec*2-1) '];'...
            '3dTstat -mean -prefix me_cest_s0_TE4.nii.gz me_cest_all_s0_moco.nii.gz[' num2str(Ns0*3+Nspec*3) '..' num2str(Ns0*4+Nspec*3-1) '];'...
            'bet2 me_cest_s0_TE1.nii.gz cest -m -n -f 0.2;'...
            'rm me_cest_all_s0_moco.nii.gz; rm me_cest_all_s0.nii.gz;']);
        warning('on'); if st, warning(cou); end; warning('off')
    end
    
    if ~isfile('MTRasym_IDE.nii.gz')
        nii = load_untouch_nii('me_cest_s0_TE1.nii.gz');
        % IDE B0 correction
        [B0map, MTRasym, ~, ~, ~] = B0_IDE('me_cest', 0, offsets);
        postfix = 'IDE';
        nii.img = B0map; save_untouch_nii(nii, ['B0_' postfix '.nii.gz']);
        nii.img = MTRasym; save_untouch_nii(nii, ['MTRasym_' postfix '.nii.gz']);
    end
end