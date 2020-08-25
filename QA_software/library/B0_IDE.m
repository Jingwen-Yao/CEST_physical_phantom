function [B0map, MTRasym, time_3d, Z, ppm_n] = B0_IDE(flag, ABTC_flag, offsets, remove_offsetInd)

% Load data and parameters
if nargin ~= 4
    remove_offsetInd = [];
end
if nargin >= 3
    ppm = offsets;
    ppm_n = ppm; ppm_n(ppm == 0) = 0.00001;
    ppm_ex = ppm(1):0.002:ppm(end);
    n = length(ppm);
elseif strcmp(flag,'me_cest')
    if ABTC_flag
        ppm = [-3.5:0.1:-2.5 -0.3:0.1:0.3 2.5:0.1:3.5];
        ppm_n = [-3.5:0.1:-2.5 -0.3 -0.2 -0.1 0.00001 0.1 0.2 0.3 2.5:0.1:3.5];
    else
        ppm = [-3.5:0.1:-2.5 -0.3:0.1:0.3 2.5:0.1:3.5];
        ppm_n = [-3.5:0.1:-2.5 -0.3 -0.2 -0.1 0.00001 0.1 0.2 0.3 2.5:0.1:3.5];
    end
    ppm_ex = -3.5:0.002:3.5;
    n = length(ppm);
elseif strcmp(flag,'cest')
    if ABTC_flag
        ppm = [-3.5:0.1:-2.5 -0.3:0.1:0.3 2.5:0.1:3.5];
        ppm_n = [-3.5:0.1:-2.5 -0.3 -0.2 -0.1 0.00001 0.1 0.2 0.3 2.5:0.1:3.5];
    else
        ppm = [-3.4:0.1:-2.6 -0.3:0.1:0.3 2.6:0.1:3.4];
        ppm_n = [-3.4:0.1:-2.6 -0.3 -0.2 -0.1 0.00001 0.1 0.2 0.3 2.6:0.1:3.4];
    end
    ppm_ex = -3.5:0.002:3.5;
    n = length(ppm);
else
    warning('Invalid flag.');
end

% Load images
nii = load_untouch_nii('cest_mask.nii.gz');
mask = double(squeeze(nii.img)); mask = imfill(mask,'holes'); mask(mask == 0) = nan;
if strcmp(flag,'me_cest')
    nii = load_untouch_nii('me_cest_TE1_moco.nii.gz');
    CEST(:,:,:,1:n) = double(squeeze(nii.img)).*repmat(mask,[1 1 1 n]);
    nii = load_untouch_nii('me_cest_TE2_moco.nii.gz');
    CEST(:,:,:,1+n:2*n) = double(squeeze(nii.img)).*repmat(mask,[1 1 1 n]);
    nii = load_untouch_nii('me_cest_s0_TE1.nii.gz');
    CEST_s0_TE(:,:,:,1) = double(squeeze(nii.img));
    nii = load_untouch_nii('me_cest_s0_TE2.nii.gz');
    CEST_s0_TE(:,:,:,2) = double(squeeze(nii.img));
elseif strcmp(flag,'cest')
    nii = load_untouch_nii('cest_moco.nii.gz');
    CEST = double(squeeze(nii.img)).*repmat(mask,[1 1 1 n]);
    nii = load_untouch_nii('cest_s0_moco.nii.gz');
    CEST_s0_TE = double(squeeze(nii.img));
else
    warning('Invalid flag.');
end

%% Normalize
if strcmp(flag,'me_cest')
    for i = 1:2
        Z(:,:,:,(i-1)*n+1:i*n) = CEST(:,:,:,(i-1)*n+1:i*n)./repmat(CEST_s0_TE(:,:,:,i),[1 1 1 n]);
    end
    
    % denoise
%     t = 3;
%     f = 2;
%     h1 = 1;
%     h2 = 50;
%     selfsim = 0;
%     for slc = 1:size(Z,3)
%         CEST_2D = squeeze(Z(:,:,slc,:));
%         Z(:,:,slc,:) = simple_nlm2Dspec(CEST_2D,t,f,h1,h2,selfsim);
%     end
    
    Z = (Z(:,:,:,1:n)+Z(:,:,:,1+n:n+n))/2;
elseif strcmp(flag,'cest')
    Z = CEST./repmat(CEST_s0_TE,[1 1 1 n]);
else
    warning('Invalid flag.');
end

if ~isempty(remove_offsetInd)
    Z = Z(:,:,:,~remove_offsetInd);
    ppm_n = ppm_n(~remove_offsetInd);
    ppm = ppm(~remove_offsetInd);
end
matSize = size(Z);

%% down sample
tic;
B0_DS = zeros([matSize(1:3) 7]);
MSE = zeros(matSize(3),7);

fun = fittype('c - c./((x-a).^2+b)'); % - f./((x-a-3).^2+e)');
b0 = [0 3 1]; % 0.5 0.01];

weight = ones(size(ppm_n)); weight(abs(ppm_n) < 0.5) = 10;

%% Loop through slices
h = waitbar(0,'JingwenSOCpipeline::B0-IDE Loop through slices...');
for slc = 1:matSize(3)
    
    Zslice_kmeans = squeeze(Z(:,:,slc,:));
    mask_slice = squeeze(mask(:,:,slc));
    
    for m = 1:7
        
        bmin = [-1/m 0 0];
        bmax = [1/m 5 2];
        B0 = zeros(matSize(1)*matSize(2),1);
        X = reshape(Zslice_kmeans,matSize(1)*matSize(2),[]);
        nnc = sum(~isnan(X(:,1)));
        
        ds = 3^m;
        if ds > nnc
            ds = nnc;
        end
        if ds ~= 0
            idx = kmeans(X,ds,'Distance','correlation');
        else
            break
        end
        
        spectra = zeros(ds,matSize(4));
        for i = 1:ds
            ind_i = find(idx == i);
            Zsub = X(ind_i,:);
            spectra(i,:) = nanmean(Zsub,1);
            if sum(isnan(spectra(i,:))) == 0 && sum(isinf(spectra(i,:))) == 0
                % spectra(i,abs(ppm_n) < 0.5) = medfilt1(spectra(i,abs(ppm_n) < 0.5),5);
                [f,gof] = fit(ppm_n',squeeze(spectra(i,:))'...
                    ,fun,'StartPoint',b0,'Weights',weight,'Lower',bmin,'Upper',bmax);
                if gof.rsquare > 0.7 && abs(f.a) < 0.6
                    off = f.a;
                else
                    f = interp1(ppm_n,squeeze(spectra(i,:)),ppm_ex,'makima'); % makima
                    ind = find(f == min(f(abs(ppm_ex) < 1)));
                    off = ppm_ex(ind(1));
                end
                B0(ind_i,1) = off;
            end
        end
        
        B0 = reshape(B0,size(mask_slice)).*mask_slice;
        
        % update Zslice
        for i = 1:matSize(1)
            for j = 1:matSize(2)
                if isnan(mask_slice(i,j)) || sum(isnan(Zslice_kmeans(i,j,:))) ~= 0
                    Zslice_kmeans(i,j,:) = nan;
                else
                    ppm_data = ppm_n - B0(i,j);
                    z_data = squeeze(Zslice_kmeans(i,j,:))';
                    Zslice_kmeans(i,j,:) = interp1(ppm_data,z_data,ppm_n,'makima','extrap');
                end
            end
        end
        
        B0_DS(:,:,slc,m) = B0;
        MSE(slc,m) = nansum(nansum(B0.^2))./nansum(mask_slice(:));
        if MSE(slc,m) < sum(MSE(slc,:))/10 && m >= 3
            break
        end
    end
    
    MSE(slc,:);
    waitbar(slc/matSize(3));
end
close(h);
time_3d = toc;

B0map = sum(B0_DS,4);
B0map(isnan(B0map)) = 0;
B0map = imgaussfilt3(B0map,1);

%% MTRasym
Z_B0 = zeros(matSize);
for k = 1:matSize(3)
    for i = 1:matSize(1)
        for j = 1:matSize(2)
            if mask(i,j,k) ~= 1 || sum(isnan(Z(i,j,k,:))) ~= 0
                Z_B0(i,j,k,:) = nan;
            else
                z_data = [squeeze(Z(i,j,k,:))' 0];
                Z_B0(i,j,k,:) = interp1([ppm_n - B0map(i,j,k), 0],z_data,ppm_n,'makima','extrap');
            end
        end
    end
end

neg_offset = (ppm >= -3.2) & (ppm <= -2.8);
pos_offset = (ppm >= 2.8) & (ppm <= 3.2);
MTRasym = squeeze(mean(Z_B0(:,:,:,neg_offset),4) - mean(Z_B0(:,:,:,pos_offset),4))*100;
