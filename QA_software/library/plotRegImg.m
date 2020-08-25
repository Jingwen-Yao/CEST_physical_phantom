function [ssimval] = plotRegImg(t2_reg,t2_targ,lev1, lev2,lev3)

t2_reg = t2_reg./nanmedian(t2_reg(:));
t2_targ = t2_targ./nanmedian(t2_targ(:));

midRow = ceil(size(t2_reg,1)/2);
midCol = ceil(size(t2_reg,2)/2);
midSli = ceil(size(t2_reg,3)/2);

t2_targ_mid = [standardizeImg(t2_targ(:,:,midSli)) ...
    standardizeImg(t2_targ(midRow,:,:)) ...
    standardizeImg(t2_targ(:,midCol,:))];
t2_targ_mid(t2_targ_mid>prctile(t2_targ_mid(:),99)) = ...
    prctile(t2_targ_mid(:),99);

t2_reg_mid = [standardizeImg(t2_reg(:,:,midSli)) ...
    standardizeImg(t2_reg(midRow,:,:)) ...
    standardizeImg(t2_reg(:,midCol,:))];
t2_reg_mid(t2_reg_mid>prctile(t2_reg_mid(:),99)) = ...
    prctile(t2_reg_mid(:),99);

t2_targ_sli = [standardizeImg(t2_targ(:,:,lev1)) ...
    standardizeImg(t2_targ(:,:,lev2)) ...
    standardizeImg(t2_targ(:,:,lev3))];
t2_targ_sli(t2_targ_sli>prctile(t2_targ_sli(:),99)) = ...
    prctile(t2_targ_sli(:),99);

t2_reg_sli = [standardizeImg(t2_reg(:,:,lev1)) ...
    standardizeImg(t2_reg(:,:,lev2)) ...
    standardizeImg(t2_reg(:,:,lev3))];
t2_reg_sli(t2_reg_sli>prctile(t2_reg_sli(:),99)) = ...
    prctile(t2_reg_sli(:),99);

[~,padNumber] = standardizeImg(t2_reg(midRow,:,:));

figure;
imshowpair(t2_targ_mid,t2_reg_mid); hold on;
line(512+padNumber(3)+[lev1 lev1],[1 256],...
    'Color','y','LineWidth',2,'LineStyle','--');
line(512+padNumber(3)+[lev2 lev2],[1 256],...
    'Color','c','LineWidth',2,'LineStyle','--');
line(512+padNumber(3)+[lev3 lev3],[1 256],...
    'Color','r','LineWidth',2,'LineStyle','--');
export_fig('Registered_maps/imgReg_mid.png','-m4');
close;

figure;
imshowpair(t2_targ_sli,t2_reg_sli); hold on;
rectInd = [16 16; 16 240; 240 16; 240 240];
drawRectangle(rectInd,'y');
rectInd(:,1) = rectInd(:,1) + 256;
drawRectangle(rectInd,'c');
rectInd(:,1) = rectInd(:,1) + 256;
drawRectangle(rectInd,'r');
axis off; 
export_fig('Registered_maps/imgReg_sli.png','-m4');
close;

% calculate SSIM
t2_targ_BW = double(edge(t2_targ_sli,'Canny'));
t2_reg_BW = double(edge(t2_reg_sli,'Canny'));

ssimval = ssim(t2_targ_BW,t2_reg_BW);

end