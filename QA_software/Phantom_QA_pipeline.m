clc; clear;

pathFolder = '/Users/yaojingwen/Box/300_SDP_2020_share_folder/CEST_phantom_QA';
% pathFolder = '/Users/Jingwen/Box/300_SDP_2020_share_folder/CEST_phantom_QA';
addpath(genpath([pathFolder '/QA_software/library']));

routineFLAG = 0; % calibration scan or routine QA scan
% ROIind = [3 4 5 6 7 8 1 2 9 12 13 14 15 16 17 10 11 18];

%% convert DICOMs to NIFTIs
pathToData = [pathFolder '/../Phantom_scan_data/Cest_Phantom_SecondBatch'];
standard_dicom2nifti(pathToData);
NIFTI_folder = [pathToData '/NIFTI'];
DICOM_folder = [pathToData '/DICOM'];
mkdir([pathToData '/Report']);
Report_folder = [pathToData '/Report'];
copyfile([pathFolder '/QA_software/Report/*'],Report_folder);

%% calculate T1/T2/T2star maps
cd(NIFTI_folder);

if ~isfile('T2map.nii.gz')
    TE = T2_mapping('MSE.nii.gz','T2');
end

if ~isfile('T2starmap.nii.gz')
    TEstar = T2_mapping('MGE.nii.gz','T2star');
end

if ~isfile('T1map.nii.gz')
    TI = T1_mapping('T1mapping');
end

%% calculate MTRasym map
% offsets may change
offsets = [-10 -7 -5 -4.5 -4:0.2:-2 -1.75:0.25:1.75 2:0.2:4 4.5 5 7 10];
process_cest_arbOffsets('me_cest', offsets);

%% register images to t1 (best resolution)
cd(NIFTI_folder);

% t2
autoreg_tk('t1.nii.gz', 't2.nii.gz', 't2.reg.nii.gz', 1, 't22t1.dat');

% T2map
autoreg_tk('t1.nii.gz', 'MSE.nii.gz', 'temp.nii.gz', 1, 'mse2t1.dat');
autoreg_tk('t1.nii.gz', 'T2map.nii.gz', 'T2map.reg.nii.gz', 0, 'mse2t1.dat');
system('rm temp.nii.gz');

% T2starmap
autoreg_tk('t1.nii.gz', 'MGE.nii.gz', 'temp.nii.gz', 1, 'mge2t1.dat');
autoreg_tk('t1.nii.gz', 'T2starmap.nii.gz', 'T2starmap.reg.nii.gz', 0, 'mge2t1.dat');
system('rm temp.nii.gz');

% T1map
autoreg_tk('t1.nii.gz', 'T1map.nii.gz', 'T1map.reg.nii.gz', 1, 'ti2t1.dat');

% CESTmap
autoreg_tk('t1.nii.gz', 'me_cest_s0_TE1.nii.gz', 'temp.nii.gz', 1, 'cest2t1.dat');
autoreg_tk('t1.nii.gz', 'MTRasym_IDE.nii.gz', 'MTRasym.reg.nii.gz', 0, 'cest2t1.dat');
autoreg_tk('t1.nii.gz', 'B0_IDE.nii.gz', 'B0.reg.nii.gz', 0, 'cest2t1.dat');

%% use T2 as template to create vial VOIs

% createROItemplate('../NIFTI_template');

if ~isfolder([pathToData '/NIFTI_template'])
    copyfile([pathFolder '/QA_software/NIFTI_template'],[pathToData '/NIFTI_template']);
end

%% register images to template
cd(NIFTI_folder);
mkdir('Registered_maps');

if ~isfile('Registered_maps/t1.nii.gz')
    % register T1
    copyfile('../NIFTI_template/t1_template.nii.gz','Registered_maps/t1_template.nii.gz');
    copyfile('../NIFTI_template/t2_template.nii.gz','Registered_maps/t2_template.nii.gz');
    cmd = ['flirt -in t1.nii.gz -ref Registered_maps/t1_template.nii.gz'...
        ' -out Registered_maps/t1.nii.gz -omat Registered_maps/regT1.mat'...
        ' -cost normcorr -searchrx -180 180 -searchry -180 180 -searchrz -180 180'];
    system(cmd);
    
    % register others with the same transformation matrix
    cmd = ['flirt -in t2.reg.nii.gz -ref ../NIFTI_template/t1_template.nii.gz'...
        ' -out Registered_maps/t2.nii.gz -applyxfm -init Registered_maps/regT1.mat;'...
        'flirt -in T2map.reg.nii.gz -ref ../NIFTI_template/t1_template.nii.gz'...
        ' -out Registered_maps/T2map.nii.gz -applyxfm -init Registered_maps/regT1.mat;'...
        'flirt -in T2starmap.reg.nii.gz -ref ../NIFTI_template/t1_template.nii.gz'...
        ' -out Registered_maps/T2starmap.nii.gz -applyxfm -init Registered_maps/regT1.mat;'...
        'flirt -in T1map.reg.nii.gz -ref ../NIFTI_template/t1_template.nii.gz'...
        ' -out Registered_maps/T1map.nii.gz -applyxfm -init Registered_maps/regT1.mat;'...
        'flirt -in MTRasym.reg.nii.gz -ref ../NIFTI_template/t1_template.nii.gz'...
        ' -out Registered_maps/MTRasym.nii.gz -applyxfm -init Registered_maps/regT1.mat;'...
        'flirt -in B0.reg.nii.gz -ref ../NIFTI_template/t1_template.nii.gz'...
        ' -out Registered_maps/B0.nii.gz -applyxfm -init Registered_maps/regT1.mat;'];
    system(cmd);
    
    cmd = ['3dcalc -a Registered_maps/T1map.nii.gz -expr 1000/a -prefix '...
        'Registered_maps/R1map.nii.gz;'...
        '3dcalc -a Registered_maps/T2map.nii.gz -expr 1000/a -prefix '...
        'Registered_maps/R2map.nii.gz'];
    system(cmd);
    
    copyfile('../NIFTI_template/ROI.nii.gz','Registered_maps/ROI.nii.gz');
    
end

% create screenshots for reports
nii = load_untouch_nii('Registered_maps/t2.nii.gz');
t2_reg = nii.img;
nii = load_untouch_nii('../NIFTI_template/t2_template.nii.gz');
t2_targ = nii.img;

lev1 = 16; lev2 = 84; lev3 = 52;
ssimval = plotRegImg(t2_reg,t2_targ,lev1,lev2,lev3);

% copy registration images to Report
copyfile([NIFTI_folder '/Registered_maps/imgReg_mid.png'],[Report_folder '/imgReg_mid.png']);
copyfile([NIFTI_folder '/Registered_maps/imgReg_sli.png'],[Report_folder '/imgReg_sli.png']);

%% extract map values for each vial
cd([NIFTI_folder '/Registered_maps']);

T2_stats = ROIstats('ROI.nii.gz','T2map.nii.gz');
% T2star_stats = ROIstats('ROI.nii.gz','T2starmap.nii.gz');
T1_stats = ROIstats('ROI.nii.gz','T1map.nii.gz');
R2_stats = ROIstats('ROI.nii.gz','R2map.nii.gz');
R1_stats = ROIstats('ROI.nii.gz','R1map.nii.gz');
MTRasym_stats = ROIstats('ROI.nii.gz','MTRasym.nii.gz');
B0_stats = ROIstats('ROI.nii.gz','B0.nii.gz');

if exist('ROIind','var')
    T2_stats = T2_stats(ROIind,:);
    T1_stats = T1_stats(ROIind,:);
    R2_stats = R2_stats(ROIind,:);
    R1_stats = R1_stats(ROIind,:);
    MTRasym_stats = MTRasym_stats(ROIind,:);
    B0_stats = B0_stats(ROIind,:);
end

%% read excel of phantom preparation
xls = [pathFolder '/Phantom_solutions.xlsx'];
[num,~,~] = xlsread(xls);

pH = num(:,2);
Pconc = num(:,3); % mM
GDconc = num(:,4); % mM
GLYconc = num(:,5); % mM
GDflag = num(:,6); % true or false
GLYflag = num(:,7); % true or false

% standard values
copyfile([pathFolder '/QA_software/Report/Result_summary.xlsx'],...
    [Report_folder '/Result_summary.xlsx']);
xls = [Report_folder '/Result_summary.xlsx'];
[num,~,~] = xlsread(xls,2);

R1AVE = num(:,6);
R2AVE = num(:,7);
CESTAVE = num(:,8);
R1STD = num(:,9);
R2STD = num(:,10);
CESTSTD = num(:,11);

% z scores of deviations from standard value
zT1 = (R1_stats(:,1) - R1AVE)./R1STD;
zT2 = (R2_stats(:,1) - R2AVE)./R2STD;
zCEST = (MTRasym_stats(:,1) - CESTAVE)./CESTSTD;

zcutoff = icdf('T',0.99,5);
zT1_Nflag = sum(zT1 > zcutoff) + sum(zT1 < -zcutoff);
zT2_Nflag = sum(zT2 > zcutoff) + sum(zT2 < -zcutoff);
zCEST_Nflag = sum(zCEST > zcutoff) + sum(zCEST < -zcutoff);
z_total = length(zT1);

%% create table of metrics - write to excel

vialIndex = cell(length(pH),1);
R1meanSD = cell(length(pH),1);
R2meanSD = cell(length(pH),1);
CESTmeanSD = cell(length(pH),1);
for i = 1:length(pH)
    vialIndex{i} = ['#' num2str(i,'%d')];
    R1meanSD{i} = [num2str(R1_stats(i,1),'%.2f') ' ± ' ...
        num2str(R1_stats(i,2),'%.2f')];
    R2meanSD{i} = [num2str(R2_stats(i,1),'%.2f') ' ± ' ...
        num2str(R2_stats(i,2),'%.2f')];
    CESTmeanSD{i} = [num2str(MTRasym_stats(i,1),'%.2f') ' ± ' ...
        num2str(MTRasym_stats(i,2),'%.2f')];
end

T = table(vialIndex,[num2str(pH,'%.2f')],[num2str(GDconc,'%.2f')],...
    [num2str(Pconc,'%.2f')],[num2str(GLYconc,'%.2f')],...
    R1meanSD,R2meanSD,CESTmeanSD,...
    T1_stats(:,1),T2_stats(:,1),MTRasym_stats(:,1),...
    T1_stats(:,2),T2_stats(:,2),MTRasym_stats(:,2));
T.Properties.VariableNames = ...
    {'Sample #','pH','Gd (mM)','Phosphate (mM)',...
    'Glycine (mM)','R1 (s-1)','R2 (s-1)','MTRasym (%)',...
    'T1mean','T2mean','CESTmean','T1std','T2std','CESTstd'};
Tplot = [{'Sample #','pH','Gd (mM)','Phosphate (mM)',...
    'Glycine (mM)','$R_1$ ($s^{-1}$)','$R_2$ ($s^{-1}$)','$MTR_{asym}$ (%)'};...
    table2cell(T(:,1:8))];

xlsFile = [Report_folder '/Result_summary.xlsx'];
writetable(T,xlsFile,'Sheet','Output','Range','A3','WriteVariableNames',false);

%% create colored table for latex

[Tplot_color] = addColorToCellArray(zT1, 6, Tplot, zcutoff);
[Tplot_color] = addColorToCellArray(zT2, 7, Tplot_color, zcutoff);
[Tplot_color] = addColorToCellArray(zCEST, 8, Tplot_color, zcutoff);

input.tableBorders = 1;
input.booktabs = 1;
[latex] = cellArray2latex(Tplot_color, input, 'c');

% save as txt
fid = fopen( [Report_folder '/MeasurementTable.txt'], 'wt');
fprintf(fid, '%s\n',char(latex)');
fclose(fid);

Tplot_routine = [{'Sample #','pH','Glycine (mM)','$MTR_{asym}$ (%)'};...
    table2cell(T(:,[1 2 5 8]))];
Tplot_routine_color = addColorToCellArray(zCEST, 4, Tplot_routine, zcutoff);
[latex] = cellArray2latex(Tplot_routine_color, input, 'c');

% save as txt
fid = fopen( [Report_folder '/MeasurementTable_routine.txt'], 'wt');
fprintf(fid, '%s\n',char(latex)');
fclose(fid);

%% create figures
IndSlice = [26 98];

nii = load_untouch_nii([NIFTI_folder '/Registered_maps/t2.nii.gz']);
T2img = int16(nii.img(:,:,IndSlice));

nii = load_untouch_nii([NIFTI_folder '/Registered_maps/ROI.nii.gz']);
ROImask = nii.img(:,:,IndSlice) > 0;
ROI_img = int16(nii.img(:,:,IndSlice));

T1mapFile = [NIFTI_folder '/Registered_maps/T1map.nii.gz'];
label_str = 'T1 (ms)';
plotMaps(T2img, T1mapFile, IndSlice, ROImask, ROI_img, [0 1000], 'w', label_str);
export_fig([Report_folder '/T1map.png'],'-m2','-transparent');
close;

T2mapFile = [NIFTI_folder '/Registered_maps/T2map.nii.gz'];
label_str = 'T2 (ms)';
plotMaps(T2img, T2mapFile, IndSlice, ROImask, ROI_img, [0 500], 'w', label_str);
export_fig([Report_folder '/T2map.png'],'-m2','-transparent');
close;

% T2starmapFile = [NIFTI_folder '/Registered_maps/T2starmap.nii.gz'];
% plotMaps(T2img, T2starmapFile, IndSlice, ROImask, ROI_img, [0 300]);

MTRmapFile = [NIFTI_folder '/Registered_maps/MTRasym.nii.gz'];
label_str = 'MTRasym at 3.0 ppm (%)';
plotMaps(T2img, MTRmapFile, IndSlice, ROImask, ROI_img, [-5 10], 'w', label_str);
export_fig([Report_folder '/MTRasymmap.png'],'-m2','-transparent');
close;

plotMaps(T2img, MTRmapFile, IndSlice, ROImask, ROI_img, [-5 10], 'w', label_str, 1);
export_fig([Report_folder '/MTRasymmap_vert.png'],'-m2','-transparent');
close;

%% create plots - relaxation
set(0,'DefaultAxesFontSize',16);
set(0,'DefaultLineLineWidth',2);

GDflag = GDflag > 0;

% T1
axisInfo.xlim = [-0.1 1.1];
axisInfo.ylabel = 'T_1 (s)';
axisInfo.ylabel2 = 'R_1 (s^{-1})';
axisInfo.imgFile = 'T1plot.png';
plotRelaxation(GDconc,GDflag,T1_stats,axisInfo,Report_folder);

% T2
axisInfo.xlim = [-0.1 1.1];
axisInfo.ylabel = 'T_2 (s)';
axisInfo.ylabel2 = 'R_2 (s^{-1})';
axisInfo.imgFile = 'T2plot.png';
plotRelaxation(GDconc,GDflag,T2_stats,axisInfo,Report_folder);

%% create plots - MTRasym

% pH dependency
GLYflag = GLYflag > 0;
GLYarray = [10 20 50];
axisInfo.xlim = [5.9 7.5];
axisInfo.xlimGD = [-0.1 1.1];
axisInfo.legend = {'GLY 10mM','GLY 20mM','GLY 50mM'};

plotCEST(pH,GLYconc,GLYflag,GDconc,GDflag,MTRasym_stats,GLYarray,...
    axisInfo,Report_folder)

%% extract scan information for report
[ScanInfo] = extractScanInfo(DICOM_folder);

% File for report
fid = fopen( [Report_folder '/ScanInformation.txt'], 'wt');
fprintf(fid, '%s\n', char(ScanInfo.txt)');
fclose(fid);

%% extract sequence information for report
cd(NIFTI_folder);

% t1w
[SeqInfo] = extractSeqInfo(DICOM_folder,'t1');
fid = fopen( [Report_folder '/T1w.txt'], 'wt');
fprintf(fid, '%s\n', char(SeqInfo.txt)');
fclose(fid);

% t2w
[SeqInfo] = extractSeqInfo(DICOM_folder,'t2');
fid = fopen( [Report_folder '/T2w.txt'], 'wt');
fprintf(fid, '%s\n', char(SeqInfo.txt)');
fclose(fid);

% t2map
[SeqInfo] = extractSeqInfo(DICOM_folder,'t2map');
fid = fopen( [Report_folder '/T2map.txt'], 'wt');
fprintf(fid, '%s\n', char(SeqInfo.txt)');
fclose(fid);

% t1map
[SeqInfo] = extractSeqInfo(DICOM_folder,'t1map');
fid = fopen( [Report_folder '/T1map.txt'], 'wt');
fprintf(fid, '%s\n', char(SeqInfo.txt)');
fclose(fid);

% cest
[SeqInfo] = extractSeqInfo(DICOM_folder,'cest');
fid = fopen( [Report_folder '/CEST.txt'], 'wt');
fprintf(fid, '%s\n', char(SeqInfo.txt)');
fclose(fid);

%% quality check
latexCheckBox = '$\text{\rlap{$\checkmark$}}\square$';
latexNoBox = '$\boxtimes$';

input.tableBorders = 0;
input.booktabs = 0;

QAtxt = [];
QAtxt_routine = [];
% registration
if ssimval >= 0.75
    check_temp = latexCheckBox;
    QA_temp = ['SSIM score = ' num2str(ssimval,'%.3f') ' $(>= 0.75)$'];
else
    check_temp = latexNoBox;
    QA_temp = ['SSIM score = ' num2str(ssimval,'%.3f') ' $(< 0.75)$. Please check'];
end
QAtxt = [QAtxt;
    [{'Registration to template'} {check_temp} {QA_temp}]];
QAtxt_routine = [QAtxt_routine;
    [{'Registration to template'} {check_temp} {QA_temp}]];

% relaxation rate measurements
if zT1_Nflag == 0
    check_temp = latexCheckBox;
    QA_temp = 'All within 99% confidence interval of standard values.';
else
    check_temp = latexNoBox;
    QA_temp = [num2str(zT1_Nflag,'%d') ' out of ' num2str(z_total,'%d')...
        ' outside the 99% confidence interval of standard values. Please check.'];
end
QAtxt = [QAtxt;
    [{'$T_1$ relaxation measurements'} {check_temp} {QA_temp}]];

if zT2_Nflag == 0
    check_temp = latexCheckBox;
    QA_temp = 'All within 99% confidence interval of standard values.';
else
    check_temp = latexNoBox;
    QA_temp = [num2str(zT2_Nflag,'%d') ' out of ' num2str(z_total,'%d')...
        ' outside the 99% confidence interval of standard values. Please check.'];
end
QAtxt = [QAtxt;
    [{'$T_2$ relaxation measurements'} {check_temp} {QA_temp}]];

% CEST measurements
if zCEST_Nflag == 0
    check_temp = latexCheckBox;
    QA_temp = 'All within 99% confidence interval of standard values.';
else
    check_temp = latexNoBox;
    QA_temp = [num2str(zCEST_Nflag,'%d') ' out of ' num2str(z_total,'%d')...
        ' outside the 99% confidence interval of standard values. Please check.'];
end
QAtxt = [QAtxt;
    [{'CEST measurements ($MTR_{asym}$)'} {check_temp} {QA_temp}]];
QAtxt_routine = [QAtxt_routine;
    [{'CEST measurements ($MTR_{asym}$)'} {check_temp} {QA_temp}]];

fixWidth = '{p{0.3\textwidth}p{0.05\textwidth}p{0.6\textwidth}}';
[latex] = cellArray2latex(QAtxt, input, 'l', fixWidth);
[latex_routine] = cellArray2latex(QAtxt_routine, input, 'l', fixWidth);

% save as txt
fid = fopen( [Report_folder '/QA.txt'], 'wt');
fprintf(fid, '%s\n',char(latex)');
fclose(fid);

% save as txt
fid = fopen( [Report_folder '/QA_routine.txt'], 'wt');
fprintf(fid, '%s\n',char(latex_routine)');
fclose(fid);

%% generate report
copyfile([pathFolder '/QA_software/Report/*'],Report_folder);

date = ScanInfo.Date;
date = strrep(date,' ','_');

% if ~routineFLAG
cd(Report_folder);
cmd = 'pdflatex Report.tex;';
system(cmd);

cmd = ['mv Report.pdf Report_' date '.pdf']; system(cmd);
copyfile(['Report_' date '.pdf'],[pathToData '/Report_' date '.pdf']);
cd(pathToData);
% else
cd(Report_folder);
cmd = 'pdflatex Report_Routine.tex;';
system(cmd);

cmd = ['mv Report_Routine.pdf Report_routine_' date '.pdf']; system(cmd);
copyfile(['Report_routine_' date '.pdf'],[pathToData '/Report_routine_' date '.pdf']);
cd(pathToData);
% end