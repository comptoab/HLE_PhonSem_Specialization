%% Second level analysis
%One sample t test and regression analysis written initally by Jin Wang 3/19/2019
%Edited for Home Literacy Environment Analyses - Anna Banaszkiewicz and Alisha B. Compton 

addpath(genpath('/gpfs51/dors2/gpc/JamesBooth/JBooth-Lab/BDL/AniaB/ELP_SemPhon_HLE/scripts/4secondlevel')); % the path of your scripts
spm_path='/gpfs51/dors2/gpc/JamesBooth/JBooth-Lab/BDL/LabCode/typical_data_analysis/spm12_elp'; %the path of spm
addpath(genpath(spm_path));

%define your data path
data=struct();
root='/gpfs51/dors2/gpc/JamesBooth/JBooth-Lab/BDL/AniaB/ELP_SemPhon_HLE/preprocessed';  %your project path
subjects={}; %manually put subject numbers say 'sub-5004' 'sub-5009', or leave this empty if you have an excel with data_info.
data_info='/gpfs51/dors2/gpc/JamesBooth/JBooth-Lab/BDL/AniaB/ELP_SemPhon_HLE/sesfive_subjects_N33_final.xlsx'; %your subject list excel, if you have regressors, it's better there too.
if isempty(subjects)
    M=readtable(data_info);
    subjects=M.subjects;
end

out_dir='/gpfs51/dors2/gpc/JamesBooth/JBooth-Lab/BDL/AniaB/ELP_SemPhon_HLE/secondlevel_results_ses-5/FINAL/noSES_regression_Phon_ChildAlone'; % your output folder of secondlevel analysis results

% find the data. If you follow the data structure and naming convention in
% the walkthrough, you won't need to change this part
global CCN;
CCN.session='ses-5';
CCN.func_pattern='sub*';
analysis_folder='analysis_ses-5';
model_deweight='deweight';

% choose your analysis method
test=2; %1 one-sample t test, 2 mutiple regression analysis

% if you have covariates in your second level modeling
cov=1; % 1 if you have covariates or 0 if you do not have covariates

% define your covariates of control for your one-sample t test. Or define your covariates of interest for your multiple regression
% analysis if you have covariates.
cov_num=3; %number of your covariates
%you can define as many covariates as you want by adding name and values.
if cov==1
    name=[];
    name{1}='ChildAlone'; %covariate of interest
    %name{2}='Mother_Edu'; % This should be your column header in your excel, should be exactly to be the same as in your excel, otherwise it won't read in these covariates.
    name{2}='Acc_PHON_mean';% This should be your column header in your excel
    name{3}='Age_M_mean';% This should be your column header in your excel
    val=[];
    for v=1:length(name)
        val{v}=M{:,name{v}};
    end
end

%%%%%%%%%%%%%%%%%%%%%%should not edit below %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize
%addpath(spm_path);
spm('defaults','fmri');
spm_jobman('initcfg');
spm_figure('Create','Graphics','Graphics');

% Dependency and sanity checks
if verLessThan('matlab','R2013a')
    error('Matlab version is %s but R2013a or higher is required',version)
end

req_spm_ver = 'SPM12 (6225)';
spm_ver = spm('version');
if ~strcmp( spm_ver,req_spm_ver )
    error('SPM version is %s but %s is required',spm_ver,req_spm_ver)
end

%Start to analyze the data from here

%load the contrast file path for each subject
scan=[];
contrast=[];
for i=1:length(subjects)
    deweight_spm=[root '/' subjects{i} '/' analysis_folder '/' model_deweight '/SPM.mat'];
    deweight_p=fileparts(deweight_spm);
    load(deweight_spm);
    contrast_names=[];
    scan_files=[];
    for ii=1:length(SPM.xCon)
        contrast_names{ii,1}=SPM.xCon(ii).name;
        scan_files{ii,1}=[deweight_p '/' SPM.xCon(ii).Vcon.fname];
    end
    contrast{i}=contrast_names;
    scan{i}=scan_files;
end

allscans=[];
for i=1:length(scan{1})
    for j=1:length(subjects)
        allscans{i}{j,1}=[scan{j}{i} ',1'];
    end
end

%make output folder for each contrast
if ~exist(out_dir)
    mkdir(out_dir);
end
cd(out_dir);
for ii=1:length(contrast{1})
    out_dirs{ii}=[out_dir '/' contrast{1}{ii}];
    if ~exist(out_dirs{ii})
        mkdir(out_dirs{ii});
    end
end

%covariates
%pass the covariates to a struct
if cov==1
    covariates.name=name;
    for i=1:cov_num
        values{i}=transpose(val{i});
    end
    covariates.values=values;
else
    covariates={};
end

if test==1 % one-sample t test
    onesample_t(out_dirs,allscans,covariates);
    
elseif test==2 %multiple regression analysis
    multiple_regression(out_dirs,allscans,covariates);
    
end



