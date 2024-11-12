%Count Repaired
%%This script counts the number of repaired volumes in each run.
%%Written by Jin Wang 5/20/2019 
%%It will print out the name of the run (the first column) the number of volumes being
%%replaced (the second column) and how many chunks of more than 6 consecutive volumes being
%%replaced (the third column) for each subject.

%Orginially written by Jin Wang 5/20/2019 and modified by Neelima Wagley 12/12/22
%Edited for Home Literacy Environment Analyses - Anna Banaszkiewicz and Alisha B. Compton 

%% Define parameters

% path to script folder
addpath(genpath('/gpfs51/dors2/gpc/JamesBooth/JBooth-Lab/BDL/Alisha/ELP_SemPhon_HLE/scripts/2preprocessing/'));

% list of participants
subjects={}; %manually put subject numbers say 'sub-5004' 'sub-5009', or leave this empty if you have an excel with data_info.

%or define with xls file, one col with header 'subjects'
data_info='/gpfs51/dors2/gpc/JamesBooth/JBooth-Lab/BDL/Alisha/sesfive_subjects_sem.xlsx'; %In this excel, there should be a column of subjects with the header (subjects). The subjects should all be sub plus numbers (sub-5002).

if isempty(subjects)
   M=readtable(data_info);
   subjects=M.subjects;
end

% path of preprocessed data files
root='/gpfs51/dors2/gpc/JamesBooth/JBooth-Lab/BDL/Alisha/ELP_preproc/';

global CCN;
CCN.func_n='*task*'; %define the functional data
CCN.ses='ses-5';

n=6; %number of consecutive volumes being replaced

%% %%%%%%%%%%%%%%should not edit below%%%%%%%%%%%%%%%%%%%%%%
    cd(root);
    writefile='count_repaired_ses5_ACsem.txt';
    fid=fopen(writefile, 'w');
    hdr= 'subjects run_name num_repaird chunks';
    fprintf(fid,'%s', hdr);
    fprintf(fid, '\n');
for i=1:length(subjects)
    func_p=[root '/' subjects{i}];
    func_f=expand_path([func_p '/[ses]/func/[func_n]/']);
    for j=1:length(func_f)
        run_n=func_f{j}(1:end-1);
        [run_p, run_name]=fileparts(run_n);
        cd(run_n);
        fileid=fopen('art_repaired.txt');
        m=fscanf(fileid, '%f');
        [num_repaired, col]=size(m);     
        x=diff(m')==1;
        ii=strfind([0 x 0], [0 1]);
        jj=strfind([0 x 0], [1 0]);
        idx=max(jj-ii);
        out=(idx>=n);
        if out==0
            chunks=0;
        else
            [chunks, col]=size(out);
        end
        fprintf(fid, '%s %s %d %d \n', subjects{i}, run_name, num_repaired, chunks);
    end
end

 cd(root);
