function acuteNormalize
% This script normalizes FLAIR scans and lesions
% Adjusted for the specified folder structure

baseDir = 'C:\Users\lrm22005\OneDrive - University of Connecticut\Research\Stroke_MRI_data\Stroke_MRI_Segmentation_Project\data\raw';
subjectDirs = dir(fullfile(baseDir, 'sub-*')); % Adjust to list subject directories
outdir = 'R:\hfp14002\PosadaLab\Luis\Research\StrokeOutcomeOptimizationProjectDemo\norm';

for i = 1 : length(subjectDirs)
    subjectDir = fullfile(subjectDirs(i).folder, subjectDirs(i).name);
    anatDir = fullfile(subjectDir, 'anat');
    dwiDir = fullfile(subjectDir, 'dwi');
    fprintf('%d/%d %s\n', i, length(subjectDirs), subjectDir);

    flairFiles = dir(fullfile(anatDir, '*_FLAIR.nii.gz'));
    
    for j = 1:length(flairFiles)
        fl = fullfile(flairFiles(j).folder, flairFiles(j).name);
        [p,n,x] = fileparts(fl);
        subj = split(n,'_');
        subj = subj{1};
        
        lesMaskDir = fullfile(baseDir, 'derivatives', 'lesion_masks', subjectDirs(i).name, 'dwi');
        lesMasks = dir(fullfile(lesMaskDir, '*_desc-lesion_mask.nii.gz'));
        
        if ~isempty(lesMasks)
            les = fullfile(lesMasks(1).folder, lesMasks(1).name);
            tr = fullfile(dwiDir, [subj, '_rec-TRACE_dwi.nii.gz']);
            if ~exist(tr, 'file')
                error('Unable to find %s\n', tr);
            end
            
            les2 = fullfile(outdir, [subj, '_lesion.nii.gz']);
            tr2 = fullfile(outdir, [subj, '_TRACE_dwi.nii.gz']);
            copyfile(les, les2);
            copyfile(tr, tr2);
        else
            les2 = [];
            tr2 = [];
        end

        fl2 = fullfile(outdir, [subj, '_FLAIR.nii.gz']); 
        copyfile(fl, fl2);
        normalizeFlair(fl2, les2, tr2);
        deletenii(fl2);
        deletenii(les2);
        deletenii(tr2);
    end
end
%end acuteNormalize()

function deletenii(fnm)
if ~isempty(fnm) && exist(fnm, 'file')
    delete(fnm)
end
%end deletenii()

function normalizeFlair(fl, les, tr)
% Attempt to address file path issues and ensure .nii format

% Check if the files are in .nii.gz format and unzip if necessary
fl_unzipped = checkAndUnzip(fl);
les_unzipped = checkAndUnzip(les);
tr_unzipped = checkAndUnzip(tr);

% Assuming matlabbatch structure is correct for your SPM version
matlabbatch{1}.spm.tools.MRI.MRnorm.anat = {fl_unzipped};
matlabbatch{1}.spm.tools.MRI.MRnorm.les = {les_unzipped};
matlabbatch{1}.spm.tools.MRI.MRnorm.t2 = {tr_unzipped};
matlabbatch{1}.spm.tools.MRI.MRnorm.modality = 3;
matlabbatch{1}.spm.tools.MRI.MRnorm.brainmask = 0;
matlabbatch{1}.spm.tools.MRI.MRnorm.bb = [NaN NaN NaN; NaN NaN NaN];
matlabbatch{1}.spm.tools.MRI.MRnorm.vox = [1 1 1];
matlabbatch{1}.spm.tools.MRI.MRnorm.DelIntermediate = 1;
matlabbatch{1}.spm.tools.MRI.MRnorm.AutoSetOrigin = 0;

spm_jobman('initcfg'); % Initialize the job configuration
spm_jobman('run',matlabbatch);

% Cleanup unzipped files if originally zipped
cleanupUnzipped(fl_unzipped, fl);
cleanupUnzipped(les_unzipped, les);
cleanupUnzipped(tr_unzipped, tr);
%end normalizeFlair()

function outPath = checkAndUnzip(inPath)
% Check if the input path is a .nii.gz file and unzip it if necessary
if endsWith(inPath, '.nii.gz')
    % Correctly remove '.gz' to avoid duplicated extension
    outPath = inPath(1:end-3); % Removes the '.gz', keeping '.nii'
    if ~exist(outPath, 'file')
        gunzip(inPath);
    end
else
    outPath = inPath;
end
if ~exist(outPath, 'file')
    error(['File does not exist: ', outPath]);
end
%end checkAndUnzip()

function cleanupUnzipped(unzippedPath, originalPath)
% Cleanup unzipped files if they were originally .nii.gz
if endsWith(originalPath, '.nii.gz') && exist(unzippedPath, 'file')
    delete(unzippedPath);
end
%end cleanupUnzipped()
