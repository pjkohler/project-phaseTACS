clear all
close all

dataFolder = '/Volumes/Denali_4D2/kohler/fMRI_EXP/phaseTACS/FIRST';
subDirs = subfolders(sprintf('%s/2017*',dataFolder),1);
subDirs = subDirs(1:end-1);
numConds = 4;

trDur = 2;
stimDur = 8;
subTractTR = 2;
subTractSecs = subTractTR * trDur;

for s = 1:length(subDirs)
    RTfiles = subfiles(sprintf('%s/stim_data/Exp_MATL_PD1010_1_Cz/RT*',subDirs{s}),1);
    runCount = 0;
    for z=1:length(RTfiles)    
        tData = load(RTfiles{z});
        if ~isempty(tData.TimeLine)
            runCount = runCount + 1;
            tempTimes = cat(1,tData.TimeLine(:).segTimeSec);
            if min(tempTimes) < subTractSecs
                error('negative timing values after subtraction');
            else
                tempTimes = tempTimes - subTractSecs;
            end
            tempConds = cat(1,tData.TimeLine(:).cndNmb);
            if ~isequal(unique(tempConds)', 1:numConds)
                error('unexpected number of conditions');
            else
            end
            for c = 1:numConds
                glm_design{s}{runCount}{c} = tempTimes(tempConds==c);
            end
        else
        end
    end
    dataFiles = subfiles(sprintf('%s/run*sc.dt.nii.gz',subDirs{s}),1);
    if length(dataFiles) ~= length(glm_design{s})
        error('number of data files different from stim files');
    else
    end
    for z=1:length(dataFiles)
        glm_struct{s}{z} = NIfTI.Read(dataFiles{z});
        % subtract TRs
        glm_data{s}{z}   = glm_struct{s}{z}.data(:,:,:,subTractTR+1:end);
        if any(isnan(glm_data{s}{s}(z)))
            error('glm data cannot have NaNs!');
        else
        end
        % remove data field from struct, but save struct for later
        glm_struct{s}{z} = rmfield(glm_struct{s}{z},'data');
    end
    [glm_results{s},glm_denoiseddata{s}] = GLMdenoisedata(glm_design{s},glm_data{s},stimDur,trDur,[],[],[],'figures');
    outStruct = glm_struct{s}{z};
    outStruct.data = cat(4,glm_results{s}.modelmd{2},glm_results{s}.SNR); % get betas and SNR
    outName = sprintf('%s/glmDenoise_out.nii.gz',subDirs{s});
    NIfTI.Write(outStruct,outName);
end

