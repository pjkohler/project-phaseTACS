clear all
close all
codeFolder = '/Users/kohler/code';
addpath(genpath([codeFolder,'/git/mrC']));
addpath(genpath([codeFolder,'/git/schlegel/matlab_lib']));
addpath(genpath([codeFolder,'/git/kay']));

dataFolder = '/Volumes/Denali_4D2/kohler/fMRI_EXP/phaseTACS/FIRST';
subDirs = subfolders(sprintf('%s/2017*',dataFolder),1);

overWrite(1) = false;
overWrite(2) = false;
createNiftis = false;
doRand = false;

includeMotion = false;
numConds = 4;
numReps = 5;
trDur = 2;
stimDur = 8;
subTractTR = 2;
subTractSecs = subTractTR * trDur;

%% get design matrix
lateStim = 0;
for s = 1:length(subDirs)
    RTfiles = subfiles(sprintf('%s/stim_data/Exp_MATL_PD1010_1_Cz/RT*',subDirs{s}),1);
    runCount = 0;
    for z=1:length(RTfiles)    
        tData = load(RTfiles{z});
        if length(tData.TimeLine) == (numConds * numReps)
            runCount = runCount + 1;
            tempTimes = cat(1,tData.TimeLine(:).segTimeSec);
            if max(tempTimes) > lateStim;
                % record latest stimulus
                lateStim = max(tempTimes);
            else
            end
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
                randTimes(:,runCount) = tempTimes;
                randConds(:,runCount) = tempConds;
            end
        else
        end
    end
    for z=1:size(randTimes,1)
        % shuffle times among runs
        randTimes(z,:) = randTimes(z,randperm(runCount));
        if z <= runCount
            % shuffle conditions within run
            randConds(:,z) = randConds(randperm(numReps*numConds),z);
        else
        end
    end
    for r = 1:runCount
        for c = 1:numConds
            glm_rand{s}{r}{c} = randTimes(randConds(:,r)==c,r);
        end
    end
    clear randTimes randConds;
end

maxTR = ceil(lateStim/trDur)+(stimDur/trDur)*2; % onset of latest block + blockduration (in TR) x 2

%% run glm
hrf = getcanonicalhrf(stimDur,trDur)';

for s = 1:length(subDirs)
    outName = sprintf('%s/glmDenoise',subDirs{s});
    figPath = sprintf('%s_figs',outName);
    matPath = sprintf('%s.mat',outName);
    if exist(matPath,'file')
        if ~overWrite(1)
            % if file exists, skip subject 
            fprintf('skipping: %s\n',subDirs{s});
            continue;
        else
            fprintf('running: %s\n',subDirs{s});
        end
    else
        fprintf('running: %s\n',subDirs{s});
    end
    if exist(figPath,'dir')
        % if figure directory exists, but not file, delete figure directory
        rmdir(figPath, 's');
    else
    end

    dataFiles = subfiles(sprintf('%s/run*vr.nii.gz',subDirs{s}),1);
    if length(dataFiles) ~= length(glm_design{s})
        error('number of data files different from stim files');
    else
    end
    motFiles = subfiles(sprintf('%s/motparam.run*vr.1D',subDirs{s}),1);
    if length(dataFiles) ~= length(motFiles)
        error('number of data files different from motion regressors');
    else
    end
    for z=1:length(dataFiles)
        glm_struct = NIfTI.Read(dataFiles{z});
        % subtract TRs
        numTR = size(glm_struct.data,4);
        glm_data{z}   = glm_struct.data(:,:,:,subTractTR+1:min(numTR,maxTR));
        mot_data{z} = load(motFiles{z});
        mot_data{z} = mot_data{z}(subTractTR+1:min(numTR,maxTR),:);
        if any(isnan(glm_data{z}(:)))
            error('glm data cannot have NaNs!');
        else
        end
        clear glm_struct;
    end
    mkdir(figPath);
    if includeMotion
        opt.extraregressors = mot_data;
    else
    end
    opt.wantparametric = true;
    [glm_results,denoised_data] = GLMdenoisedata(glm_design{s},glm_data,stimDur,trDur,'assume',hrf,opt,figPath);
    save(matPath,'-struct','glm_results','-v7.3');
    save(matPath,'-append','denoised_data');
    clear glm_data;
end

%% DO RANDOM GLM

if doRand;
    for s = 1
        outName = sprintf('%s/glmDenoiseRand',subDirs{s});
        figPath = sprintf('%s_figs',outName);
        matPath = sprintf('%s.mat',outName);
        dataFiles = subfiles(sprintf('%s/run*vr.nii.gz',subDirs{s}),1);
        for z=1:length(dataFiles)
            glm_struct = NIfTI.Read(dataFiles{z});
            % subtract TRs
            numTR = size(glm_struct.data,4);
            glm_data{z}   = glm_struct.data(:,:,:,subTractTR+1:min(numTR,maxTR));
            if any(isnan(glm_data{z}(:)))
                error('glm data cannot have NaNs!');
            else
            end
            clear glm_struct;
        end
        if exist(figPath,'dir')
            % if figure directory exists, but not file, delete figure directory
            rmdir(figPath, 's');
        else
        end
        mkdir(figPath)
        opt.wantparametric = true;
        [glm_results,denoised_data] = GLMdenoisedata(glm_rand{s},glm_data,stimDur,trDur,'assume',hrf,opt,figPath);
        save(matPath,'-struct','glm_results','-v7.3');
        save(matPath,'-append','denoised_data');
    end
else
end

%% INDIVIDUAL RUN GLMs
% For demonstration purposes, we will adopt the strategy of assuming a fixed HRF
% in the analysis of each run.  Here we compute a canonical HRF.

for s = 1:length(subDirs)
    dataFiles = subfiles(sprintf('%s/run*vr.nii.gz',subDirs{s}),1);    
    for z=1:length(dataFiles)
        outName = sprintf('glmIndividual_run%0.0f',z);
        outPath = sprintf('%s/%s',subDirs{s},outName);
        matPath = sprintf('%s.mat',outPath);
        if exist(matPath,'file')
            if ~overWrite(2)
                % if file exists, skip subject 
                fprintf('skipping glm run %0.0f: %s\n',z,subDirs{s});
                continue;
            else
                fprintf('running glm run %0.0f: %s\n',z,subDirs{s});
            end
        else
            fprintf('running glm run %0.0f: %s\n',z,subDirs{s});
        end
        glm_struct = NIfTI.Read(dataFiles{z});
        % subtract TRs
        numTR = size(glm_struct.data,4);
        glm_data = glm_struct.data(:,:,:,subTractTR+1:min(numTR,maxTR));
        if any(isnan(glm_data(:)))
            error('glm data cannot have NaNs!');
        else
        end
        clear glm_struct;
        glm_results = GLMestimatemodel(glm_design{s}(z),glm_data,stimDur,trDur,'assume',hrf,0);
        save(matPath,'-struct','glm_results','-v7.3');
        clear glm_data;
    end
end

%% SPLIT RUN GLMs
for s = 1:length(subDirs)
    fprintf('running splits %s\n',subDirs{s});
    dataFiles = subfiles(sprintf('%s/run*vr.nii.gz',subDirs{s}),1);    
    % figure out splits
    set1{1} = 1:2; %combnk(1:length(dataFiles),2);
    set2 = cellfun(@(x) find(~ismember(1:length(dataFiles),x)),set1,'uni',false)';
    for z=1:length(dataFiles)
        glm_struct = NIfTI.Read(dataFiles{z});
        % subtract TRs
        numTR = size(glm_struct.data,4);
        glm_data{z}   = glm_struct.data(:,:,:,subTractTR+1:min(numTR,maxTR));
        if any(isnan(glm_data{z}(:)))
            error('glm data cannot have NaNs!');
        else
        end
        clear glm_struct;
    end
    for z=1:length(set1)
        splits = {set1{z} set2{z}};
        for p = 1:length(splits)
            outName = sprintf('glmSplit_run%s',num2str(splits{p},'%0.0f'));
            outPath = sprintf('%s/%s',subDirs{s},outName);
            matPath = sprintf('%s.mat',outPath);
            figPath = sprintf('%s_figs',outPath);
             if exist(figPath,'dir')
            % if figure directory exists, but not file, delete figure directory
                rmdir(figPath, 's');
             else
            end
            mkdir(figPath);
            opt.wantparametric = true;
            ix = splits{p};
            [glm_results,denoised_data] = GLMdenoisedata(glm_design{s}(ix),glm_data(ix),stimDur,trDur,'assume',hrf,opt,figPath);
            save(matPath,'-struct','glm_results','-v7.3');
            save(matPath,'-append','denoised_data');
        end
    end
    clear glm_data;
end

%% CREATE NIFTIS
if createNiftis
    for s = 1:length(subDirs)
        fprintf('creating niftis %s\n',subDirs{s});
        dataFiles = subfiles(sprintf('%s/run*vr.nii.gz',subDirs{s}),1);
        glm_struct = NIfTI.Read(dataFiles{1});
        for z=1:(length(dataFiles)+1)
            outStruct.hdr = glm_struct.hdr;
            outStruct.method = glm_struct.method;
            if z > length(dataFiles) 
                % do all runs
                outName = 'glmDenoise';
            else
                % do individual runs
                outName = sprintf('glmIndividual_run%0.0f',z);
            end
            outPath = sprintf('%s/%s',subDirs{s},outName);
            niiPath = sprintf('%s.nii.gz',outPath);
            matPath = sprintf('%s.mat',outPath);
            glm_results = load(matPath,'model*','R2','SNR');

            betas = glm_results.modelmd{2};
            snr = betas./glm_results.modelse{2}; % divide by voxel-wise standard error to get SNR
            R2 = glm_results.R2./100; % variance explained, divide by 100
            tempData(:,:,:,1:2:size(betas,4)*2) = betas;
            tempData(:,:,:,2:2:size(betas,4)*2) = snr;
            outStruct.data = cat(4,tempData,R2,glm_results.SNR); % betas + beta SNR, R2 and overall SNR
            outStruct.hdr.dim(5) = size(outStruct.data,4);
            if exist(niiPath,'file')
                system(sprintf('rm %s',niiPath));
            else
            end
            NIfTI.Write(outStruct,niiPath);
            % assign new labels to sub-bricks
            afniCommand{1} = sprintf(...
                'source ~/.bashrc; cd %s; 3drefit -fbuc -redo_bstat -relabel_all_str ''beta0 snr0 beta90 snr90 beta180 snr180 beta270 snr270 modelR2 modelSNR'' %s',subDirs{s}, niiPath);
            % go from volume to surface
            tempSub = split(subDirs{s},'_');
            subID = tempSub{end};
            if strcmp(subID,'nl-0034');
                subID = 'skeri0004';
            else
            end
            afniCommand{2} = sprintf(...
                'source ~/.bashrc; cd %s; rm *h.%s.niml.dset; mriVol2Surf.py %s %s.nii.gz --surfvol %s_fs4_SurfVol_ns_Alnd_Exp+orig',subDirs{s},outName,subID,outName,subID);
            status = system(afniCommand{1});
            if status == 127
                error('afni 3drefit command failed!');
            else
            end
            if z > length(dataFiles) 
                status = system(afniCommand{2});
                if status == 127
                    error('afni mriVol2Surf command failed!');
                else
                end
            else
            end
            clear tempData;
            clear outStruct;
        end
    end
else
    fprintf('skipping nifti creation \n');
end

%% DO STATS
for s = 1:length(subDirs)
    outName = sprintf('%s/glmDenoise',subDirs{s});
    matPath = sprintf('%s.mat',outName);
    glm_results = load(matPath,'model*','R2*','SNR','parametric');
    snrIdx = find(glm_results.SNR>2);
    betas = glm_results.modelmd{2};
    snr = betas./glm_results.modelse{2}; 
    for z=1:4; 
        temp = snr(:,:,:,z); 
        temp = temp(snrIdx);
        meanSnr(z,s) = mean(temp); 
        temp = betas(:,:,:,z); 
        temp = temp(snrIdx);
        %temp = temp(temp>0);
        meanBeta(z,s) = mean(temp); 
    end
    designOut{s} = glm_results.parametric.designmatrix;
    % R2 stuff
    runColor = {'k','b','r','g','ma','cy'};
    numRuns = size(glm_results.R2run,4);
    figure;
    for z = 1:(numRuns+1)
        if z>numRuns
            R2vals = glm_results.R2;
        else
            R2vals = glm_results.R2run(:,:,:,z);
        end
        n = histc(R2vals(:),0:0.5:5); 
        subplot(numRuns+1,1,z);
        hold on;
        barH = bar(0:0.5:5,n,'histc');
        set(barH,'facecolor','none','edgecolor',runColor{z});
        ylim([0,3*10^5]);
        hold off
    end
    set(gcf,'Units','centimeters');
    curPos = get(gcf,'position');
    curPos(3) = 10;
    curPos(4) = numRuns * 10;
    set(gcf,'position',curPos);
    clear glm_results;
end

allMean = mean(meanBeta,2);
allStd = std(meanBeta,0,2)./size(meanBeta,2);
plot(1:4,allMean);
hold on;
errorb(1:4,allMean,allStd);

%% plot design matrices
figure;
for s = 1:length(subDirs)
    designImg = designOut{s}(:,1:4);
    subplot(2,5,s);
    tempSub = split(subDirs{s},'_');
    subID = tempSub{end};
    if strcmp(subID,'nl-0034');
        subID = 'skeri0004';
    else
    end
    set(gcf,'Units','points','Position',[100 100 350 500]);
    imagesc(designImg);
    colormap(gray);
    colorbar;
    xlabel('Conditions');
    ylabel('Time points');
    title(sprintf('Design matrix %s \n (%d time points, %d conditions)',subID,size(designImg,1),size(designImg,2)));
end


