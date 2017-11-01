clear all;
close all;
%% ADD PATHS
codeFolder = '/Users/kohler/code';
addpath(genpath([codeFolder,'/git/schlegel/matlab_lib']));
addpath(genpath([codeFolder,'/git/mrC']));
addpath(genpath([codeFolder,'/git/MRI/matlab']));

wangROIs = {'V1','V2' 'V3','hV4' 'VO1' 'VO2' 'PHC1' 'PHC2' ...
            'TO1' 'LO2' 'LO1' 'V3B' 'V3A' 'IPS0' 'IPS1' 'IPS2' 'IPS3' 'SPL1' };
        % note truncated version with only 18 ROIs

TOP_DIR='/Volumes/Denali_4D2/kohler/fMRI_EXP/phaseTACS/FIRST';

SUB_DIR = subfolders(sprintf('%s/2017*',TOP_DIR),1);
sesNames = subfolders(sprintf('%s/2017*',TOP_DIR));
sesNames{1} = '20170517_skeri0004';
subNames = cellfun(@(x) split_string(x,'_',2), sesNames,'uni',false);

montageObj = containers.Map({'20170517_skeri0004','20170720_nl-0037','20170829_nl-0032'},...
                            repmat({'P07-OZ-P08'},1,3));
montageObj = [montageObj;containers.Map({'20170731_nl-0032','20170801_nl-0037','20170801_skeri0004','20170818_nl-0011'},...
                            repmat({'CZ-OZ'},1,4))];
montageObj = [montageObj;containers.Map({'20170803_skeri0004','20170817_nl-0037','20170822_nl-0032'},...
                            repmat({'T7-OZ-T8'},1,3))];                           

roiLabel = [cellfun(@(x) [x,'-L'],wangROIs,'uni',false),cellfun(@(x) [x,'-R'],wangROIs,'uni',false)];

for s=1:length(SUB_DIR);
    disp(['Reading ',SUB_DIR{s},' ...']);
    curROI = subfiles(sprintf('%s/bl.wang_atlas_tacs.rs.nii.gz',SUB_DIR{s}),1);
    curDATA = subfiles(sprintf('%s/glmDenoise.nii.gz',SUB_DIR{s}),1);
    tmpData = mriRoiExtract(curDATA,curROI);
    roiData(:,:,1,s) = reshape(tmpData,length(wangROIs),2);
    runDATA = subfiles(sprintf('%s/glmIndividual_run*.nii.gz',SUB_DIR{s}),1);
    numRuns(s) = length(runDATA);
    for z=1:numRuns(s)
        tmpData = mriRoiExtract(runDATA{z},curROI);
        roiData(:,:,z+1,s) = reshape(tmpData,length(wangROIs),2);
    end
end

roiMean = cellfun(@(x) nanmean(x), roiData,'uni',false);

save('tacsROIs.mat');

%% MAKE FIGURES
close all;
dataIdx = 1:2:8;
fontSize = 12;
condLabels = {'0','90','180','270'};
gcaOpts = {'tickdir','out','box','off','fontsize',fontSize,'fontname','Helvetica','linewidth',1};
for s=3%1:length(SUB_DIR);
    figure;
    for r = 1:3
        subplot(1,3,r)
        hold on
        pH(1) = plot(1:4,roiMean{r,1,1,s}(:,dataIdx),'-b','linewidth',2);
        arrayfun(@(x) plot(1:4,roiMean{r,1,x,s}(:,dataIdx),'-b','linewidth',1),2:numRuns(s)+1);
        hold on;
        pH(2) = plot(1:4,roiMean{r,2,1,s}(:,dataIdx),'-r','linewidth',2);
        arrayfun(@(x) plot(1:4,roiMean{r,2,x,s}(:,dataIdx),'-r','linewidth',1),2:numRuns(s)+1);
        title(sprintf('%s: %s %s',subNames{s},wangROIs{r},montageObj(sesNames{s})),'fontsize',12,'fontname','helvetica');
        xlim([0.5,4.5]);
        %ylim([0,2]);
        if r == 1
            legend(pH(1:2),{'left','right'},'location','southwest');
        else
        end
        set(gca,gcaOpts{:},'xtick',1:4,'xticklabel',condLabels);
    end
end