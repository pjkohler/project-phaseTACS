function script05_roiStats(useSplits)
    %% ADD PATHS
    codeFolder = '/Users/kohler/code';
    addpath(genpath([codeFolder,'/git/schlegel/matlab_lib']));
    addpath(genpath([codeFolder,'/git/mrC']));
    addpath(genpath([codeFolder,'/git/MRI/matlab']));
    addpath(genpath([codeFolder,'/git/libsvm']));
    addpath(genpath([codeFolder,'/git/export_fig']));

    if nargin < 1
        useSplits = true;
    else
    end

    wangROIs = {'V1','V2' 'V3','hV4' 'VO1' 'VO2' 'PHC1' 'PHC2' ...
                'TO1' 'LO2' 'LO1' 'V3B' 'V3A' 'IPS0' 'IPS1' 'IPS2' 'IPS3' 'SPL1' };
            % note truncated version with only 18 ROIs
    wangROIs{end+1} = 'all';

    TOP_DIR='/Volumes/Denali_4D2/kohler/fMRI_EXP/phaseTACS/FIRST';
    FIG_DIR = sprintf('%s/FIGURES',TOP_DIR);
    
    if ~exist(FIG_DIR,'dir')
        mkdir(FIG_DIR);
    else
    end
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
        if useSplits
            splitData = subfiles(sprintf('%s/glmSplit*.nii.gz',SUB_DIR{s}),1);
        else
            splitData = subfiles(sprintf('%s/glmIndividual_run*.nii.gz',SUB_DIR{s}),1);
        end
        numSplits(s) = length(splitData);
        for z=1:numSplits(s)
            tmpData = mriRoiExtract(splitData{z},curROI);
            roiData(:,:,z+1,s) = reshape(tmpData,length(wangROIs),2);
        end
    end
    
    if useSplits
        save('~/Desktop/roiData_splits.mat');
    else
         save('~/Desktop/roiData_runs.mat');
    end
    
    load('/Users/kohler/Desktop/roiData_splits.mat','roiData','SUB_DIR','numSplits');
    FIG_DIR = '/Users/kohler/Desktop/FIGURES';

    roiFields = size(roiData);
    roiDataTemp = cell(size(roiData,1)+1,size(roiData,2),size(roiData,3),size(roiData,4));
    for z=1:prod(roiFields(2:end));
        roiDataTemp{19,z} = cat(1,roiData{1:3,z});
    end
    roiData(19,:,:,:) = roiDataTemp(19,:,:,:);
    roiDataSelect = cell(size(roiData));
    for s=1:length(SUB_DIR);
        % pick voxels with top 25% SNR (based on full model)
        featSelect = cellfun(@(x) find(x(:,10) > prctile(x(:,10),75) ),roiData(:,:,1,s),'uni',0);
        featSelect = repmat(featSelect,[1,1,numSplits(s)+1]);
        roiDataSelect(:,:,1:numSplits(s)+1,s) = cellfun(@(x,y) x(y,:),roiData(:,:,1:numSplits(s)+1,s),featSelect,'uni',0);
    end

    roiMean = cellfun(@(x) nanmean(x), roiData,'uni',false);
    roiMeanSelect = cellfun(@(x) nanmean(x), roiDataSelect,'uni',false);
    
    

    %% MAKE FIGURES
    close all;
    dataIdx = 1:2:8;
    subLabels = {'A','B','C','D'};
    montLabels = {'CZ-OZ','P07-OZ-P08','T7-OZ-T8'};

    colors = {'b','r'};
    fontSize = 12;
    condLabels = {'-90','0','90','180'};
    lWidth = 2;
    gcaOpts = {'tickdir','out','box','off','fontsize',fontSize,'fontname','Helvetica','linewidth',lWidth,'ticklength',[.03,.03]};
    roisToPlot = [1,2,3];
    for s=1:length(SUB_DIR);
        switch split_string(sesNames{s},'_',2)
            case 'skeri0004'
                subIdx(s) = 1;
            case 'nl-0037'
                subIdx(s) = 2;
            case 'nl-0032'
                subIdx(s) = 3;
            case 'nl-0011'
                subIdx(s) = 4;
            otherwise
                msg = sprintf('\n unknown subject %s\n',split_string(sesNames{s},'_',2));
                error(msg);
        end
        montIdx(s) = find(cell2mat(cellfun(@(x) strcmp(x,montageObj(sesNames{s})),montLabels,'uni',false)));
        figure;
        curVals = double(cat(1,roiMeanSelect{roisToPlot,:,:,s}));
        curVals = curVals(:,dataIdx);
        yMax = ceil(max(curVals(:)));
        yMin = floor(min(curVals(:)));
        for r = 1:length(roisToPlot)
            for h = 1:2
                subplot(length(roisToPlot),2,h+(r-1)*2);
                hold on
                pH(h) = plot(1:4,roiMeanSelect{r,h,1,s}(:,dataIdx),sprintf('-%so',colors{h}),'linewidth',2);
                arrayfun(@(x) plot(1:4,roiMeanSelect{r,h,x,s}(:,dataIdx),sprintf('--%so',colors{h}),'linewidth',1),2:numSplits(s)+1);
                hold off;
                xlim([0.5,4.5]);
                ylim([yMin,yMax]);
                if h ==1
                    text(min(get(gca,'xlim'))*1.15,max(get(gca,'ylim'))*.95,sprintf('LH %s',wangROIs{r}),'fontsize',fontSize,'fontname','helvetica');
                    if r == 1
                        text(max(get(gca,'xlim'))+diff(get(gca,'xlim'))*.15,max(get(gca,'ylim'))+diff(get(gca,'ylim'))*.15,sprintf('subject %s: %s',subLabels{subIdx(s)},montageObj(sesNames{s})),...
                            'fontsize',fontSize,'fontname','helvetica','horizontalalignment','center');
                    elseif r == 2
                        ylabel('Beta weight (% sign. change)','fontsize',fontSize,'fontname','helvetica');
                    end
                else
                    text(min(get(gca,'xlim'))*1.15,max(get(gca,'ylim'))*.95,sprintf('RH %s',wangROIs{r}),'fontsize',fontSize,'fontname','helvetica');
                end
                set(gca,gcaOpts{:},'xtick',1:4,'xticklabel',condLabels);
            end
        end
        xlabel({'Phase delay'; '(relative to visual stimulus)'},'fontsize',fontSize,'fontname','helvetica');
        set(gcf,'Units','centimeters');
        curPos = get(gcf,'position');
        curPos(3) = 10;
        curPos(4) = 16;
        set(gcf,'position',curPos);
        figName = sprintf('%s/%s_sub%s_beta_splits.pdf',FIG_DIR,montageObj(sesNames{s}),subLabels{subIdx(s)});
        export_fig(figName,'-pdf','-transparent',gcf);
    end


%     %% CORR PLOTS
%     figure;
%     condColor = {'bo','ro','go','mao'};
%     gcaOpts = {'tickdir','out','box','off','fontsize',fontSize,'fontname','Helvetica','linewidth',1,'ticklength',[.03,.03]};
% 
%     for c = 1:4
%         subplot(4,2,1+(c-1)*2);
%         plot(leftData(c,:),leftData(c+4,:),condColor{c});
%         xlim([-5,5]);
%         ylim([-5,5]);
%         set(gca,gcaOpts{:});
%         subplot(4,2,2+(c-1)*2);
%         plot(rightData(c,:),rightData(c+4,:),condColor{c});
%         xlim([-5,5]);
%         ylim([-5,5]);
%         set(gca,gcaOpts{:});
%     end
%     set(gcf,'Units','centimeters');
%     curPos = get(gcf,'position');
%     curPos(3) = 20;
%     curPos(4) = 40;
%     set(gcf,'position',curPos);



    %% CLASSIFICATION
    roisToPlot = [1,2,3];
    dataIdx = 1:2:8;
    colors = {'-b','-r'};
    fontSize = 12;
    condLabels = {'-90','0','90','180'};
    lWidth = 2;
    gcaOpts = {'tickdir','out','box','off','fontsize',fontSize,'fontname','Helvetica','linewidth',lWidth,'ticklength',[.03,.03]};
    
    % prepare data
    close all
    offSet = [-.1 .1];
    for s=1:length(SUB_DIR);
        
        if useSplits
            figName = sprintf('%s/%s_sub%s_splits.pdf',FIG_DIR,montageObj(sesNames{s}),curSub);
            leaveOut = 1;
        else
            figName = sprintf('%s/%s_sub%s.pdf',FIG_DIR,montageObj(sesNames{s}),curSub);
            leaveOut = 2;
        end
        figure;
        for r = 1:length(wangROIs)
            leftData = [];
            rightData = [];
            targetArray = [];
            chunkArray = [];
            for z=1:numSplits(s)
                tmpData = roiDataSelect{r,1,z+1,s}(:,dataIdx);
                leftData = [leftData; permute(tmpData,[2,1])];
                tmpData = roiDataSelect{r,2,z+1,s}(:,dataIdx);
                rightData = [rightData; permute(tmpData,[2,1])];
                targetArray = [targetArray; (1:4)'];
                chunkArray = [chunkArray; ones(4,1)*z];
            end

            leftData = double(leftData);
            rightData = double(rightData);
            leftClf{r,s} = MVPA.CrossValidation(leftData,targetArray,chunkArray,'zscore','data','partitioner',leaveOut);
            rightClf{r,s} = MVPA.CrossValidation(rightData,targetArray,chunkArray,'zscore','data','partitioner',leaveOut);
            fourWayAcc(:,1,r,s) = 100.*diag(leftClf{r,s}.confusion)./sum(leftClf{r,s}.confusion,2);
            fourWayAcc(:,2,r,s) = 100.*diag(rightClf{r,s}.confusion)./sum(rightClf{r,s}.confusion,2);
            
            if ismember(r,roisToPlot)
                subplot(length(roisToPlot),1,r)
           
                hold on
                plot(0:100,ones(1,101)*25,'-','color',[.75 .75 .75],'linewidth',lWidth);
                 pH(1) = plot((1:4)+offSet(1),fourWayAcc(:,1,r,s),'-ro','linewidth',lWidth); hold on
                pH(2) = plot((1:4)+offSet(2),fourWayAcc(:,2,r,s),'-bo','linewidth',lWidth);
                ylim([0,100]);
                xlim([0.75,4.25]);
                if r == 1
                    title(sprintf('subject %s: %s',curSub,montageObj(sesNames{s})),'fontsize',fontSize,'fontname','helvetica');
                elseif r == 2
                    ylabel('SVM accuracy (%, 1 vs all)','fontsize',fontSize,'fontname','helvetica');
                elseif r == 3
                    if s == 1
                        legend(pH(1:2),{'left','right'},'location','northeast','box','off');
                    else
                    end
                    xlabel({'Phase delay'; '(relative to visual stimulus)'},'fontsize',fontSize,'fontname','helvetica');
                end
                set(gca,gcaOpts{:},'ytick',0:25:100,'xtick',1:4,'xticklabel',condLabels);

                hold off
            else
            end
        end
        set(gcf,'Units','centimeters');
        curPos = get(gcf,'position');
        curPos(3) = 8;
        curPos(4) = 16;
        set(gcf,'position',curPos);
        export_fig(figName,'-pdf','-transparent',gcf);
        %close(gcf);
    end
   %%
   roisToPlot = [1,2,3,19];
   subColors = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; 0.4940 0.1840 0.5560];
   subSymbols = {'x','sq','o','d'};
   close all
   figure;
   overallClf(:,:,1) = cell2mat(cellfun(@(x) x.mean,leftClf(1:3,:),'uni',false)');
   overallClf(:,:,2) = cell2mat(cellfun(@(x) x.mean,rightClf(1:3,:),'uni',false)');
   for m = 1:length(montageList);
       subplot(1,3,m);
       hold on
       title(montageList(m),'fontsize',fontSize,'fontname','helvetica');
       plot(0:100,ones(1,101)*25,'-','color',[.75 .75 .75],'linewidth',lWidth);
       for s=1:max(subIdx)
           meanToPlot = squeeze(mean(fourWayAcc(:,:,roisToPlot,montIdx==m & subIdx ==s)));
           if ~isempty(meanToPlot)
               plot((1:length(roisToPlot))+offSet(1),meanToPlot(1,:),subSymbols{s},'color',subColors(s,:),'markersize',10);
               pH(s) = plot((1:length(roisToPlot))+offSet(2),meanToPlot(2,:),subSymbols{s},'color',subColors(s,:),'markersize',10);
           else
           end
       end
       xlim([0.75,length(roisToPlot)+.25]);
       ylim([0,100]);
       set(gca,gcaOpts{:},'ytick',0:25:100,'xtick',1:length(roisToPlot),'xticklabel',wangROIs(roisToPlot));
       if m == 1
           ylabel('SVM accuracy (%, 4-way)','fontsize',fontSize,'fontname','helvetica');
       elseif m == 3
         lH = legend(pH,{'subA','subB','subC','subD'},'fontsize',fontSize,'fontname','helvetica','box','off');
         lPos = get(lH, 'pos');
         lPos(1) = lPos(1) * 1.1;
         lPos(2) = lPos(2) * 1;
         set(lH, 'pos',lPos);
       end
       hold off
       aboveChanceRatio(m,1) = length(find(overallClf(montIdx==m,:,1)>.25));%./numel(overallClf(montIdx==m,:,1));
       aboveChanceRatio(m,2) = length(find(overallClf(montIdx==m,:,2)>.25));%./numel(overallClf(montIdx==m,:,2));
   end
   set(gcf,'Units','centimeters');
   curPos = get(gcf,'position');
   curPos(3) = 20;
   curPos(4) = 6;
   set(gcf,'position',curPos);
   figName = sprintf('%s/allclf_splits.pdf',FIG_DIR);
   export_fig(figName,'-pdf','-transparent',gcf);

   

end

