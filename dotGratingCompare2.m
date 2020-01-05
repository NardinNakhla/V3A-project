% get firing rate matrices
% compute d prime
% compute d prime difference per unit? then scatter plot for both V3 and
% MT?
% dotsnames={'slu048b','slu047b','slu046c','slu045c','slu044b','slu022b','slu017b',...
%     'ytu326b','ytu331a','ytu332b','ytu334c','ytu336c','slu053k'};
dotsnames={'ytu326b','ytu331a','ytu332b','ytu334c','ytu336c','slu017b','slu022b',...
    'slu023a','slu045c','slu044b','slu046c','slu047b',...
    'slu048b','slu050c','slu055a','slu060a','slu062c'};
% dotsnames={'ytu326b','ytu331a','ytu332b','ytu334c','ytu336c','slu017b','slu022c',...
%     'slu023a','slu045d','slu044c','slu046d','slu047c',...
%     'slu048c','slu050d','slu055a','slu060b','slu062d'};
gratnames={'slu048d','slu047e','slu046e','slu045e','slu044d','slu022d','slu017c',...
    'ytu326a','ytu331b','ytu332d','ytu334b','ytu336a','ytu328a','slu062e',...
    'slu060d', 'slu058d','slu053i','slu050e'};
%'ytu328a','slu062e','slu060d', 'slu058d','slu053i,j','slu050e,f'
% gratnames={'slu017c'};
% dotsnames={'slu017b'};
count=1;
dots=1;
grat=0;
comparison=0;
contrast=100;%8%40;%100
coherence=100;%10;%55;%100
Fs=10000;
timeWin=(0.05*Fs:0.4*Fs);%0.4*Fs);
DIThresh=0:0.05:1.5;

pathD='C:\research\data\SuperTuneSpkTrains\';
pathG='C:\research\data\PlaidSpkTrains\';
V3GratdP=[]; MTGratdP=[]; V3DotsdP=[]; MTDotsdP=[];
if dots
    numExp=length(dotsnames);
end
if grat
    numExp=length(gratnames);
end
   
for j=1:numExp
    if grat
    if exist(['C:\research\V3 things\V3 categorized2\',gratnames{1,j}(1:end-1),'_V3categ2.mat'],'file')==2
        load(['C:\research\V3 things\V3 categorized2\',gratnames{1,j}(1:end-1),'_V3categ2.mat']);
    else
        load(['C:\research\V3 things\V3 categorized2\',gratnames{1,j}(1:end),'_V3categ2.mat']); 
    end
    end
    if dots
    if exist(['C:\research\V3 things\V3 categorized2\',dotsnames{1,j}(1:end-1),'_V3categ2.mat'],'file')==2
        load(['C:\research\V3 things\V3 categorized2\',dotsnames{1,j}(1:end-1),'_V3categ2.mat']);
    else
        load(['C:\research\V3 things\V3 categorized2\',dotsnames{1,j}(1:end),'_V3categ2.mat']); 
    end
    end
%v3categ=sortrows(v3categ2(v3categ2(:,4)<1.5,:));
v3categ=sortrows(v3categ2);

V3units=v3categ((v3categ(:,3)<=4),1:2);%|v3categ(:,3)==4
MTunits=v3categ(v3categ(:,3)==5,1:2);
if dots
dotsparams=load(['C:\research\data\RFiles\',dotsnames{1,j},'_TrialStructure.mat']);
end
if grat
gratparams=load(['C:\research\data\RFiles\',gratnames{1,j},'_TrialStructure.mat']);
sizes(j)=gratparams.file.taskDialogValues.minPatchRadius;
end
%%
if dots
mincoh=dotsparams.file.taskDialogValues.maxCoherence-(dotsparams.file.taskDialogValues.numberOfCoherences-1)*...
    (dotsparams.file.taskDialogValues.coherenceStep);
coherences=mincoh:dotsparams.file.taskDialogValues.coherenceStep:dotsparams.file.taskDialogValues.maxCoherence;
coherenceidx=find(coherences==coherence);
end
if grat
contrasts=gratparams.file.taskDialogValues.contrastArray;
contrastidx=find(contrasts>=contrast-15 & contrasts<=contrast+15);
end
%%
for ci=1:size(V3units,1)+size(MTunits,1)
    if ci<=size(V3units,1)
        ch=V3units(ci,1);
        u=V3units(ci,2);
    else
        ch=MTunits(ci-size(V3units,1),1);
        u=MTunits(ci-size(V3units,1),2);
    end
if dots
% dotfiring=load(['C:\research\data\SuperTuneFiringMatrix\',dotsnames{1,j},num2str(ch),num2str(u),'firingMat']);
% dotfiring1=squeeze(dotfiring.firing(:,1,:));
dotspktrain=load([pathD,dotsnames{1,j},num2str(ch),num2str(u),'spktrain.mat']);
dotfiringall=sum(dotspktrain.spktrain(timeWin,:,:,:,:,:,:),1)*Fs/length(timeWin);%size(gratspktrain.spktrain,1);
dotfiring1=sum(dotspktrain.spktrain(timeWin,:,1,:,:,:,:),1)*Fs/length(timeWin);%size(gratspktrain.spktrain,1);
    
dotfiring=mean(dotfiring1,5);
sizedotstim=size(dotfiring);%direction,motiontype,pos
speed=dotsparams.file.taskDialogValues.speedStepDegPerSec^dotsparams.file.taskDialogValues.minSpeedDegPerSec;

numsizdot=dotsparams.file.taskDialogValues.numberOfRadii;
numgridposdot=dotsparams.file.taskDialogValues.superTuneColumns*dotsparams.file.taskDialogValues.superTuneRows;

dotspktrainbl=load([pathD,dotsnames{1,j},num2str(ch),num2str(u),'spktrain_bl.mat']);
dotbaseline=squeeze(sum(dotspktrainbl.spktrain_bl,1))*Fs/size(dotspktrainbl.spktrain_bl,1);
blDots=mean(dotbaseline(:));

[maxfirD,ID]=max(dotfiring(:));
[dtim,ddir,motTyp,dpos,dtrial,dsiz,dcoh]= ind2sub(size(dotfiring),ID);
dirfirdot=squeeze(dotfiring(:,:,motTyp,dpos,dtrial,dsiz,coherenceidx));
dirfirdot_best= squeeze(dotfiring(:,:,motTyp,dpos,dtrial,dsiz,dcoh));
muPrefdirD_best=dirfirdot_best(ddir);
if ddir<=dotsparams.file.taskDialogValues.numberOfDirections/2
    muNulldirD_best=dirfirdot_best(ddir+(dotsparams.file.taskDialogValues.numberOfDirections/2));
else
    muNulldirD_best=dirfirdot_best(ddir-(dotsparams.file.taskDialogValues.numberOfDirections/2));
end
%SI_best=(muPrefdirD_best- muNulldirD_best)/((muPrefdirD_best-blDots)+ (muNulldirD_best-blDots));
SI_best=1-((muNulldirD_best-blDots)/(muPrefdirD_best-blDots));

%  if numgridposdot>1
%      [ddir,dpos]= ind2sub(size(dotfiring),ID);
%      dirfirdot=squeeze(dotfiring(:,dpos));
%      dirdotspktrain=squeeze(sum(dotspktrain.spktrain(:,:,1,dpos,:),1))*Fs/size(dotspktrain.spktrain,1);
%  elseif numgridposdot==1
%      [ddir]= ID;
%      dirfirdot=squeeze(dotfiring);
%      dirdotspktrain=squeeze(sum(dotspktrain.spktrain(:,:,1,:),1))*Fs/size(dotspktrain.spktrain,1);
%  end
% time,directions,numMotion,rows*columns,trialsPerFeature,sizes,coherences

 [h,p] = ttest(dotbaseline(:),dotfiringall(:));
 if coherences(coherenceidx)==100
    keepCriteriaD=(p<=0.05)&(SI_best>0.3)&(max(dotfiring1(:))>0);%>2
 else
    keepCriteriaD=(p<=0.05)& (SI_best>0.3)&(max(dotfiring1(:))>0);
 end
% keepCriteriaD= maxfirD>2*blDots;%mean(dotfiring(:));
end


if grat
%     gratingfiring=load(['C:\research\data\PlaidFiringMatrix\',gratnames{1,j},num2str(ch),num2str(u),'firingMat.mat']);
%     gratingfiring=gratingfiring.firing;
    gratspktrain=load([pathG,gratnames{1,j},num2str(ch),num2str(u),'spktrain.mat']);
    gratingfiring1=sum(gratspktrain.spktrain(timeWin,:,:,:,:,:,:),1)*Fs/length(timeWin);%size(gratspktrain.spktrain,1);
    gratingfiring=mean(gratingfiring1,7);
    sizegratstim=size(gratingfiring); %time,directions,speed,pos,size,contrast, tria#
    gratspktrainbl=load([pathG,gratnames{1,j},num2str(ch),num2str(u),'spktrain_bl.mat']);
    baseline2=squeeze(sum(gratspktrainbl.spktrain_bl,1))*Fs/size(gratspktrainbl.spktrain_bl,1);
    blGrat=mean(baseline2(:));
    speedsgrat=zeros(sizegratstim(3),1);
    for i=1:sizegratstim(3)
        speedsgrat(i)=(gratparams.file.taskDialogValues.TFbase^...
            (gratparams.file.taskDialogValues.cyclesPerSecond+i-1))/...
            (2*gratparams.file.taskDialogValues.spatialFrequency1);
        %temporal freq/spatial freq
        %spatial frequesncy multiplied by 2 because it's multiplied by 2 in
        %monkeylab code (some crazy person programmed it that way)
    end
    
    [maxfirG,IG]=max(gratingfiring(:));
     [gtime,gdir,gspd,gpos,gsiz,gcont]= ind2sub(size(gratingfiring),IG); 
    
    if comparison==1
        spdIdx=find(speedsgrat==speed);
        if isempty(spdIdx)
            display('stimulus speeds do not match')
            sorry
        end
    else
        spdIdx=gspd;%speedsgrat(bestSpdGrat);
    end

    dirfirgrat=squeeze(gratingfiring(:,:,gspd,gpos,gsiz,contrastidx));
    dirfirgrat_best= squeeze(gratingfiring(:,:,gspd,gpos,gsiz,gcont));
    muPrefdirG_best=dirfirgrat_best(gdir);
    if gdir<=gratparams.file.taskDialogValues.numberOfDirections/2
        muNulldirG_best=dirfirgrat_best(gdir+(gratparams.file.taskDialogValues.numberOfDirections/2));
    else
        muNulldirG_best=dirfirgrat_best(gdir-(gratparams.file.taskDialogValues.numberOfDirections/2));
    end
    %SI_bestG=(muPrefdirG_best- muNulldirG_best)/((muPrefdirG_best-blGrat)+ (muNulldirG_best-blGrat));
    SI_bestG=1-(( muNulldirG_best-blGrat)/(muPrefdirG_best-blGrat));


    
    %keepCriteriaG= maxfirG>2*blGrat;%mean(gratingfiring(:));
    [h,p] = ttest(baseline2(:),gratingfiring1(:));
    if gcont==contrastidx
        keepCriteriaG=(p<=0.05)&& SI_bestG>0.3;
    else
        keepCriteriaG=(p<=0.05)&& SI_bestG>0.3;
    end


end    

if comparison
    keepCriteria = ddir==gdir;
elseif dots && ~grat
    keepCriteria=keepCriteriaD;    
elseif grat && ~dots
    keepCriteria=keepCriteriaG;    
elseif grat && dots
    keepCriteria=keepCriteriaG||keepCriteriaD;
end

 if keepCriteria
    goodch(count)=ch;
    goodunit(count)=u;
    if dots && ~isempty(dirfirdot)
    muPrefdirD=dirfirdot(ddir);
    %VarPrefdirD=var(dirdotspktrain(ddir,:))?
    if ddir<=dotsparams.file.taskDialogValues.numberOfDirections/2
        muNulldirD=dirfirdot(ddir+(dotsparams.file.taskDialogValues.numberOfDirections/2));
       % VarNulldirD=var(dirdotspktrain(ddir+(dotsparams.file.taskDialogValues.numberOfDirections/2),:));
    else
        muNulldirD=dirfirdot(ddir-(dotsparams.file.taskDialogValues.numberOfDirections/2));
        %VarNulldirD=var(dirdotspktrain(ddir-(dotsparams.file.taskDialogValues.numberOfDirections/2),:));
    end
%     dprimeDots(count)=(muPrefdirD- muNulldirD)/sqrt((VarPrefdirD+VarNulldirD)/2);
    if ci<=size(V3units,1)
        %dprimeDotsV3=(muPrefdirD- muNulldirD)/((muPrefdirD-blDots)+ (muNulldirD-blDots));
        %dprimeDotsV3=1-(muNulldirD/muPrefdirD);
        dprimeDotsV3=1-((muNulldirD-blDots)/(muPrefdirD-blDots));
    else
        %dprimeDotsMT=(muPrefdirD- muNulldirD)/((muPrefdirD-blDots)+ (muNulldirD-blDots));
        %dprimeDotsMT=1-(muNulldirD/muPrefdirD);
        dprimeDotsMT=1-((muNulldirD-blDots)/(muPrefdirD-blDots));
    end
    else
        if ci<=size(V3units,1)
        dprimeDotsV3=[];
        else
        dprimeDotsMT=[];
         end
    end
    
    if grat && ~isempty(dirfirgrat)
    muPrefdirG=dirfirgrat(gdir);
    %VarPrefdirG=var(dirgratspktrain(gdir,:))?
    if gdir<=gratparams.file.taskDialogValues.numberOfDirections/2
        muNulldirG=dirfirgrat(gdir+(gratparams.file.taskDialogValues.numberOfDirections/2));
       %VarNulldirG=var(dirgratspktrain(gdir+(gratparams.file.taskDialogValues.numberOfDirections/2),:));
    else
        muNulldirG=dirfirgrat(gdir-(gratparams.file.taskDialogValues.numberOfDirections/2));
       %VarNulldirG=var(dirgratspktrain(gdir-(gratparams.file.taskDialogValues.numberOfDirections/2),:));
    end
%     dprimeGrat(count)=(muPrefdirG- muNulldirG)/sqrt((VarPrefdirG+VarNulldirG)/2);
    if ci<=size(V3units,1)
        %dprimeGratV3=(muPrefdirG- muNulldirG)/((muPrefdirG-blGrat)+ (muNulldirG-blGrat));
        dprimeGratV3=1-(( muNulldirG-blGrat)/(muPrefdirG-blGrat));
    else
       % dprimeGratMT=(muPrefdirG- muNulldirG)/((muPrefdirG-blGrat)+ (muNulldirG-blGrat));
        dprimeGratMT=1-(( muNulldirG-blGrat)/(muPrefdirG-blGrat));
    end
    else
        if ci<=size(V3units,1)
        dprimeGratV3=[];
        else
        dprimeGratMT=[];
         end
    end
    if ci>size(V3units,1)
if dprimeDotsMT>10
    display('r')
end
    end
    count=count+1;
    if ci<=size(V3units,1)
        if grat
            V3GratdP=[V3GratdP,dprimeGratV3];
        end
        if dots
            V3DotsdP=[V3DotsdP,dprimeDotsV3];
        end
    else
        if grat
            MTGratdP=[MTGratdP,dprimeGratMT];
        end
        if dots
            MTDotsdP=[MTDotsdP,dprimeDotsMT];
        end
    end

        
 end

clear dprimeGratV3 dprimeGratMT dprimeDotsV3 dprimeDotsMT
 end
end

for di=1:length(DIThresh)
    thresh=DIThresh(di);
    nDmt(di)=length(MTDotsdP(MTDotsdP<=thresh))/length(MTDotsdP);
    nGmt(di)=length(MTGratdP(MTGratdP<=thresh))/length(MTGratdP);
    nDv3(di)=length(V3DotsdP(V3DotsdP<=thresh))/length(V3DotsdP);
    nGv3(di)=length(V3GratdP(V3GratdP<=thresh))/length(V3GratdP);
end
%[p,h,stats] = ranksum(MTDotsdP,V3DotsdP);
%[h,p]=ttest(MTDotsdP,V3DotsdP(1:64));
% Two-sample Kolmogorov-Smirnov test
%[h,p] = kstest2(MTDotsdP(MTDotsdP>0),V3DotsdP(V3DotsdP>0));
if dots
[h,p] = kstest2(MTDotsdP(1:min([length(V3DotsdP),length(MTDotsdP)])),V3DotsdP(1:min([length(V3DotsdP),length(MTDotsdP)])));
figure
plot(DIThresh,nDmt,'LineWidth',2)
hold on
plot(DIThresh,nDv3,'LineWidth',2)
legend('MT','v3')
title(['Dots - coherence ',num2str(coherence), ', p = ',num2str(p)])
xlabel('Selectivity Index','LineWidth',8)
ylabel('percentile','LineWidth',8)



figure
histogram(MTDotsdP(MTDotsdP<=10 & MTDotsdP>=-1),7)
hold on
histogram(V3DotsdP(V3DotsdP<=10 & V3DotsdP>=-1),10)
legend('MT','V3')
title(['Dots, coherence = ',num2str(coherence)])
xlabel('selectivity index [(mui_p-mui_n)/(mui_p+mui_n-2*bl)]')
save(['C:\research\2019 updates\figures data\dots_',num2str(coherence),'.mat'],...
    'V3DotsdP','MTDotsdP','DIThresh','nDmt','nDv3')
end

if grat
%[p,h,stats] = ranksum(MTGratdP,V3GratdP);
%[h,p] = kstest2(MTGratdP(MTGratdP>0),V3GratdP(V3GratdP>0));
[h,p] = kstest2(MTGratdP,V3GratdP);
figure
plot(DIThresh,nGmt,'LineWidth',2)
hold on
plot(DIThresh,nGv3,'LineWidth',2)
legend('MT','v3')
title(['Grating - contrast ',num2str(contrast) ', p = ',num2str(p)])
xlabel('Selectivity Index')
ylabel('percentile')

figure
histogram(MTGratdP(MTGratdP<=10 & MTGratdP>=-1),11)
hold on
histogram(V3GratdP(V3GratdP<=10& V3GratdP>=-1),10)
legend('MT','V3')
title(['grat, contrast = ',num2str(contrast)])
xlabel('selectivity index [(mui_p-mui_n)/(mui_p+mui_n-2*bl)]')
save(['C:\research\2019 updates\figures data\grat_',num2str(contrast),'.mat'],...
    'V3GratdP','MTGratdP','DIThresh','nGmt','nGv3')
end
% figure
% histogram(MTDotsdP(MTDotsdP>DIThresh & MTDotsdP<=10),50)
% hold on
% histogram(V3DotsdP(V3DotsdP>DIThresh & V3DotsdP<=10),50)
% legend(['MT , ',num2str(length(MTDotsdP(MTDotsdP>DIThresh & MTDotsdP<=10))*100/length(MTDotsdP(MTDotsdP<=10))),'% > thresh']...
%     ,['V3 ,',num2str(length(V3DotsdP(V3DotsdP>DIThresh & V3DotsdP<=10))*100/length(V3DotsdP(V3DotsdP<=10))),'% > thresh'])
% title(['Dots, coherence = ',num2str(coherence), ', thresh = ', num2str(DIThresh)])
% xlabel('selectivity index [(mui_p-mui_n)/(mui_p+mui_n-2*bl)]')
% figure
% histogram(MTGratdP(MTGratdP>DIThresh & MTGratdP<=10),50)
% hold on
% histogram(V3GratdP(V3GratdP>DIThresh & V3GratdP<=10),50)
% legend(['MT , ',num2str(length(MTGratdP(MTGratdP>DIThresh & MTGratdP<=10))*100/length(MTGratdP(MTGratdP<=10))),'% > thresh']...
%     ,['V3 ,',num2str(length(V3GratdP(V3GratdP>DIThresh & V3GratdP<=10))*100/length(V3GratdP(V3GratdP<=10))),'% > thresh'])
% title(['grat, contrast = ',num2str(contrast), ', thresh = ', num2str(DIThresh)])
% xlabel('selectivity index [(mui_p-mui_n)/(mui_p+mui_n-2*bl)]')
% STATS=mwwtest(MTDotsdP,V3DotsdP)
% STATS2=mwwtest(MTGratdP,V3GratdP)
% figure
% ROC_data = roc_curve(V3DotsdP,MTDotsdP);
% annotation('textbox',[.83 .5 .1 .2],'String','Can responses to dots distinguish V3 from MT?','EdgeColor','none')
% figure
% ROC_data2 = roc_curve(MTGratdP,V3GratdP);
% annotation('textbox',[0.83 .5 .1 .2],'String','Can responses to grating distinguish V3 from MT?','EdgeColor','none')
% figure
% ROC_data3 = roc_curve(V3DotsdP,V3GratdP);
% annotation('textbox',[.83 .5 .1 .2],'String','Can V3 responses distinguish dots stim from grating stim?','EdgeColor','none')
% figure
% ROC_data4 = roc_curve(MTGratdP,MTDotsdP);
% annotation('textbox',[.83 .5 .1 .2],'String','Can MT responses distinguish dots stim from grating stim?','EdgeColor','none')
