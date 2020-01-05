
names={'slu062a','slu058c','slu055a','slu050a','slu048a','slu047a','slu046b','slu045b','slu044a','slu023a','slu022a','slu017b',...
'ytu326b','ytu331a','ytu332b','ytu334c','ytu336c'};

path='C:\research\data\SuperTuneSpkTrains\';
Fs=10000;
SIV3=[];
SIMT=[];
PV3=[];
PMT=[];
count=1;
nameV3=[];
nameMT=[];
timeWin=(0.05*Fs:0.35*Fs);%0.4*Fs);

for j=1:length(names)
params=load(['C:\research\data\RFiles\',names{1,j},'_TrialStructure.mat']);
if exist(['C:\research\Synchrony things\chunum\',names{1,j}(1:end-1),'Chunum.mat'],'file')==2
    load(['C:\research\Synchrony things\chunum\',names{1,j}(1:end-1),'Chunum.mat'])
else
    load(['C:\research\Synchrony things\chunum\',names{1,j}(1:end),'Chunum.mat'])
end
load(['C:\research\V3 things\V3 categorized2\',names{1,j}(1:end-1),'_V3categ2.mat']);
v3categ=sortrows(v3categ2);
%v3categ=sortrows(v3categ2(v3categ2(:,4)<1.5,:));
V3units=v3categ((v3categ(:,3)<=4),1:2);%|v3categ(:,3)==4
MTunits=v3categ(v3categ(:,3)==5,1:2);
tuningW=zeros(length(v3categ(:,1)),3);
for ci=1:size(V3units,1)+size(MTunits,1)
    if ci<=size(V3units,1)
        ch=V3units(ci,1);
        u=V3units(ci,2);
    else
        ch=MTunits(ci-size(V3units,1),1);
        u=MTunits(ci-size(V3units,1),2);
    end
firing=load(['C:\research\data\SuperTuneFiringMatrix\',names{1,j},num2str(ch),num2str(u),'firingMat']);
numMot=size(firing.firing,2);

spktrain=load([path,names{1,j},num2str(ch),num2str(u),'spktrain.mat']);
% time,directions,numMotion,rows*columns,trialsPerFeature,sizes,coherences
spktrainbl=load([path,names{1,j},num2str(ch),num2str(u),'spktrain_bl.mat']);
baseline=sum(spktrainbl.spktrain_bl,1)*Fs/size(spktrainbl.spktrain_bl,1);
allstimfir=sum(spktrain.spktrain(timeWin,:,:,:,:,:,:),1)*Fs/length(timeWin);

 [h,p] = ttest(baseline(:),allstimfir(:));
 keepCriteria=(p<=0.05);
%
bl=mean(baseline(:));

 if keepCriteria && numMot==3
    goodch(count)=ch;
    goodunit(count)=u;
    for c = 1:numMot
        firing1=mean(sum(spktrain.spktrain(timeWin,:,c,:,:,:,:),1)*Fs/length(timeWin),5);
        %firing1=firing.firing(:,c,:,:,:);
        [maxfir,I]=max(firing1(:));
        P(c)=maxfir;
        [dtm,ddir,typ,dpos,dsiz,dcoh]= ind2sub(size(firing1),I);
        dirfir=squeeze(firing1(:,:,typ,dpos,dsiz,dcoh));
        muPrefdirD=dirfir(ddir);
        if ddir<=params.file.taskDialogValues.numberOfDirections/2
            muNulldirD=dirfir(ddir+(params.file.taskDialogValues.numberOfDirections/2));
        else
            muNulldirD=dirfir(ddir-(params.file.taskDialogValues.numberOfDirections/2));
        end
       % SI(c)=(muPrefdirD- muNulldirD)/((muPrefdirD-bl)+ (muNulldirD-bl));
        SI(c)=1-((muNulldirD-bl)/(muPrefdirD-bl));
    end
    if ci<=size(V3units,1)
        SIV3=[SIV3;SI];
        PV3=[PV3;P];
        nameV3=[nameV3,{[names{1,j},num2str(ch),num2str(u)]}];
        clear SI
    else
        SIMT=[SIMT;SI];
        PMT=[PMT;P];
        nameMT=[nameMT,{[names{1,j},num2str(ch),num2str(u)]}];
        clear SI
    end
    count=count+1;
 end


end

end
%%
% [xx,yy]=find(SIMT>3);
% SIMT(xx,:)=[];
% [xxx,yyy]=find(SIV3>3);
% SIV3(xxx,:)=[];
% figure;scatter3(SIV3(:,1),SIV3(:,2),SIV3(:,3))
% xlabel('translation')
% ylabel('spiral')
% zlabel('deformation')
% figure;scatter3(SIMT(:,1),SIMT(:,2),SIMT(:,3))
% xlabel('translation')
% ylabel('spiral')
% zlabel('deformation')
% 
% figure;scatter3(PV3(:,1),PV3(:,2),PV3(:,3))
% xlabel('translation')
% ylabel('spiral')
% zlabel('deformation')
% title('V3')
% figure;scatter3(PMT(:,1),PMT(:,2),PMT(:,3))
% xlabel('translation')
% ylabel('spiral')
% zlabel('deformation')
% title('MT')
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1
%find index
for kkk=1:length(nameV3)
    if strcmp(nameV3{1,kkk},'slu017b71')
        indexC=kkk;
    end
    if strcmp(nameV3{1,kkk},'slu050a72')
        indexS=kkk;
    end
end
% figure
% %-- Plot the axis system
% [h,hg,htick]=terplot;
% %-- Plot the data ...
% hter=ternaryc(SIV3(:,1),SIV3(:,2),SIV3(:,3));
% %-- ... and modify the symbol:
% set(hter,'marker','o','markerfacecolor','none','markersize',4)
% hlabels=terlabel('translation','spirals','deformation');
% title('V3 SI')
% hold on
% hter=ternaryc(SIV3([indexS,indexC],1),SIV3([indexS,indexC],2),SIV3([indexS,indexC],3));
% %-- ... and modify the symbol:
% set(hter,'marker','o','markerfacecolor','r','markersize',8)
% 
% 
% figure
% %-- Plot the axis system
% [h,hg,htick]=terplot;
% %-- Plot the data ...
% hter=ternaryc(SIMT(:,1),SIMT(:,2),SIMT(:,3));
% %-- ... and modify the symbol:
% set(hter,'marker','o','markerfacecolor','none','markersize',4)
% hlabels=terlabel('translation','spirals','deformation');
% title('MT SI')
% 
% figure
% %-- Plot the axis system
% [h,hg,htick]=terplot;
% %-- Plot the data ...
% hter=ternaryc(PV3(:,1),PV3(:,2),PV3(:,3));
% %-- ... and modify the symbol:
% set(hter,'marker','o','markerfacecolor','none','markersize',4)
% hlabels=terlabel('translation','spirals','deformation');
% title('V3 maximum firing')
% 
% figure
% %-- Plot the axis system
% [h,hg,htick]=terplot;
% %-- Plot the data ...
% hter=ternaryc(PMT(:,1),PMT(:,2),PMT(:,3));
% %-- ... and modify the symbol:
% set(hter,'marker','o','markerfacecolor','none','markersize',4)
% hlabels=terlabel('translation','spirals','deformation');
% title('MT maximum firing')
% 
% 
% %%
% % [xx,yy]=find(SIMT>3);
% % SIMT(xx,:)=[];
% % [xxx,yyy]=find(SIV3>3);
% % SIV3(xxx,:)=[];
% %%
% figure;scatter(SIV3(:,1),SIV3(:,2))%,SIV3(:,3))
% V3trsp=sum(SIV3(:,1)>SIV3(:,2))/length(SIV3(:,1))*100;
% title(['V3,  SI_t > SI_s % = ',num2str(V3trsp)])
% xlabel('translation')
% ylabel('spiral')
% hold on
% plot(SIV3(:,1),SIV3(:,1))
% % xlim([0 max([SIV3(:);SIMT(:)])])
% % ylim([0 max([SIV3(:);SIMT(:)])])
% xlim([0 2])
% ylim([0 2])
% 
% [h,p,ci,stats] = ttest(SIV3(:,1)-SIV3(:,2));
% 
% figure
% histogram(SIV3(:,1)-SIV3(:,2))
% title(['tr-sp V3 p value for t test = ', num2str(p)])
% %% 2
% figure;scatter(SIMT(:,1),SIMT(:,2),'b')%,SIMT(:,3))
% xlabel('translation')
% ylabel('spiral')
% MTtrsp=sum(SIMT(:,1)>SIMT(:,2))/length(SIMT(:,1))*100;
% title(['MT,  SI_t > SI_s % = ',num2str(MTtrsp)])
% hold on
% plot(SIV3(:,1),SIV3(:,1))
% hold on
% scatter(SIV3(:,1),SIV3(:,2),'r')%,SIMT(:,3))
% % xlim([0 max([SIV3(:);SIMT(:)])])
% % ylim([0 max([SIV3(:);SIMT(:)])])
% xlim([0 2])
% ylim([0 2])
% 
% [h,p,ci,stats] = ttest(SIMT(:,1)-SIMT(:,2));
% 
% figure
% histogram(SIMT(:,1)-SIMT(:,2))
% title(['tr-sp MT p value for t test = ', num2str(p)])
% %% 3
% figure;scatter(SIV3(:,1),SIV3(:,3))%,SIV3(:,3))
% xlabel('translation')
% ylabel('deformation')
% V3trdf=sum(SIV3(:,1)>SIV3(:,3))/length(SIV3(:,1))*100;
% title(['V3,  SI_t > SI_d % = ',num2str(V3trdf)])
% hold on
% plot(SIV3(:,1),SIV3(:,1))
% % xlim([0 max([SIV3(:);SIMT(:)])])
% % ylim([0 max([SIV3(:);SIMT(:)])])
% xlim([0 2])
% ylim([0 2])
% 
% [h,p,ci,stats] = ttest(SIV3(:,1)-SIV3(:,3));
% 
% figure
% histogram(SIV3(:,1)-SIV3(:,3))
% title(['tr-df V3 p value for t test = ', num2str(p)])
% %% 4
% figure;scatter(SIMT(:,1),SIMT(:,3))%,SIMT(:,3))
% xlabel('translation')
% ylabel('deformation')
% MTtrdf=sum(SIMT(:,1)>SIMT(:,3))/length(SIMT(:,1))*100;
% title(['MT,  SI_t > SI_d % = ',num2str(MTtrdf)])
% hold on
% plot(SIV3(:,1),SIV3(:,1))
% % xlim([0 max([SIV3(:);SIMT(:)])])
% % ylim([0 max([SIV3(:);SIMT(:)])])
% xlim([0 2])
% ylim([0 2])
% 
% [h,p,ci,stats] = ttest(SIMT(:,1)-SIMT(:,3));
% 
% figure
% histogram(SIMT(:,1)-SIMT(:,3))
% title(['tr-df MT p value for t test = ', num2str(p)])
% %% 5
% figure;scatter(SIV3(:,2),SIV3(:,3))%,SIV3(:,3))
% xlabel('spiral')
% ylabel('deformation')
% V3spdf=sum(SIV3(:,2)>SIV3(:,3))/length(SIV3(:,2))*100;
% title(['V3,  SI_s > SI_d % = ',num2str(V3spdf)])
% hold on
% plot(SIV3(:,1),SIV3(:,1))
% % xlim([0 max([SIV3(:);SIMT(:)])])
% % ylim([0 max([SIV3(:);SIMT(:)])])
% xlim([0 2])
% ylim([0 2])
% 
% [h,p,ci,stats] = ttest(SIV3(:,2)-SIV3(:,3));
% 
% figure
% histogram(SIV3(:,2)-SIV3(:,3))
% title(['sp-df V3 p value for t test = ', num2str(p)])
% %% 6
% figure;scatter(SIMT(:,2),SIMT(:,3))%,SIMT(:,3))
% xlabel('spiral')
% ylabel('deformation')
% MTspdf=sum(SIMT(:,2)>SIMT(:,3))/length(SIMT(:,1))*100;
% title(['MT,  SI_s > SI_d % = ',num2str(MTspdf)])
% hold on
% plot(SIV3(:,1),SIV3(:,1))
% % xlim([0 max([SIV3(:);SIMT(:)])])
% % ylim([0 max([SIV3(:);SIMT(:)])])
% xlim([0 2])
% ylim([0 2])
% 
% [h,p,ci,stats] = ttest(SIMT(:,2)-SIMT(:,3));
% 
% figure
% histogram(SIMT(:,2)-SIMT(:,3))
% title(['sp-df MT p value for t test = ', num2str(p)])
% 
% data=[SIV3(1:175,:);SIMT];
% %[p,tbl] = anova2(data,175)
% 
pt = ranksum(SIV3(:,1),SIMT(:,1))
ps = ranksum(SIV3(:,2),SIMT(:,2))
pd = ranksum(SIV3(:,3),SIMT(:,3))
% %%
% figure;bubbleplot(SIV3(:,1),SIV3(:,2),SIV3(:,3),6,[],'^');
% grid on
% box on
% xlabel('Translation DSI')
% ylabel('Spirals DSI')
% zlabel('Deformation DSI')
% xlim([0 2])
% ylim([0 2])
% zlim([0 2])
% title('V3A')
% hold on
% hihi=scatter3(SIV3([indexS,indexC],1),SIV3([indexS,indexC],2),...
%     SIV3([indexS,indexC],3),80,'r','^');
% %-- ... and modify the symbol:
% set(hihi,'markerfacecolor','r');
% 
% figure;bubbleplot(SIMT(:,1),SIMT(:,2),SIMT(:,3),6,[],'o');
% grid on
% box on
% xlabel('Translation DSI')
% ylabel('Spirals DSI')
% xlim([0 2])
% ylim([0 2])
% zlim([0 2])
% title('MT')
% zlabel('Deformation DSI')


save(['C:\research\2019 updates\figures data\BubblePlotsDSI.mat'],...
    'SIMT','SIV3','indexS','indexC','pt','ps','pd')
