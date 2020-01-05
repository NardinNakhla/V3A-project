%getTrialFiles('ytu280','g','132.216.57.96')
%clear variables
%function [responses1,firing2]=supertunenana(name,ch,unit)
%name='ytu279g';
% firingAve=zeros(32,2);
% firingGp=zeros(32,9);
  Fs=10000;
 name='slu053i';%slu050e
 %'slu048d'151,'slu047e','slu046e','slu045e','slu044d'42,'slu022d','slu017c',...
   % 'ytu326a','ytu331b','ytu332d'141,'ytu334b','ytu336a'
%load(['C:\research\Synchrony things\chunum\',name(1:end),'Chunum.mat']);%(1:end-1)
if exist(['C:\research\Synchrony things\chunum\',name(1:end-1),'Chunum.mat'],'file')==2
    load(['C:\research\Synchrony things\chunum\',name(1:end-1),'Chunum.mat'])
else
    load(['C:\research\Synchrony things\chunum\',name(1:end),'Chunum.mat'])
end
% ci2=1;
 chunum=chunum;%unit3 3 5 6 14 15 18 19 20 26 27 29 unit4 18 20 27
 %chunum=chunum(goodcells,:);
load(['C:\research\data\RFiles\',name,'_TrialStructure.mat'])
rows=file.taskDialogValues.rows;
columns=file.taskDialogValues.columns;

for ci=1:length(chunum)
    counting=1;
ch=chunum(ci,1);
unit=chunum(ci,2);
[spikeMat,stimLength,baseLineLength]=plx2NFile_plaidnew(name,...
    ch,unit,1);
Contrasts=unique(spikeMat(:,9));
%spikeMat=SpikeMat(SpikeMat(:,9)==Contrast,:);
  t1=0;
  t2=mean(stimLength);
 time=floor((t2-t1)*Fs);
 time_bl=floor(mean(baseLineLength)*Fs);
%in spikemat column 2 trial number, column 3 direction of motion (8
%directions), column 4 type of motion (3 types), column 5 speed of
%translation, 
%spike times are aligned from stim on. 0 means stim on
featureIdx=3; %the 8 directions
xaxis=unique(spikeMat(:,featureIdx));
spds=unique(spikeMat(:,5));
numSpd=length(spds);
count=1:length(xaxis);
vec=zeros(length(unique(spikeMat(:,featureIdx))),1);
gridIndeces=unique(spikeMat(:,11));
sizes=unique(spikeMat(:,7));
startstim=0;
endstim=mean(stimLength)*1000;
%trialsPerFeature=length(spikeMat(spikeMat(:,1)==-1000 & spikeMat(:,4)==1 
%& spikeMat(:,7)==1,1))/length(unique(spikeMat(:,3)));

firing=zeros(length(xaxis),numSpd,length(gridIndeces),length(sizes),length(Contrasts));
firing_bl=firing;
responses1=zeros(length(xaxis),numSpd,rows,columns);
%firing2=zeros(8,3,9,trialsPerFeature);
% spikes_bl=spikeMat(spikeMat(:,1)<=0 ,:);
% firing_bl=length(spikes_bl(:,1))/(numTrials*(mean(baseLineLength)));
for m=1:length(Contrasts)
    contrast=Contrasts(m);
for j=1:length(sizes)
    siz=sizes(j);
for l=1:numSpd
        spd=spds(l);
        figure
        a=zeros(max(gridIndeces),4);
for k=1:max(gridIndeces)
    gridIndex=k;

spikes=spikeMat(spikeMat(:,1)>startstim & spikeMat(:,5)==spd & ...
    spikeMat(:,11)==gridIndex & spikeMat(:,7)==siz & spikeMat(:,9)==contrast,:);
trialsPerFeature=floor(length(unique(spikeMat(spikeMat(:,1)==-1000 & ...
     spikeMat(:,5)==spd & spikeMat(:,11)==gridIndex & spikeMat(:,7)==siz & spikeMat(:,9)==contrast,2)))...
     /length(unique(spikeMat(:,3))));
if counting==1
spktrain=zeros(time,length(xaxis),numSpd,length(gridIndeces),length(sizes),length(Contrasts),trialsPerFeature);
spktrain_bl=zeros(time_bl,length(xaxis),numSpd,length(gridIndeces),length(sizes),length(Contrasts),trialsPerFeature);
counting=2;
end
 for i=count
% hist(spikeMat(spikeMat(:,3)==i,1));
vec(i)=length(spikes(spikes(:,featureIdx)==xaxis(i),1));
vec_rep_sep=spikes(spikes(:,featureIdx)==xaxis(i),:);
spikes_bl=spikeMat(spikeMat(:,1)<=0 & spikeMat(:,1)>=(0-time_bl)*1000/Fs& spikeMat(:,5)==spd &...
    spikeMat(:,11)==gridIndex & spikeMat(:,featureIdx)==xaxis(i)& spikeMat(:,7)==siz & spikeMat(:,9)==contrast,:);
firing_bl(i,l,k,j,m)=length(spikes_bl(:,1))/...
    (trialsPerFeature*(mean(baseLineLength)));
for o=1:trialsPerFeature
    trialIdx=unique(spikeMat(spikeMat(:,1)==-1000 & spikeMat(:,5)==spd &...
    spikeMat(:,11)==gridIndex & spikeMat(:,7)==siz & spikeMat(:,9)==contrast & spikeMat(:,featureIdx)==xaxis(i),2));
    spktrain1=zeros(time,1);
indices=round(vec_rep_sep(vec_rep_sep(:,2)==trialIdx(o),1)*(Fs/1000));
indices=indices(indices<=time);
spktrain1(indices)=1;
    spktrain(:,i,l,k,j,m,o)=spktrain1;
    spktrain2=zeros(time_bl,1);
    if trialIdx(o)> length(baseLineLength)
            indices_bl=round((spikes_bl(spikes_bl(:,2)==trialIdx(o),1)+(mean(baseLineLength)*1000))*(Fs/1000));
    else
    indices_bl=round((spikes_bl(spikes_bl(:,2)==trialIdx(o),1)+(baseLineLength(trialIdx(o))*1000))*(Fs/1000));
    end
    indices_bl=indices_bl((indices_bl>0)&(indices_bl<=time_bl));
    spktrain2(indices_bl)=1;
    spktrain_bl(:,i,l,k,j,m,o)=spktrain2;
end
firing(i,l,k,j,m)=vec(i)/(trialsPerFeature*(endstim-startstim)/1000);
%firing_bl(i)=length(spikes_bl(spikes_bl(:,featureIdx)==xaxis(i),1))/(trialsPerFeature*startstim/1000);
%ylim([0 100])
 end

 
%%
nana(k)=subplot(rows,columns,k); %was 3 l k
plot(nana(k),xaxis,squeeze(firing(:,l,k,j,m)));
 hold on
 plot(xaxis,firing_bl(:,l,k,j,m),'r')
%polar(xaxis,vec)
a(k,:)=axis();

end
 [RFcenter,I]=max(firing(:));
 [i1,i2,RFcenterIdx]= ind2sub(size(firing),I);
axis([nana],[0 330 0 max(a(:,4))])

axes('Units','Normal');
h = title(['speed',num2str(spd),', size',num2str(siz),', contrast',num2str(contrast),', ',name,num2str(ch),num2str(unit)]);
set(gca,'visible','off')
set(h,'visible','on')

end
%responses1=reshape(squeeze(firing(:,:,:,j)), [length(xaxis) numSpd columns rows]);
for jjj=1:rows
    responses1(:,:,jjj,:)=squeeze(firing(:,:,(jjj-1)*columns+1:jjj*columns,j,m));
end
firing_bl2=firing_bl(:,:,:,j,m);
figure
 ShowSuperTuneAlt4(responses1,rows, columns,mean(firing_bl2(:)),numSpd);
%title([name num2str(ch),' ', num2str(unit)])
title(['speed',num2str(spd),', size',num2str(siz),', contrast',num2str(contrast),', ',name,num2str(ch),num2str(unit)]);
%legend(['baseline = ',num2str(mean(firing_bl2(:)))])
end
end
% save(['C:\research\data\SuperTuneSpkTrains\',name,num2str(chunum(ci,1)),num2str(chunum(ci,2)),'spktrain_bl.mat'],'spktrain_bl')
% save(['C:\research\data\SuperTuneSpkTrains\',name,num2str(chunum(ci,1)),num2str(chunum(ci,2)),'spktrain.mat'],'spktrain')

save(['C:\research\data\PlaidSpkTrains\',name,num2str(ch),num2str(unit),'spktrain_bl.mat'],'spktrain_bl','Fs')
save(['C:\research\data\PlaidSpkTrains\',name,num2str(ch),num2str(unit),'spktrain.mat'],'spktrain','Fs','RFcenterIdx')
save(['C:\research\data\PlaidFiringMatrix\',name,num2str(ch),num2str(unit),'firingMat.mat'],'firing')

%%
% metadata=0;
% %firing_bl=0.1;
% %responses=reshape(firing,[8 3 3 3]);
% % if strcmp(name,'ytu196a')
% % responses1(:,:,1,:)=firing(1:4,:,1:3);
% % responses1(:,:,2,:)=firing(1:4,:,4:6);
% % responses1(:,:,3,:)=firing(1:4,:,7:9);
% % else
% 
%     responses1(:,:,1,:)=firing(:,:,1:3);
% responses1(:,:,2,:)=firing(:,:,4:6);
% responses1(:,:,3,:)=firing(:,:,7:9);
% responses1(:,:,4,:)=firing(:,:,10:12);
% firingAve(ch,:)=[mean(firing(:));mean(firing_bl(:))];
% firingGp(ch,:)=mean(squeeze(mean(firing,1)),1);
% %end
% if ~isempty(spikes)|| ~isempty(spikes_bl)
% figure
 %ShowSuperTuneAlt4(responses1,metadata,mean(firing_bl(:)),numSpd);
% title([name num2str(ch),' ', num2str(unit)])
% legend(['baseline = ',num2str(mean(firing_bl(:)))])
% end
%end
 clearvars -except firingAve channel firingGp chunum name Fs rows columns
end
% mana=[firingAve,zeros(32,1),firingGp];
% figure
% scatter(firingAve(:,1),firingAve(:,2))
% xlabel('stim on')
% ylabel('baseline')
%clear variables
%end
%savefig([name num2str(ch) num2str(unit) '_TuningWheels.fig'])
%end
%end
%axis([nana],[min(min(a(:,1))) max(max(a(:,2))) min(min(a(:,3))) max(max(a(:,4)))])
%axis([min(a(:,1)) max(a(:,2)) min(a(:,3)) max(a(:,4))])