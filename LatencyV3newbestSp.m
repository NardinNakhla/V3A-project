
plaid=0;

path='C:\research\data\SuperTuneSpkTrains\';
% names={'slu062a','slu058c','slu055a','slu050a','slu048a','slu047a','slu046b',...
%     'slu045b','slu044a','slu023a','slu022a','slu017b',...
%     'ytu326c','ytu331a','ytu332b','ytu334c','ytu336c'};%dots tasks
names={'slu017b','slu022b','slu023a','slu044b','slu045c',...
    'slu046c','slu047b','slu048b','slu055a',...
    'slu050c','slu058a','slu060a','slu062c',...
    'ytu326b','ytu331a','ytu332b','ytu334c','ytu336c'};
moreSpds={'slu022c','slu044c','slu045d','slu046d','slu047c','slu048c',...
    'slu050d','slu058b','slu060b','slu062d'};
Fs=10000;
bin_ms=1;
bin=(bin_ms/1000)*Fs;%10 ms (10/1000)*Fs
win_ms=15;
win=(win_ms/1000)*Fs;
latV3=[];
latMT=[];

for j=1:length(names)
    load(['C:\research\V3 things\V3 categorized2\',names{1,j}(1:end-1),'_V3categ2.mat']);
    v3categ=sortrows(v3categ2);
    V3units=v3categ((v3categ(:,3)<=4),1:2);%|v3categ(:,3)==4
    MTunits=v3categ(v3categ(:,3)==5,1:2);
    spd2=~cellfun(@isempty, strfind(moreSpds,names{1,j}(1:end-1)));
    spd2idx=find(spd2,1);
    speed=[];
    if ~isempty(spd2idx)
        params1=load(['C:\research\data\RFiles\',names{1,j},'_TrialStructure.mat']);
        speed(1)=params1.file.taskDialogValues.minSpeedDegPerSec;
        params2=load(['C:\research\data\RFiles\',moreSpds{1,spd2idx},'_TrialStructure.mat']);
        speed(2)=params2.file.taskDialogValues.minSpeedDegPerSec;
    else
        params1=load(['C:\research\data\RFiles\',names{1,j},'_TrialStructure.mat']);
        speed(1)=params1.file.taskDialogValues.minSpeedDegPerSec;
    end
    for ci=1:size(V3units,1)+size(MTunits,1)%:length(chunum(:,1))%[1,3:8,11];
        if ci<=size(V3units,1)
            ch=V3units(ci,1);
            unit=V3units(ci,2);
        else
            ch=MTunits(ci-size(V3units,1),1);
            unit=MTunits(ci-size(V3units,1),2);
        end

        if ~isempty(spd2idx)            
            load([path,names{1,j}(1:end),num2str(ch),num2str(unit),'spktrain.mat']);
            load([path,names{1,j}(1:end),num2str(ch),num2str(unit),'spktrain_bl.mat']);
            spktrain_bl1=spktrain_bl;
            spktrain1=spktrain;
            firing1=load(['C:\research\data\SuperTuneFiringMatrix\',names{1,j},num2str(ch),num2str(unit),'firingMat']);
            
            load([path,moreSpds{1,spd2idx}(1:end),num2str(ch),num2str(unit),'spktrain.mat']);
            load([path,moreSpds{1,spd2idx}(1:end),num2str(ch),num2str(unit),'spktrain_bl.mat']);
            spktrain_bl2=spktrain_bl;
            spktrain2=spktrain;            
            firing2=load(['C:\research\data\SuperTuneFiringMatrix\',moreSpds{1,spd2idx},num2str(ch),num2str(unit),'firingMat']);
            firing=cat(length(size(firing1.firing))+1,firing1.firing,firing2.firing);
            spktrain=cat(length(size(spktrain1))+1,spktrain1(1:4000,:,:,:,:,:,:),spktrain2(1:4000,:,:,:,:,:,:));
            spktrain_bl=cat(length(size(spktrain_bl1))+1,spktrain_bl1(1:4000,:,:,:,:,:,:),spktrain_bl2(1:4000,:,:,:,:,:,:));            
        else
            load([path,names{1,j}(1:end),num2str(ch),num2str(unit),'spktrain.mat']);
            load([path,names{1,j}(1:end),num2str(ch),num2str(unit),'spktrain_bl.mat']);
            load(['C:\research\data\SuperTuneFiringMatrix\',names{1,j},num2str(ch),num2str(unit),'firingMat']);
        end
        baseline=squeeze(sum(spktrain_bl,1))*Fs/size(spktrain_bl,1);
        allstimfir=squeeze(sum(spktrain,1))*Fs/size(spktrain,1);
        [h,p] = ttest(baseline(:),allstimfir(:));
        keepCriteria=(p<=0.05);
        if keepCriteria
            if length(size(firing))>=5
                firing1=firing(:,1,:,:,3,:);
                [maxfir,I]=max(firing1(:));
                [ddir,typ,dpos,dsiz,dcoh,dspd]= ind2sub(size(firing1),I);
                spktrain=spktrain(:,ddir,1,dpos,:,dsiz,3,dspd);
            elseif length(size(firing))>=3

                firing1=firing(:,1,:,:);
                [maxfir,I]=max(firing1(:));
                if length(size(firing))==4 && isempty(spd2idx)
                    [ddir,typ,dpos,dsiz,dspd]= ind2sub(size(firing1),I);
                    spktrain=spktrain(:,ddir,1,dpos,:,dsiz,dspd);
                else
                    [ddir,typ,dpos,dspd]= ind2sub(size(firing1),I);
                    spktrain=spktrain(:,ddir,1,dpos,:,dspd);
                end
            end
            soso=reshape(spktrain,[size(spktrain,1) size(spktrain,2)*size(spktrain,3)*...
                size(spktrain,4)*size(spktrain,5)*size(spktrain,6)*size(spktrain,7)*size(spktrain,8)]);
            spktrain_bl=spktrain_bl(:,:,1,:,:,:,:);
            soso_bl=reshape(spktrain_bl,[size(spktrain_bl,1) size(spktrain_bl,2)*size(spktrain_bl,3)*...
                size(spktrain_bl,4)*size(spktrain_bl,5)*size(spktrain_bl,6)*size(spktrain_bl,7)*size(spktrain_bl,8)]);
            count=1;
            for k=1:bin:size(soso,1)-bin
                if k+win<=size(soso,1)
                    lala=sum(soso((k:k+win),:),1);
                    psthst(count)=mean(lala);
                    
                    lala_bl=sum(soso_bl((k:k+win),:),1);
                    psthst_bl(count)=mean(lala_bl);
                    count=count+1;
                else
                    lala=sum(soso((k:end),:),1);
                    psthst(count)=mean(lala);
                    
                    lala_bl=sum(soso_bl((k:end),:),1);
                    psthst_bl(count)=mean(lala_bl);
                    count=count+1;
                end
            end
            if dspd==1
                params=params1;
            elseif dspd==2
                params=params2;
            end
            xaxis=[(1:length(psthst))*params.file.taskDialogValues.interval/length(psthst)];
%             figure
%             plot(xaxis,psthst)
%             hold on
%             plot(xaxis,psthst_bl)
            psth.condition1=[psthst_bl psthst];
            parameters.nconditions=1;
            parameters.trial_duration=(params.file.taskDialogValues.interval+...
                params.file.taskDialogValues.holdFixationGrace)/1000*Fs/bin;% /1000 to convert to sec
            parameters.baseline_duration=params.file.taskDialogValues.holdFixationGrace/1000*Fs/bin;
            parameters.stim_duration=params.file.taskDialogValues.interval/1000*Fs/bin;
            parameters.analysis_period=250/1000*Fs/bin;
            parameters.analysis_startTime=40/1000*Fs/bin;
            [lat]=latency(psth,parameters,'percentile')*bin*1000/Fs;
            if ci<=size(V3units,1)
                latV3=[latV3; lat speed(dspd)];
            else
                latMT=[latMT; lat speed(dspd)];
            end
        end
    end
end
latV3=latV3(~isnan(latV3(:,1)),:);
latMT=latMT(~isnan(latMT(:,1)),:);

figure
histogram(latV3(latV3(:,2)==3,1))
hold on
histogram(latMT(latMT(:,2)==3,1))


figure
histogram(latV3(latV3(:,2)==4,1),40)
hold on
histogram(latMT(latMT(:,2)==4,1),30)


figure
histogram(latV3(latV3(:,2)==5,1))
hold on
histogram(latMT(latMT(:,2)==5,1))
%%
mean1V3=mean(latV3(latV3(:,2)==3,1));
stderr1V3=std(latV3(latV3(:,2)==3,1))/sqrt(length(latV3(latV3(:,2)==3,1)));

mean2V3=mean(latV3(latV3(:,2)==4,1));
stderr2V3=std(latV3(latV3(:,2)==4,1))/sqrt(length(latV3(latV3(:,2)==4,1)));

mean3V3=mean(latV3(latV3(:,2)==5,1));
stderr3V3=std(latV3(latV3(:,2)==5,1))/sqrt(length(latV3(latV3(:,2)==5,1)));

x=[2^3 2^4 2^5];
yV3=[mean1V3 mean2V3 mean3V3];
errorV3=[stderr1V3 stderr2V3 stderr3V3];

mean1MT=mean(latMT(latMT(:,2)==3,1));
stderr1MT=std(latMT(latMT(:,2)==3,1))/sqrt(length(latMT(latMT(:,2)==3,1)));

mean2MT=mean(latMT(latMT(:,2)==4,1));
stderr2MT=std(latMT(latMT(:,2)==4,1))/sqrt(length(latMT(latMT(:,2)==4,1)));

mean3MT=mean(latMT(latMT(:,2)==5,1));
stderr3MT=std(latMT(latMT(:,2)==5,1))/sqrt(length(latMT(latMT(:,2)==5,1)));

%x=[3 4 5];
yMT=[mean1MT mean2MT mean3MT];
errorMT=[stderr1MT stderr2MT stderr3MT];
y=[yMT; yV3]

figure
b1=bar(x,y')%,'FaceAlpha',0.5)
%b1.FaceAlpha={0.5 0.5}
% hold on
% b2=bar(x,yV3,'FaceColor',[0,0.7,0.7])
% b2.FaceAlpha=0.5
legend('MT','V3A')
hold on
er=errorbar(x-1.2,yMT,errorMT,'LineWidth',2);
er.Color = [0 0 0];
er.LineStyle = 'none';
hold on
er=errorbar(x+1.2,yV3,errorV3,'LineWidth',2);
er.Color = [0 0 0];
er.LineStyle = 'none';
xlabel('stimulus speed [deg/s]')
ylabel('mean latency [ms]')
[h,p]=ttest2(latMT(:,1),latV3(:,1))
%[p]=ranksum(latMT(:,1),latV3(:,1))
medianMT=median(latMT(:,1))
medianV3=median(latV3(:,1))
V332=[latV3(latV3(:,2)==5,1)];
V316=[latV3(latV3(:,2)==4,1)];
V38=[latV3(latV3(:,2)==3,1)];

MT32=[latMT(latMT(:,2)==5,1)];
MT16=[latMT(latMT(:,2)==4,1)];
MT8=[latMT(latMT(:,2)==3,1)];
p32=ranksum( V332 , MT32 );
meanMT=mean(latMT(:,1));
meanV3=mean(latV3(:,1));
p = ranksum( latMT(:,1) , latV3(:,1) )
p16 = ranksum( V316 , MT16 )
p8 = ranksum( V38 , MT8 )
figure
histogram(latMT(:,1),40)
hold on
histogram(latV3(:,1),40)
%%
errorMT=[std(latMT(:,1))/sqrt(length(latMT(:,1)))];
errorV3=[std(latV3(:,1))/sqrt(length(latV3(:,1)))];
y=[meanMT; meanV3]

figure
b1=bar([1 1.7],y',0.5,'FaceAlpha',0.5)

hold on
er=errorbar(1,y(1),errorMT,'LineWidth',2);
er.Color = [0 0 0];
er.LineStyle = 'none';
hold on
er=errorbar(1.7,y(2),errorV3,'LineWidth',2);
er.Color = [0 0 0];
er.LineStyle = 'none';
ylabel('mean latency [ms]')
 [h,p]=ttest2(latMT(:,1),latV3(:,1))
save(['C:\research\2019 updates\figures data\latency.mat'],...
    'latMT','latV3','errorMT','errorV3','y','p')