names={'slu060b','slu062a','slu058c','slu055a','slu050a','slu048a','slu047a',...
    'slu046b','slu045b','slu044a','slu023a','slu022a','slu017b','ytu326b',...
    'ytu331a','ytu332b','ytu334c','ytu336c','ytu328a','slu053l'};
% {'ytu326b','ytu331a','ytu332b','ytu334c','ytu336c','slu017b','slu022b',...
%     'slu023a','slu045c','slu044b','slu046c','slu047b',...
%     'slu048b','slu050c','slu055a','slu060a','slu062c'}
% slu053??
Fs=10000;
dirSelCountv3=0;
transSelCountv3=0;
totalv3=0;
cmplxtotv3=0;
dirSelCountMT=0;
transSelCountMT=0;
totalMT=0;
cmplxtotMT=0;
count=1;
siV3=[];
siMT=[];
si2V3=[];si2V3S=[];si2V3Y=[];
si2MT=[];
count2=0;mt=0; v3=0;
timeWin=(0.05*Fs:0.35*Fs);%0.4*Fs);

%% loop over all units all days to create structs containing pref dir motion
for j=1:length(names)
    %params= load(['C:\research\data\RFiles\',names{1,j},'_TrialStructure.mat']);
    
    if strcmp(names{1,j},'ytu328a')
        path='C:\research\data\PlaidSpkTrains\';
        firingpath='C:\research\data\PlaidFiringMatrix\';
    else
        path='C:\research\data\SuperTuneSpkTrains\';
        firingpath='C:\research\data\SuperTuneFiringMatrix\';
    end
    if exist(['C:\research\Synchrony things\chunum\',names{1,j}(1:end-1),'Chunum.mat'],'file')==2
        load(['C:\research\Synchrony things\chunum\',names{1,j}(1:end-1),'Chunum.mat'])
    else
        load(['C:\research\Synchrony things\chunum\',names{1,j}(1:end),'Chunum.mat'])
    end
    if exist(['C:\research\V3 things\V3 categorized2\',names{1,j}(1:end-1),'_V3categ2.mat'],'file')==2
        load(['C:\research\V3 things\V3 categorized2\',names{1,j}(1:end-1),'_V3categ2.mat']);
    else
        load(['C:\research\V3 things\V3 categorized2\',names{1,j}(1:end),'_V3categ2.mat']); 
    end
    v3categ=sortrows(v3categ2);
   %v3categ=sortrows(v3categ2(v3categ2(:,4)<1.5,:));
    V3units=v3categ((v3categ(:,3)<=4),1:2);
    MTunits=v3categ(v3categ(:,3)==5,1:2);
    
    for ci=1:size(V3units,1)+size(MTunits,1)
        if ci<=size(V3units,1)
            ch=V3units(ci,1);
            u=V3units(ci,2);
        else
            ch=MTunits(ci-size(V3units,1),1);
            u=MTunits(ci-size(V3units,1),2);
        end
        
        spktrain=load([path,names{1,j},num2str(ch),num2str(u),'spktrain.mat']);
        % time,directions,numMotion,rows*columns,trialsPerFeature,sizes,coherences
        spktrainbl=load([path,names{1,j},num2str(ch),num2str(u),'spktrain_bl.mat']);
        baseline=squeeze(sum(spktrainbl.spktrain_bl,1))*Fs/size(spktrainbl.spktrain_bl,1);
        allstimfir=squeeze(sum(spktrain.spktrain,1))*Fs/size(spktrain.spktrain,1);
        [h,p] = ttest(baseline(:),allstimfir(:));
        keepCriteria=(p<=0.05)&(max(allstimfir(:))>0);%2
        if (max(allstimfir(:))<=0)
            count2=count2+1;
        end
        if keepCriteria
            load([firingpath,names{1,j},num2str(ch),num2str(u),'firingMat']);
            dirfir=zeros(size(firing,2),size(firing,1));
            for i=1:size(firing,2)
                if size(firing,5)>1
                     firing1=mean(sum(spktrain.spktrain(timeWin,:,i,:,:,:,3),1)*Fs/length(timeWin),5);

                    %firing1=firing(:,i,:,:,3);
                else
                %firing1=firing(:,i,:,:,:);
                        firing1=mean(sum(spktrain.spktrain(timeWin,:,i,:,:,:,:),1)*Fs/length(timeWin),5);

                end
                [maxfir,I]=max(firing1(:));
                [time,ddir(i),typ,dpos,dsiz,dcoh]= ind2sub(size(firing1),I);
                dirfir(i,:)=squeeze(firing1(:,:,typ,dpos,dsiz,dcoh));
                muPrefdir(i)=dirfir(i,ddir(i));
                if ddir(i)<=size(firing,1)/2
                    muNulldir=dirfir(i,ddir(i)+(size(firing,1)/2));
                else
                    muNulldir=dirfir(i,ddir(i)-(size(firing,1)/2));
                end
                SI(i)=(muPrefdir(i)- muNulldir)/((muPrefdir(i)-mean(baseline(:)))+ (muNulldir-mean(baseline(:))));
                SI2(i)=1-( (muNulldir-mean(baseline(:)))/(muPrefdir(i)-mean(baseline(:))));
                diffi(count)=SI2(i)-SI(i);
                count=count+1;
            end
            
            if ci<=size(V3units,1)

                if SI2(1)>=0.5
                    transSelCountv3=transSelCountv3+1;
                end
                if size(firing,2)>1
                    cmplxtotv3=cmplxtotv3+1;
                    if max(SI2)>=0.5
                        dirSelCountv3=dirSelCountv3+1;
                    end
                end
                totalv3=totalv3+1;
                siV3=[siV3 SI(1)];
                if strcmp(names{1,j}(1),'s')
                si2V3S=[si2V3S SI2(1)];
                else
                si2V3Y=[si2V3Y SI2(1)];
                end
                si2V3=[si2V3 SI2(1)];

            end
            
            if ci>size(V3units,1)

                if SI2(1)>=0.5
                    transSelCountMT=transSelCountMT+1;
                end
                if size(firing,2)>1
                    cmplxtotMT=cmplxtotMT+1;
                    if max(SI2)>=0.5
                        dirSelCountMT=dirSelCountMT+1;
                    end
                end
                totalMT=totalMT+1;
                siMT=[siMT SI(1)];
                si2MT=[si2MT SI2(1)];
            end
            
        end
         if ci<size(V3units,1)
             v3=v3+1;
         else
             mt=mt+1;
         end
        clear SI SI2 ddir
    end
end
transSelCountv3
dirSelCountv3
totalv3
cmplxtotv3
transSelCountMT
dirSelCountMT
totalMT
cmplxtotMT
% figure
% histogram(siV3,30)
% hold on
% histogram(siMT,30)
% figure
% histogram(si2V3,100)
% mean(si2V3)
% mean(siV3)
% hold on
% histogram(si2MT,100)
% mean(si2MT)
% mean(siMT)