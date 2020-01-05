names={'slu048d','slu047e','slu046e','slu045e','slu044d','slu022d','slu017c',...
    'ytu326a','ytu331b','ytu332d','ytu334b','ytu336a','ytu328a','slu060d','slu062e',...
     'slu058d','slu053i','slu050e'};

% slu053??
Fs=10000;
dirSelCountV3=0;
dirSelCountMT=0;

transSelCount=0;
total=0;
cmplxtot=0;
v3total=0;
MTtotal=0;
timeWin=(0.05*Fs:0.35*Fs);%0.4*Fs);

%% loop over all units all days to create structs containing pref dir motion
for j=1:length(names)
    %params= load(['C:\research\data\RFiles\',names{1,j},'_TrialStructure.mat']);
    count=1;
        path='C:\research\data\PlaidSpkTrains\';
        firingpath='C:\research\data\PlaidFiringMatrix\';


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
        keepCriteria=(p<=0.05);
        if keepCriteria
            load([firingpath,names{1,j},num2str(ch),num2str(u),'firingMat']);
            dirfir=zeros(1,size(firing,1));
                if length(size(firing))>4
                %firing1=firing(:,:,:,:,3);
                    firing1=mean(sum(spktrain.spktrain(timeWin,:,:,:,:,3,:),1)*Fs/length(timeWin),7);
                else
                  % firing1=firing(:,:,:,:);
                   firing1=mean(sum(spktrain.spktrain(timeWin,:,:,:,:,:),1)*Fs/length(timeWin),7);
                end
                [maxfir,I]=max(firing1(:));
                [dtime,ddir,spd,dpos,dsiz,dcon]= ind2sub(size(firing1),I);
                dirfir=squeeze(firing1(:,:,spd,dpos,dsiz,dcon));
                muPrefdir=dirfir(ddir);
                if ddir<=size(firing,1)/2
                    muNulldir=dirfir(ddir+(size(firing,1)/2));
                else
                    muNulldir=dirfir(ddir-(size(firing,1)/2));
                end
                %SI=(muPrefdir- muNulldir)/((muPrefdir-mean(baseline(:)))+ (muNulldir-mean(baseline(:))));
                SI=1-((muNulldir-mean(baseline(:)))/(muPrefdir-mean(baseline(:))));
            
            if ci<=size(V3units,1)
                if SI>=0.5
                    dirSelCountV3=dirSelCountV3+1;
                end
                 v3total=v3total+1;
            else
                 if SI>=0.5
                    dirSelCountMT=dirSelCountMT+1;
                 end
                 MTtotal=MTtotal+1;
            end
            
        end
        total=total+1;
    end
end

dirSelCountMT
MTtotal
dirSelCountV3
v3total

