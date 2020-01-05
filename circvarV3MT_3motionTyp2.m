names={'slu062a','slu058c','slu055a','slu050a','slu048a','slu047a',...
    'slu046b','slu045b','slu044a','slu023a','slu022a','slu017b','ytu326b',...
    'ytu331a','ytu332b','ytu334c','ytu336c'};

path='C:\research\data\SuperTuneSpkTrains\';
Fs=10000;
SIV32=[];
SIMT2=[];
TunWMT2=[];
TunWV32=[];
circvarV32=[];
circvarMT2=[];
SIV31=[];
SIMT1=[];
TunWMT1=[];
TunWV31=[];
circvarV31=[];
circvarMT1=[];
SIV33=[];
SIMT3=[];
TunWMT3=[];
TunWV33=[];
circvarV33=[];
circvarMT3=[];
RSQV33=[];RSQV32=[];RSQV31=[];RSQMT3=[];RSQMT2=[];RSQMT1=[];
blV33=[];blV32=[]; blV31=[]; blMT1=[]; blMT2=[];blMT3=[];
timeWin=(0.05*Fs:0.35*Fs);%0.4*Fs);
newcriteria=1
for j=1:length(names)
    params=load(['C:\research\data\RFiles\',names{1,j},'_TrialStructure.mat']);
    if exist(['C:\research\Synchrony things\chunum\',names{1,j}(1:end-1),'Chunum.mat'],'file')==2
        load(['C:\research\Synchrony things\chunum\',names{1,j}(1:end-1),'Chunum.mat'])
    else
        load(['C:\research\Synchrony things\chunum\',names{1,j}(1:end),'Chunum.mat'])
    end
    load(['C:\research\V3 things\V3 categorized2\',names{1,j}(1:end-1),'_V3categ2.mat']);
    v3categ=sortrows(v3categ2);
  % v3categ=sortrows(v3categ2(v3categ2(:,4)<1.5,:));
    V3units=v3categ((v3categ(:,3)<=4),1:2);%|v3categ(:,3)==4
    MTunits=v3categ(v3categ(:,3)==5,1:2);
    for ci=1:size(V3units,1)+size(MTunits,1)
        if ci<=size(V3units,1)
            ch=V3units(ci,1);
            u=V3units(ci,2);
        else
            ch=MTunits(ci-size(V3units,1),1);
            u=MTunits(ci-size(V3units,1),2);
        end
        firing=load(['C:\research\data\SuperTuneFiringMatrix\',names{1,j},num2str(ch),num2str(u),'firingMat']);
        if size(firing.firing,2)==3
            for motTyp=1:3
%                 if motTyp==1
%                     firing1=firing.firing(:,1,:,:,:);
%                 elseif motTyp==2
%                     firing1=firing.firing(:,2,:,:,:);
%                 elseif motTyp==3
%                     firing1=firing.firing(:,3,:,:,:);
%                 end
                spktrain=load([path,names{1,j},num2str(ch),num2str(u),'spktrain.mat']);
                % time,directions,numMotion,rows*columns,trialsPerFeature,sizes,coherences
                spktrainbl=load([path,names{1,j},num2str(ch),num2str(u),'spktrain_bl.mat']);
                baseline=squeeze(sum(spktrainbl.spktrain_bl,1))*Fs/size(spktrainbl.spktrain_bl,1);
                allstimfir=sum(spktrain.spktrain(timeWin,:,:,:,:,:,:),1)*Fs/length(timeWin);
                if motTyp==1
                    firing2=mean(sum(spktrain.spktrain(timeWin,:,1,:,:,:,:),1)*Fs/length(timeWin),5);
                elseif motTyp==2
                    firing2=mean(sum(spktrain.spktrain(timeWin,:,2,:,:,:,:),1)*Fs/length(timeWin),5);
                elseif motTyp==3
                    firing2=mean(sum(spktrain.spktrain(timeWin,:,3,:,:,:,:),1)*Fs/length(timeWin),5);
                end
                
                
                [h,p] = ttest(baseline(:),allstimfir(:));
                keepCriteria=(p<=0.05)&(max(allstimfir(:))>0);
                %
                bl=mean(baseline(:));
                [maxfir,I]=max(firing2(:));
                [time,ddir,typ,dpos,dsiz,dcoh]= ind2sub(size(firing2),I);
                dirfir=squeeze(firing2(:,:,typ,dpos,dsiz,dcoh));
                
                if keepCriteria
                    muPrefdirD=dirfir(ddir);
                    if ddir<=params.file.taskDialogValues.numberOfDirections/2
                        muNulldirD=dirfir(ddir+(params.file.taskDialogValues.numberOfDirections/2));
                    else
                        muNulldirD=dirfir(ddir-(params.file.taskDialogValues.numberOfDirections/2));
                    end
                    %dprimeDots=(muPrefdirD- muNulldirD)/((muPrefdirD-bl)+ (muNulldirD-bl));
                    dprimeDots=1-((muNulldirD-bl)/(muPrefdirD-bl));
                    
                   % a=mtfit(dirfir);
                    xaxis=0:360/length(dirfir):(length(dirfir)-1)*(360/length(dirfir));
%                     fitted=vonMises(a',xaxis'*pi/180);
%                     rsq=sum((fitted-mean(dirfir)).^2)/sum((dirfir-mean(dirfir)).^2);
%                     peak=find(fitted==max(fitted));
%                     k=floor(length(dirfir)/2)-peak;
%                     if k>0
%                         fittedshift=circshift(fitted,k);
%                     else
%                         fittedshift=circshift(fitted,8+k);
%                     end
%                     [pks,locs,TunWidth,~] = findpeaks(fittedshift,xaxis);
                    circvar=circ_var(xaxis'*pi/180,squeeze(dirfir'));
%                     if circvar==0
%                         display('r')
%                     end
                    if motTyp==1
                        if ci<=size(V3units,1)
                            SIV31=[SIV31,dprimeDots];
                           % TunWV31=[TunWV31,TunWidth];
                            blV31=[blV31,bl];
                            circvarV31=[circvarV31,circvar];
                           % RSQV31=[RSQV31 rsq];
                        else
                            SIMT1=[SIMT1,dprimeDots];
                            %TunWMT1=[TunWMT1,TunWidth];
                            blMT1=[blMT1,bl];
                            circvarMT1=[circvarMT1,circvar];
                            %RSQMT1=[RSQMT1 rsq];
                        end
                    elseif motTyp==2
                        if ci<=size(V3units,1)
                            SIV32=[SIV32,dprimeDots];
                           % TunWV32=[TunWV32,TunWidth];
                            blV32=[blV32,bl];
                            circvarV32=[circvarV32,circvar];
                            %RSQV32=[RSQV32 rsq];
                        else
                            SIMT2=[SIMT2,dprimeDots];
                            %TunWMT2=[TunWMT2,TunWidth];
                            blMT2=[blMT2,bl];
                            circvarMT2=[circvarMT2,circvar];
                            %RSQMT2=[RSQMT2 rsq];
                        end
                    elseif motTyp==3
                        if ci<=size(V3units,1)
                            SIV33=[SIV33,dprimeDots];
                           % TunWV33=[TunWV33,TunWidth];
                            blV33=[blV33,bl];
                            circvarV33=[circvarV33,circvar];
                           % RSQV33=[RSQV33 rsq];
                        else
                            SIMT3=[SIMT3,dprimeDots];
                           % TunWMT3=[TunWMT3,TunWidth];
                            blMT3=[blMT3,bl];
                            circvarMT3=[circvarMT3,circvar];
                           % RSQMT3=[RSQMT3 rsq];
                        end
                        
                    end
                    
                end
            end
        end
        
        
    end
end
%% new criteria
maxSIMT=max([SIMT1;SIMT2;SIMT3]);
maxSIV3=max([SIV31;SIV32;SIV33]);

% if newcriteria
%     circvarMT1h=circvarMT1(maxSIMT>0.3&circvarMT1>0);
%     circvarMT2h=circvarMT2(maxSIMT>0.3&circvarMT2>0);
%     circvarMT3h=circvarMT3(maxSIMT>0.3&circvarMT3>0);
%     circvarV31h=circvarV31(maxSIV3>0.3&circvarV31>0);
%     circvarV32h=circvarV32(maxSIV3>0.3&circvarV32>0);
%     circvarV33h=circvarV33(maxSIV3>0.3&circvarV33>0);
% end

if newcriteria
    circvarMT1h=circvarMT1(maxSIMT>0.3);
    circvarMT2h=circvarMT2(maxSIMT>0.3);
    circvarMT3h=circvarMT3(maxSIMT>0.3);
    circvarV31h=circvarV31(maxSIV3>0.3);
    circvarV32h=circvarV32(maxSIV3>0.3);
    circvarV33h=circvarV33(maxSIV3>0.3);
end
% if newcriteria2
%     circvarMT1h=circvarMT1(SIMT1>0.3);
%     circvarMT2h=circvarMT2(SIMT2>0.3);
%     circvarMT3h=circvarMT3(SIMT3>0.3);
%     circvarV31h=circvarV31(SIV31>0.3);
%     circvarV32h=circvarV32(SIV32>0.3);
%     circvarV33h=circvarV33(SIV33>0.3);
% end
%%
varMT3=nanmedian(circvarMT3h);
varV33=nanmedian(circvarV33h);
varMT2=nanmedian(circvarMT2h);
varV32=nanmedian(circvarV32h);
varMT1=nanmedian(circvarMT1h);
varV31=nanmedian(circvarV31h);
varMTall=nanmedian([circvarMT1h,circvarMT2h,circvarMT3h]);
varV3all=nanmedian([circvarV31h,circvarV32h,circvarV32h]);

[p1,h,stats] = ranksum(circvarMT1h,circvarV31h)
figure
histogram(circvarMT1h,40)
hold on
histogram(circvarV31h,40)
%legend(['MT tuning width'],['V3 tuning width'])
legend(['MT circular var, mean =',num2str(varMT1)],['V3 circular var, mean =',num2str(varV31)])
title(['Circular Variance Translation motion, p value ',num2str(p1)])


[p2,h,stats] = ranksum(circvarMT2h,circvarV32h)
figure
histogram(circvarMT2h,40)
hold on
histogram(circvarV32h,40)
%legend(['MT tuning width'],['V3 tuning width'])
legend(['MT circular var, mean =',num2str(varMT2)],['V3 circular var, mean =',num2str(varV32)])
title(['Direction circular var for spiral motion, p value ',num2str(p2)])

[p3,h,stats] = ranksum(circvarMT3h,circvarV33h)
figure
histogram(circvarMT3h,40)
hold on
histogram(circvarV33h,40)
%legend(['MT tuning width'],['V3 tuning width'])
legend(['MT circular var, mean =',num2str(varMT3)],['V3 circular var, mean =',num2str(varV33)])
title(['Direction circular var for deformation motion, p value ',num2str(p3)])
tabledata=[varMT1,varV31,p1;varMT2,varV32,p2;varMT3,varV33,p3]

save(['C:\research\2019 updates\figures data\circVar.mat'],...
    'circvarMT1h','circvarV31h','circvarMT2h','circvarV32h','circvarMT3h','circvarV33h','tabledata')

[pMT,h,stats] = ranksum(circvarMT1h,[circvarMT2h,circvarMT3h])
[pV3,h,stats] = ranksum(circvarV31h,[circvarV32h,circvarV33h])


%%
% [p,h,stats] = ranksum(circvarMT1,circvarMT)
% figure
% histogram(circvarMT1,40)
% hold on
% histogram(circvarMT,40)
% legend(['translational motion'],['complex motion'])
% title(['Direction circular var for MT, p value ',num2str(p)])
% 
% % [p,h,stats] = ranksum(circvarV31,circvarV3)
% % figure
% % histogram(circvarV31,40)
% % hold on
% % histogram(circvarV3,40)
% % legend(['translational motion'],['complex motion'])
% % title(['Direction circular var for V3, p value ',num2str(p)])
% % 
% % [p,h,stats] = ranksum([circvarMT1,circvarMT],[circvarV31,circvarV3])
% % figure
% % histogram([circvarMT1,circvarMT],40)
% % hold on
% % histogram([circvarV31,circvarV3],40)
% % legend(['MT circular var, median =',num2str(varMTall)],['V3 circular var, median =',num2str(varV3all)])
% % title(['Direction circular var, all motion types (VonMises), p value ',num2str(p)])
% %%
% allcirvar=[circvarMT1,circvarMT,circvarV31,circvarV3];
% area=blanks(length(allcirvar));
% motion=area;
% area(1:length([circvarMT1,circvarMT]))='M';
% area(length([circvarMT1,circvarMT])+1:end)='V';
% motion(1:length(circvarMT1))='S';
% motion(length(circvarMT1)+1:length(circvarMT)+length(circvarMT1))='C';
% motion(length(circvarMT1)*2+1:length(circvarMT1)*2+length(circvarV31))='S';
% motion(length(circvarMT1)*2+length(circvarV31)+1:end)='C';
% 
% p2 = anovan(allcirvar,{area motion},'model','interaction','varnames',{'Cortical Area',...
%     'Motion Type'});
% %[p3,tbl] = friedman(data,192);
% 
% %% new criteria
% if newcriteria
%     TunWMT=TunWMT(RSQMT>=0.5);
%     TunWMT1=TunWMT1(RSQMT1>=0.5);
%     TunWV31=TunWV31(RSQV31>=0.5);
%     TunWV3=TunWV3(RSQV3>0.3);
% end
% %% other tests
% 
% meanMT=median(TunWMT);
% meanV3=median(TunWV3);
% meanMT1=median(TunWMT1);
% meanV31=median(TunWV31);
% meanMTall=median([TunWMT1,TunWMT]);
% meanV3all=median([TunWV31,TunWV3]);
% [p,h,stats] = ranksum(TunWMT,TunWV3)
% figure
% histogram(TunWMT,40)
% hold on
% histogram(TunWV3,40)
% legend(['MT tuning width, mean =',num2str(meanMT)],['V3 tuning width, mean =',num2str(meanV3)])
% title(['Direction Tuning width for complex motion (VonMises), p value ',num2str(p)])
% %
% [p,h,stats] = ranksum(TunWMT1,TunWMT)
% figure
% histogram(TunWMT1,40)
% hold on
% histogram(TunWMT,40)
% legend(['translational motion'],['complex motion'])
% title(['Direction Tuning width for MT (VonMises), p value ',num2str(p)])
% %
% [p,h,stats] = ranksum(TunWV31,TunWV3)
% figure
% histogram(TunWV31,40)
% hold on
% histogram(TunWV3,40)
% legend(['translational motion'],['complex motion'])
% title(['Direction Tuning width for V3 (VonMises), p value ',num2str(p)])
% %
% [p,h,stats] = ranksum(TunWMT1,TunWV31)
% figure
% histogram(TunWMT1,40)
% hold on
% histogram(TunWV31,40)
% legend(['MT tuning width, median =',num2str(meanMT1)],['V3 tuning width, median =',num2str(meanV31)])
% title(['Direction Tuning width for simple motion (VonMises), p value ',num2str(p)])
% 
% [p,h,stats] = ranksum([TunWMT1,TunWMT],[TunWV31,TunWV3])
% figure
% histogram([TunWMT1,TunWMT],40)
% hold on
% histogram([TunWV31,TunWV3],40)
% legend(['MT tuning width, median =',num2str(meanMTall)],['V3 tuning width, median =',num2str(meanV3all)])
% title(['Direction Tuning width, all motion types (VonMises), p value ',num2str(p)])
% median([TunWMT1,TunWMT])
% median([TunWV31,TunWV3])
% 
