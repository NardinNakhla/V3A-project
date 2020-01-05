function  [spikeMat,stimLength,baseLineLength] = plx2NFile_plaidnew(fname,chnum,spnum,plxflag)

% creates NFile out of MonkeyLab data and Plexon data
%
% fname: plexon file name 'fname.plx'
% chnum: channel number
% spnum: spike unit number 1,2,3,4 corresponding a,b,c,d
% NFile: name of out file, which will be saved as fname11N
%
% 
if nargin < 5
    NFile = fname;
end

%close(ftp);
%ftpObj = ftp(ip,'monkeylab','monkeylab123');

% MonkeyLab files
% DIRPATH = 'C:/research/data/RFiles';
% DIRPATH = 'C:/Users/ywcui/Documents/MT_UprobData/data/PFiles';
% plexon files directory
% PLXDIR = 'C:/research/data/plexon files/';
% PLXDIR = 'C:/Users/ywcui/Documents/MT_UprobData/data/plexonfiles/';

%global DIRPATH PLXDIR RAWDPATH TASKDIR;
DIRPATH = 'C:\research\data\RFiles';
PLXDIR = 'C:\research\data\plexon files\';

NFile = [DIRPATH '\' NFile num2str(chnum) num2str(spnum) 'N'];
%cd(ftpObj, ['/Users/roberto/DATA/' fname(1) '/' fname(2) '/' fname]);
%try
%    mget(ftpObj,[fname '_TrialStructure.mat'],DIRPATH);
%catch
%   disp('trial super structure not generated yet. Please generate on MonkeyLab computer first');
%   disp('Exiting plx2NFile');
%   return;
%end
if exist([NFile '.mat' ],'file') 
   disp('WARNING: NFile already exists'); 
   %reply = input('Delete and Regenerate? Y/N [Y]: ', 's');
%    if isempty(reply)
%        reply = 'N';
%    end
%    if strcmp(reply,'Y') || strcmp(reply,'y')
%        delete([NFile '.mat' ]);
%    else
       disp(['loading ' NFile]);
       baseLineLength=[];
       load(NFile);
       numspikes=[];
       return;
  % end
end
load([DIRPATH '/' fname '_TrialStructure.mat'])
%delete([DIRPATH '/' fname '_TrialStructure.mat']);

[a,numtrials]=size(taskTrials);
plxfname = [PLXDIR fname '-01.plx'];
disp(['Reading Plexon File ' plxfname]);


if plxflag == 1     % 1 is plx, 2 is pl2    
    plxfname = [PLXDIR fname '-01.plx'];
    
    [n, evts9, sv] = plx_event_ts(plxfname, 9); % end of trial
    [n, evts8, sv] = plx_event_ts(plxfname, 8); % start collecting data
    [n, evts7, sv] = plx_event_ts(plxfname, 7); % stop stimulus
    [n, evts6, sv] = plx_event_ts(plxfname, 6); % stimulus on
    [n, evts5, sv] = plx_event_ts(plxfname, 5); % fixating 
    [n, evts4, sv] = plx_event_ts(plxfname, 4); % fix On
    [n, evts3, sv] = plx_event_ts(plxfname, 3); % start trial
else    
    plxfname = [PLXDIR fname '-01.pl2'];
    
    [evts] = PL2EventTs(plxfname, 9); % end of trial
    evts9=evts.Ts;
    [evts] = PL2EventTs(plxfname, 8); % start collecting data
    evts8=evts.Ts;
    [evts] = PL2EventTs(plxfname, 7); % stop stimulus
    evts7=evts.Ts;
    [evts] = PL2EventTs(plxfname, 6); % stimulus on
    evts6=evts.Ts;
    [evts] = PL2EventTs(plxfname, 5); % fixating 
    evts5=evts.Ts;
    [evts] = PL2EventTs(plxfname, 4); % fix On
    evts4=evts.Ts;
    [evts] = PL2EventTs(plxfname, 3); % start trial
    evts3=evts.Ts;
end

EndTrialT = evts9;
StartTrialT = PopOutGoodTimeStamp(evts3, EndTrialT);
FixDotT = PopOutGoodTimeStamp(evts4, EndTrialT);
FixOnT = PopOutGoodTimeStamp(evts5, EndTrialT);
StimOnT = PopOutGoodTimeStamp(evts6, EndTrialT);
StimOffT = PopOutGoodTimeStamp(evts7, EndTrialT);
TimeStamps = struct('StartTrialT',StartTrialT,'FixDotT',FixDotT,'FixOnT',FixOnT,...
    'StimOnT',StimOnT,'StimOffT',StimOffT,'EndTrialT',EndTrialT);

fprintf(' total # of trials %d \n', length(evts3));
fprintf(' # of trials with fixation on %d \n', length(EndTrialT));

% remove trials from the plexon data (evts9) that either
% missing the end of trial event or are not correct trials in the monkey
% lab data
% badTrials = [];
% for j=1:numtrials
%     if  ~isempty(trials(j).endOfTrial) & strcmp(trials(j).endOfTrial.eot,'correct')
%     else
%        badTrials= [badTrials j];
%     end;
% end;
% badTrials
% badTrials(find(badTrials > length(evts9))) = [];
% evts9(badTrials) = [];


% find goodtrials in fixation acquired
% goodtrials5 = 1:numtrials;
% goodtrials5(badTrials) = [];


goodtrialsll = [];
badTrials = [];
m=1;
for j=1:numtrials
    if  ~isempty(taskTrials(j).endOfTrial) && strcmp(taskTrials(j).endOfTrial.eot,'correct')      %trials(j).trialEnd==0 | trials(j).trialEnd==5   % trialEnd is 0 or 5 (code) when it is a "good trial" (in lablib)
        %if  ~isempty(taskTrials(j).endOfTrial) & strcmp(taskTrials(j).endOfTrial.eot,'computer')      %trials(j).trialEnd==0 | trials(j).trialEnd==5   % trialEnd is 0 or 5 (code) when it is a "good trial" (in lablib)
        goodtrialsll(m)=j;                            % goodtrials11 is the thing that actually counts - the fixation/stimulus parameters were all "good"
        stimLength(m)=taskTrials(j).endTrialState-taskTrials(j).startStimulus;
        baseLineLength(m)=taskTrials(j).startStimulus-taskTrials(j).fixationOn;
        m=m+1;
    else
        badTrials= [badTrials j];
    end;
end;
fprintf(' # of trials with correct outcome %d \n', length(goodtrialsll));


startStim = [];
fixationOn = [];
for j=1:length(goodtrialsll)
    startStim(end+1) = taskTrials(goodtrialsll(j)).startStimulus;
    fixationOn(end+1) = taskTrials(goodtrialsll(j)).fixating.time;
end;
% figure, plot(startStim-startStim(1)); hold on; plot(StimOnT-StimOnT(1),'r*')

% save badTrials
% evts9(badTrials) = [];


%     error('Different numbers of trials in Lablib and Plexon data');
%     return;
% end;
if (m~=length(evts9))
    disp('Different numbers of trials in Monkeylab and Plexon data');
    m = length(evts9);
end;


% create a column for each parameter in annular dots that contains the bin
% number so for direction min=45,step=45,num=8 you have an 8 element array
% with zero corresponding to 45 and 7 corresponding to 360 and the column
% in spike mat will indicate which bin each direction is


% go through all the trials --> same as length(evts9) and number of correct
% trials

NgoodTrial = length(goodtrialsll);
disp(['Number of good trials ' num2str(NgoodTrial)]);

% single out time stamp for good trials

LFPMat = cell(NgoodTrial,1);

StimCondAll = zeros(NgoodTrial,10);
for i = 1:length(goodtrialsll)
    if ~isstruct(taskTrials(goodtrialsll(i)).trialTaskValues)
        disp(['missing trialTaskValues at taskTrials(' num2str(goodtrialsll(i)) ')']);
        continue;
    else
        StimCond(1) = taskTrials(goodtrialsll(i)).trialTaskValues.direction;
        StimCond(2) = taskTrials(goodtrialsll(i)).trialTaskValues.directionBin;
        StimCond(3) = taskTrials(goodtrialsll(i)).trialTaskValues.speed;
        StimCond(4) = taskTrials(goodtrialsll(i)).trialTaskValues.speedBin;
        StimCond(5) = taskTrials(goodtrialsll(i)).trialTaskValues.size;
        StimCond(6) = taskTrials(goodtrialsll(i)).trialTaskValues.sizeBin;
        StimCond(7) = taskTrials(goodtrialsll(i)).trialTaskValues.contrast;
        StimCond(8) = taskTrials(goodtrialsll(i)).trialTaskValues.contrastBin;
        StimCond(9) = taskTrials(goodtrialsll(i)).trialTaskValues.currentGridIndex;
        StimCond(10) = taskTrials(goodtrialsll(i)).trialTaskValues.currentGridIndexCond;
    end
    StimCondAll(i,:) = StimCond;
end
Nval = zeros(10,1);
for i=1:10
    Nval(i) = length(unique(StimCondAll(:,i)));
end
parcol = find(Nval>1,1);
% if Nval(4)>1
    disp('Direction Tuning - Plaid');
    ParList = StimCondAll(:,parcol);
% end

ParVal = unique(ParList);
NParVal = length(ParVal);

% if length(EndTrialT)==NgoodTrial+1
%     % for whatever reason, the last trial was not saved by monkeylab
%     % we thus inferred the param for the last trial by checking which
%     % parameter has a missing trial
%     Nrepeats = length(EndTrialT)\nParVal;
%     Completion = zeros(NParVal,1);
%     for i=1:NgoodTrial
%         parindex = find(ParList(i)==ParVal);
%         Completion(parindex)=Completion(parindex)+1;
%     end
%     
%     LastTrial = find(Completion<Nrepeats);
%     ii = find(ParList==ParVal(LastTrial),1);
%     StimCondLastTrial = StimCondAll(ii,:);
%     StimCondAll = [StimCondAll; StimCondLastTrial];
%     NgoodTrial = NgoodTrial+1;
% end

for chi = 1:length(chnum)
    ch = chnum(chi);
%     LFPFile = [DIRPATH '/' fname num2str(ch) 'LFP'];
    for spi = 1:length(spnum)
        sp = spnum(spi);
        
        NFile = [DIRPATH '/' fname num2str(ch) num2str(sp) 'N'];
%         
%         RawDataFile = [RAWDPATH '/' fname num2str(ch) num2str(sp) 'NRawData'];
%         
%         if exist([NFile '.mat' ],'file')
%             disp('WARNING: NFile already exists, replacing');
%             %    reply = input('Delete and Regenerate? Y\n [Y]: ', 's');
%             %    if isempty(reply)
%             %        reply = 'N';
%             %    end
%             %    if strcmp(reply,'Y') || strcmp(reply,'y')
%             %        delete([NFile '.mat' ]);
%             %    else
%             %        disp(['loading ' NFile]);
%             %        load(NFile);
%             %        return;
%             %    end
%         end
        
%         % read neural data
         [nspk, spts] = plx_ts(plxfname, ch, sp);
%         if nspk==0
%             disp(' Spike signal does not exist ');
%         end
%         [numwaves, npw, wavetimes, spkwave] = plx_waves_v(plxfname, ch, sp);
%         
%         % analog signal (LFP)
%         [adfreq, nLFP, LFPstartT, fn, LFP] = plx_ad_v(plxfname, ch-1);
%         
%         if nLFP==0
%             disp(' Analog signal does not exist ');
%         end
%         
%         if nLFP==0 && nspk==0 % no data exist for the selected unit
%             continue;
%         end
%         % spike time relative to LFP recording
%         % sptinLFP = spts-LFPstartT;
%         % sptinLFP = round(sptinLFP*adfreq)+1-npw/4;
%         % LFP = DeSpikeLFP(LFP,spkwave,sptinLFP);
%         tLFP = LFPstartT+(0:length(LFP)-1)/adfreq;
%         
%         
%         % keyboard
%         RawLFP = struct('t',tLFP,'LFP',LFP);
%         RawData = struct('RawLFP',RawLFP,'RawSpike',spts,'SpkWave',spkwave,'TimeStamps',TimeStamps,...
%             'StimCond',StimCondAll);
%         save(RawDataFile, 'RawData');
%         disp(['saved raw data to ' RawDataFile])
%         
        
        % spki = randi(length(spts),200,1);
        % figure(1);clf;
        % for i=1:200
        %     plot((1:32)/10,spkwave(spki(i),:));hold on;
        % end
        %  figure(1);clf;
        % for i=1:200
        %     lfpsnip = LFP(sptinLFP(spki(i))+(-20:20));
        %     plot((-20:20)*1000/adfreq,lfpsnip-mean(lfpsnip));hold on;
        % end
        
        % count number of valid spikes
        count = 0;
        for i = 1:NgoodTrial
            % find all spike times which are between fix acquired and correct
            k=find(spts>FixOnT(i)&spts<EndTrialT(i));
            count =count+length(k);
        end
        
        spikeMat = zeros(count,12);
        c=1;
        for i = 1:NgoodTrial
            % find all spike times which are between fix acquired and correct
            k=find(spts>FixOnT(i)&spts<EndTrialT(i));
            for j = 0:length(k)
                if 0==j
                    % fake spike at -1000 insures that trial is represented; first entry will always be -1000 (in col 1)
                    spikeMat(c,1) = -1000;
                else
                    % time of spike is when the spike occurred subtrated by the
                    % "stimulus on" - these times represent time after stimulus
                    spikeMat(c,1) = floor((spts(k(j)) - StimOnT(i))*1000);
                end
                
                spikeMat(c,2) = i;
                spikeMat(c,3:12) = StimCondAll(i,:);
                c = c + 1;
            end
            
%             k = find( tLFP > FixOnT(i) & tLFP < StimOffT(i));
%             % align with stimulus onset
%             LFPMat{i} = struct('t',tLFP(k)-StimOnT(i),'LFP',LFP(k));
        end

        disp(['saved NFile to ' NFile])
        save(NFile,'spikeMat','stimLength','baseLineLength')
    end
%     disp(['saved LFPFile to ' LFPFile])
%     save(LFPFile,'LFPMat')
end
plx_close(plxfname);
return;
end