function [psds,freqs,T,chanData] = sleepSpectrogram(chans,refs, startTime,endTime)
% This function asks for an EDF file to load, extracts the channels of interest
% along with the references suggested and creates separate spectrograms for
% each of them. You can specify whether you want to create a single subplot
% or two figures. 
% INPUTS:
% chans as a cell array e.g. {'CZ','FZ'}
% refs as a cell array e.g. {'CII', 'PZ'}
% startTime: value between 0 and length of EDF in hours specifying where to start reading the file
% endTime: value between 0 and length of EDF in hours (and greater than startTime) specifying where to end reading the file
%OUTPUTS:
% A figure window will pop up with crosshairs to allow clicking on the
% figure to get time points, this process repeats until 'stop' is typed
% into the command window. After each click it is necessary to go to the
% command window and hit return to allow for the user to click on the
% figure again. The outputted time will (on command window) will be
% relative to absolute time i.e. relative to start of the full file. 
% psds: spectrograms for all the channels in chans organized as T x F x
% chans
% freqs: frequencies for the spectrograms
% T: time for the spectrograms
% chanData: the extracted chan data so that it can be saved out if desired
%%%%%% 
% EDITED: 10/21/2021: added line traces
% EDITED: 10/25/2021: have added ability to select segment of time to
% extract
% EDITED: 10/27/2021: The EEG traces are made in same plot as spectrograms,
% they are time limited to 5 minutes. A second plot with full spectrogram
% is also plotted. Changed normalization and added way to jump to a
% specific time point.
% LAST AUTHOR: Anirudh W

if length(chans) ~= length(refs)
    disp('Please enter a reference for every channel or the first ref will be used for all')
    for i = 2:length(chans)
        refs{i} = refs{1};
    end
end

if ~(endTime > startTime)
    disp('Please ensure end time is later than start time')
    disp('setting end time to 1')
    endTime = 1;
end

addpath(genpath('../sleepstaging'))

[file,path] = uigetfile('*.edf');
filename = [path,file];
hdr = read_edf(filename);
Fs = hdr.Fs;
startTime = (startTime*3600*Fs)/hdr.nSamples;
endTime = (endTime*3600*Fs)/hdr.nSamples;

if startTime <0 | startTime >1 | endTime<0 | endTime >1
    disp('Please enter appropriate start and end times')
    disp('Assuming that the start time is 0')
    startTime =0;
    endTime = 1;
end

chanInds = []; 
% extract channel/ref locations
totalChans = length(chans);
chans(totalChans+1:totalChans+6) = {'F3', 'F4', 'C3', 'C4', 'O1', 'O2'}; % F3 F4 C3 C4 O1 O2
for i = chans
    tmp = find(strcmp(hdr.label,i));
    if isempty(tmp)
        disp(['channel name ', i, ' is incorrect, please correct and continue']) 
        return
    end
    chanInds = [chanInds, tmp];
end
totalRefs = length(refs);
refs(totalRefs+1:totalRefs+2) = {'T1','T2'};
for i = refs
    tmp = find(strcmp(hdr.label,i));
    if isempty(tmp)
        disp(['channel name ', i, ' is incorrect, please correct and continue']) 
        return
    end
    chanInds = [chanInds, tmp];
end
fprintf('\n <strong> Extracting data... </strong>\n (could take upto 10 minutes) ')
parfor i = 1:length(chanInds)
    chanData(i,:)= read_edf(filename, hdr, floor(startTime*hdr.nSamples) + 1,floor(endTime*hdr.nSamples),chanInds(i));
end
disp('Data extracted!')
%% Generate the spectrograms and filter data
ref =[];
Fs = 1024;
N = 1;
numSamps = Fs*N;
lastSamp =numSamps * floor(length(chanData)/numSamps);
chanData = chanData(:,1:lastSamp);

refData = chanData(size(chanData,1)-totalRefs-1:size(chanData,1)-totalRefs,:)';
chanData_red = chanData(1:size(chanData,1)-totalRefs-8,:)';

wo = 60/(Fs/2);  
bw = wo/35;
[b,a] = iirnotch(wo,bw);
dataFilt = filtfilt(b,a,chanData_red-refData);

[b, a] = butter(2, (.1)/(Fs/2), 'high');
dataFilt = filtfilt(b,a,dataFilt);

params.tapers = [.5,4,1]; % assuming 2 second segments this means our bin size is 1 Hz
params.Fs = 1024;
params.pad = -1;
params.fpass = [.5,50];

[psds,T, freqs]  = mtspecgramc(dataFilt,[4,2],params);

[b, a] = butter(2, (50)/(Fs/2), 'low');
dataFilt = filtfilt(b,a,dataFilt);
dataFilt = dataFilt(1:10:end,:);
%% Plotting full spectrograms
figure1 = figure('WindowState','maximized')

for i = 1:(totalChans)
    psdsCurr = squeeze(psds(:,:,i));
    g1(i) = subplot((totalChans),1,i)
    imagesc(T,freqs,10*log10(psdsCurr)');
    set(gca,'ydir','norm')
    set(gca, 'Layer','top');
    caxis([prctile(10*log10(psdsCurr(:)),25),prctile(10*log10(psdsCurr(:)),99)]);
    title([chans{i},'-', refs{i}])
    xlabel('Time (s)')
    set(gca,'Fontsize', 14)
    ylabel('Frequency')
    h = colorbar;
    h.Position(1) = h.Position(1) + .1;
    h.Position(3) = 0.01;
    h.Position(4) = 0.36;
    ylabel(h,'Power (dB)')
    colormap(flipud(brewermap([],'spectral')))
    ylim([0 50]);
end
linkaxes(g1,'x')
%% Plotting spectrograms, but delimited in time
figure2 = figure('WindowState','maximized')

if (totalChans)>1
    for i = 1:(totalChans)
        psdsCurr = squeeze(psds(:,:,i));
        g(i) = subplot((totalChans+1),1,i)
        imagesc(T,freqs,10*log10(psdsCurr)');
        set(gca,'ydir','norm')
        set(gca, 'Layer','top');
        caxis([prctile(10*log10(psdsCurr(:)),25),prctile(10*log10(psdsCurr(:)),99)]);

        hold on
        plot((10/Fs:10/Fs:length(dataFilt)*(10/Fs)),1.25*normalize(dataFilt(:,i),'zscore','robust')+40,'Linewidth', 1.2, 'Color', [.96,.68,.78,.75])
        
        title([chans{i},'-', refs{i}])
        xlabel('Time (s)')
        set(gca,'Fontsize', 14)
        ylabel('Frequency')
        h = colorbar;
        h.Position(1) = h.Position(1) + .1;
        h.Position(3) = 0.02;
        h.Position(4) = 0.2;
        ylabel(h,'Power (dB)')
        colormap(flipud(brewermap([],'spectral')))
        ylim([0 50]);
    end
    startChan = size(chanData,1)-6-totalRefs-2;
    endChan = size(chanData,1)-totalRefs-2;
    g(i+1) = subplot(totalChans+1,1,i+1);
    plotEEGtraces(chanData(startChan+1:endChan,:)', mean(chanData(endChan+totalRefs+1:endChan+totalRefs+2,:),1));
    linkaxes(g,'x')
else  
    g(1) = subplot(211);
    image(T,freqs,10*log10(psds)','CDataMapping','scaled');
    set(gca,'ydir','norm')
    set(gca, 'Layer','top');
    caxis([prctile(10*log10(psds(:)),25),prctile(10*log10(psds(:)),99)]);
    hold on
    plot((10/Fs:10/Fs:length(dataFilt)*10/Fs),1.25*normalize(dataFilt(:,i),'zscore','robust')+40,'Linewidth', 1.2, 'Color', [.96,.68,.78 .75])
    
    title([chans{1},'-', refs{1}])
    xlabel('Time (s)')
    set(gca,'Fontsize', 14)
    ylabel('Frequency')
    h = colorbar;
    h.Position(1) = h.Position(1) + .1;
    h.Position(3) = 0.01;
    h.Position(4) = 0.36;

    ylabel(h,'Power (dB)')
    colormap(flipud(brewermap([],'spectral')))
    xlim([0 max(T)]);
    ylim([0 50]);
    
    g(2) = subplot(212);
    startChan = size(chanData,1)-6-totalRefs-2;
    endChan = size(chanData,1)-totalRefs-2;
    h = plotEEGtraces(chanData(startChan+1:endChan,:)', mean(chanData(endChan+totalRefs+1:endChan+totalRefs+2,:),1));
    linkaxes(g,'x')
end
xlim([0,300]); % 5 mins at a time

%% Defining sleep stages

flagEnd =0;
while ~flagEnd
    [x,y] = ginput(1);
    x = x/3600 + floor(startTime*hdr.nSamples)/(3600*Fs);
    disp(['Time Point: ', num2str(x*3600), 's and ', num2str(x), ' hours']);
    hr_tmp = hours(x);
    [h,m,s] = hms(hr_tmp);
    disp(['Time Point: ',num2str(h),':', num2str(m),':', num2str(s)]);
    z = input('Hit enter to continue identifying time points, type stop to stop \n Or type T to enter a time point to jump to','s');
    if strcmp(z,'stop')
        flagEnd = 1;
    elseif strcmp(z,'T')
        tp = input('Enter time point in seconds to jump to');
        set(g(1),'xlim',[tp,tp+300]);
    end
end

