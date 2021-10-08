function [psds,freqs,T,chanData] = sleepSpectrogram(chans,refs)
% This function asks for an EDF file to load, extracts the channels of interest
% along with the references suggested and creates separate spectrograms for
% each of them. You can specify whether you want to create a single subplot
% or two figures. 
% INPUTS:
% chans as a cell array e.g. {'CZ','FZ'}
% refs as a cell array e.g. {'CII', 'PZ'}
%OUTPUTS:
% A figure window will pop up with crosshairs to allow clicking on the
% figure to get time points, this process repeats until 'stop' is typed
% into the command window. After each click it is necessary to go to the
% command window and hit return to allow for the user to click on the
% figure again. 
% psds: spectrograms for all the channels in chans organized as T x F x
% chans
% freqs: frequencies for the spectrograms
% T: time for the spectrograms
% chanData: the extracted chan data so that it can be saved out if desired
%%%%%% 
% LAST EDITED: 10/8/2021
% LAST AUTHOR: Anirudh W

if length(chans) ~= length(refs)
    disp('Please enter a reference for every channel or the first ref will be used for all')
    for i = 2:length(chans)
        refs{i} = refs{1};
    end
end

[file,path] = uigetfile('*.edf');
filename = [path,file];
hdr = read_edf(filename);
chanInds = []; 
% extract channel/ref locations
for i = chans
    tmp = find(strcmp(hdr.label,i));
    if isempty(tmp)
        disp(['channel name ', i, ' is incorrect, please correct and continue']) 
        return
    end
    chanInds = [chanInds, tmp];
end
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
    chanData(i,:)= read_edf(filename, hdr, 1,hdr.nSamples,chanInds(i));
end
disp('Data extracted!')
%%
ref =[];
Fs = 1024;
N = 1;
numSamps = Fs*N;
lastSamp =numSamps * floor(length(chanData)/numSamps);
chanData = chanData(:,1:lastSamp);

refData = chanData(length(chans)+1:end,:)';
chanData = chanData(1:length(chans),:)';

wo = 60/(Fs/2);  
bw = wo/35;
[b,a] = iirnotch(wo,bw);
dataFilt = filtfilt(b,a,chanData-refData);

[b, a] = butter(2, (.1)/(Fs/2), 'high');
dataFilt = filtfilt(b,a,dataFilt);

params.tapers = [.5,10,1]; % assuming 2 second segments this means our bin size is 1 Hz
params.Fs = 1024;
params.pad = -1;
params.fpass = [.5,50];

[psds,T, freqs]  = mtspecgramc(dataFilt,[10,2],params);
%% Plotting spectrograms
figure1 = figure('WindowState','maximized')

if length(chans)>1
    for i = 1:length(chans)
        psdsCurr = squeeze(psds(:,:,i));
        g(i) = subplot(length(chans),1,i)
        image(T,freqs,10*log10(psdsCurr)','CDataMapping','scaled');
        set(gca,'ydir','norm')
        set(gca, 'Layer','top');
        caxis([prctile(10*log10(psdsCurr(:)),50),prctile(10*log10(psdsCurr(:)),90)]);

        title([chans{i},'-', refs{i}])
        xlabel('Time (s)')
        set(gca,'Fontsize', 14)
        ylabel('Frequency')
        h = colorbar;
        ylabel(h,'Power (dB)')
        colormap(flipud(brewermap([],'spectral')))
        ylim([0 50]);
    end
    linkaxes(g,'x')
else  
    image(T,freqs,10*log10(psds)','CDataMapping','scaled');
    set(gca,'ydir','norm')
    set(gca, 'Layer','top');
    caxis([prctile(10*log10(psds(:)),60),prctile(10*log10(psds(:)),90)]);

    title([chans{1},'-', refs{1}])
    xlabel('Time (s)')
    set(gca,'Fontsize', 14)
    ylabel('Frequency')
    h = colorbar;
    ylabel(h,'Power (dB)')
    colormap(flipud(brewermap([],'spectral')))
    xlim([0 max(T)]);
    ylim([0 50]);
end

%% Defining sleep stages

flagEnd =0;
while ~flagEnd
    [x,y] = ginput(1);
    disp(['Time Point: ', num2str(x)]);
    z = input('Hit enter to continue identifying time points, type stop to stop','s');
    if strcmp(z,'stop')
        flagEnd = 1;
    end
end

