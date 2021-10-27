function h = plotEEGtraces(chanTraces, refTraces)
    Fs = 1024;
    chanTraces = bsxfun(@minus, chanTraces, refTraces');
    
    [b, a] = butter(2, (.1)/(Fs/2), 'high');
    chanTraces = filtfilt(b,a,chanTraces);
    
    [b, a] = butter(2, (40)/(Fs/2), 'low');
    chanTraces = filtfilt(b,a,chanTraces);
    
    chanTraces = chanTraces(1:10:end,:);
    chanTraces = bsxfun(@plus, normalize(chanTraces,'zscore','robust'), 3*[1:size(chanTraces,2)]);
    
    h = plot((10/Fs:10/Fs:length(chanTraces)*10/Fs), chanTraces,'Linewidth', 1.25);
    grid on
    xlabel('Time (s)')
    ylabel('Signal')
    set(gca,'Fontsize', 14)
    legend({'F3', 'F4', 'C3', 'C4', 'O1', 'O2'}, 'Orientation', 'horizontal')
    ylim([0,23])
end