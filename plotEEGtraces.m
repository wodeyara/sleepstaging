function [h,chanTraces] = plotEEGtraces(chanTraces, refTraces)
    Fs = 1024;
    chanTraces = bsxfun(@minus, chanTraces, refTraces');
    
    [b, a] = butter(2, (.1)/(Fs/2), 'high');
    chanTraces = filtfilt(b,a,chanTraces);
    
    [b, a] = butter(2, (50)/(Fs/2), 'low');
    chanTraces = filtfilt(b,a,chanTraces);
    
    chanTraces = chanTraces(1:10:end,:);
    chanTraces = bsxfun(@plus, (2/3) * normalize(chanTraces,'zscore','robust'), 5*[1:size(chanTraces,2)]);
    
    h = plot((10/Fs:10/Fs:length(chanTraces)*10/Fs)/60, chanTraces,'Linewidth', 1.25);
    grid on
    xlabel('Time (mins)')
    ylabel('Signal')
    set(gca,'Fontsize', 14)
    legend({'F3', 'F4', 'C3', 'C4', 'O1', 'O2'}, 'Orientation', 'horizontal')
    ylim([2,32])
end