function lbls = printChanLabels()
% this function returns the labels for the channels
% lbls can then be printed out 

[file,path] = uigetfile('*.edf');
filename = [path,file];
hdr = read_edf(filename);
lbls = hdr.label;