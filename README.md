# sleepstaging
The main file to use is the sleepSpectrogram function. To this function two inputs need to be passed: channels and references. See the function help for examples of how this can be done. If there is a need to know what the channels are that exist in the dataset please use: printChanLabels without an input, it can also give an output if desired. 

Currently the file assumes there are the following EEG channels in the data to add as line traces: {'C4', 'C3', 'F4', 'F3', 'O1', 'O2'} These are first plotted as line traces referenced to the mastoid. Spectrograms of the channels requested are then generated and a crosshairs appears, click anywhere the first time. The spectrogram also has an overlaid line trace of the time series itself. 

The actual staging could proceed as follows: 
1. Crosshairs appear for the first time, click anywhere to get rid of them. The program then waits for the user to zoom/move the time series around. 
2. Hit Enter on the command line, this will bring back the crosshairs on the spectrogram, this will allow the user to click where desired to get a time estimate of that location. This is pasted on the command line. 
3. Repeat process across the spectrogram. 

Analysis details: 
The plotted time series have been downsampled 10 times, this shouldn't matter as they have also been filtered below 50 Hz
The spectrogram happens in 4s intervals with 2s overlaps. 
