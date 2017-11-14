These scripts are wrapper and analysis scripts for FlySongClusterSegment.

The basic pipeline for analysis is as follows:

1) Generate a training set for your species or clade of interest, using many examples of each potential type of pulse song.
CHOOSE EITHER 1a OR 1b:

1a) If your data have low levels of background noise:
Generate a training set based on some subset of your data (~10 files).
As input to the wrapCreateTraining, you may:
-use the full path to a folder containing .wav files
-make a .txt file with each line representing the full path to a .wav file
-For each species or line in your dataset, create a .txt file listing all recordings from that species.
-...and save a list of these .txt files in the .txt file 'folderorlist'.
Run:
wrapCreateTraining(folderorlist);

OR
1b) If your data have more background noise:
Estimate noise levels for all recordings in 'songlist', saving those with maximum of the PSD below 'maxmaxP'.
Run:
GetSongNoiseLevels(songlist,maxmaxP,fscsdir);
Use the resulting _withmaxP.csv file to identify clean recordings (e.g., sdNoise < .05, maxP < 10).
Manually extract example pulse trains from these clean recordings and save as .wav files to 'wavfolder'.
Create a training set by randomizing these trains and saving as 'outname' with sampling frequency 'outfs'.
Run:
randomizeWavForTraining(wavfolder, outname, outfs);

THE FOLLOWING STEP IS NOT NECESSARY FOR MOST USERS, BUT ONLY FOR TESTING PERFORMANCE OF FSCS:
1c) (optional) Generate a validation set for optimizing FSCS parameters.
Manually extract example pulse trains from multiple recordings per species and save as .wav files to 'signalfolder'.
Manually extract example recording subsets with no pulse song and save to 'noisefolder'.
Create a validation set by randomizing these trains and saving as 'outname' with sampling frequency 'outfs'
Run:
randomizeAndCombineSignalNoiseWav(signalfolder, noisefolder, outname, outfs);


2) Generate templates using your training set and (optionally) validate
Determine the value to use for 'minnoise', or the minimum noise threshold for a recording.
Recommendations should have been provided when you ran wrapCreateTraining (1a) or GetSongNoiseLevels (1b).
You can also use the minimum value of 'convthreshold' from the .csv files resulting from those scripts.
Set 'trainfile' to the _trainingset.wav file generated in step 1.
If using only a training set, set 'valfile' to the same as 'trainfile' OR use an example .wav file for your species.
The default is set for automated specification of which templates are signal or noise, 
based on how many initial calls of each template fall below the noise threshold (ampthresh).
Run:
CreateAndAssignTemplatesOpts(trainfile,valfile,fscspath,minnoise); 

2b) Check the accuracy of template assignments, and modify options if necessary.
First, see whether you can tell which templates represent signal and which represent noise
based on inspection of the template histograms. Try openfig('examples/example_templates1.fig')
to see what shapes signal and noise templates usually take.
Next (optionally), examine the _SignalTemplatesAssigned.fig figure to determine whether true pulses have been called accurately.
If you observe a high false positive rate, re-run creation and assignment with a higher 'noisepost'.
If you observe a high false negative rate, re-run creation and assignment with a lower 'noisepost'.
If you observe an extremely high false positive rate, determine which of the templates is/are causing it; 
templates are plotted in this order: r--k--g--m--c--r-.k-.g-.m-.c-.r-*k-*g-*m--c-*r-^k-^g-^m-^c-^
If one of the templates appears to be mis-classified as signal, rather than noise, either
-load the '_outputCreateTemplates.mat' object and replace 'outputData.isNoise' with the correct categorization (0 for signal, 1 for noise).
-or re-run CreateAndAssignTemplatesOpts and categorize the noise template(s) with 'n' rather than 's' 
 (open the '_TemplateHistograms.fig' to see what the offending template(s) look(s) like).
If using automated specification of signal and noise templates, check the fraction '< noiseThreshold' in the '_TemplateHistograms.fig' figure.
Then re-run creation and assignment with an 'ampthresh' that would exclude this value but retain the signal templates.
If necessary, run:
CreateAndAssignTemplatesOpts(trainfile,valfile,fscspath,minnoise,...
    ampthresh,noisepost); 

2c) (very optional) Use your validation file to estimate false positive and false negative rates for a wide range of parameters.
Add a column called 'numpeaks' to the '.order.csv' file output by randomizeAndCombineSignalNoiseWav (this file will be 'categorytimes'). 
Manually count the number of pulses in each interval of the validation .wav file and enter these into the new column.
Choose a parameter (like noise_posterior_threshold or smoothingLength_noise) from 'makeParameterStructure' to vary.
Set the name of this parameter as 'varname' and give it a prefix to include in the output filenames ('varpref')
Set 'filelist' to a file where you would like to save the names of individual output files.
Set 'tempfile' to the '_outputCreateTemplates.mat' file produced by creating templates.
Set 'valsignaldir' and 'valnoisedir' to the directories where you have the example .wav files that make up the validation set.
Run:
for i=1:length(valstotest)
    SeparateSongsAssignTemplatesAndGetFPFN(tempfile,valsignaldir,...
        valnoisedir,categorytimes,valstotest(i),filelist,varname,varpref);
end
An example of this is in 'testMultipleNoisePostsSeparate.m'
Use the output of this to create a table and a plot of precision and recall across all tested values.
Set 'groupname' to any label you would like for this species group, and leave other variables as before.
Run:
function makePRcurve(filelist,varname,varpref,groupname);

2d) (optional) Use clustering to create sets of templates for different song types.
Re-run creation and assignment to group signal templates into 'numsongtypes' clusters.
Run:
CreateAndAssignTemplatesOpts(trainfile,valfile,fscspath,minnoise,...
    [],[],[],[],clusttemp,numsongtypes);
Check the clustering by looking at the provided plots.
Compare '_SignalTemplatesAssigned.fig' with '_clustered_SignalTemplatesAssigned.fig'.
Look at the '_TemplateHistograms_group{n}.fig' figures, particularly to check whether noise templates and signal templates cluster together.
If necessary, load the '_clustered_outputCreateTemplates.mat' file, and modify 'templateGroupings' and/or 'isNoiseTemplateGrouping'.


3) Assign templates and estimate song parameters for all recordings of interest (can be done in parallel).
Set 'songfile' to a .wav file of a recording and 'tempfile' to the '_outputCreateTemplates.mat' file
Set 'Species' to whatever species/line/group the recording comes from.
Set 'outdir' to the directory in which you want to save summary files and plots.
Set 'plotassign' to 'true' if you want to see the assignment of pulses for each recording.
Optional (otherwise the script uses the default values):
Read the comments in lines 46 - 51 of 'oneSongAssignAndSummarize' and choose options for 'summaxIPI','minno', and 'minpulse'
Run: 
oneSongAssignAndSummarize(songfile, tempfile, Species, ...
    outdir,[],[],[],plotassign);
OR
oneSongAssignAndSummarize(songfile, tempfile, Species, ...
    outdir, summmaxIPI, minno, minpulse, plotassign);
To run in parallel, create a list of recordings for a single species ('specieslist').
Run:
spID = fopen(specieslist);
spFiles = textscan(spID,'%s');
fclose(spID);
for s = 1:length(spFiles)
    oneSongAssignAndSummarize(spFiles(s), tempfile, Species, ...
    outdir, summmaxIPI, minno, minpulse, plotassign);
end

4) Summarize and plot song parameters for a set of recordings.
Make a list ('filelist') of all the output files for your song parameter of interest (e.g., ls outdir/*ipi*summarystats.csv > filelist).
Set 'plotlabel' to your preferred y-axis label.
Run:
summarizeSongStatDist(filelist, plotlabel, repopath);

