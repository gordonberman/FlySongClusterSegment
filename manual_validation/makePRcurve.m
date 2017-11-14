function makePRcurve(filelist,varlabel,specpref,groupname)
    %Make a precision-recall curve plot
    %Also make plots of the values of noisepost (or other variable) against
    %precision and recall.
    %Input:
    %filelist (string): a list of .mat files containing output from
    %EstimateFPFN
    %varlabel (string): the name of the variable to plot
    %specpref (string): the string preceding the value of this variable in
    %the filename (filename must look like this:
    %something_specpref{value}_something)
    %groupname (string): the clade of flies analyzed, for title/filename
    %Output:
    %fullfile(outdir,[basename 'PRPlot.fig/.jpeg']): a figure with the PR 
    %plot
    %fullfile(outdir,[basename 'PRby' varlabel '.fig/.jpeg']): a figure 
    %with the plot of variable versus precision and recall
    %fullfile(outdir,[basename '.withprecisionrecall.' varlabel '.csv']): 
    %a table with variable value, precision, and recall
    %To do: Sort by recall in order to avoid self-crossing on plot
    %To do: Enable plotting of multiple results files
    
    %replace _ with \_ for labels
    if ~isempty(strfind(varlabel,'_'))
        texvarlabel = strrep(varlabel,'_','\_');
    else
        texvarlabel = varlabel;
    end

    %Read in file list
    folfileID=fopen(filelist);
    folders=textscan(folfileID,'%s');
    fclose(folfileID);
    allfiles=folders{1}; %scans as a cell array of cell arrays
    
    %Make variable value, precision, and recall cell arrays to fill in
    varval = zeros(length(allfiles),1);
    precision = zeros(length(allfiles),1);
    recall = zeros(length(allfiles),1);
    fscore = zeros(length(allfiles),1);
    
    %Loop over the files, and extract the variable value
    for f = 1:length(allfiles)
        [~, fpfnname, ~] = fileparts(allfiles{f});
        fpfnsplit = strsplit(fpfnname,'_');
        IndexVar = strfind(fpfnsplit, specpref);
        Index = find(not(cellfun('isempty',IndexVar)));
        try
            varval(f) = str2num(strrep(fpfnsplit{Index},specpref,''));
        catch MEinputs
            fprintf(1,'Specpref found multiple places in filename:\n');
            fpfnsplit{Index}
            rethrow(MEinputs);
        end
        %Get fn and fd for recall and precision
        m = matfile(allfiles{f});
        precision(f) = 1 - m.fd;
        recall(f) = 1 - m.fn;
        fscore(f) = (2*precision(f)*recall(f))/(precision(f)+recall(f));
    end
    
    %Print varval(s) at maximum fscore
    indmaxf = find(fscore == max(fscore));
    fprintf(['The following values of noise_posterior_threshold yield\n' ...
        'maximum F-score (%.3f):\n'],max(fscore));
    disp(varval(indmaxf));
    
    %Combine variable value, recall, and precision in a table
    prtable = table(varval,precision,recall,fscore);
    prtable.Properties.VariableNames{'varval'} = varlabel;
    
    %Sort the table by variable value
    %New: For PR plot, sort the table by RECALL, since it makes the plot prettier
    prtable = sortrows(prtable,'recall','ascend');
    
    %Make the PR plot
    plot(prtable.recall,prtable.precision);
    title(['Precision-Recall curve for FSCS in ' groupname ', varying ' texvarlabel]);
    xlabel('Recall (TPR)');
    ylabel('Precision (PPV)');
    %New: make limits [0 1]... or not?
    %xlim([0 1]);
    %ylim([0 1]);
    [outdir,basename,~] = fileparts(filelist);
    savefig(fullfile(outdir,[basename 'PRPlot.fig']));
    saveas(gcf,fullfile(outdir,[basename 'PRPlot.jpeg']));
    close;
    
    %Make plots of variable value versus precision and recall
    %Sort the table by variable value for this plot
    prtable = sortrows(prtable,varlabel,'ascend');
    
    plot(prtable.(1),prtable.precision,prtable.(1),prtable.recall);
    title(['Precision/recall by ' texvarlabel ' for FSCS in ' groupname]);
    xlabel(texvarlabel);
    ylabel('Precision/Recall');
    %New: make limits [0 1]
    xlim([min(prtable.(1)) max(prtable.(1))]);
    ylim([0 1]);
    legend('Precision (PPV)', 'Recall (TPR)');
    savefig(fullfile(outdir,[basename 'PRby' varlabel '.fig']));
    saveas(gcf,fullfile(outdir,[basename 'PRby' varlabel '.jpeg']));
    close;
    
    %Save the table
    writetable(prtable,fullfile(outdir,[basename '.withprecisionrecall.' varlabel '.csv']));
end