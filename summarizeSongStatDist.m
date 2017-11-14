function summarizeSongStatDist(filelist, plotlabel, repopath, whichsumm, ...
    litfile, dosort)
    %Take a list of the output from doStatHist and summarize across
    %all files (need 'species' column to determine group to plot)
    %Input:
    %filelist (character): a file listing the output files from doStatHist
    %to be summarized
    %plotlabel (character): y-axis label for plots
    %whichsumm (character): column name of the summary statistic to be used
    %for each fly
    %litfile (character, optional): a .csv file with literature values, of
    %the same format as the files in 'filelist'
    %repopath: path to directory containing scripts
    %To do: add more sub-directories for organization
    %To do: make some of the plotting and saving optional
    %To do: Can the mean + SD be added without the points using errorbarjitter? 
    %This may be a way to incorporate literature values.
    %To do: restore builtin histogram function
    if nargin < 3 || isempty(repopath)
        repopath = pwd; %assume scripts are in current directory or subdirectory
    end 
    addpath(genpath(repopath));
    
    %Set default for whichsumm
    if nargin < 4 || isempty(whichsumm)
        whichsumm = 'StatMaxDensity';
    end
    
    %Sort by name or leave in order of file listing
    if nargin < 6 || isempty(dosort)
        dosort = true;
    end
    
    %Save files with basename from filelist
    [outdir,outbase,~] = fileparts(filelist);
    
    %Check whether input file has been generated
    if exist(fullfile(outdir,strjoin({outbase, [whichsumm '.csv']},'_')),...
            'file') == 2
        fullstattable = readtable(fullfile(outdir,strjoin({outbase, ...
        [whichsumm '.csv']},'_')),'Delimiter',',');
    else
        fprintf(1,'Generating %s\n',strjoin({outbase, [whichsumm '.csv']},'_'));
        %Read in files
        folfileID=fopen(filelist);
        folders=textscan(folfileID,'%s');
        fclose(folfileID);
        allfiles=folders{1}; %scans as a cell array of cell arrays

        %Join stats in one table
        for f=1:length(allfiles)
            [~,songfile,~] = fileparts(allfiles{f});
            newstattable = readtable(allfiles{f},'Delimiter',',',...
                    'ReadVariableNames',true);
            if size(newstattable,2) < 10
                fprintf(1,'Skipping %s:\nNot enough columns\n',songfile);
                continue
            end
            if size(newstattable,1) > 1
                if ~ismember('WavFile',newstattable.Properties.VariableNames)
                    fnunknown = cell(size(newstattable,1),1);
                    fnunknown(:) = {'unknown'};
                    newstattable.WavFile = fnunknown;
                end
            else
                if ~ismember('WavFile',newstattable.Properties.VariableNames)
                    newstattable.WavFile = {songfile};
                end
            end
            %convert TemplateGroup column to cell array of strings
            if isnumeric(newstattable.TemplateGroup)
                newstattable.TemplateGroup = ...
                    cellstr(num2str(newstattable.TemplateGroup));
            end
            if exist('fullstattable','var') ~= 1
                fullstattable = newstattable;
            else
                try
                    fullstattable = vertcat(fullstattable, newstattable);
                catch MEnv
                    size(fullstattable)
                    size(newstattable)
                    fullstattable.Properties.VariableNames
                    newstattable.Properties.VariableNames
                    class(fullstattable.TemplateGroup)
                    class(newstattable.TemplateGroup)
                    allfiles{f}
                    rethrow(MEnv)
                end
            end
        end

        %Write table
        writetable(fullstattable,fullfile(outdir,strjoin({outbase, ...
            [whichsumm '.csv']},'_')));
    end
    
    %Eliminate outliers
    otorem = sum(fullstattable.(whichsumm) < 0);
    if otorem > 0
        fprintf(1,'   %.0f outliers removed for %s < 0\n', otorem, whichsumm);
        fullstattable = fullstattable(fullstattable.(whichsumm) > 0,:);
    end
    
    %Add literature values from litfile, if available
    %Read in already written file, if it exists
    if nargin > 4 && ~isempty(litfile) && exist(litfile,'file')==2
        fprintf(1,'Adding values from literature: %s\n',litfile);
        outbase = [outbase '_withlit'];
        withlittabname = fullfile(outdir,strjoin({outbase, [whichsumm '.csv']},'_'));
        if exist(withlittabname,'file') ~= 2
            littable = readtable(litfile,'Delimiter',',');
            %If a 'TemplateGroup' column exists, convert TemplateGroup (numeric
            %in original table, cell in literature table to identify references)
            if ismember('TemplateGroup', fullstattable.Properties.VariableNames)
                if isnumeric(fullstattable.TemplateGroup)
                    fullstattable.TemplateGroup = cellstr(num2str(fullstattable.TemplateGroup));
                end
            end
            fullstattable = vertcat(fullstattable,littable);
            writetable(fullstattable,withlittabname);
        else
            fullstattable = readtable(withlittabname,'Delimiter',',');
        end
    end
    
    %Make violin plot
    %For distributionPlot, need to use species column AND group column as 'groups'
    %Paste together separately for each species, if more than one group
    if ~ismember('speciesgroup', fullstattable.Properties.VariableNames)
        fullstattable.speciesgroup = fullstattable.Species;
        if ismember('TemplateGroup', fullstattable.Properties.VariableNames)
            if isnumeric(fullstattable.TemplateGroup)
                fullstattable.TemplateGroup = cellstr(num2str(fullstattable.TemplateGroup));
            end
            fullstattable.TemplateGroup = strtrim(fullstattable.TemplateGroup);
            uspeciesorig = unique(fullstattable.Species);
            for u=1:length(uspeciesorig)
                indexus = strfind(fullstattable.speciesgroup,uspeciesorig(u));
                wrows = find(not(cellfun('isempty', indexus)));
                if length(unique(fullstattable{wrows,'TemplateGroup'})) > 1
                    fullstattable(wrows,'speciesgroup') = ...
                        strcat(fullstattable{wrows,'Species'},{'\_'},...
                        fullstattable{wrows,'TemplateGroup'});
                end
            end
        end
    end
    %sort by speciesgroup
    %add option to customize sort order
    if dosort
        fullstattable = sortrows(fullstattable,'speciesgroup','ascend');
    end
    fprintf(1,'Plotting output by species (saving to\n%s)\n',...
        strjoin({outbase, whichsumm, 'violin.fig'},'_'));
    %figure
    distributionPlot(fullstattable.(whichsumm),'groups',fullstattable.speciesgroup);
    ylabel(plotlabel);
    fprintf(1,'Saving figure\n');
    savefig(gcf,fullfile(outdir,strjoin({outbase, whichsumm, 'violin.fig'},'_')),'compact');
    close;
    
    %Make error bar jitter plot
    %First convert the table to the format required for errorbarjitter
    %Need to make a table (or possibly an array) with the same number of 
    %columns as options in 'speciesgroup'
    forebjname = fullfile(outdir,strjoin({outbase, [whichsumm 'forjitter.csv']},'_'));
    if exist(forebjname,'file') ~= 2
        if dosort
            uspecies = unique(fullstattable.speciesgroup)
        else
            uspecies = unique(fullstattable.speciesgroup, 'stable')
        end
        newdataarray = NaN(height(fullstattable),length(uspecies));
        for u=1:length(uspecies)
            %find rows containing this species
            %wlarray = fullstattable.speciesgroup == uspecies(u);
            %wrows = find(wlarray(:,width(fullstattable)+1)); %This may work?
            indexus = strfind(fullstattable.speciesgroup,uspecies(u));
            wrows = find(not(cellfun('isempty', indexus)));
            %assign these values to the appropriate column in newdataarray
            newdataarray(wrows,u) = fullstattable.(whichsumm)(wrows);
        end
        uspeciesvn = strrep(uspecies, '\_', '_'); %replace '\_' in names
        ndaastable = array2table(newdataarray,'VariableNames',uspeciesvn);
        writetable(ndaastable, forebjname);
    else
        ndaastable = readtable(forebjname,'Delimiter',',','ReadVariableNames',true);
        newdataarray = table2array(ndaastable);
        uspeciesvn = ndaastable.Properties.VariableNames;
        uspecies = strrep(uspeciesvn, '_', '\_'); %reinsert '\_' into names
    end
    
    %Try errorbarjitter function
    %Add more options for this plot?
    %figure
    errorbarjitter_mod(newdataarray,gcf,'colheaders',uspecies,'barends','yes');
    ylabel(plotlabel);
    savefig(gcf,fullfile(outdir,strjoin({outbase, whichsumm, 'jitter.fig'},'_')),'compact');
    close;
    
end
