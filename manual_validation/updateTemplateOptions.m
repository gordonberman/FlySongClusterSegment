function updateTemplateOptions(tempfile,varargin)
    %allow modification of 'options' saved when creating templates
    %print out whether modifying or creating a new variable
    %throw a warning if the class or size of the variable is different from
    %the old one
    %varargin should be like so:
    %'min_noise_threshold',3e-3,'noise_posterior_threshold',0.2, etc.
    %Get info on all variables in addition to options, to save them back
    m = matfile(tempfile);
    options = m.options;
    for i = 1:(length(varargin)/2)
        if isfield(options,varargin{i*2-1})
            fprintf(1,' Modifying options variable %s\n',varargin{i*2-1});
            if class(options.(varargin{i*2-1})) ~= class(varargin{i*2})
                error([' Attempt to change variable class from %s to ' ...
                    '%s; aborting\n'],class(options.(varargin{i*2-1})),...
                    class(varargin{i*2}))
            end
            if sum(size(options.(varargin{i*2-1})) == ...
                    size(varargin{i*2-1})) ~= length(size(options.(varargin{i*2-1})))
                fprintf(1,' Warning: Size of %s is changing from...\n',...
                    varargin{i*2-1});
                disp(size(options.(varargin{i*2-1})))
                fprintf(1,' to...\n');
                disp(size(varargin{i*2-1}))
            end
        end
        options.(varargin{i*2-1}) = varargin{i*2};
    end
    %write out the new tempfile with a modified filename
    [tempdir,tempbase,~] = fileparts(tempfile);
	newtempfile = fullfile(tempdir,[tempbase '_mod' datestr(now,'mmddyy') '.mat']);
    foundversion = false;
    modnum = 1;
    while(foundversion == false)
        if exist(newtempfile,'file') ~= 2
            copyfile(tempfile, newtempfile);
            foundversion = true;
        else
            modnum = modnum + 1;
            newtempfile = fullfile(tempdir,[tempbase '_mod' num2str(modnum) ...
                '_' datestr(now,'mmddyy') '.mat']);
        end
    end
    fprintf(1,'Saving new options to %s\n',newtempfile);
    mnew = matfile(newtempfile,'Writable',true);
    mnew.options = options;
end