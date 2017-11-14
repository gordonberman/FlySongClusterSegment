function [out]=tybutter(in,f,g,type)

% function [out]=tybutter(in,f,g,type)
%
% Description:	A wrapper for the Matlab Signal Processing toolbox
%		Butterworth filter.  This function creates and applies a 4th
%		order, zero lag digital Butterworth filter using the "butter" and
%		"filtfilt" commands.  See the help files on those commands for
%		additional information.
%		
%	Inputs:
%		in = the data to be filtered (assumes a columnwise data matrix)
%		f = the cutoff frequency
%		g = the recording frequency
%		type = 'low', 'high', 'stop' or 'bandpass' for a low-pass, high-pass, band
%			stop or band pass filter (default = low).  If you wish to use the 
%     bandstop or bandpass filters filter you need to give f as [low high].
%
% Example:
%   filteredpts=tybutter(pts,40,250,'low');
%
% Author:	Ty Hedrick
% Changelog:	6/27/00 - initial version
% Changelog:	1/11/02 - added 'type' parameter
% Changelog:	8/14/03 - added 'stop' instructions
% Changelog:  2/27/06 - added 'pass' instructions

if exist('type')~=1
	type='low';
end

poles=4; % 4 pole filter

out=in;
[b,a]=butter(poles,f./(g/2),type); % run the Matlab "butter" command
out(:,:)=filtfilt(b,a,in(:,:)); % apply the filter
end