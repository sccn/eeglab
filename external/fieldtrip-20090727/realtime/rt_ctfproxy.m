function rt_ctfproxy(cfg)

% RT_CTFPROXY reads continuous data from shared memory on a CTF
% acquisition system and writes it to a FieldTrip buffer.
%
% The FieldTrip buffer is a network transparent server that allows the
% acquisition client to stream data to it. An analysis client can connect
% to read the data upon request. Multiple clients can connect simultaneously,
% each analyzing a specific aspect of the data concurrently.
%
% Since the CTF shared memory interface is only available on the
% acquisition machine itself, this function must run on the acquisition
% machine. The buffer to which the data is streamed is available through
% the network, so the actual analysis can be done elsewhere.
%
% Use as
%   rt_ctfproxy(cfg)
%
% The target to write the data to is configured as
%   cfg.target.datafile      = string, target destination for the data (default = 'buffer://localhost:1972')
%   cfg.target.dataformat    = string, default is determined automatic
%
% To stop this realtime function, you have to press Ctrl-C

% Copyright (C) 2008, Robert Oostenveld
%
% $Log: not supported by cvs2svn $
% Revision 1.2  2009/01/14 21:16:52  marvger
% changes related to realtime processing
%
% Revision 1.1  2008/10/24 08:51:38  roboos
% new implementation
%

% set the defaults
if ~isfield(cfg, 'target'),             cfg.target = [];                                  end
if ~isfield(cfg.target, 'datafile'),    cfg.target.datafile = 'buffer://localhost:1972';  end
if ~isfield(cfg.target, 'dataformat'),  cfg.target.dataformat = [];                       end % default is to use autodetection of the output format

% use another function to do the actual work
cfg.source.datafile   = 'shm://';
cfg.source.dataformat = 'ctf_shm';
rt_fileproxy(cfg);

