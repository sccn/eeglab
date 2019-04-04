% icademo() - a sample ICA analysis script using the ICA/ERP package 
%             of Matlab functions distributed via
%             http://www.sccn.ucsd.edu/eeglab
%
% Notes:
%   Reads ascii ERP data file:            pnas.adt
%   Reads example ascii ERP data files:   chan.nam   chan14.locs (=chan_locs)
%
% Author: Scott Makeig, SCCN/INC/UCSD, La Jolla, 5-18-97 
%
% See also: runica(), icavar(), icaproj(), icaact()

% Copyright (C) 5-18-97 Scott Makeig, SCCN/INC/UCSD, scott@sccn.ucsd.edu
%
% This file is part of EEGLAB, see http://www.eeglab.org
% for the documentation and details.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
% THE POSSIBILITY OF SUCH DAMAGE.

% This version tested on package version
% Added envproj(), used maxmap in compplot() 6-18-97 -sm
% Reworked to use actual sample ERP data 6-27-97 -sm
% Pruned again 7/18/97 -sm
% Adjusted envproj, plotproj, icaproj to remove datamean, distributing 
%     mean offset (if any) among component projections 7/23/97 -sm
% Fixed plotproj() call - it can't label arbitrry channel lists -sm 
% Adjusted envproj etc. 12/08/97 -sm
% Added timtopo() 1-16-98 -sm
% Added more formatting, fixed compplot and eegmovie plotting 10-9-99 -sm
% 12-19-00 revised icaproj args -sm
% 01-05-01 advertised the web tutorial, fixed chanproj(), improved eegmovie() -sm
% 01-12-01 removed 'sphere' arg from plotproj(), icaact() -sm
% 01-25-02 reformated help & license, added links -ad 

   chans  = 14;  % data channels (rows in data matrix)
   frames = 312; % frames per data epoch (columns in data submatrix)
   epochs = 2;   % epochs in simulated data
   srate  = 312.5; % data sampling rate in Hz
   limits = [0 995 0 0]; % epoch time limits in msec, y limits unspecified
   off    = [25 -25 0 0];   % successive figure offset in pixels
   icadefs  % run icadefs to see whether the user has changed the XXX reference.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Begin
%
fprintf('\n *************************************************************\n');
  fprintf(' * Demonstration script showing several ICA package analysis *\n')
  fprintf(' * and plotting routines applied to sample lowpass ERP data  *\n');
  fprintf(' *                                                           *\n');
  fprintf(' * Note: Use of these and other toolbox functions is now     *\n');
  fprintf(' *           documented in an online tutorial at:            *\n');
  fprintf(' *         http://www.cnl.salk.edu/~scott/tutorial/          *\n');
fprintf(' *************************************************************\n\n');

v=version;
if v(1) < '5'
  fprintf('Not all segments may work for Matlab version 4.\n')
end

% add path to access 
if ~exist('pnas.flt')
    p = which('eeglab');
    p = p(1:findstr(p,'eeglab.m')-1);
    addpath([ p 'sample_data' ] );
end

% name of channel locations file
chan_locs  = 'pnas_chan14.locs';
chan_locs2 = 'pnas_chan.locs';


figure
set(gcf,'Color',BACKCOLOR);
plot([-4 -4 4 4 -4],[-4.5 3 3 -4.5 -4.5],'w','Linewidth',3);hold on;
axis([-4 4 -5 5]);
t=text(0,2,'icademo()');
set(t,'FontSize',14,'HorizontalAlignment','center');
t=text(0,1,' Demonstration script showing several ICA toolbox ');
set(t,'FontSize',14,'HorizontalAlignment','center');
t=text(0,0,'analysis and plotting routines'); 
set(t,'FontSize',14,'HorizontalAlignment','center');
t=text(0,-1,'applied to sample lowpass ERP data');
set(t,'FontSize',14,'HorizontalAlignment','center');
t=text(0,-2.5,'See http://www.cnl.salk.edu/~scott/tutorial/');
set(t,'FontSize',14,'HorizontalAlignment','center');
t=text(0,-3.5,'for more details.');
set(t,'FontSize',14,'HorizontalAlignment','center');
axis('off')

if ~exist('icademoauto')
    fprintf('\n****> Hit any key to continue: '); pause; fprintf('\n\n'); %%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Load an ascii ERP data matrix 
%
fprintf('Loading sample data...\n');
   load pnas.adt -ascii
fprintf('Reshaping sample data...\n');
   data = [pnas(15:28,:) pnas(1:14,:)]; % lapses (will plot red), hits (blue)
   clear pnas
% Scale to uV
data = data/24.2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Lowpass filter (if available)
%
lowpass = 30; % filter band lowpass freq in Hz
if exist('firls')  % if signal processing toolbox is installed . . .
  fprintf('Lowpass filtering the data using eegfilt()...\n');
  % data = eegfilt(data,srate,lowpassHz,highpassHz,frames); 
   [data,filtwts] = eegfilt(data,srate,0,lowpass,frames); 
                                        % lowpass filter below 30 Hz
else                                    % filter each epoch similarly
  fprintf('Script eegfilt() requires the Matlab signal processing toolbox.\n');
end
pos = get(gcf,'Position'); % record figure position to use as a base for later figures

if ~exist('icademoauto')
    fprintf('\n****> Hit any key to continue: '); pause; fprintf('\n\n'); %%%
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Overplot the 2 312-point data epochs at each channel using plotdata.m
% Assume that the epoch timebase is 0 to 995 msec.
% Read 4-character channel identifiers from the file 'chan.nam'
% Use the default color order.
%
fprintf('\nOverplotting the two sample auditory ERP target responses\n')
fprintf('        of %d chans by %d frames using plotdata()...\n\n', ...
                               size(data,1),size(data,2));
figure('Position',pos+off/2); % #2a

plotdata(data,frames,[0 995 -10 10],'ERP Data', chan_locs2 ,0,'(2 conditions)');

if ~exist('icademoauto')
    fprintf('\n****> Hit any key to continue: '); pause; fprintf('\n\n'); %%%
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
fprintf('\nOverplotting the two sample auditory ERP target responses\n')
fprintf('        of %d chans by %d frames using plottopo()...\n\n', ...
                               size(data,1),size(data,2));
%
figure('Position',pos+off); % #2b

% plotttopo(data,chan_locs,frames,limits,title,axsize,colors) 

plottopo(data,chan_locs,frames,[0 995 -10 10],'ERP Data');

if ~exist('icademoauto')  
    fprintf('\n****> Hit any key to continue: '); pause; fprintf('\n\n'); %%%
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% run timtopo.m
%
fprintf('\nPlotting the Hit responses with topo maps at specified times ...\n')
figure('Position',pos+off*3/2); % #2c

% timtopo(data,chan_locs,[limits],[plottimes]','title',[plotchans],[voffsets]);

timtopo(data(:,1:frames),chan_locs,[0 995 -10 10],[250 320 390 500],'Target Hits');

if ~exist('icademoauto')
    fprintf('\n****> Hit any key to continue: '); pause; fprintf('\n\n'); %%%
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% run eegplot.m
%
fprintf(...
'Plotting the same data as two consecutive 1-sec epochs using eegplot().\n');
fprintf('Try using the on-screen and menu control elements...\n')

%  >> eegplot(data,srate,spacing,eloc_file,windowlength,title)

     %eegplotold(data,srate,0,chan_locs,1,'Two data epochs using eegplot()')
     eegplot(data,'srate',srate, 'spacing', 0,'eloc_file',readlocs(chan_locs), ...
				'winlength',1,'title','Two data epochs using eegplot()')
     %set(gcf,'Position',pos+2*off); % #3 - eegplot() makes its own figure

if ~exist('icademoauto')
    fprintf('\n****> Hit any key to continue: '); pause; fprintf('\n\n'); %%%
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
fprintf(...
'Plotting the same data as two consecutive 1-sec epochs using eegplotold()...\n');
fprintf('Thus, eegplotold() can be used to embed plots in larger figures.\n')
%
% >> eegplot('noui',data,srate,spacing,eloc_file,startsec,color)
%
figure('Position',pos+2.5*off); 
subplot(1,2,1)
  %eegplot('noui',data(:,1:312),'srate',srate,'spacing', 0,'eloc_file', chan_locs, 'color', {{'b'}} );
  eegplotold('noui',data(:,1:312),srate,0,chan_locs,0,'b');
  title('Lapses')
subplot(1,2,2)
  %eegplot('noui',data(:,1:313:624),'srate',srate,'spacing', 0,'eloc_file', chan_locs, 'color', {{'b'}} );
  eegplotold('noui',data(:,313:624),srate,0,chan_locs,0,'r');
  title('Hits')

if ~exist('icademoauto')
    fprintf('\n****> Hit any key to continue: '); pause; fprintf('\n\n'); %%%
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% run the ICA decomposition using runica.m
%
% Decompose the data using the ICA algorithm
% Specify the baseline vector to be sure that the baselines are zeroed.
%         data       = input data (size(data) = [chans,frames*epochs])
%         frames     = number of data points in each data epoch 
%         basevector = data point vector to use in zeroing baseline 
%
fprintf('Now decompose both epochs at once using the ICA algorithm, runica() ...\n\n');

if ~exist('icademoauto')
    fprintf('\n****> Hit any key to continue: '); pause; fprintf('\n\n'); %%%
end

[weights,sphere,activations,bias,signs,lrates] = runica(data);

fprintf('\nDone!\n');

if ~exist('icademoauto')
    fprintf('\n****> Hit any key to continue: '); pause; fprintf('\n\n'); %%%
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Compute the activation waveforms for the resulting ICA components using icaact.m
%
  datamean = mean(data')';
  activations = icaact(data,weights*sphere,datamean);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Reorder the first 10 components by increasing in-epoch latency to the
% point of maximum variance accounted for. Do not take into account 
% channel 14 in the variance calculation; only use channels 1:13.
%
fprintf('Re-sorting the resulting components by LATENCY ')
fprintf('          using compsort() ...\n\n');
fprintf('Note: runica() orders components by size - comport() not commonly needed!!\n')

% [windex,maxvar,maxframe,maxepoch,maxmap] ...
%                = compsort(data,weights,sphere,datamean, ...
%                                         frames,nlargest,chanlist);

  [windex,maxvar,maxframe,maxepoch,maxmap] = ...
                   compsort(data,weights,sphere,1,...
                                          frames,10,[1:13]);

fprintf('\nDone!\n');

if ~exist('icademoauto')
    fprintf('\n****> Hit any key to continue: '); pause; fprintf('\n\n'); %%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Plot the activation waveforms of all the components 
%        as reordered by compsort()
% Identify the components using the original component numbers.
%
fprintf('Plotting the component activation waveforms using plotdata()')
fprintf('using plotdata() ....\n\n');
figure('Position',pos+3*off); 

  plotdata(activations(windex,:),frames,limits,...
                 'ICA Component Activations', ...
                          windex,0,'(2 conditions)');

if ~exist('icademoauto')
    fprintf('\n****> Hit any key to continue: '); pause; fprintf('\n\n'); %%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Plot the contributions to the first epoch of data by the first reordered 
% 10 ICA components. Superimpose their projections on the original data.
%
fprintf('Plotting contributions of the 1st 10 ICA components ');
fprintf('using plotproj()...\n\n');
figure('Position',pos+4*off); 

%   [projdata] = plotproj(data,weights,compnums, ...
%                            title,limits,chanlist,channames,colors);
 [projdata] = plotproj(data(:,1:frames),weights*sphere,windex(1:10), ...
               'First epoch (comps. 1-10)',limits,[1:10], chan_locs2 );

if ~exist('icademoauto')
    fprintf('\n****> Hit any key to continue: '); pause; fprintf('\n\n'); %%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Plot the contributions to the first epoch of data by the 5 largest
% ICA components. Superimpose their projections on the original data.
%
fprintf('Plotting contributions of the 5 largest ICA components ');
fprintf('(windex order) using projtopo()...\n\n');
figure('Position',pos+4*off); 

%       >> [projdata] = projtopo(data,weights,[compnums],chan_locs,...
%                                 'title',[limits],colors,chans);

 [projdata] = projtopo(data(:,1:frames),weights*sphere,1:5,chan_locs,...
               'First epoch',limits);

if ~exist('icademoauto')
    fprintf('\n****> Hit any key to continue: '); pause; fprintf('\n\n'); %%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Plot the envelopes of the data epoch 1 and the first 10 components
%
fprintf('Plotting the [min,max] "envelopes" of the data ');
fprintf('and largest 6 components using envproj().\n(Data envelope is in white).\n\n');
figure('Position',pos+5*off); 

% >> [envdata] = envproj(data,weights,compnums, ...
%                        title,limits,chanlist,compnames,colors);
  [envdata] = envproj(data(:,1:floor(0.6*frames)),weights*sphere,1:6, ...
                  'Envelopes of data and largest 6 components (cond. 1)',...
                     [0 995*0.6 0 0]);

if ~exist('icademoauto')
    fprintf('\n****> Hit any key to continue: '); pause; fprintf('\n\n'); %%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
fprintf('Now plotting envelopes PLUS scalp maps using envtopo().\n');
figure('Position',pos+5.5*off); 
%
%  >> envtopo(data,weights,chan_locs,[limits],[compnums],...
%                                   'title',[plotchans],[voffsets]);

  envtopo(data(:,1:floor(0.6*frames)),weights*sphere,'chanlocs',readlocs(chan_locs),...
         'limits', [0 995*0.6],'compnums',[2 3 4 6],'title', 'Largest Components');

if ~exist('icademoauto')
    fprintf('\n****> Hit any key to continue: '); pause; fprintf('\n\n'); %%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Plot a close-up of the decomposition of epoch 2, channel 2.
% Focus on the time interval 250 to 400 msec.
%
  startframe = fix(0.250*srate)+1; % start 250 msec after epoch start
  endframe   = fix(0.400*srate);   % end 400 msec after epoch start

fprintf('Plotting contributions of the 1st 5 ICA components ');
fprintf('at Cz using chanproj().\nData is in white.\n\n');
figure('Position',pos+6*off); 

% [chandata] = chanproj(projdata,chan,ncomps,framelist,limits,titl,clrs);

  [chandata] = chanproj(projdata,1,5,startframe:endframe,...
                                 [0 995 0 0],'Epoch 1, channel Fz',0);
fprintf('\n Note that two components make up most of the response peak.\n');

if ~exist('icademoauto')
    fprintf('\n****> Hit any key to continue: '); pause; fprintf('\n\n'); %%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Plot the scalp maps of the first four reordered components 
% Read the electrode location file chan_locs
%
fprintf('Plotting 1st 4 ICA component scalp maps using compmap()...\n');
figure('Position',pos+7*off); 

% compmap (winv,eloc_file,index,title,pagesize,labels,printflag)
  compmap(maxmap,chan_locs,1:4,'Maps of 1st 4 ICA Components',...
                                                   [2 2],windex(1:4));

if ~exist('icademoauto')
    fprintf('\n****> Hit any key to continue: '); pause; fprintf('\n\n'); %%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Display the projected epoch-1 data for reordered ICA component windex(7) 
% together with a scalp map at its point of max variance
%
fprintf('Plotting the scalp maps and 1st-epoch projection of\n');
fprintf('   the button press component using compplot()...\n\n');
figure('Position',pos+8*off); 

% compplot(data,plotframe,chan_file,xstart,srate,title);

  frame = maxframe(7);
  projp3 = icaproj(data(:,frames+1:2*frames),weights*sphere,windex(7));

  compplot(projp3,frame,chan_locs,0,srate,'ICA Button Press Component');

if ~exist('icademoauto')
    fprintf('\n****> Hit any key to continue: '); pause; fprintf('\n\n'); %%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Make and display a brief movie of a small segment of the first data epoch
%
fprintf('Making a brief movie of an N200 portion of the 1st data epoch\n')
fprintf('            using eegmovie() & seemovie()...\n\n');
figure('Position',pos+9*off); 

% [Movie,Colormap] = eegmovie(data,srate,elec_locs,title,movieframes,minmax);

try,
  [Movie,Colormap] = eegmovie(data(:,313:412),srate,chan_locs,...
                                              'Demo of eegmovie()',75:85,0);
% Display the movie slowly, then 5 times forward/back
%
  seemovie(Movie,-5,Colormap); 
catch,
    disp('Problem generating movie');
end

if ~exist('icademoauto')
    fprintf('\n****> Hit any key to continue: '); pause; fprintf('\n\n'); %%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Run testica script using optimum parameters 
%    (nsources = ncomponents to separate)
% Sources are super-Gaussian 
%
fprintf('Last, running testica() - a simulation script for estimating the accuracy\n')
fprintf('         of ICA separation for a given data size and source distribution ...\n\n');

fprintf('Perfect separation would produce an identity matrix.\n');
fprintf('This is not possible here, since 14 sources are mixed into only 6 channels.\n');
fprintf('The figure will show that the 6 largest sources are still recovered clearly ...\n');

if ~exist('icademoauto')
    fprintf('\n****> Hit any key to continue: '); pause; fprintf('\n\n'); %%%
    testica(chans,epochs*frames,14,-0.05,1.4);
else 
    [testresult] = testica(chans,epochs*frames,14,-0.05,1.4);
end

% model 14 sources with a wide range of sizes & positive kurtosis

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(gcf,'Position',pos+10*off); 

if ~exist('icademoauto')
    fprintf('\n****> Hit any key to continue: '); pause; fprintf('\n\n'); %%%    
end; 
fprintf('\n *************************************************************\n');
fprintf(' *               For a tutorial on the toolbox:                *\n');
fprintf(' *                                                             *\n');
fprintf(' *          http://www.cnl.salk.edu/~scott/tutorial/           *\n');
fprintf(' *                                                             *\n');
fprintf(' *      For an overview of available scripts  >> help ica      *\n');
fprintf(' *      For information about using ICA, see  icafaq.html      *\n');
fprintf(' *      For further details and examples, see icabib.html      *\n');
fprintf(' *************************************************************\n\n')

return
