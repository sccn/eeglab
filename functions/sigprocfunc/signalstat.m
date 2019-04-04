% signalstat()  -  Computes and plots statistical characteristics of a signal,
%                  including the data histogram, a fitted normal distribution,
%                  a normal ditribution fitted on trimmed data, a boxplot, and
%                  the QQ-diagram. The estimates value are printed in a panel and
%                  can be read as output. Optionally, a topographic map (see TOPOPLOT)
%                  can be plotted.
%                  The boxplot and the Kolmogorov-Smirnov test require the 
%                  MATLAB Statistics Toolbox.
%
% Usage:
%   >>  signalstat( data )
%   >>  signalstat( data, plotlab, dlabel, percent );
%   >>  [M,SD,sk,k,med,zlow,zhi,tM,tSD,tndx,ksh] = ...
%                 signalstat( data, plotlab, dlabel, percent, dlabel2, map, chan_locs );
%
% Inputs:
%   data        - data vector
%
% Optional inputs:
%   plotlab     - 1: default->plot  |  0: ->no plot
%   dlabel      - A label for the data ([]: default->'Potential [V]')
%   percent     - percentage of data to exclude for trimmed mean & SD ([]:default->5)
%                 Excluded is 'percent'/2 high % and 'percent'/2 low %
%   dlabel2     - A title label for the statistics table
%   map         - Data vector to be displayed as topographic map. If a single integer,
%                 only the corresponding electrode location is displayed
%   chan_locs   - name of an EEG electrode position file (See >> topoplot example for format).
%                 Can also be a structure (see >> help pop_editset)
%
% Outputs:
%   M,SD        - mean and standard deviation
%   sk,k        - skewness and excess kurtosis
%   med         - median
%   zlow,zhi    - low and high 'percent/2'-Percentile ('percent/2'/100-Quantile)
%   tM,tSD      - trimmed mean and SD, removing data<zlow and data>zhigh
%   tndx        - index of the data retained after trimming
%   ksh         - output flag of the Kolmogorov-Smirnov test at level p=0.05 
%                 0: data could be normally distributed; 1: data are not normally distributed 
%                 -1: test could not be executed 
%
% Author: Luca Finelli, CNL / Salk Institute - SCCN, 2 August 2002
%
% See also: 
%   pop_signalstat(), qqdiagram(), eeglab() 

% Copyright (C) 2002 Luca Finelli, Salk/SCCN, La Jolla, CA

% Note: 
% QQDIAGRAM IS EQUIVALENT TO PERCENTILE/PERCENTILE PLOT
% X = EEG.data(5,:); % data
% Y = randn(1, 1000); % gaussan random distribution
% figure; qqdiagram(X, Y,  2);
% figure; plot(prctile(X,2), prctile(Y,2));

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

function [M,SD,sk,k,med,zlow,zhi,tM,tSD,tndx,ksh] = signalstat( data, plotlab, dlabel, percent, dlabel2, map, chan_locs);

M=[]; SD=[]; sk=[]; k=[]; med=[]; zlow=[]; zhi=[]; tM=[]; tSD=[]; tndx=[]; ksh=[]; 

istats=1;	
%hs = help('stats');
% toolbx = ver;
% if isempty(hs) || all([~any(strcmpi({toolbx.Name},'statistics toolbox')), ~any(strcmpi({toolbx.Name},'statistics and machine learning toolbox'))])
%     disp('signalstat() note: the boxplot (not shown) requires the MATLAB Statistics Toolbox or Statistics and Machine Learning Toolbox');
%     istats=0;
% end

if (nargin<8 && nargin>5) && min(size(map))~=1
		error('signalstat(): the map input must be a vector')
end

if nargin<7 && nargin>5
	disp('signalstat(): no location file for the topographic map')
	help signalstat;
	return
end

if nargin < 6
	map = [];
end

if nargin < 5
	dlabel2 = '';
end

if nargin>3 
	if isempty(percent)
		percent=5;
	end
	if any(percent > 100) || any(percent < 0)
		error('signalstat(): percent must be between 0 and 100');
	end
end

if nargin < 4
	percent = 5;
end

if (nargin < 3 || isempty(dlabel))
	dlabel='Potential [V]';
end
		
if nargin < 2
  plotlab=1;
end;	
	
if ~isnumeric(plotlab)
	error('signalstat(): plotlab must be numeric');
end

if plotlab ~= 0 && plotlab ~= 1
		error('signalstat(): plotlab must be 0 or 1');
end
if nargin < 1
	help signalstat;
	return;
end;	

if ndims(data)>2
	error('signalstat(): data must be a vector (1-dim signal)')
end

if ~isreal(data)
	error('signalstat(): data cannot be complex')
end

fprintf('signalstat(): computing statistics...\n');

% Statistical characteristics
%----------------------------
pnts=length(data);   % number of data points
rg=max(data)-min(data);

M=mean(data);        % mean
med=median(data);    % median

vr=var(data);        % variance (N-1 normalized)
SD=std(data);        % standard deviation

if istats
	sk=skewness(data,0); % skewness (third central moment divided by
                         % the cube of the standard deviation)
    k=kurtosis(data,0)-3;  % kurtosis (fourth central  moment divided by 
                         % fourth power of the standard deviation)
else
	sk=NaN;
	k=kurt(data)-3;
end

% Checks on skewness and kurtosis
%--------------------------------
sklab='Distribution is symmetric';
if sk>0.01
	sklab='Distribution is right-skewed';
elseif sk < -0.01
	sklab='Distribution is left-skewed';
end

klab='';
if k>0.01
	klab='Distribution is super-Gaussian'; % i.e. kurtosis bigger then Gaussian
elseif k < -0.01
	klab='Distribution is sub-Gaussian';
end

% Estimates without the highest and lowest 'percent'/2 % of data
%---------------------------------------------------------------
pc=percent/100;

zlow = quantile(data,(pc / 2));   % low  quantile
zhi  = quantile(data,1 - pc / 2); % high quantile
tndx = find((data >= zlow & data <= zhi & ~isnan(data)));

tM=mean(data(tndx)); % mean with excluded pc/2*100% of highest and lowest values
tSD=std(data(tndx)); % trimmed SD

% Selected central tendency estimator
%------------------------------------
cte=M;

% Normal fit
%-----------
if istats
	alpha=0.05;          % 1-alpha confidence interval
	[muhat,sigmahat,muci,sigmaci] = normfit(data,alpha);
end

nbins=max(50,round(pnts/100));
[nel,binpos]=hist(data,nbins);
dx=binpos(2)-binpos(1);               % bin width

datafit=normpdf(binpos,cte,SD);       % estimated pdf
datafit=datafit*pnts*dx;

tdatafit=normpdf(binpos,cte,tSD);     % estimated pdf with trimmed SD
tdatafit=tdatafit*pnts*dx;

%datarnd=normrnd(cte,sigmahat,1,pnts); % synthetic data

if istats
	% Goodness-of-fit hypothesis test
	%--------------------------------
	kstail = 0; % 0 = 2-sided test
	
	CDF=normcdf(data,cte,sigmahat);        % estimated cdf
	
	[ksh,ksp,ksstat,kscv] = kstest(data,[data', CDF'],alpha,kstail); % Kolmogorov-Smirnov test
	
	kstestlab='Kolmogorov-Smirnov test: verified Gaussian';
	kscol=[0.2 1 0.2];
	if ksh
		kstestlab='Kolmogorov-Smirnov test: not Gaussian';
		kscol=[.7 .3 .3];
	end 
end

% Graphics
%-------------------------------------

if plotlab
  figure
  try, icadefs; set(gcf, 'color', BACKCOLOR); catch, end
  COLOR = [0.56 .66 .9];
  set(gcf,'NumberTitle','off','Name','Signal statistics -- signalstat()')
  fwidth=800;  % figure size in pixels
  fheight=600;
    
  funits=get(0,'Units');
  set(0,'Units','Pixel')
  scnsize=get(0,'ScreenSize');
  fpos=[round((scnsize(3)-fwidth)/2),round((scnsize(4)-fheight)/2),fwidth,fheight];
  set(0,'Units',funits)
   
  set(gcf,'Position',fpos)
 
  % Plotting the histogram
  %------------------------
  subplot(2,2,1)
  
  hist(data,nbins);  % 'XLim',[-125 125]
  uapos = get(gca,'Position');
  xlim  = get(gca,'XLim');
  set(gca,'FontSize',14)
  xlabel(dlabel)
  title('Data Histogram and Fitted Normal PDF')
 
  % HM=pnts*normpdf(M,muhat,sigmahat)/2; % FWHM height
  % plot([cte-sd cte+sd],[HM HM],'r--','LineWidth',2)

  % Overplotting a normal distribution
  %-----------------------------------
  hold on
  h1=plot(binpos,datafit,'c','LineWidth',2);
  set(gca,'XLim',xlim)
  
  % Overplotting a normal distribution from trimmed SD
  %---------------------------------------------------
  h2=plot(binpos,tdatafit,'y');
  ymin=get(gca,'YLim');
  plot([zlow zlow],[0 ymin(2)/20],'y','LineWidth',2) % low  percentile
  plot([zhi  zhi], [0 ymin(2)/20],'y','LineWidth',2) % high percentile
  set(gca,'XLim',xlim)

  % Overplotting a mean and zero line
  %----------------------------------
 % xmin=get(gca,'XLim');
 % ymin=get(gca,'YLim');
  plot([0 0],ymin,'k')
  h3=plot([cte cte],ymin,'r--','LineWidth',2);
  set(gca,'Color',COLOR,'XMinorTick','on','XLim',xlim)
  
  if istats
	  set(gca,'XTick',[])	  
  elseif ~istats && strcmp(dlabel,'Potential [V]')
	  set(gca,'XTick',[-125, -75, -25, 0, 25, 75,  125],...
			  'XTickLabel',['-125' ; ' -75' ; ' -25' ; '  0 ' ; ' 25 ' ; ' 75 ' ; ' 125'])
  end
  
  if strcmp(dlabel,'Potential [V]')
	  set(gca,'XLim',[-125 125])
  end

  set(gca,'FontSize',10)
  H=[h1 h2 h3];
  legend(H,'Gaussian fit','Trimmed G.fit','Mean') 
  legend boxoff
  
  zoom off

  % Boxplot
  %--------------------
  if istats
	  subplot(2,2,3)
	  
	  boxplot(data,1,'+',0,1.5)
	  lapos=get(gca,'Position');
	  set(gca,'Position', [uapos(1) uapos(2)-uapos(4)/3 uapos(3) uapos(4)/3])
	  nlapos=get(gca,'Position');
	  hold on
	  ymin2=get(gca,'YLim');
	  plot([0 0],[0 ymin2(2)],'k')
	  plot([cte cte],[0 ymin(2)],'r--','LineWidth',2)
	  set(gca,'FontSize',14,'XMinorTick','on') 
      set(gca,'XLim',xlim)
	  
      if strcmp(dlabel,'Potential [V]')
		  set(gca,'XTick',[-125 -75 -25 0 25 75  125],...
				  'XTickLabel',['-125' ; ' -75' ; ' -25' ; '  0 ' ; ' 25 ' ; ' 75 ' ; ' 125'],...
				  'XLim',[-125 125])
	  end
	 
	  xlabel(dlabel)
	  ylabel('')
	  zoom off
  end
  
  % QQ plot
  %--------
  subplot(2,2,2)
  
  qqdiagram(data)
  apos=get(gca,'Position');
  set(gca,'Position', [1-uapos(1)-uapos(3) uapos(2)-uapos(4)/3 (uapos(3)) (uapos(4)+uapos(4)/3)])
  set(gca,'XTick',[-4 -2 0 2 4])
  xmin=get(gca,'XLim');
  ymin=get(gca,'YLim');
  hold on
  plot([xmin(1) xmin(1)+diff(xmin)/20],[zlow zlow],'y-','LineWidth',2)
  plot([xmin(1) xmin(1)+diff(xmin)/20],[zhi  zhi] ,'y-','LineWidth',2)
  set(gca,'XLim',xmin);
  %plot([0 0],ymin,'k--')
  set(gca,'FontSize',14)
  xlabel('Standard Normal Quantiles [Std.Dev.]')
  if strcmp(dlabel,'Potential [V]')
	  ylabel('Ordered Observations [V]')
  elseif strcmp(dlabel,'Component Activity')
	  ylabel('Ordered Observations [rel. V]')
  else
	  ylabel('Ordered Observations')
  end
  
  title('QQ Plot (Data vs Standard Normal)')
  set(gca,'Color',COLOR)
  
  % TOPO plot   
  %---------
  if (~isempty(map))
	  sbplot(7,9,6)
	  % th=axes('Position',[]);
	  % subplot('Position',[.10 .86 .20 .14]); 
	  fprintf('signalstat(): plotting a topographic map...\n');
	  if length(map) == 1
		  topoplot(map,chan_locs,'electrodes','off', ...
				   'style', 'blank', 'emarkersize1chan', 10);
	  else
		  topoplot(map,chan_locs,'electrodes','off');
	  end
	  axis('square')
  end 

  % Color schemes
  %--------------
  nero    = [0 0 0];
  bordeau = [0.7 0.3 0.3];
  rosso   = [1 0 0];
  roschi1 = [1 .3 0.4];
  giallo1 = [1 .9 0];
  arancio = [1 .5 .3];
  verchi1 = [0.2 1 0.2];
  verchi2 = [0.7 1 0.7];
  verchi3 = [0.4 1 0.4];
  verscu1 = [0.1 0.7 0.2];
  bluchi1 = [0.4 0.9 1];
  bluchi2 = [.2 .5 .7];
  grichia = [.95 .95 .95];
  
  % Data axis
  %------------
  bgcolor=COLOR;
  
  dah = axes('Position',[uapos(1) .1 1-2*uapos(1) .25]);
  set(dah,'Box','on','Color',bgcolor,'XTick',[],'YTick',[],'FontName','Courier','FontSize',12,'FontWeight','demi')
  title(dlabel2,'FontWeight','bold','FontSize',14,'FontName','Arial');
  text(0.05,0.9,['Mean:          ' num2str(M,3)]  ,'Color',bordeau,'FontName','Courier','FontSize',12,'FontWeight','demi')
  text(0.05,0.8,['Trimmed mean:  ' num2str(tM,3)] ,'Color',giallo1,'FontName','Courier','FontSize',12,'FontWeight','demi')

  text(0.05,0.5,['Standard dev.: ' num2str(SD,4)] ,'Color',bordeau,'FontName','Courier','FontSize',12,'FontWeight','demi')
  text(0.05,0.4,['Trimmed st.d.: ' num2str(tSD,4)],'Color',giallo1,'FontName','Courier','FontSize',12,'FontWeight','demi')
  text(0.05,0.3,['Variance:      ' num2str(vr,4)] ,'Color',giallo1,'FontName','Courier','FontSize',12,'FontWeight','demi')
  text(0.05,0.2,['Range:         ' num2str(rg,4)] ,'Color',giallo1,'FontName','Courier','FontSize',12,'FontWeight','demi') 
  text(0.05,0.1,['Data points:   ' num2str(pnts)] ,'Color',bordeau,'FontName','Courier','FontSize',12,'FontWeight','demi')
   
  text(0.4,0.9,[num2str(percent/2/100,'%1.3f') '-quantile: ' num2str(zlow,3)] ,'Color',verchi2,...
	   'FontName','Courier','FontSize',12,'FontWeight','demi')
  text(0.4,0.8,['0.5  -quantile: ',num2str(med,3),'  (median)'],'Color',verchi1,'FontName','Courier','FontSize',12,'FontWeight','demi')
  text(0.4,0.7,[num2str((100-percent/2)/100,'%1.3f') '-quantile:  ' num2str(zhi,3)] ,'Color',verchi2,...
	   'FontName','Courier','FontSize',12,'FontWeight','demi')
  
  text(0.4,0.3,['Excess kurtosis: ' num2str(k, 3) ' (near 0 if Gaussian)'] ,'Color',verchi1,...
				'FontName','Courier','FontSize',12,'FontWeight','demi')
  text(0.4,0.2,klab,'Color',verchi2,'FontName','Courier','FontSize',12,'FontWeight','demi')
 
  if istats
	  text(0.4,0.5,['Skewness: ' num2str(sk,3) ' (near 0 if Gaussian)'] ,'Color',verchi1,...
		   			'FontName','Courier','FontSize',12,'FontWeight','demi')
	  text(0.4,0.4,sklab,'Color',verchi2,'FontName','Courier','FontSize',12,'FontWeight','demi')
      text(0.4,0.1,kstestlab,'Color',kscol,'FontName','Courier','FontSize',12,'FontWeight','demi')
  end
  axcopy;
end

%--------------------------------------------------
% clone of the normpdf function of the stat toolbox
function fitvals = normpdf(myvals,mymean,mystd)
if nargin < 3,
    mystd = 1;
end
if nargin < 2;
    mymean = 0;
end
if length(mymean) < length(myvals)
	tmpmean = mymean;
	mymean = zeros(size(myvals));
	mymean(:) = tmpmean;
end
if length(mystd) < length(myvals)
	tmpmean = mystd;
	mystd = zeros(size(myvals));
	mystd(:) = tmpmean;
end
mymean(1:10);
mystd(1:10);

fitvals = zeros(size(myvals));
tmp = find(mystd > 0);
if any(tmp)
    myvalsn = (myvals(tmp) - mymean(tmp)) ./ mystd(tmp);
    fitvals(tmp) = exp(-0.5 * myvalsn .^2) ./ (sqrt(2*pi) .* mystd(tmp));
end
tmp1 = find(mystd <= 0);
if any(tmp1)
    tmp2   = NaN;
    fitvals(tmp1) = tmp2(ones(size(tmp1)));
end
