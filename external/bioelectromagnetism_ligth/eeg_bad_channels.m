function [p] = eeg_bad_channels(p,plot)

% eeg_bad_channels - Calculate various indices of bad voltage channels
%
% Useage: [p] = eeg_bad_channels(p,plot)
% 
% Returns OK channel index in p.elec.data.OK (Nelec,1)
% 
% Inputs are the eeg_toolbox p struct and a boolean 'plot'
% command to generate various plots.
% 
% Some of the calculations require the variance of the
% ERP potential.
% 
% 

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:50 $

% Licence:  GNU GPL, no express or implied warranties
% History:  07/2000, Darren.Weber_at_radiology.ucsf.edu
%           04/2001, Darren.Weber_at_radiology.ucsf.edu
%           -   modified inputs to p struct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% default initialisations
if ~exist('plot', 'var'), plot = 0; end

if ~exist('p','var'),
   [p] = eeg_toolbox_defaults('create');
end
if isempty(p.elec.data),
   [p] = elec_open(p);
end
if isempty(p.volt.data),
   [p] = eeg_open(p);
end

% Initialise bad_channels array
bad_channels = zeros(1,size(p.elec.data,1));

elec = str2num(char(p.elec.data.label))';
X = p.elec.data.x;
Y = p.elec.data.y;
Z = p.elec.data.z;
n_elec = length(X);

% Open output results text file
[path,name,ext] = fileparts(strcat(p.volt.path, filesep, p.volt.file));
outputfile = fullfile(path,strcat(name, '.', 'volt.bad'));
FID = fopen(outputfile,'w');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VARIANCE
% report max variance > 40 and mean variance >< (mean +/- SD mean_var)

if ~isempty(p.volt.var),
	
	max_var = max(p.volt.var);
	extreme_maxvar = or((max_var < 1),(max_var > 40));
	
	meanvar = mean(p.volt.var);
	Xmeanvar = ones(size(meanvar)) .* mean(meanvar);
	SDmeanvar = ones(size(meanvar)) .* std(meanvar);
	Critmeanvar = (2 .* SDmeanvar) + Xmeanvar;
	extreme_meanvar = meanvar > Critmeanvar;
	
	fprintf ('\nExtreme max variance (< 1 or > 40) or mean variance (Mean + 2*StDev):\n\n');
	fprintf (FID,'\nExtreme max variance (< 1 or > 40) or mean variance (Mean + 2*StDev):\n\n');
	
	extreme_var = or( extreme_maxvar, extreme_meanvar );
	bad_channels = bad_channels + extreme_var;
	
	out = [elec(extreme_var); maxvar(extreme_var); meanvar(extreme_var)];
	fprintf ('elec     max    mean\n',out);
	fprintf ('%4d  %6.2f  %6.2f\n',out);
	fprintf (FID,'elec     max    mean\n',out);
	fprintf (FID,'%4d  %6.2f  %6.2f\n', out);
	
	if plot,
		figure('name','Extreme Mean Variance'); hold on
		plot(meanvar,'g');
		plot(Xmeanvar,'b');
		plot(Critmeanvar,'r');
	end
else
    fprintf('EEG_BAD_CHANNELS: p.volt.var is empty.\n');
    extreme_var = zeros(1,n_elec);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DRIFT
% linear fit >> 0.0   &mean_volt_outlyers = mean(volt_outlyers');&   >= ???? xSD from mean linear fit of all elecs

p = eeg_linfit(p);

slope = p.volt.fitslope(1,:);

mean_slope = mean(slope);
std_slope = std(slope);
excess_slope = or( (slope > (mean_slope + 3 * std_slope)),(slope < (mean_slope - 3 * std_slope)));

%mean_int = mean(intercept);
%std_int = std(intercept);
%excess_intercept = or( (intercept > (mean_int + 3 * std_int)),(intercept < (mean_int - 3 * std_int)) );
%poor_fit = or(excess_slope, excess_intercept);

poor_fit = excess_slope;
bad_channels = bad_channels + poor_fit;

fprintf ('\nChannels with poor linear fit:\n\n');
fprintf (FID,'\nChannels with poor linear fit:\n\n');

fprintf ('mean slope = %10.4f\n\n', mean_slope);
fprintf (FID,'mean slope = %10.4f\n\n', mean_slope);
%fprintf ('mean intercept = %6.2f\n', mean_int);
%fprintf (FID,'mean intercept = %6.2f\n', mean_int);

out = [elec(poor_fit);slope(poor_fit)]; %intercept(poor_fit)];
fprintf ('elec   slope\n',out);
fprintf ('%4d  %10.4f\n',out);
fprintf (FID,'elec   slope\n',out);
fprintf (FID,'%4d  %10.4f\n', out);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CORRELATIONS & VOLTAGE VARIANCE
% calculates mean correlation of elec with surrounding elecs
% excludes elecs with extreme variance in surround elecs correlations

% identity matrix of electrode (column) correlations
COR = corrcoef(p.volt.data);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find nearest neighbour electrodes, excluding PAN and centroid/ref
% Calculate inter-electrode distances & nearest neighbour electrodes
% Would be nice to use the dsearchn function, but it only returns
% one nearest neighbour, rather than N neighbours.

% Check if distance file already calculated (speedy option)
[path,name,ext] = fileparts(strcat(p.elec.path, filesep, p.elec.file));
distances_file = fullfile(path,strcat(name,'.dist.mat'));
if exist(distances_file) == 2,
    load(distances_file);
else
    [dist_e, dist_d] = elec_distance(elec',X,Y,Z);
    clear X Y Z;
    save(distances_file,'dist_e','dist_d');
end

% obtain nearest neighbour electrode labels (numbers), in NN_e
[NN_e, NN_d] = elec_distance_nn(dist_e, dist_d);
clear dist_e dist_d NN_d;

% remove first column of NN_e because each electrodes nearest neighbour is itself.
NN_e = NN_e(:,2:end);

correlation = zeros(6,n_elec); % initialise nearest neighbour correlations

for e = 1:n_elec
    nn = 0;
    good_nn = 0;
    while good_nn < 6     % Get 6 good nearest neighbour correlations
        nn = nn + 1;
        nn_elec = NN_e(e,nn);              % get nn elec number
        [t,r] = strtok(char(nn_elec),','); % parse nn_elec and get
        nn_elec = str2num(r(2:end));       % value of nn_elec as number
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % check nn_elec is "good" so far => avoid contaminating good
        % channels with correlations from bad channels
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if extreme_var(nn_elec) == 0
            good_nn = good_nn + 1;
            correlation(good_nn,e) = COR(nn_elec, e);
        end
    end
end

fprintf ('\nChannels with correlations < 0.65 :\n\n');
fprintf (FID,'\nChannels with correlations < 0.65 :\n\n');

mean_corr = mean(abs(correlation));

neg_cor = mean_corr < 0;
small_cor = mean_corr < 0.65;       % 0.65 is reasonable correlation
bad_cor = or(neg_cor, small_cor);

bad_channels = bad_channels + bad_cor;

out = [ elec(bad_cor); mean_corr(bad_cor) ];
fprintf ('elec  correlation\n');
fprintf ('%4d  %11.2f\n', out);
fprintf (FID,'elec  correlation\n');
fprintf (FID,'%4d  %11.2f\n', out);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ALL BAD CHANNELS

p.elec.data.OK = (bad_channels == 0)';

fprintf ('\n\nAll Bad Channels :\n\n');
fprintf (FID,'\n\nAll Bad Channels :\n\n');

fprintf('%4d', elec(bad_channels > 0));
fprintf(FID,'%4d', elec(bad_channels > 0));

fprintf ('\n\n');
fprintf (FID,'\n\n');

% Plot voltage/variance data
if plot,
    time = p.volt.timeArray(:,1);
    figure('name','All Voltage');
    plot(time,p.volt.data);
    figure('name','Bad Voltage');
    plot(time,p.volt.data(:,(bad_channels > 0)) );
    figure('name','OK Voltage');
    plot(time,p.volt.data(:,(bad_channels == 0)));
    
    if ~isempty(p.volt.var),
        figure('name','All Variance');
        plot(time,p.volt.var);
        figure('name','Bad Variance');
        plot(time,p.volt.var(:, (bad_channels > 0)) );
        figure('name','OK Variance');
        plot(time,p.volt.var(:, (bad_channels == 0)));
    end
end

fclose('all');
return
