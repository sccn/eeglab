% eeg_load_data_script - load ERP data
% Load data (680 rows x 124 cols)
eegpath = 'D:\data_emse\ptsdpet\link14hz\'; cd(eegpath);
ext = '_link14hz';
groups = {'c' 'p'};
cond = {'tac' 'oac'};
%cond = {'tac_oac'};
subs = [1:10];

% Define sample rate & epoch parameters
sample_rate = 2.5;
epoch_start = -200;
epoch_end = 1500;
points = 681;
timeArray = meshgrid(epoch_start:sample_rate:epoch_end,1)';
timeNonZero = find(timeArray);
timeZero = find(timeArray == 0);

for g = 1:length(groups),
    for i=1:length(subs),
        for c=1:length(cond),
            
            % Load data
            file = sprintf('%s%s%02d%s%s.txt',eegpath,groups{g},subs(i),cond{c},ext);
            data = sprintf(  '%s%02d%s%s',            groups{g},subs(i),cond{c},ext);
            if ~exist(data,'var')
                fprintf('loading %s\n',file);
                load(file);
                
                % Interpolate the Zero value
                volt = eval(data);
                InterpZero = interp1( timeArray(timeNonZero), volt, 0, 'cubic' );
                volt = [volt(1:timeZero-1,:); InterpZero; volt(timeZero:end,:)];
                eval(strcat(data,' = volt;'))
                
            else
                %dat = eval(data); dat = dat';
                %file = sprintf('%s%s%02d%s%s.dat',eegpath,groups{g},subs(i),cond{c},ext);
                %eeg_write_ascii(file,dat,'%12.4f');
            end
        end
    end
end

clear g i c file data timeNonZero volt InterpZero
