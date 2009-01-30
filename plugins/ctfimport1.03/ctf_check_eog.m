
ctf = ctf_read_res4('/data/meg/WA_singleTrials/WA_HiN_alltrials_1250msec.ds')
ctf = ctf_read_markerfile([],ctf)
ctf.class = ctf_read_classfile(ctf)

% -----------------------------
% This section is based on prior knowledge of the artifact names used to
% identify VEOG and HEOG artifacts during ctf preprocessing

classNames = {ctf.class.data(:).name};
EYEBLINKaIndex = strmatch('EYEBLINKa',classNames);
EYEBLINKaTrials = ctf.class.data(EYEBLINKaIndex).trials;

EYEBLINKdIndex = strmatch('EYEBLINKd',classNames);
EYEBLINKdTrials = ctf.class.data(EYEBLINKdIndex).trials;

EYEMVMTaIndex = strmatch('EYEMVMTa',classNames);
EYEMVMTaTrials = ctf.class.data(EYEMVMTaIndex).trials;

EYEMVMTdIndex = strmatch('EYEMVMTd',classNames);
EYEMVMTdTrials = ctf.class.data(EYEMVMTdIndex).trials;

% find all unique artifact trials
EYEBLINKTrials = unique([EYEBLINKaTrials, EYEBLINKdTrials]);
EYEMVMTTrials  = unique([EYEMVMTaTrials, EYEMVMTdTrials]);
EOGBADTrials = unique([EYEBLINKaTrials, EYEBLINKdTrials,...
                    EYEMVMTaTrials, EYEMVMTdTrials]);

allTrials = 1:ctf.setup.number_trials;

thresholdDetect.veog.amp = 50e-6; % 50 uV
thresholdDetect.veog.der = 25e-3; % 25 mV/sec
thresholdDetect.heog.amp = 25e-6; % 25 uV
thresholdDetect.heog.der = 25e-3; % 25 mV/sec

% ------------------------------

cd(ctf.folder)

% -----------------------------
% This section here is based on prior extraction of the STIM, EEG121 and
% EEG122 channels, it could be replaced with the ctf_read_meg4 function
% to do this directly

stimeegFile = fullfile(ctf.folder,'WA_HiN_alltrials_1250msec_stimeeg.mat');
stimeegData = load(stimeegFile);
ctf.data = stimeegData.data;
ctf.sensor.label = stimeegData.labels
% -----------------------------

stim = strmatch('STIM',ctf.sensor.label);
veog = strmatch('EEG121',ctf.sensor.label);
heog = strmatch('EEG122',ctf.sensor.label);

blPoints = ctf.setup.pretrigger_samples;
timeSec = ctf.setup.time_sec;
epochs = size(ctf.data,3);


% For continuous data, use arbitrary epoch lengths for viewing
%i = 1;
%epochs(i,:) = 1:10000;
%timeSec = epochs(i,:) * ctf.setup.sample_sec;
%while (i * 10000) < size(ctf.data,1),
%    epochs(i,:) = (i * 10000) + epochs(1,:);
%    i = i + 1;
%end

fig = figure;

for e = allTrials,
  
  set(fig,'Name',sprintf('Trial %04d',e))
  
  veogData = ctf.data(:,veog,e);
  heogData = ctf.data(:,heog,e);

  veogBaseline = mean(veogData(1:blPoints));
  heogBaseline = mean(heogData(1:blPoints));

  veogDataBaselined = veogData - veogBaseline;
  heogDataBaselined = heogData - heogBaseline;

  %plot(timeSec,[veogDataBaselined,heogDataBaselined]);
  %legend('VEOG','HEOG')

  % EEG differential must be normalized by the sample rate (sec)
  veogDataDiff = [0; diff(veogDataBaselined)] / ctf.setup.sample_sec;
  heogDataDiff = [0; diff(heogDataBaselined)] / ctf.setup.sample_sec;
  
  % check if this is an artifact trial
  EYEBLINKaTrial = find(EYEBLINKaTrials == e);
  EYEBLINKdTrial = find(EYEBLINKdTrials == e);
  EYEMVMTaTrial  = find(EYEMVMTaTrials == e);
  EYEMVMTdTrial  = find(EYEMVMTdTrials == e);
  
  % -------------------------------
  % VEOG check

  
  
  % amplitude threshold
  t = thresholdDetect.veog.amp; % 80 uV
  % detect values above threshold
  eyeblinkA = find(abs(veogDataBaselined) > t);
  eyeblinkWave = veogDataBaselined * 0;
  eyeblinkWave(eyeblinkA) = max(abs(veogDataBaselined));

  veogAxisAmp = subplot(2,2,1);
  plot(timeSec, [veogDataBaselined, eyeblinkWave]); hold on
  legend('VEOG', 'VEOG Blink (Amp)')
  plot(timeSec, [veogData*0, (veogData*0)+t, (veogData*0)-t],'k:'); hold off
  ylabel('Volt')
  xlabel('Time (sec)')
  set(veogAxisAmp,'XLim',[ctf.setup.start_sec, ctf.setup.end_sec]);
  ampLimit = thresholdDetect.veog.amp + (thresholdDetect.veog.amp * 0.5);
  set(veogAxisAmp,'YLim',[-ampLimit, ampLimit])
  set(veogAxisAmp,'YTick',-ampLimit:10e-6:ampLimit)
  if EYEBLINKaTrial,
    set(veogAxisAmp,'color',[1 0.75 0.75])
  else
    set(veogAxisAmp,'color',[1 1 1])
  end
  
  
  
  % detect derivatives above 20 mV/sec
  t = thresholdDetect.veog.der;
  eyeblinkD = find(abs(veogDataDiff) > t);
  eyeblinkWave = veogDataDiff * 0;
  eyeblinkWave(eyeblinkD) = max(abs(veogDataDiff));
  
  veogAxisDer = subplot(2,2,3);
  plot(timeSec,[veogDataDiff, eyeblinkWave]); hold on
  legend('diff(VEOG)/sec', 'VEOG Blink (Der)')
  plot(timeSec, [veogData*0, (veogData*0)+t, (veogData*0)-t],'k:'); hold off
  ylabel('Volt / sec')
  xlabel('Time (sec)')
  set(veogAxisDer,'XLim',[ctf.setup.start_sec, ctf.setup.end_sec]);
  ampLimit = thresholdDetect.veog.der + (thresholdDetect.veog.der * 0.5);
  set(veogAxisDer,'YLim',[-ampLimit, ampLimit])
  %set(veogAxisDer,'YTick',-ampLimit:10e-6:ampLimit)
  if EYEBLINKdTrial,
    set(veogAxisDer,'color',[1 0.75 0.75])
  else
    set(veogAxisDer,'color',[1 1 1])
  end
  
  % -------------------------------
  % HEOG check

  
  
  % amplitude threshold
  t = thresholdDetect.heog.amp;
  % detect values above threshold
  eyeMoveA = find(abs(heogDataBaselined) > t);
  eyeMoveWave = heogDataBaselined * 0;
  eyeMoveWave(eyeMoveA) = max(abs(heogDataBaselined));
  
  heogAxisAmp = subplot(2,2,2);
  plot(timeSec, [heogDataBaselined, eyeMoveWave]); hold on
  legend('HEOG', 'HEOG Move (Amp)')
  plot(timeSec, [heogData*0, (heogData*0)+t, (heogData*0)-t],'k:'); hold off
  ylabel('Volt')
  xlabel('Time (sec)')
  set(heogAxisAmp,'XLim',[ctf.setup.start_sec, ctf.setup.end_sec]);
  ampLimit = thresholdDetect.heog.amp + (thresholdDetect.heog.amp * 0.5);
  set(heogAxisAmp,'YLim',[-ampLimit, ampLimit])
  set(heogAxisAmp,'YTick',-ampLimit:10e-6:ampLimit)
  if EYEMVMTaTrial,
    set(heogAxisAmp,'color',[1 0.75 0.75])
  else
    set(heogAxisAmp,'color',[1 1 1])
  end


  
  % HEOG movements may accelerate at about 20 mV/sec
  t = thresholdDetect.heog.der;
  eyeMoveD = find(abs(heogDataDiff) > t);
  eyeMoveWave = heogDataDiff * 0;
  eyeMoveWave(eyeMoveD) = max(abs(heogDataDiff));

  heogAxisDer = subplot(2,2,4);
  plot(timeSec, [heogDataDiff, eyeMoveWave]); hold on
  legend('diff(HEOG)/sec', 'HEOG Move (Der)')
  plot(timeSec, [heogData*0, (heogData*0)+t, (heogData*0)-t],'k:'); hold off
  ylabel('Volt / sec')
  xlabel('Time (sec)')
  set(heogAxisDer,'XLim',[ctf.setup.start_sec, ctf.setup.end_sec]);
  ampLimit = thresholdDetect.heog.der + (thresholdDetect.heog.der * 0.5);
  set(heogAxisDer,'YLim',[-ampLimit, ampLimit])
  %set(heogAxisDer,'YTick',-ampLimit:10e-6:ampLimit)
  if EYEMVMTdTrial,
    set(heogAxisDer,'color',[1 0.75 0.75])
  else
    set(heogAxisDer,'color',[1 1 1])
  end
  
  pause(10)

end

return
