function makeDatFromCTF(filesearch,datfilename,typego,typebutton)

% makeDatFromCTF
% filesearch is a string that if used with dir, can identify all files
% necssary to include. Make sure you're in the right directory. It should
% lead to files with columns of continuous STIM channel output
% typego lists the types of stim that one should expect button responses
% to.
% typebutton is the button response one should get. 0 if not buttonpress

if (nargin < 3)
    typego = [6 7 8 9];
    typebutton = [2 0 0 0];
end

respvals = [33554432 65536 262144 327680 393216 458752 524288 589824];
resptype = [2 1 4 5 6 7 8 9];
maxtypeval = 589824; % above this should be button presses only.
RESPLENGTH = 50;

fs = 1200;
buttonwait = fs*1.5;
numHeaderLines = 19;
files = dir(filesearch);

allStim = [];
for i = 1:length(files)
    fprintf('Getting Data from %s ... \n',files(i).name);
    [time,stim] = textread(files(i).name,'%f %f','delimiter','\t','headerlines',2);
    allStim = [allStim;stim];
    clear time stim;
end

fid = fopen(datfilename,'W');

for i = 1:numHeaderLines
    fprintf(fid,'\n');
end

fprintf(fid,'Trial\tResp\tType\tCorrect\tLatency\tStim/Resp\n');
fprintf(fid,'-----\t----\t----\t-------\t-------\t---------\n');
fprintf(fid,'1\t0\t1\t1\t1000\tStim\n');

i = 1;
trialnum = 2;
while i <= length(allStim)
    nextJump = 1;
    tmp = allStim(i:min(i+RESPLENGTH,length(allStim)));
    if (tmp(1) > 10 && ~isempty(find(respvals == max(tmp))) && max(tmp) <= maxtypeval)
        type = resptype(find(respvals == max(allStim(i:i+RESPLENGTH))));

        correct = 1;
        latency = 0;
        buttonPressed = 0;

        if (~isempty(find(typego == type))) % s2 occured, look for button press
            validButton = typebutton(find(typego == type));
            j = i + RESPLENGTH;
            correct = 0;

            while j < i + buttonwait % see if a button was pressed during the wait period
                tmp = allStim(j:min(j+RESPLENGTH,length(allStim)));
                if (~isempty(tmp) && tmp(1) > 10 && ~isempty(find(respvals == max(tmp))))
                    latency = round((j - i) * 1000 / fs);
                    buttonPressed = resptype(find(respvals == max(tmp)));
                    j = i + buttonwait;

                end
                j = j + 1;
            end
            
            if (buttonPressed == validButton)
                correct = 1;
            end

        end

        if (~buttonPressed) % get Latency if not button pressed
            j = i + RESPLENGTH;
            cont = 1;
            while (cont && j < length(allStim))
                j = j + 1;
                tmp = allStim(j:min(j+RESPLENGTH,length(allStim)));
                if (isempty(tmp))
                    cont = 0;
                elseif (tmp(1) > 10 && ~isempty(find(respvals == max(tmp))) && max(tmp) <= maxtypeval)
                    cont = 0;
                end
            end
            latency = round((j - i) * 1000 / fs);
            nextJump = (j - i) - 1;
        else
            nextJump = round(latency * fs / 1000) + RESPLENGTH;
        end
        fprintf(fid,'%d\t%d\t%d\t%d\t%d\tStim\n',trialnum,buttonPressed,type,correct,latency);
        trialnum = trialnum+1;
    end
    i = i + nextJump;
end

fclose(fid);

fprintf('***********************************************\n');
fprintf('%s created successfully: %d events found\n\n',datfilename,trialnum);
