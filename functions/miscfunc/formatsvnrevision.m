function formatsvnrevision(v1, v2)

disp('Getting revision information...');
system(sprintf('svn log --verbose -r%d:%d > tmprev.txt', v1, v2));

oldpwd = pwd;
newpwd = fileparts(which('eeglab.m'));
cd(newpwd);

fid = fopen('tmprev.txt', 'r');

state = 1;
while ~feof(fid)
    txt = fgetl(fid);
    if length(txt) > 4 && strcmpi(txt(1:4), '----')
        state = 1;
        if feof(fid), return; end;
    end;
    if state == 1,
        % get rev number
        txt  = fgetl(fid);
        rev = strtok(txt);
        revnum = rev(2:end);
        
        % get function name
        txt = fgetl(fid);
        txt = fgetl(fid);
        [tmp file] = strtok(txt);
        file = deblank(file);
        [tmp file ext] = fileparts(file);
        file = [ file ext ];
        state = 2;
    end;
    if isempty(txt)
        state = 3;
    end;
    if state == 3
        revinfo = fgetl(fid);
        revinfo = revinfo(1:end);
        state = 1;
        fprintf('** %s, %s (SVN %s - Arno)\n', file, revinfo, revnum);
    end;
end;
cd(oldpwd);
