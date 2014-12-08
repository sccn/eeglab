function std_plotmat(usrdat,des)

%retreiving subject
listsubj = unique(usrdat.subjects);
indtmp =  inputgui( 'geometry', { [1 1]}, 'uilist', { ...
                   { 'style', 'text', 'string', 'Select subject' }, ...
                   { 'style', 'popupmenu', 'string', listsubj}} );
subj = listsubj{cell2mat(indtmp)};

des
