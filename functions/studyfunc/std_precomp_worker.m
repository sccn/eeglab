function std_precomp_worker(filename, varargin);

g = struct(varargin{2:end});

% load dataset
% ------------
[STUDY ALLEEG] = pop_loadstudy('filename', filename);
if ~isfield(g, 'design'), g.design = STUDY.currentdesign; end;

% run std_precomp (THIS IS THE PART WE WANT TO PARALELIZE)
% ---------------
% for index = 1:length(STUDY.design(g.design).cell)
%     [STUDY ALLEEG] = std_precomp(STUDY, ALLEEG, varargin{:}, 'cell', index);
% end;
[STUDY ALLEEG] = std_precomp(STUDY, ALLEEG, varargin{:});
