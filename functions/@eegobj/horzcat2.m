function this = horzcat(varargin);
    this = varargin{1};
    for index = 2:length(varargin)
        this.EEG(index) = varargin{index}.EEG;
    end;
