function obj2 = fft( obj1, varargin )

obj2 = unitaryopp(@fft, obj1, varargin{:});
