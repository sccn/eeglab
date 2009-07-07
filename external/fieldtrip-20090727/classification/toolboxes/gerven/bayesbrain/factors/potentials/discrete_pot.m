classdef discrete_pot < pot
% DISCRETE_POT discrete potential class 
%
%   SEE ALSO:
%       multinomial_pot
%
%   Copyright (C) 2008  Marcel van Gerven
%
%   $Log: not supported by cvs2svn $
%

    methods
        function obj = discrete_pot(ddomain,dsize)
            % constructor

            obj = obj@pot([],ddomain,dsize);
        end
    end
end