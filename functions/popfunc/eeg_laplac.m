% eeg_laplac() - gives the laplacian for the data contained in EEG
% for the values 1:stepl:end with electrodes position given by EEG.chanlocs. 
% Step must be a number between 1 and length(EEG.data) 
% The method used for computig the laplacian is the described by Perrin et al.
% 1989 using a spline interpolation between the electrodes. This function 
% requires the functions g2() and gm_1() to be in the path or in the same dir. 
% Usage:
%        >>  [laplac time] = eeg_laplac(EEG, step);
%        >>  laplac = eeg_laplac(EEG, step);
% Inputs:
%   EEG     - A EEG structure containing the EEG.data an the EEG.chanlocs
%   step    - The step for subsampling the data
%    
% Outputs:
%   laplac     - A matrix containing the value of laplacian for the data
%   time       - Time spent in calculation  
%
% See also: del2map()
%
% Author: Juan Sebastian Gonzalez, DFI / Universidad de Chile
% 
% Message from Jrgen Kayser, one of the main author of the CSD (Current
% Source Density) toolbox: "I do not find any differences between his Laplacian 
% solution and our CSD transform (i.e., apart from the lack of a smoothing
% constant, the choice of the m constant, and the specific spherical locations 
% assigned to any given electrode). Even more important would be a warning that
% the flexibility constant m is fixed at 3 (2 = max flexible, 3 = medium, 
% 4=medium rigid, 5+=increasingly less flexible splines), which is a debatable 
% choice and not recommended by us. The choice of the m constant (and the 
% regulization or smoothing constant) should depend on the objective of the
% study."

% Copyright (C)2006  Juan Sebastian Gonzalez, DFI / Universidad de Chile
%                                          
% This file is part of EEGLAB, see http://www.eeglab.org
% for the documentation and details.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
% THE POSSIBILITY OF SUCH DAMAGE.

function [laplac, time] = eeg_laplac(EEG, paso)

if nargin < 2
   help eeg_laplac;
   return;
end

tic
n = length(EEG.data(:,1));
a = zeros(n,3);
g = zeros(n);
for i=1:n
    a(i,1)=EEG.chanlocs(i).sph_phi;
    a(i,2)=EEG.chanlocs(i).sph_theta;
end
V = zeros(1,n+1);
aux = a*pi/180;
matriz_g = zeros(n+1);
for i=1:n
    matriz_g(i,1:n) = g2(aux(i,1)*ones(n,1),aux(i,2)*ones(n,1),aux(:,1),aux(:,2));
    g(i,:) = gm_1(aux(i,1)*ones(n,1),aux(i,2)*ones(n,1),aux(:,1),aux(:,2));
end
matriz_g(:,n+1) = 1;
matriz_g(n+1,:) = 1;
inv_m = inv(matriz_g);
k = 1;
largo = length(EEG.data(1,:));
laplac = zeros(n,round((largo-1)/paso));
for j=1:paso:largo
    % ciclo para sacar los anglos phi y theta en grados    
    for i=1:n
        V(i) = EEG.data(i,j);
    end
    % calculo la constante a partir de la matriz invertida 
    constant = inv_m*V';
    % una vez que tengo las constantes, calculo el laplaciano segun la formula
    % de perrin 89 corregida
    for i=1:n
        laplac(i,k) = g(i,:)*constant(1:n);
    end
    k = k + 1;
end
time = toc;
%---------------------------------------------------------------------
% internal functions g2 gm_1
function valor=g2(r1,r2,r3,r4)
%g entrega la sumatoria hasta 7 con m=3 de la ec. (2) de la tecnical note del
%Electrical Geodesics, salvo que esta esta bien

x=cos(r1).*cos(r2).*cos(r3).*cos(r4) + ...
    cos(r1).*sin(r2).*cos(r3).*sin(r4) + ...
    sin(r1).*sin(r3);
L = length(x);
for i=1:L
    if x(i) > 1
        x(i) = 1;
    
    elseif x(i) < -1
        x(i) = -1;
        
    end
end

P = zeros(25,length(x));
for n=1:25
    aux = legendre(n,x); % calculo los polinomios de legendre de orden n %
    P(n,:) = aux(1,:); % tomo el que me sirve %
end

n = 1:25;

for i=1:25
    aux2(i,:) = ((2*n(i)+1)./(n(i).*(n(i)+1)).^3)*ones(1,L);
end
    valor = 1/4/pi*sum(aux2.*P);
 %-----------------------------------------------------   
function valor=gm_1(r1,r2,r3,r4)
%g entrega la sumatoria hasta 7 con m=3 de la ec. (2) de la tecnical note del
%Electrical Geodesics

x=cos(r1).*cos(r2).*cos(r3).*cos(r4) + ...
    cos(r1).*sin(r2).*cos(r3).*sin(r4) + ...
    sin(r1).*sin(r3);
L = length(x);
for i=1:L
    if x(i) > 1
        x(i) = 1;
    
    elseif x(i) < -1
        x(i) = -1;
        
    end
end
P = zeros(25,length(x));
for n=1:25
    aux=legendre(n,x);
    P(n,:)=aux(1,:);
end

n=1:25;
%el problema es calcular esto para cada columna...sin hacer un for
%respuesta, hacer for de 7 vueltas jajaj
for i=1:25
    aux2(i,:)=(2*n(i)+1)/(n(i)*(n(i)+1))^2*ones(1,L);
end

    valor=1/4/pi*sum(aux2.*P);
    


