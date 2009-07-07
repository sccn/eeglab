% this script tests classification procedures; 
% add procedure below for benchmarking. 
%
% benchmarking saved as preprocessors.log
% 
%   Copyright (c) 2009, Marcel van Gerven
%
%   $Log: not supported by cvs2svn $
%

fclose('all');
close all
clear all

%% specify procedure here

procedures = {
  { standardizer() kernelmethod() } ... 
  { preprocessor('prefun',@(x)(log10(x))) standardizer() kernelmethod() } ...
  { preprocessor('prefun',@(x)(log10(x))) standardizer() rbmhierarchy('nbatches',10,'rbms',{rbmlayer('numhid',100,'maxepoch',3) rbmlayer('numhid',10,'maxepoch',3)}) kernelmethod() } ...
  { preprocessor('prefun',@(x)(log10(x))) standardizer() pcanalyzer('proportion',0.8) kernelmethod() } ...
  { preprocessor('prefun',@(x)(log10(x))) standardizer() tsner('dims',2,'initial_dims',30) kernelmethod() } ...
};

descriptions = { ...
  'standardizer' ...
  'log transform + standardizer' ...
  'restricted boltzmann machine' ...
  'principal component analysis' ...
  't-SNE algorithm' ...  
};

%% start analysis

% iterate over all specified classification procedures
fid = fopen('preprocessors.log','w+');
for c=1:length(procedures)

  tic;
  [acc,p] = test_procedure(procedures{c});
   
  fprintf(fid,'%f\t%f\t%f\t%s\n',acc,p,toc,descriptions{c});

end
fclose(fid);
