% This make.m is used under Windows

if exist('OCTAVE_VERSION','builtin')
        % Octave 2.9 or higher supports mex-interface
        mex   svmtrain.c svm.cpp svm_model_matlab.c
        % mex svmpredict.c svm.cpp svm_model_matlab.c  % is currently (Octave 2.9.8+) not working 
else
        mex -O -c svm.cpp
        mex -O -c svm_model_matlab.c
        mex -O svmtrain.c svm.obj svm_model_matlab.obj
        mex -O svmpredict.c svm.obj svm_model_matlab.obj
end;