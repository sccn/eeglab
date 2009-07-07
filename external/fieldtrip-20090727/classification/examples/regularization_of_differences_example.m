fclose('all');
close all
clear all
clc
format long

timedivs=12;%6%[1,2,3:3:90]
    
for filenum=4:10
   
    clear cfg
    
    data = load(sprintf('/scratch/ali/outputf-1_3_3_90_filterer/%d/freq%d',filenum,timedivs));
    design=data.labels;
    ddata=data.data;
    clear data
    
    data{1}=ddata;
    data{2}=ddata;
    data{3}=ddata;
    clear ddata;
    
    data{1}.powspctrm=data{1}.powspctrm(design==1,13:15,7:13,:);%12,8:9
    data{2}.powspctrm=data{2}.powspctrm(design==2,13:15,7:13,:);
    data{3}.powspctrm=data{3}.powspctrm(design==3,13:15,7:13,:);

    rand('twister',1); randn('state',1);

    for c=1:length(data)
        data{c}.time = 1;
        data{c}.dimord = 'rpt_chan_freq';
       
        if c==1
            ddata=data{c}.powspctrm;
        else
            ddata=cat(1,ddata,data{c}.powspctrm);
        end
    end
    data=ddata;
    clear ddata
    
    rand('twister',1); randn('state',1);
    myproc = clfproc({ ...
        standardizer() ...
        regdif('MaxIter',obj.MaxIter_Intern,'method','ROV','lambda',10^-3,'divnum',obj.divnum)...
        });
     
    cv = crossvalidator('procedure',myproc,'cvfolds',10,'randomize',true);
    cv = cv.validate(data,design);
    
    accuracy(filenum)=evaluate(cv.post,cv.design,'metric','accuracy')
    
    cvg{filenum}=cv;
end