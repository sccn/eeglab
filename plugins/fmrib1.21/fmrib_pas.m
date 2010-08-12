% fmrib_pas() - Remove pulse
%   artifacts from EEG data collected inside an MRI machine.
%
%   This program removes pulse artifacts from EEG data collected
%   collected inside the MRI scanner. A choice of methods is available.
%
%   The first choice is Optimal Basis Set [Niazy06].
%   This method aligns all the pulse artifacts, in each EEG channel
%   separately, in a matrix and performs a Principal Component Analysis 
%   (PCA) on the data.  The first N PCs (the Optimal Basis Set) are then 
%   fitted to each artifact in that channel. The process is repeated for
%   each channel.  The other three methods are based on [Allen98] with 
%   improvements to better capture and subtract the artifacts. 
%   Basically, in these methods a statistical measurement is used to find a 
%   template for the artifact  at each heart beat.  A window of 30 artifacts 
%   centred around the artifact being processed is used.  The 30 artifacts 
%   are processed by taking the mean ('mean' method), the median 
%   ('median' method) or a Gaussian-weighted
%   ('gmean' method to emphasise shape of current artifact) mean.  It is
%   recommended to use the first ('obs') method as it generally better
%   fits the true artifact.
%
% Usage:
%    >> EEGOUT= fmrib_pas(EEG,qrsevents,method)
% or >> EEGOUT= fmrib_pas(EEG,qrsevents,method,npc)
%
% Inputs:
%   EEG: EEGLAB data structure
%   qrsevents: vector of QRS event locations specified in samples.
%   method: 'obs' for artifact principal components.  You need to specify
%               'npc', which is the number of PC to use in fitting the
%               artifacts.  If unsure, use 4.
%           'mean' for simple mean averaging.
%           'gmean' for Gaussian weighted mean.
%           'median' for median filter.
%
%
% [Niazy06] R.K. Niazy, C.F. Beckmann, G.D. Iannetti, J.M. Brady, and 
%   S.M. Smith (2005) Removal of FMRI environment artifacts from EEG data 
%   using optimal basis sets. NeuroImage 28 (3), pages 720-737.
%
%
% [Allen98] Allen et.al., 1998, Identification of EEG events in the MR
%   scanner: the problem of pulse artifact and a method for its
%   subtraction. NeuroImage8,229-239.
%
%
%
% Also See pop_fmrib_pas
%
%   Author:  Rami K. Niazy
%   
%   Copyright (c) 2006 University of Oxford

% Copyright (C) 2006 University of Oxford
% Author:   Rami K. Niazy, FMRIB Centre
%           rami@fmrib.ox.ac.uk
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% SEP 20, 2005
% Template range now based on median instead of mean

% AUG 16, 2005
% Fixed bug dealing with single precision data.

% JUN 3, 2005
% Released

% APR 20, 2005
% in OBS mode artifact matrix now 'constant'  demeaned instead
%  of linear trend removal

% DEC 23, 2004
% updated (c)

% DEC 22, 2004
% fixed end artifact bug
% fixed copyright

% Dec16, 2004
% Added RAPCO
% re-written help

% Nov 4, 2004
% Updated brog bar
% for extra steps

% Oct 26, 2004
% Changed input to accept
% vector of event location
% instead of event type to
% allow for scripting.





function EEG = fmrib_pas(EEG,QRSevents,method,npc)

nargchk(3,4,nargin);

switch method
    
case 'obs'
    
    if nargin < 4
        error('Incorrect number of input arguments'); 
    end
     %init
    %-----
    fs=EEG.srate; 
    [channels samples]=size(EEG.data);
    pcs=npc-1;

    %standard delay between QRS peak and artifact (allen,1998)  
    delay=round(0.21*fs);

    Gwindow=20; 
    GHW=floor(Gwindow/2);
    rcount=0;
    firstplot=1;

    % memory allocation
    %------------------
    avgdata=zeros(Gwindow,round(fs*2));
    drift=zeros(1,round(fs*2)*Gwindow);
    fitted_art=zeros(channels,samples);
    mPA=zeros(1,round(2*fs));
    peakplot=zeros(1,samples);
    mEye=eye(channels);

    %Extract QRS events
    %------------------

    peakplot(QRSevents)=1;
    sh=zeros(1,delay);
    np=length(peakplot);
    peakplot=[sh peakplot(1:np-delay)];
    peakI=find(peakplot>0);
    peakcount=length(peakI);
    pa=peakcount;
    
    %make filter
    %-----------
    a=[0 0 1 1];
    f=[0 0.4/(fs/2) 0.9/(fs/2) 1];
    ord=round(3*fs/0.5);
    fwts=firls(ord,f,a);

    
   
    %Artifact Subtraction
    %---------------------

    % for the first window/2 points use arthemitic mean for averageing.
    % findg mean QRS peak-to-peak (R-to-R) interval
    for ch=1:channels
        if ch==1
         %init waitbar
        %------------
            barth=5;
            barth_step=barth;
            Flag25=0;
            Flag50=0;
            Flag75=0;
            fprintf('\nPulse artifact subtraction in progress...Please wait!\n');
        end
        
        RR=diff(peakI);
        mRR=median(RR);
        sRR=std(RR);
        %PArange=round(1.25*(mRR+sRR)/2); % half RR
        PArange=round(1.5*mRR/2); % half RR
        midP=PArange+1;
        
%         if ch==1
%             while ((QRSevents(pa-1)+2*PArange-1) > samples)
%                 pa=pa-1;
%             end
%             steps=channels*pa;
%         end
        
        if ch==1
            while (peakI(pa)+PArange > samples)
                pa=pa-1;
            end
            steps=channels*pa;
            peakcount=pa;
        end
        
        eegchan=filtfilt(fwts,1,double(EEG.data(ch,:)));
        pcamat=zeros(pa-1,2*PArange+1);
        dpcamat=pcamat;
        for p=2:pa
            pcamat(p-1,:)=eegchan(peakI(p)-PArange:peakI(p)+PArange);
        end
        pcamat=detrend(pcamat','constant')';
         meaneffect=mean(pcamat);
         dpcamat=detrend(pcamat,'constant');
        [apc,ascore,asvar]=pca_calc(dpcamat');

        papc=[meaneffect' ascore(:,1:pcs)];
   
        
       try
           pad_fit=double(papc)*(double(papc)\...
                double(detrend(EEG.data(ch,peakI(1)-PArange:...
                peakI(1)+PArange)','constant')));
                      
           fitted_art(ch,peakI(1)-PArange:peakI(1)+PArange)=...
                pad_fit';
       catch
       end
       
        

            

        %-----------------------------------------------------
        for p=2:GHW+1

           pad_fit=double(papc)*(double(papc)\...
                double(detrend(EEG.data(ch,peakI(p)-PArange:...
                peakI(p)+PArange)','constant')));
                      
           fitted_art(ch,peakI(p)-PArange:peakI(p)+PArange)=...
                pad_fit';

            %update bar
            %----------               
            percentdone=floor(((ch-1)*pa+p)*100/steps);      
            if floor(percentdone)>=barth
                if percentdone>=25 & Flag25==0
                    fprintf('25%% ')
                    Flag25=1;
                elseif percentdone>=50 & Flag50==0
                    fprintf('50%% ')
                    Flag50=1;
                elseif percentdone>=75 & Flag75==0
                    fprintf('75%% ')
                    Flag75=1;
                elseif percentdone==100
                    fprintf('100%%\n')
                else
                    fprintf('.')
                end

                while barth<=percentdone
                    barth=barth+barth_step;
                end
                if barth>100
                    barth=100;
                end
            end 
        end
       
        %---------------- Processing of central data ---------------------
        %cycle through peak artifacts identified by peakplot
        rcount=GHW;
        for p=GHW+2:peakcount-GHW-1 
            PreP=ceil((peakI(p)-peakI(p-1))/2);
            PostP=ceil((peakI(p+1)-peakI(p))/2);
            if PreP > PArange
                PreP=PArange;
            end
            if PostP > PArange
                PostP=PArange;
            end
            
           pad_fit=double(papc(midP-PArange:midP+PArange,:))*...
               (double(papc(midP-PArange:midP+PArange,:))\...
                double(detrend(EEG.data(ch,peakI(p)-PArange:...
                peakI(p)+PArange)','constant')));
                      
           fitted_art(ch,peakI(p)-PreP:peakI(p)+PostP)=...
                pad_fit(midP-PreP:midP+PostP)';
        
                          

            %update bar
            %----------   
            percentdone=floor(((ch-1)*pa+p)*100/steps);
            if floor(percentdone)>=barth
                if percentdone>=25 & Flag25==0
                    fprintf('25%% ')
                    Flag25=1;
                elseif percentdone>=50 & Flag50==0
                    fprintf('50%% ')
                    Flag50=1;
                elseif percentdone>=75 & Flag75==0
                    fprintf('75%% ')
                    Flag75=1;
                elseif percentdone==100
                    fprintf('100%%\n')
                else
                    fprintf('.')
                end

                while barth<=percentdone
                    barth=barth+barth_step;
                end
                if barth>100
                    barth=100;
                end
            end 
        end
    
        %---------------- Processing of last GHW peaks------------------
        sectionPoints=samples-(peakI(peakcount-GHW)-PreP)+1;

        
        for p=peakcount-GHW:peakcount
            try
               
               pad_fit=double(papc(midP-PArange:midP+PArange,:))*...
                   (double(papc(midP-PArange:midP+PArange,:))\...
                    double(detrend(EEG.data(ch,peakI(p)-PArange:...
                    peakI(p)+PArange)','constant')));

               fitted_art(ch,peakI(p)-PreP:peakI(p)+PostP)=...
                    pad_fit(midP-PreP:midP+PostP)';
                
             
            catch
            end


            %update bar
            %----------      
            percentdone=floor(((ch-1)*pa+p)*100/steps);
            if floor(percentdone)>=barth
                if percentdone>=25 & Flag25==0
                    fprintf('25%% ')
                    Flag25=1;
                elseif percentdone>=50 & Flag50==0
                    fprintf('50%% ')
                    Flag50=1;
                elseif percentdone>=75 & Flag75==0
                    fprintf('75%% ')
                    Flag75=1;
                elseif percentdone==100
                    fprintf('100%%\n')
                else
                    fprintf('.')
                end

                while barth<=percentdone
                    barth=barth+barth_step;
                end
                if barth>100
                    barth=100;
                end
            end        
        end
        
    end

    EEG.data=EEG.data-fitted_art;
      

otherwise

    %init
    %-----
    fs=EEG.srate; 
    [channels samples]=size(EEG.data);

    %standard delay between QRS peak and artifact (allen,1998)  
    delay=round(0.21*fs);

    Gwindow=20; 
    GHW=floor(Gwindow/2);
    rcount=0;
    firstplot=1;

    % memory allocation
    %------------------
    avgdata=zeros(channels,round(fs*2),Gwindow);
    drift=zeros(channels,round(fs*2)*Gwindow);
    clean=zeros(channels,samples);
    mPA=zeros(channels,round(2*fs));
    mMag=zeros(channels,Gwindow);
    peakplot=zeros(1,samples);
    mEye=eye(channels);

    %Extract QRS events
    %------------------

    peakplot(QRSevents)=1;
    sh=zeros(1,delay);
    np=length(peakplot);
    peakplot=[sh peakplot(1:np-delay)];
    peakI=find(peakplot>0);
    peakcount=length(peakI);

    %Artifact Subtraction
    %---------------------

    % for the first window/2 points use arthemitic mean for averageing.
    % findg mean QRS peak-to-peak (R-to-R) interval
    for p=2:GHW+1
        RR(p-1)=peakI(p)-peakI(p-1);
    end
    mRR=mean(RR);
    sRR=std(RR);
    PArange=round((mRR+sRR)/2); % half RR

    t=(0:peakI(GHW+1)+PArange-1)/fs;

    % subtraction of low freq trend; could use detrend.m function instead
    for n=1:channels
        pdrift=polyfit(t,EEG.data(n,1:peakI(GHW+1)+PArange),1);
        drift(n,1:peakI(GHW+1)+PArange)=polyval(pdrift,t);
    end

    EEG.data(:,1:peakI(GHW+1)+PArange)=...
        EEG.data(:,1:peakI(GHW+1)+PArange)-...
        drift(:,1:peakI(GHW+1)+PArange);


    for p=2:GHW+1
        avgdata(:,1:2*PArange+1,p-1)=...
            EEG.data(:,peakI(p)-PArange:peakI(p)+PArange);
    end

    for chan =1:channels
        avgdata(chan,:,:)=detrend(squeeze(avgdata(chan,:,:)),'constant');
    end

    mPA(:,1:2*PArange+1)=mean(avgdata(:,1:1+2*PArange,1:GHW),3);

    if peakI(1) > PArange

        alphaV=...
            sum(detrend(EEG.data(:,peakI(1)-PArange:peakI(1)+PArange)','constant')'...
            .*mPA(:,1:2*PArange+1),2)./...
            sum(mPA(:,1:2*PArange+1).*mPA(:,1:2*PArange+1),2);

        alphaD=diag(alphaV);

    else

        alphaV=sum(detrend(EEG.data(:,1:peakI(1)+PArange)','constant')'...
            .*mPA(:,PArange-peakI(1)+2:2*PArange+1),2)./...
            sum(mPA(:,PArange-peakI(1)+2:2*PArange+1).*...
            mPA(:,PArange-peakI(1)+2:2*PArange+1),2);

        alphaD=diag(alphaV);

    end


    EEG.data(:,1:peakI(GHW+1)+PArange)=...
        EEG.data(:,1:peakI(GHW+1)+PArange)+...
        drift(:,1:peakI(GHW+1)+PArange);


    %init waitbar
    %------------
    barth=5;
    barth_step=barth;
    Flag25=0;
    Flag50=0;
    Flag75=0;
    fprintf('\nPulse artifact subtraction in progress...Please wait!\n');


    if peakI(1) > PArange


        clean(:,peakI(1)-PArange:peakI(1)+PArange)=...
            EEG.data(:,peakI(1)-PArange:peakI(1)+PArange)-alphaD*mPA(:,1:2*PArange+1);
        startpoints=peakI(1)-PArange-1;
        clean(:,1:startpoints)=EEG.data(:,1:startpoints);
    else

        clean(:,1:peakI(1)+PArange)=...
            EEG.data(:,1:peakI(1)+PArange)-alphaD*mPA(:,PArange-peakI(1)+2:2*PArange+1);
    end

    %update bar
    %----------
    percentdone=floor(1*100/peakcount);      
    if floor(percentdone)>=barth
        if percentdone>=25 & Flag25==0
            fprintf('25%% ')
            Flag25=1;
        elseif percentdone>=50 & Flag50==0
            fprintf('50%% ')
            Flag50=1;
        elseif percentdone>=75 & Flag75==0
            fprintf('75%% ')
            Flag75=1;
        elseif percentdone==100
            fprintf('100%%\n')
        else
            fprintf('.')
        end

        while barth<=percentdone
            barth=barth+barth_step;
        end
        if barth>100
            barth=100;
        end
    end 


    %-----------------------------------------------------
    for p=2:GHW+1

        alphaV=sum(EEG.data(:,peakI(p)-PArange:peakI(p)+PArange).*...
            mPA(:,1:2*PArange+1),2)./...
            sum(mPA(:,1:2*PArange+1).*mPA(:,1:2*PArange+1),2);

        alphaD=diag(alphaV);

        clean(:,peakI(p)-PArange:peakI(p)+PArange)=...
            EEG.data(:,peakI(p)-PArange:peakI(p)+PArange)-...
            alphaD*mPA(:,1:2*PArange+1);

        %update bar
        %----------               
        percentdone=floor(p*100/(peakcount-1));      
        if floor(percentdone)>=barth
            if percentdone>=25 & Flag25==0
                fprintf('25%% ')
                Flag25=1;
            elseif percentdone>=50 & Flag50==0
                fprintf('50%% ')
                Flag50=1;
            elseif percentdone>=75 & Flag75==0
                fprintf('75%% ')
                Flag75=1;
            elseif percentdone==100
                fprintf('100%%\n')
            else
                fprintf('.')
            end

            while barth<=percentdone
                barth=barth+barth_step;
            end
            if barth>100
                barth=100;
            end
        end 
    end

    %---------------- Processing of central data ---------------------
    %cycle through peak artifacts identified by peakplot
    for p=GHW+2:peakcount-GHW-1 
        PreP=ceil((peakI(p)-peakI(p-1))/2); % interval before peak artifact
        PostP=ceil((peakI(p+1)-peakI(p))/2);% interval after peak artifact
        sectionPoints=peakI(p+GHW)+PostP-(peakI(p-GHW)-PreP-1);

        % subtract drift
        t=(0:sectionPoints-1)/fs;
        for n=1:channels
            pdrift=...
                polyfit(t,EEG.data(n,peakI(p-GHW)-PreP:peakI(p+GHW)+PostP),1);
            drift(n,1:sectionPoints)=polyval(pdrift,t);
        end

        EEG.data(:,peakI(p-GHW)-PreP:peakI(p+GHW)+PostP)=...
            EEG.data(:,peakI(p-GHW)-PreP:peakI(p+GHW)+PostP)-...
            drift(:,1:sectionPoints);

        rr=1;

        for r=p-GHW:p+GHW
            mMag(:,rr)=mean(abs(EEG.data(:,peakI(r)-PreP:peakI(r)+PostP)),2);
            rr=rr+1;
        end

        minMag=min(mMag,[],2);

        %only count those EEG data with no motion or line artifact by rejecting
        %excessively large variations.  Use channel 1 as ref (Cz).
        oldavgdata=avgdata;
        oldrcount=rcount;
        rcount=0;
        for r=p-GHW:p+GHW
            if (mean(abs(EEG.data(1,peakI(r)-PreP:peakI(r)+PostP))) ...
                    < 4*minMag(1) ) ...
                    & (max(abs(EEG.data(1,peakI(r)-PreP:peakI(r)+PostP)))...
                    < 2* 175)        
                avgdata(:,1:PreP+PostP+1,rcount+1)=...
                    EEG.data(:,peakI(r)-PreP:peakI(r)+PostP);
                rcount=rcount+1;
            end
        end

        for chan =1:channels
            avgdata(chan,:,:)=detrend(squeeze(avgdata(chan,:,:)),'constant');
        end

    %     rcount
        %if more than 8 good data sections found then do averaging, else use
        %previous.
        if rcount > 5
            switch method
            case 'gmean'
                if p==peakcount-GHW-1 
                    mPA=mean(avgdata(:,:,1:rcount),3);
                elseif rcount==(2*GHW+1)
                    G=gausswin(rcount);
                    G=G/sum(G);
                    for r=1:rcount
                        avgdata(:,:,r)=avgdata(:,:,r)*G(r);
                    end
                    mPA=sum(avgdata(:,:,1:rcount),3); 
                else
                    mPA=median(avgdata(:,:,1:rcount),3);
                end
            case 'mean'
                mPA=mean(avgdata(:,:,1:rcount),3);
            case 'median'
                mPA=median(avgdata(:,:,1:rcount),3);
            end
        else
            switch method
            case 'mean'
                mPA=mean(oldavgdata(:,:,1:oldrcount),3);
            case 'median'
                mPA=median(oldavgdata(:,:,1:oldrcount),3);
            case 'gmean'
                mPA=mean(oldavgdata(:,:,1:oldrcount),3);
            end
        end

        alphaV=sum(detrend(EEG.data(:,peakI(p)-PreP:peakI(p)+PostP)','constant')'.*...
            mPA(:,1:PreP+PostP+1),2)./...
            sum(mPA(:,1:PreP+PostP+1).*mPA(:,1:PreP+PostP+1),2);

        alphaD=diag(alphaV);

        EEG.data(:,peakI(p-GHW)-PreP:peakI(p+GHW)+PostP)=...
            EEG.data(:,peakI(p-GHW)-PreP:peakI(p+GHW)+PostP)+...
            drift(:,1:sectionPoints);


        clean(:,peakI(p)-PreP:peakI(p)+PostP)=...
            EEG.data(:,peakI(p)-PreP:peakI(p)+PostP)-...
            alphaD*mPA(:,1:PreP+PostP+1);

        %update bar
        %----------               
        percentdone=floor(p*100/(peakcount-1));      
        if floor(percentdone)>=barth
            if percentdone>=25 & Flag25==0
                fprintf('25%% ')
                Flag25=1;
            elseif percentdone>=50 & Flag50==0
                fprintf('50%% ')
                Flag50=1;
            elseif percentdone>=75 & Flag75==0
                fprintf('75%% ')
                Flag75=1;
            elseif percentdone==100
                fprintf('100%%\n')
            else
                fprintf('.')
            end

            while barth<=percentdone
                barth=barth+barth_step;
            end
            if barth>100
                barth=100;
            end
        end 
    end

    %---------------- Processing of last GHW peaks------------------
    sectionPoints=samples-(peakI(peakcount-GHW)-PreP)+1;

    if samples-peakI(peakcount) >= PostP

        alphaV=sum(detrend(EEG.data(:,peakI(peakcount)-...
            PreP:peakI(peakcount)+PostP)','constant')'.*...
            mPA(:,1:PreP+PostP+1),2)./...
            sum(mPA(:,1:PreP+PostP+1).*mPA(:,1:PreP+PostP+1),2);

        alphaD=diag(alphaV);

        clean(:,peakI(peakcount)-PreP:peakI(peakcount)+PostP)=...
            EEG.data(:,peakI(peakcount)-PreP:peakI(peakcount)+PostP)- ...
                            alphaD*mPA(:,1:PreP+PostP+1);

        endpoints=samples-peakI(peakcount)-PostP;

        clean(:,samples-endpoints+1:samples)=...
            EEG.data(:,samples-endpoints+1:samples);        

    else

        alphaV=sum(detrend(EEG.data(:,peakI(peakcount)-PreP:samples)','constant')'.*...
            mPA(:,1:samples-(peakI(peakcount)-PreP)+1),2)./...
            sum(mPA(:,1:samples-(peakI(peakcount)-PreP)+1).*...
            mPA(:,1:samples-(peakI(peakcount)-PreP)+1),2);

        alphaD=diag(alphaV);

        clean(:,peakI(peakcount)-PreP:samples)=...
            EEG.data(:,peakI(peakcount)-PreP:samples)-...
            alphaD*mPA(:,1:samples-(peakI(peakcount)-PreP)+1);
    end



    for p=peakcount-GHW:peakcount-1

        alphaV=sum(detrend(EEG.data(:,peakI(p)-PreP:peakI(p)+PostP)','constant')'.*...
            mPA(:,1:PreP+PostP+1),2)./...
            sum(mPA(:,1:PreP+PostP+1).*mPA(:,1:PreP+PostP+1),2);

        alphaD=diag(alphaV);

        clean(:,peakI(p)-PreP:peakI(p)+PostP)=...
            EEG.data(:,peakI(p)-PreP:peakI(p)+PostP)-...
            alphaD*mPA(:,1:PreP+PostP+1);

        %update bar
        %----------               
        percentdone=floor(p*100/(peakcount-1));
        if floor(percentdone)>=barth
            if percentdone>=25 & Flag25==0
                fprintf('25%% ')
                Flag25=1;
            elseif percentdone>=50 & Flag50==0
                fprintf('50%% ')
                Flag50=1;
            elseif percentdone>=75 & Flag75==0
                fprintf('75%% ')
                Flag75=1;
            elseif percentdone==100
                fprintf('100%%\n')
            else
                fprintf('.')
            end

            while barth<=percentdone
                barth=barth+barth_step;
            end
            if barth>100
                barth=100;
            end
        end           
    end
    EEG.data=clean;
        
end
return;
