% corrmap_plot_v2() - uses cluster info and displays summary plots for Matlab 7.5.0
%
% Usage:
%             >> corrmap_plot(CORRMAP,CORRMAP_temp,comp,chan,rep,m,in,fin,plot_ics,mediaplot,tt,flg)
%
% Inputs:
%   CORRMAP         - CORRMAP structure - main
%   CORRMAP_temp    - CORRMAP structure with temporary info from all iterations
%   comp            - cell with components (ALLEEG.icawinv)
%   chan            - channel locations info
%   rep             - auxiliary info to identify ICs from the same dataset
%   m               - index that informs whether to plot summary for 1st or 2nd step
%   in & fin        - flags to indicate which plots are going to be displayed
%   plot_ics        - indices of ICs that are going to be displayed
%   mediaplot       - average plots
%   tt              - threshold info
%   flg             - flag to inform whether it is running in "auto" or "manual" mode
%
% Outputs:
%      summary plots
%
%
%  See also:  corrmap() and pop_corrmap()
%
% Authors: Filipa Campos Viola, 25/01/2008, MRC-IHR, Southampton, UK
% (f.viola@soton.ac.uk)

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) F. Campos Viola, MRC-IHR
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

% revised by F Campos-Viola - corrmap1.01 (30/01/2009)

% revised by F. Campos-Viola make Info Box compatible with Matlab 7.5.0. (28/01/2009); 

function corrmap_plot_v2(CORRMAP,CORRMAP_temp,comp,chan,rep,m,in,fin,plot_ics,mediaplot,tt,flg)

%%%%%%%%%%% fixed parameters for plots %%%%%%%%%%%%%%%%%%

%default values for "boxes" - PLOT max 35 ICs per window
topoplotopts = { };
min=0.04;
space=0.02;
wid=0.12;
leng=wid;
delta=space+wid;
y=min:delta:1-0.3; % positions for IC plots
x=min:delta:1; % positions for IC plots

len_x=length(x); %positions for the individual IC topo plots
len_y=length(y);


%%%%%%%%%%%%%%%% array with different colors for datasets with more than one IC %%%%%%%%%%%%%%%%%%
col=[0.7,0.25,0.5,0.3,0.7,0.5,0.2,0.8,0.4,0.7,1,0,0,0.6,0.2,0.4,0,0.7,0.8,0.25,0.6,1,0,0.5,0.7,0.8,0.9,0.8,0.4,0.7,0.15,0.8,0.4,0.6,0.2,0.25,0.5,0.9,0.4,0.5,0.2,0.8];

sz=4;  %smallest font size

%%%%%%%%%
for v=in:fin;

    if flg==0
        val_th=str2num(CORRMAP.input.corr_th);
        dum=0;
    else
        val_th=CORRMAP.clust.best_th;
        dum=1;
        tres=16;%needs to be changed if the threshold range changes. Current range: 0.95:-0.01:0.80
    end

    %info message for the user
    fprintf('>');
    fprintf('Plotting ICs with a correlation value above %11.4g. \n',val_th);
    fprintf('>');
    fprintf('Selecting %11.4g ICs from each dataset. \n',CORRMAP.input.ics_sel);

    if v==1
        fprintf('>');
        fprintf('Using IC %11.4g from dataset %11.4g as template IC. \n',CORRMAP.template.ic,CORRMAP.template.index);
        fprintf('> \n');
    else
        fprintf('>');
        fprintf('Using average map from previous calculation as template IC. \n');
        fprintf('> \n');
    end

    k=1; %index for positioning head plots
    l=5; %index for positioning head plots


    if  plot_ics(v)<=35

        n_top=plot_ics(v);
        val=1; %one figure
        in=1;

    else
        n_top=35;
        val=2; %two figure
        in=1;
        fprintf('>');
        fprintf('Warning: Clustered ICs will be split in two plots. \n');
        fprintf('> \n');
    end

    %plotting
    numb=1; %while variable

    while numb<=val;
        figure()
        k=1; %index for positioning head plots
        l=5; %index for positioning head plots
        for i=in:n_top

            axes('position',[x(k)  y(l)  wid  leng])

            topoplot(comp{CORRMAP.corr.sets{v}(i)}(:,CORRMAP.corr.ics{v}(i)), chan,topoplotopts{:});

            h3=annotation('textbox',[x(k)  y(l) 0.02  0.02]);

            if rep{m}(CORRMAP.corr.sets{v}(i))>1
                set(h3,'Color',[col(CORRMAP.corr.sets{v}(i)+1) col(CORRMAP.corr.sets{v}(i)+2) col(CORRMAP.corr.sets{v}(i)+3)])
            else
                set(h3,'Color','k')
            end

            set(h3,'FitHeightToText','on')
            set(h3,'String',['r=', num2str(CORRMAP.corr.abs_values{v}(i),'%11.3g'),', set ',num2str(CORRMAP.datasetindices(CORRMAP.corr.sets{v}(i)),'%11.4g'),' IC ',num2str(CORRMAP.corr.ics{v}(i),'%11.4g')])
            set(h3,'FontSize',sz)
            set(h3,'FontUnits','normalized')
            set(h3,'LineStyle','none')
            hold on

            k=k+1;

            if k>7
                l=l-1;
                k=1;
            end

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % PLOTING: template IC (only 1st plot), average map, correlation plot
        % and automatic threshold plot (only if th=auto)


        %%Correlation plot - parameters
        l_var=length(CORRMAP.corr.abs_values{v});
        x_ax=1:l_var;

        line=val_th.*ones(1,l_var); %threshold line

        if v==1

            %template map
            axes('position',[x(2)+0.04 0.72  wid*1.8  wid*1.8])
            topoplot(comp{CORRMAP.template.index}(:,CORRMAP.template.ic), chan,topoplotopts{:});

            %average map
            axes('position',[x(3)+0.06 0.72  wid*1.8  wid*1.8])
            topoplot(reshape(mediaplot(tt,v,:),length(chan),1), chan,topoplotopts{:});


            %info template plot
            h2=annotation('textbox',[x(2)+0.04 0.70 0.05 0.05]);
            set(h2,'Color','k')
            set(h2,'FitHeightToText','on')
            set(h2,'FontSize',sz)
            set(h2,'FontUnits','normalized')
            set(h2,'LineStyle','none')
            set(h2,'HorizontalAlignment','center')
            set(h2,'String',['set ',num2str(CORRMAP.datasetindices(CORRMAP.template.index),'%11.4g'),' IC ',num2str(CORRMAP.template.ic)],'Interpreter','none')

            %info average plot
            h2=annotation('textbox',[x(3)+0.06 0.70 0.05 0.05]);
            set(h2,'FontSize',sz)
            set(h2,'FontUnits','normalized')
            set(h2,'String',['mean_r= ', num2str(CORRMAP.clust.mean_corr(v),'%6.2f \n')],'Interpreter','none')
            set(h2,'LineStyle','none')
            set(h2,'FitHeightToText','on')
            set(h2,'HorizontalAlignment','center')
            % title template
            h2=annotation('textbox',[x(2)+0.04 0.91 0.05 0.05]);
            set(h2,'String','Template IC')
            set(h2,'LineStyle','none')
            set(h2,'FontSize',sz)
            set(h2,'FontUnits','normalized')
            set(h2,'FitHeightToText','on')
            set(h2,'HorizontalAlignment','center')
            % title average
            h2=annotation('textbox',[x(3)+0.06 0.91 0.05 0.05]);
            set(h2,'String','Average map')
            set(h2,'LineStyle','none')
            set(h2,'FontSize',sz)
            set(h2,'FontUnits','normalized')
            set(h2,'FitHeightToText','on')
            set(h2,'HorizontalAlignment','center')

            %correlation plot
            axes('position',[x(5)  0.76  wid*3  0.18])
            plot(x_ax,CORRMAP.corr.abs_values{v},x_ax,line);
            xlabel('ICs','FontSize',sz,'FontUnits', 'normalized')
            ylabel('Correlation (abs)','FontSize',sz,'FontUnits', 'normalized')
            set(gca,'FontSize',sz,'FontUnits', 'normalized')%assures that numbers in axis adjust to the size of the figure
            box off
            grid on
            axis tight
            %title correlation
            h2=annotation('textbox',[0.65 0.96 0.02 0.02]);
            set(h2,'String','Correlation plot')
            set(h2,'FontSize',sz+1)
            set(h2,'FontUnits','normalized')
            set(h2,'LineStyle','none')
            set(h2,'FitHeightToText','on')
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        else

            %only plots average map
            axes('position',[x(2)+0.04 0.72  wid*1.8  wid*1.8])
            topoplot(reshape(mediaplot(tt,v,:),length(chan),1), chan,topoplotopts{:});
            %info template plot
            h2=annotation('textbox',[x(2)+0.06 0.70 0.05 0.05]);
            set(h2,'Color','k')
            set(h2,'FitHeightToText','on')
            set(h2,'FontSize',sz+1)
            set(h2,'FontUnits','normalized')
            set(h2,'LineStyle','none')
            set(h2,'HorizontalAlignment','center')
            set(h2,'String',['mean_r= ', num2str(CORRMAP.clust.mean_corr(v),'%11.4g \n')],'Interpreter','none')
            % title template
            h2=annotation('textbox',[x(2)+0.06 0.91 0.05 0.05]);
            set(h2,'String','Average Map')
            set(h2,'LineStyle','none')
            set(h2,'FontSize',sz+1)
            set(h2,'FontUnits','normalized')
            set(h2,'FitHeightToText','on')
            set(h2,'HorizontalAlignment','center')

            %plotting correlation plot
            axes('position',[x(4)+0.04 0.74  wid*1.6  wid*1.6],'FontUnits', 'normalized','Units','normalized')
            plot(x_ax,CORRMAP.corr.abs_values{v},x_ax,line);
            box off
            grid on
            %axis tight
            xlabel('ICs','FontSize',sz,'FontUnits', 'normalized')
            ylabel('Correlation (abs)','FontSize',sz,'FontUnits', 'normalized')
            set(gca,'FontSize',sz,'FontUnits', 'normalized')%assures that numbers in axis adjust to the size of the figure
            %title correlation
            h2=annotation('textbox',[x(4)+0.04 0.92 0.05 0.05]);
            set(h2,'String','Correlation Plot')
            set(h2,'LineStyle','none')
            set(h2,'FontSize',sz+1)
            set(h2,'FontUnits','normalized')
            set(h2,'FitHeightToText','on')
            set(h2,'HorizontalAlignment','center')

            if dum==1
                %plotting automatic threshold plot
                axes('position',[x(6)+0.04 0.74  wid*1.6  wid*1.6], 'FontUnits', 'normalized','Units','normalized')

                %next steps: findind min(a) - min function was not working - min
                %value is needed to ylim and to plot best threshold line
                sims=CORRMAP_temp.clust.similarity;
                c=sort(sims);
                d=c(1):0.0001:1;

                plot(1:tres,CORRMAP_temp.clust.similarity,'k')
                hold on
                plot(tt(1).*ones(1,length(d)),d,'r') %best threshold line - can be removed

                xlabel('Iterations','FontSize',sz,'FontUnits', 'normalized')
                ylabel('Similarity','FontSize',sz,'FontUnits', 'normalized')
                box off
                grid on
                % placing legend in the right side of the plot

                mid=round(length(tres)/2);
                if tt<=mid;
                    text(tt+0.01,c(2),'\leftarrow suggested threshold','HorizontalAlignment','left','color','r','FontSize',sz,'FontUnits', 'normalized') % legend for best threshold line - can be removed
                else
                    text(tt+0.01,c(2),' suggested threshold \rightarrow','HorizontalAlignment','right','color','r','FontSize',sz,'FontUnits', 'normalized') % legend for best threshold line - can be removed
                end


                if c(1)<1.00000
                    ymax=1;
                    ymin=c(1)-0.001;
                    axis([1 tres ymin ymax])
                else
                    axis tight
                end

                set(gca,'FontSize',sz,'FontUnits', 'normalized') %assures that numbers in axis adjust to the size of the figure
                %title automatic threshold
                h2=annotation('textbox',[x(6)+0.04 0.92 0.05 0.05]);
                set(h2,'String','Threshold Plot')
                set(h2,'LineStyle','none')
                set(h2,'FontSize',sz+1)
                set(h2,'FontUnits','normalized')
                set(h2,'FitHeightToText','on')
                set(h2,'HorizontalAlignment','center')

            end
        end




        %positions for info box
        gap=0.02; % gap between lines - vertical direction  %%%%%
        start=0.01;
        fline=0.94; %position to start info box
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %info box
        h3=annotation('textbox',[start fline gap gap]);
        set(h3,'Color','k')
        set(h3,'FontSize',sz)
        set(h3,'FontUnits','normalized')
        set(h3,'String','INFO:','Interpreter','none')
        set(h3,'LineStyle','none')
        set(h3,'FitBoxToText','on')
        fline=fline-gap;

        h3=annotation('textbox',[start fline gap gap]);
        set(h3,'String',['Template: ',CORRMAP.template.setname,'; Set ',num2str(CORRMAP.datasetindices(CORRMAP.template.index)),'; IC ',num2str(CORRMAP.template.ic),';'],'Interpreter','none') %fixed by initial parameters
        set(h3,'LineStyle','none')
        set(h3,'FontSize',sz)
        set(h3,'FontUnits','normalized')
        set(h3,'FitBoxToText','on')

        fline=fline-gap;

        h3=annotation('textbox',[start fline gap gap]);
        set(h3,'String',['Number of datasets: ',num2str(length(CORRMAP.datasets.setnames))],'Interpreter','none') %fixed by initial parameters
        set(h3,'LineStyle','none')
        set(h3,'FontSize',sz)
        set(h3,'FontUnits','normalized')
        set(h3,'FitBoxToText','on')
        fline=fline-gap;

        h3=annotation('textbox',[start fline gap gap]);
        set(h3,'String',['Correlation threshold: ',num2str(val_th),' (green line)'],'Interpreter','none') %fixed by initial parameters
        set(h3,'LineStyle','none')
        set(h3,'FontSize',sz)
        set(h3,'FontUnits','normalized')
        set(h3,'FitBoxToText','on')
        fline=fline-gap;

        h3=annotation('textbox',[start fline gap gap]);
        set(h3,'String',['Max ICs from each dataset: ',num2str(CORRMAP.input.ics_sel)],'Interpreter','none') %fixed by initial parameters
        set(h3,'LineStyle','none')
        set(h3,'FontSize',sz)
        set(h3,'FontUnits','normalized')
        set(h3,'FitBoxToText','on')
        fline=fline-gap;

        h3=annotation('textbox',[start fline gap gap]);
        set(h3,'String',['Cluster: ',num2str(CORRMAP.clust.ics(v)),' ICs from ',num2str(CORRMAP.clust.sets.number(v)),' sets'],'Interpreter','none')
        set(h3,'LineStyle','none')
        set(h3,'FontSize',sz)
        set(h3,'FontUnits','normalized')
        set(h3,'FitBoxToText','on')
        fline=fline-gap;


        if length(CORRMAP.clust.sets.absent{v})==1&&CORRMAP.clust.sets.absent{v}==0

            h3=annotation('textbox',[start fline gap gap]);
            set(h3,'String','All datasets contribute.')
            set(h3,'LineStyle','none')
            set(h3,'FontSize',sz)
            set(h3,'FontUnits','normalized')
            set(h3,'FitBoxToText','on')
            fline=fline-gap;
        else

            h3=annotation('textbox',[start fline gap gap]);
            set(h3,'String','Sets not contributing: ')
            set(h3,'LineStyle','none')
            set(h3,'FontSize',sz)
            set(h3,'FontUnits','normalized')
            set(h3,'FitBoxToText','on')
            valx=start; % gap in horizontal direction


            if  length(CORRMAP.clust.sets.absent{v})<=10
                fline=fline-gap;
                for z=1:length(CORRMAP.clust.sets.absent{v})
                    h3=annotation('textbox',[valx fline gap gap]);
                    set(h3,'String',['#',num2str(CORRMAP.clust.sets.absent{v}(z)),'; '],'Interpreter','none')
                    set(h3,'LineStyle','none')
                    set(h3,'FontSize',sz)
                    set(h3,'FontUnits','normalized')
                    set(h3,'FitBoxToText','on')
                    valx=valx+0.015;
                end
                fline=fline-gap;

            else
                nnn=1;
                fline=fline-gap;
                beg=1;
                val_f=round(length(CORRMAP.clust.sets.absent{v})/2);

                while nnn<3 %value is 3 because I'm assuming that I'm going to split indices in just 2 lines

                    for z=beg:val_f
                        h3=annotation('textbox',[valx fline gap gap]);
                        set(h3,'String',['#',num2str(CORRMAP.clust.sets.absent{v}(z)),'; '],'Interpreter','none')
                        set(h3,'LineStyle','none')
                        set(h3,'FontSize',sz)
                        set(h3,'FontUnits','normalized')
                        set(h3,'FitBoxToText','on')
                        valx=valx+0.015;
                    end

                    fline=fline-gap;
                    beg=round(length(CORRMAP.clust.sets.absent{v})/2)+1;
                    val_f=length(CORRMAP.clust.sets.absent{v});
                    nnn=nnn+1;
                    valx=start;

                end %close while

            end %close if - condition to check if datasets not contributing need to be in 2 lines

        end



        if CORRMAP.input.ics_sel>1

            h3=annotation('textbox',[start fline gap gap]);
            set(h3,'String',[num2str(CORRMAP.clust.sets.more_oneIC(v)),' datasets contribute with more than 1 IC.'],'Interpreter','none')
            set(h3,'LineStyle','none')
            set(h3,'FontSize',sz)
            set(h3,'FontUnits','normalized')
            set(h3,'FitBoxToText','on')
        end

        fline=fline-gap;
        if v==2

            h3=annotation('textbox',[start fline gap gap]);
            set(h3,'String',['Similarity =',num2str(CORRMAP.clust.similarity,'%11.4f')])
            set(h3,'LineStyle','none')
            set(h3,'FontSize',sz)
            set(h3,'FontUnits','normalized')
            set(h3,'FitBoxToText','on')

        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %main title
        h=annotation('textbox',[0.1 0.99 0.9 0.01]);
        set(h,'Color','k')
        set(h,'String',CORRMAP.input.title)
        set(h,'FontSize',sz+2)
        set(h,'FontUnits','normalized')
        set(h,'LineStyle','none')
        set(h,'HorizontalAlignment','center')
        %set(h1,'FitHeightToText','on')

        numb=numb+1; %updating while variable
        in=36; %updating index i in 'for' - second plot  %check this
        n_top=plot_ics(v); %updating last index i in 'for' - second plot

    end % closes while

end %closes main loop