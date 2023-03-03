function motion_correction_quantification
for ii = 1:5
    pixelSize(ii) = ei{ii}.thorExp.widthUM/ei{ii}.thorExp.pixelX;
end
%% motion correction
    ntrials = 40; %si = [Lb_T ArL_L_T Lbs_T Ab_t_T Ab_i_T Abs_t_T Abs_i_T Ar_t_D ArL_t_D Ars_t_D Ar_t_T ArL_t_T Ars_t_T Ar_i_D ArL_i_D Ars_i_D Ar_i_T ArL_i_T Ars_i_T];
    si = [Lb Ab_On Ab_Off Ar_On Ar_Off ArL_On ArL_L ArL_Off Ars_On Ars_Off Lbs Abs_On Abs_Off];
    si = [Lb Ab_On Ab_Off Ab_Offc Ar_On Ar_Off Ar_Offc ArL_On ArL_L ArL_Off ArL_Offc Ars_On Ars_Off Ars_Offc Lbs Abs_On Abs_Off Abs_Offc];
   
    Rs = o.Rs(:,si); Rs_MC = o.RsMC(:,si);
    titles = rasterNamesTxt(si);
     
     ff = makeFigureRowsCols(107,[1 0.5 4 0.5],'RowsCols',[1 18],'spaceRowsCols',[0.1 0.016],'rightUpShifts',[0.035 0.18],'widthHeightAdjustment',[-18 -305]);
    set(gcf,'color','w'); set(gcf,'Position',[5 5 6.9 0.75]);
%     [Y,E] = discretize(1:49,3);
    all_speeds = []; cTxt = {'NB-A','NB-AL','NB-A*','Ar-Off','ArL-Off','Ar*-Off'}; 
    ind = 1;
    for cn = 1:18
        mean_speed_over_trials = [];
        aThisSpeed = [];
        for an = 1:size(Rs,1)
%             thisSpeed = nanmean(Rs{an,cn}.speed);
            thisSpeed = nanmean(Rs_MC{an,cn}.sp_rasters)*pixelSize(ii);
            mean_speed_over_trials(an,:) = thisSpeed;
        end
        axes(ff.h_axes(1,cn));
        hold on;
        xs = linspace(-2,2,size(thisSpeed,2));
        N = length(xs);
        halfN = floor(N/2);
        mspeed_before = mean_speed_over_trials(:,1:halfN);
        mspeed_after = mean_speed_over_trials(:,(halfN+1):(halfN+halfN));
        mspeed = mean(mean_speed_over_trials(:,1:N)); semspeed = std(mean_speed_over_trials(:,1:N))/sqrt(5);
        
        plot(xs,mspeed);
        shadedErrorBar(xs,mspeed,semspeed);
%         changePosition(gca,[0.1 0.15 -0.05 -0.15]);
%         if cn == 1 || cn ==3
%             put_axes_labels(gca,{'',[0 0 0]},{{'Speed (cm/sec)'},[0 00 0]});
%         end
        xbTxt = [2.5 7.5 12.5]-3; ybTxt = 31;
%         text(xbTxt(1),ybTxt+5,cTxt{cn},'FontSize',5);
        ylim([0 14]);
        box off;
        format_axes(gca);
        plot([0 0],[0 12],'m','linewidth',0.25);
%         if cn == 3
%             hx = xlabel('Time (sec)');
%         end
        format_axes(gca);
        text(-1,15,titles{cn},'FontSize',6);
        if cn > 1
            set(gca,'YTick',[]);
        else 
            ylabel('Distance (\mum)');
        end
        [within,dvn,xlabels] = make_within_table({'BA','cols'},[2,halfN]);
        dataT = make_between_table({[mspeed_before mspeed_after]},dvn);
        ra = RMA(dataT,within);
        all_ras{cn} = ra;
%         if ismember(cn,[1 2 3])
%             set(gca,'XTick',[])
%         end
    end
    save_pdf(ff.hf,mData.pdf_folder,sprintf('MC_all'),600);

%% Speed Figure
while 1
    ntrials = 40; %si = [Lb_T ArL_L_T Lbs_T Ab_t_T Ab_i_T Abs_t_T Abs_i_T Ar_t_D ArL_t_D Ars_t_D Ar_t_T ArL_t_T Ars_t_T Ar_i_D ArL_i_D Ars_i_D Ar_i_T ArL_i_T Ars_i_T];
    si = [Lb Ab_On Ab_Off Ar_On Ar_Off ArL_On ArL_L ArL_Off Ars_On Ars_Off Lbs Abs_On Abs_Off];
    Rs = o.Rs(:,si); Rs_MC = o.RsMC(:,si);
    titles = rasterNamesTxt(si);
    %% speed
     ff = makeFigureRowsCols(107,[1 0.5 4 0.5],'RowsCols',[1 13],'spaceRowsCols',[0.1 0.016],'rightUpShifts',[0.035 0.18],'widthHeightAdjustment',[-18 -305]);
    set(gcf,'color','w'); set(gcf,'Position',[5 5 6.9 0.75]);
%     [Y,E] = discretize(1:49,3);
    all_speeds = []; cTxt = {'NB-A','NB-AL','NB-A*','Ar-Off','ArL-Off','Ar*-Off'}; 
    ind = 1;
    for cn = 1:13
        mean_speed_over_trials = [];
        aThisSpeed = [];
        for an = 1:size(Rs,1)
%             thisSpeed = nanmean(Rs{an,cn}.speed);
            thisSpeed = nanmean(Rs_MC{an,cn}.sp_rasters)*pixelSize(ii);
            mean_speed_over_trials(an,:) = thisSpeed;
        end
        axes(ff.h_axes(1,cn));
        hold on;
        xs = linspace(-2,2,size(thisSpeed,2));
        N = length(xs);
        halfN = floor(N/2);
        mspeed_before = mean_speed_over_trials(:,1:halfN);
        mspeed_after = mean_speed_over_trials(:,(halfN+1):(halfN+halfN));
        mspeed = mean(mean_speed_over_trials(:,1:N)); semspeed = std(mean_speed_over_trials(:,1:N))/sqrt(5);
        
        plot(xs,mspeed);
        shadedErrorBar(xs,mspeed,semspeed);
%         changePosition(gca,[0.1 0.15 -0.05 -0.15]);
%         if cn == 1 || cn ==3
%             put_axes_labels(gca,{'',[0 0 0]},{{'Speed (cm/sec)'},[0 00 0]});
%         end
        xbTxt = [2.5 7.5 12.5]-3; ybTxt = 31;
%         text(xbTxt(1),ybTxt+5,cTxt{cn},'FontSize',5);
        ylim([0 14]);
        box off;
        format_axes(gca);
        plot([0 0],[0 12],'m','linewidth',0.25);
%         if cn == 3
%             hx = xlabel('Time (sec)');
%         end
        format_axes(gca);
        text(-1,15,titles{cn},'FontSize',6);
        if cn > 1
            set(gca,'YTick',[]);
        else 
            ylabel('Distance (\mum)');
        end
        [within,dvn,xlabels] = make_within_table({'BA','cols'},[2,halfN]);
        dataT = make_between_table({[mspeed_before mspeed_after]},dvn);
        ra = RMA(dataT,within);
        all_ras{cn} = ra;
%         if ismember(cn,[1 2 3])
%             set(gca,'XTick',[])
%         end
    end
    save_pdf(ff.hf,mData.pdf_folder,sprintf('MC_all'),600);
    %%
    for cn = 1:18
        ra = all_ras{cn};
        pValsMC(cn,:) = ra.ranova.pValue([3 5 7])';
    end
    %%
    
    %% speed
     ff = makeFigureRowsCols(107,[1 0.5 4 0.5],'RowsCols',[1 13],'spaceRowsCols',[0.1 0.016],'rightUpShifts',[0.035 0.18],'widthHeightAdjustment',[-18 -305]);
    set(gcf,'color','w'); set(gcf,'Position',[5 5 6.9 0.75]);
%     [Y,E] = discretize(1:49,3);
    all_speeds = []; cTxt = {'NB-A','NB-AL','NB-A*','Ar-Off','ArL-Off','Ar*-Off'}; 
    ind = 1;
    for cn = 1:13
        mean_speed_over_trials = [];
        aThisSpeed = [];
        for an = 1:size(Rs,1)
            thisSpeed = nanmean(Rs{an,cn}.speed);
%             thisSpeed = nanmean(Rs_MC{an,cn}.sp_rasters);
            mean_speed_over_trials(an,:) = thisSpeed;
        end
        if ismember(cn,[1 2 3 11 12 13]);
            delete(ff.h_axes(1,cn));
            continue;
        end
        axes(ff.h_axes(1,cn));
        hold on;
        xs = linspace(-2,2,size(thisSpeed,2));
        N = length(xs);
        halfN = floor(N/2);
        mspeed_before = mean_speed_over_trials(:,1:halfN);
        mspeed_after = mean_speed_over_trials(:,(halfN+1):(halfN+halfN));
        mspeed = mean(mean_speed_over_trials(:,1:N)); semspeed = std(mean_speed_over_trials(:,1:N))/sqrt(5);
        
        plot(xs,mspeed);
        shadedErrorBar(xs,mspeed,semspeed);
%         changePosition(gca,[0.1 0.15 -0.05 -0.15]);
%         if cn == 1 || cn ==3
%             put_axes_labels(gca,{'',[0 0 0]},{{'Speed (cm/sec)'},[0 00 0]});
%         end
        xbTxt = [2.5 7.5 12.5]-3; ybTxt = 31;
%         text(xbTxt(1),ybTxt+5,cTxt{cn},'FontSize',5);
        ylim([0 30]);
        box off;
        format_axes(gca);
        plot([0 0],[0 25],'m','linewidth',0.25);
%         if cn == 3
%             hx = xlabel('Time (sec)');
%         end
        format_axes(gca);
        text(-1,31,titles{cn},'FontSize',6);
        if cn > 4
            set(gca,'YTick',[]);
        else 
            ylabel('Speed (cm/s)');
        end
        [within,dvn,xlabels] = make_within_table({'BA','cols'},[2,halfN]);
        dataT = make_between_table({[mspeed_before mspeed_after]},dvn);
%         ra = RMA(dataT,within);
        all_ras_S{cn} = ra;

    end
    save_pdf(ff.hf,mData.pdf_folder,sprintf('apeed_all'),600);
    %%
    for cn = 1:13
        if ismember(cn,[1 2 3 11 12 13])
            pVals(cn,:) = [1 1 1];
            continue;
        end
        ra = all_ras_S{cn};
        pVals(cn,:) = ra.ranova.pValue([3 5 7])';
    end
    %%
    break;
end
