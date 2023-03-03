for ii = 1:5
    pixelSize(ii) = ei{ii}.thorExp.widthUM/ei{ii}.thorExp.pixelX;
end
%% population vector and correlation sensory for time rasters with MC motion profile
  
    ff = makeFigureRowsCols(108,[6 3 2.1 3.2],'RowsCols',[4 2],'spaceRowsCols',[0.025 0.01],'rightUpShifts',[0.14 0.08],'widthHeightAdjustment',[10 -45]);
%     set(gcf,'color','w');    set(gcf,'Position',[10 3 3.5 3.75]);
    MY = 8; ysp = 1; mY = 0; % responsive cells
    stp = 0.25*magfac; widths = ([0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5]+0.2)*magfac; gap = 0.33*magfac;
    adjust_axes(ff,[mY MY],stp,widths,gap,{'Cell #'});
    
    an = 1;G = 'C';
    si = [Ab_On Ar_On]; headstr = {'Brake (B)','No-Brake (NB)'}; headstrpos = [0 -0.14 0 0]; xlimss = [-2 2];
    si = [Lb ArL_L]; headstr = {'Brake (B)','No-Brake (NB)'}; headstrpos = [0 -0.14 0 0]; xlimss = [-2 2];
    si = [M_On M_Off]; headstr = {'Motion Onset (MOn)','Motion Offset (MOff)'}; headstrpos = [0 -0.14 0.05 0]; xlimss = [-1.5 1.5];
    Rs = o.Rs(:,si);mR = o.mR(:,si);
    Rs_MC = o.RsMC(:,si);
    ntrials = 50;
    props1 = get_props_Rs(Rs,ntrials);
    good_FR = cell_list_op(props1,{'good_zMI','good_Gauss'});
    good_FR = cell_list_op(props1,{'vals','good_zMI'});
    good_FR = cell_list_op(props1,{'vals'});
    [CRc,aCRc,mRR] = find_population_vector_corr(Rs,mR,good_FR,0);
    mRRm = [];
    for ii = 1:size(mRR,2)
        m_mRR_ii(ii) = min(min(mRR{an,ii}));        M_mRR_ii(ii) = max(max(mRR{an,ii}));
        m_CRc_ii(ii) = min(min(CRc{an,ii}));        M_CRc_ii(ii) = max(max(CRc{an,ii}));
        m_aCRc_ii(ii) = min(min(aCRc{ii}));        M_aCRc_ii(ii) = max(max(aCRc{ii}));
    end
    m_mRR = min(m_mRR_ii); M_mRR = max(M_mRR_ii); m_CRc = min(m_CRc_ii); M_CRc = max(M_CRc_ii); m_aCRc = min(m_aCRc_ii); M_aCRc = max(M_aCRc_ii);
    cbar_p_shift = [0.01 0.09 -0.03 -0.3];
    for ii = 1:size(ff.h_axes,2)
        axes(ff.h_axes(1,ii));
        ha = ff.h_axes(1,ii);
        szln = size(Rs_MC{an,ii}.sp_rasters,2)-1;
        mean_speed_over_trials = nanmean(Rs_MC{an,ii}.sp_rasters(:,1:szln))*pixelSize(ii);
        semspeed_over_trials = std(Rs_MC{an,ii}.sp_rasters(:,1:szln))/sqrt(10);
%         mean_speed_over_trials = nanmean(Rs{an,ii}.speed)*pixelSize(ii);
%         semspeed_over_trials = std(Rs{an,ii}.speed)/sqrt(10);
        xsMC = linspace(-2,2,size(mean_speed_over_trials,2));
        xsMC = linspace(-1.5,1.5,size(mean_speed_over_trials,2));
        plot(xsMC,mean_speed_over_trials);
        shadedErrorBar(xsMC,mean_speed_over_trials,semspeed_over_trials);
        set(gca,'XTick',[],'YTick',[0 8]);
        ylim([0 8]);
        xlim(xlimss);
        if ii == 1
          ylabel('MC \mum');
        end
        format_axes(gca);
        pos = get(ha,'Position');% ha = axes; set(ha,'units','inches');set(ha,'Position',pos);
        changePosition(ha,[0 pos(4)-0.63 0 -0.42]);
        
        %******* make new axes and plot mean and gaussian fitting
        pos = get(ha,'Position'); ha = axes; set(ha,'units','inches');set(ha,'Position',pos);
        changePosition(ha,[0 pos(4)+0.13 0 0]);
        szln = size(Rs{an,ii}.speed,2)-1;
        spud = Rs{an,ii}.speed(:,1:szln)*pixelSize(ii);
        mean_speed_over_trials = nanmean(spud);
        semspeed_over_trials = std(spud)/sqrt(10);
        xsMC = linspace(-2,2,size(mean_speed_over_trials,2));
        xsMC = linspace(-1.5,1.5,size(mean_speed_over_trials,2));
        plot(xsMC,mean_speed_over_trials);
        shadedErrorBar(xsMC,mean_speed_over_trials,semspeed_over_trials);
        set(gca,'XTick',[],'YTick',[0 30]);
        ylim([0 30]);
        xlim(xlimss);
        if ii == 1
        ylabel('Sp cm/s');
        end
        format_axes(gca);
        textstr = headstr{ii};set_axes_top_text_no_line(ff.hf,gca,textstr,headstrpos);
%         if ii == 1
%           
% %           textstr = sprintf('Motion Onset (MOn)'); set_axes_top_text_no_line(ff.hf,gca,textstr,[0 -0.093 0.05 0]);
%         else
%           textstr = sprintf('No-Brake (NB)'); set_axes_top_text_no_line(ff.hf,gca,textstr,[0 -0.14 0 0]);
% %           textstr = sprintf('Motion Offset (MOff)'); set_axes_top_text_no_line(ff.hf,gca,textstr,[0 -0.093 0.05 0]);
%         end
        %********************
        tR = Rs{an,ii};
        tmRR = mRR{an,ii};
        axes(ff.h_axes(2,ii));
        if ii < 5
            xdata = [0 75 150];
        else
            xdata = [-2 0 2];
            xdata = [-1.5 0 1.5];
        end
        ydata = [1 size(tmRR,1)];
        imagesc(xdata,ydata,tmRR,[m_mRR M_mRR]); set(gca,'Ydir','normal');
        if ii == 1
            hyl = ylabel('Cell #');
            changePosition(hyl,[25 0 0]);
        end
        set(gca,'YTick',[],'XTick',[]);
        
        ylims = ylim;
        set(gca,'YTick',[1 floor(ylims(2))],'XTick',[]);
        textstr = sprintf('%d',floor(ylims(2)));
%         set_axes_top_text_no_line(ff.hf,gca,textstr,[0 -0.11 0 0]);
        
        format_axes(gca);
        mM = min(tmRR(:)); MM = max(tmRR(:)); 
        if ii == 1
        [hc,hca] = putColorBar(gca,cbar_p_shift,[mM MM],5,'eastoutside',[0.09 0.11 0.09 0.16]);
        end
        
        axes(ff.h_axes(3,ii));
        tmRR = CRc{an,ii};
        if sum(isnan(tmRR(:))) > 0
            tmRR = fillmissing(tmRR,'linear',2,'EndValues','nearest');
            if sum(isnan(tmRR(:))) > 0
                tmRR = tmRR'; tmRR = fillmissing(tmRR,'linear',2,'EndValues','nearest'); tmRR = tmRR';
            end
        end
        
            xdata = [-2 0 2];
            xdata = [-1.5 0 1.5];

        ydata = [1 size(tmRR,1)];
        imagesc(xdata,ydata,tmRR,[m_CRc M_CRc]); set(gca,'Ydir','normal');
        if ii == 1
            ylabel('Time (s)');
        end
        set(gca,'YTick',[],'XTick',[]);

        ylims = ylim;
        textstr = sprintf('%d',ceil(ylims(2)));
%         set_axes_top_text_no_line(ff.hf,gca,textstr,[0 -0.07 0 0]);
        format_axes(gca);
        if ii == 1
        mM = -0.3;min(tmRR(:)); MM = max(tmRR(:)); [hc,hca] = putColorBar(gca,cbar_p_shift,[mM MM],5,'eastoutside',[0.09 0.11 0.09 0.16]);
        end
        
        axes(ff.h_axes(4,ii));
        tmRR = aCRc{ii};

%             xdata = [0 7.5 15];

        ydata = [1 size(tmRR,1)];
        imagesc(xdata,ydata,tmRR,[m_aCRc M_aCRc]); set(gca,'Ydir','normal');
        if ii == 1
            ylabel('Time (s)');
        end
        set(gca,'YTick',[]);
 
            xlabel('Time (s)');
   
        ylims = ylim;
        textstr = sprintf('%d',ceil(ylims(2)));
%         set_axes_top_text_no_line(ff.hf,gca,textstr,[0 -0.07 0 0]);
        format_axes(gca);
        if ii == 1
        mM = -0.3;min(tmRR(:)); MM = max(tmRR(:)); [hc,hca] = putColorBar(gca,cbar_p_shift,[mM MM],5,'eastoutside',[0.09 0.11 0.09 0.16]);
        end
    end
    save_pdf(ff.hf,mData.pdf_folder,sprintf('pop_vectors_%d_%c.pdf',ntrials,G),600);

%% for running stats on before after speeds
ii = 2
mean_speed_over_trials = nanmean(Rs_MC{an,ii}.sp_rasters(:,1:szln))*pixelSize(ii);
MC = Rs_MC{an,ii}.sp_rasters(:,1:szln);
group = [ones(1,18) 2*ones(1,18)];
[h,p,st] = ttest2(mean_speed_over_trials(find(group==1)),mean_speed_over_trials(find(group==2)))
% [within,dvn,xlabels,withinD] = make_within_table({'Cond','ET','CT'},[2,3,2]); withinD3 = withinD;
% dataT = make_between_table({avar},dvn);
% ra = RMA(dataT,within,{0.05,{'lsd','hsd','bonferroni'}});
% ra.ranova
% print_for_manuscript(ra)

%% population vector and correlation sensory for time rasters
  
    ff = makeFigureRowsCols(108,[6 3 2 2.9],'RowsCols',[3 2],'spaceRowsCols',[0.03 -0.01],'rightUpShifts',[0.1 0.11],'widthHeightAdjustment',[10 -85]);
%     set(gcf,'color','w');    set(gcf,'Position',[10 3 3.5 3.75]);
    MY = 8; ysp = 1; mY = 0; % responsive cells
    stp = 0.15*magfac; widths = ([0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5]+0.2)*magfac; gap = 0.19*magfac;
    adjust_axes(ff,[mY MY],stp,widths,gap,{'Cell #'});
    
    an = 1;G = 'C';
    si = [Ab_On Ar_On];
%     si = [Lb ArL_L];
%     si = [M_On M_Off];
    Rs = o.Rs(:,si);mR = o.mR(:,si);
    ntrials = 50;
    props1 = get_props_Rs(Rs,ntrials);
    good_FR = cell_list_op(props1,{'good_zMI','good_Gauss'});
    good_FR = cell_list_op(props1,{'vals','good_zMI'});
    good_FR = cell_list_op(props1,{'vals'});
    [CRc,aCRc,mRR] = find_population_vector_corr(Rs,mR,good_FR,0);
    mRRm = [];
    for ii = 1:size(mRR,2)
        m_mRR_ii(ii) = min(min(mRR{an,ii}));        M_mRR_ii(ii) = max(max(mRR{an,ii}));
        m_CRc_ii(ii) = min(min(CRc{an,ii}));        M_CRc_ii(ii) = max(max(CRc{an,ii}));
        m_aCRc_ii(ii) = min(min(aCRc{ii}));        M_aCRc_ii(ii) = max(max(aCRc{ii}));
    end
    m_mRR = min(m_mRR_ii); M_mRR = max(M_mRR_ii); m_CRc = min(m_CRc_ii); M_CRc = max(M_CRc_ii); m_aCRc = min(m_aCRc_ii); M_aCRc = max(M_aCRc_ii);
    cbar_p_shift = [-0.007 0.09 -0.03 -0.3];
    for ii = 1:size(ff.h_axes,2)
        tR = Rs{an,ii};
        tmRR = mRR{an,ii};
        axes(ff.h_axes(1,ii));
        if ii < 5
            xdata = [0 75 150];
        else
            xdata = [-2 0 2];
        end
        ydata = [1 size(tmRR,1)];
        imagesc(xdata,ydata,tmRR,[m_mRR M_mRR]); set(gca,'Ydir','normal');
        if ii == 1
            ylabel('Cell #');
        end
        set(gca,'YTick',[],'XTick',[]);

        ylims = ylim;
        textstr = sprintf('%d',floor(ylims(2)));
        set_axes_top_text_no_line(ff.hf,gca,textstr,[0 -0.13 0 0]);
        if ii == 1
          textstr = sprintf('Brake (B)'); set_axes_top_text_no_line(ff.hf,gca,textstr,[0 -0.093 0 0]);
%           textstr = sprintf('Motion Onset (MOn)'); set_axes_top_text_no_line(ff.hf,gca,textstr,[0 -0.093 0.05 0]);
        else
          textstr = sprintf('No-Brake (NB)'); set_axes_top_text_no_line(ff.hf,gca,textstr,[0 -0.093 0 0]);
%           textstr = sprintf('Motion Offset (MOff)'); set_axes_top_text_no_line(ff.hf,gca,textstr,[0 -0.093 0.05 0]);
        end
        format_axes(gca);
        mM = min(tmRR(:)); MM = max(tmRR(:)); 
        if ii == 1
        [hc,hca] = putColorBar(gca,cbar_p_shift,[mM MM],5,'eastoutside',[0.09 0.11 0.09 0.16]);
        end
        
        axes(ff.h_axes(2,ii));
        tmRR = CRc{an,ii};
        if sum(isnan(tmRR(:))) > 0
            tmRR = fillmissing(tmRR,'linear',2,'EndValues','nearest');
            if sum(isnan(tmRR(:))) > 0
                tmRR = tmRR'; tmRR = fillmissing(tmRR,'linear',2,'EndValues','nearest'); tmRR = tmRR';
            end
        end
        
            xdata = [-2 0 2];

        ydata = [1 size(tmRR,1)];
        imagesc(xdata,ydata,tmRR,[m_CRc M_CRc]); set(gca,'Ydir','normal');
        if ii == 1
            ylabel('Time (s)');
        end
        set(gca,'YTick',[],'XTick',[]);

        ylims = ylim;
        textstr = sprintf('%d',ceil(ylims(2)));
%         set_axes_top_text_no_line(ff.hf,gca,textstr,[0 -0.07 0 0]);
        format_axes(gca);
        if ii == 1
        mM = -0.3;min(tmRR(:)); MM = max(tmRR(:)); [hc,hca] = putColorBar(gca,cbar_p_shift,[mM MM],5,'eastoutside',[0.09 0.11 0.09 0.16]);
        end
        
        axes(ff.h_axes(3,ii));
        tmRR = aCRc{ii};

%             xdata = [0 7.5 15];

        ydata = [1 size(tmRR,1)];
        imagesc(xdata,ydata,tmRR,[m_aCRc M_aCRc]); set(gca,'Ydir','normal');
        if ii == 1
            ylabel('Time (s)');
        end
        set(gca,'YTick',[]);
 
            xlabel('Time (s)');
   
        ylims = ylim;
        textstr = sprintf('%d',ceil(ylims(2)));
%         set_axes_top_text_no_line(ff.hf,gca,textstr,[0 -0.07 0 0]);
        format_axes(gca);
        if ii == 1
        mM = -0.3;min(tmRR(:)); MM = max(tmRR(:)); [hc,hca] = putColorBar(gca,cbar_p_shift,[mM MM],5,'eastoutside',[0.09 0.11 0.09 0.16]);
        end
    end
    save_pdf(ff.hf,mData.pdf_folder,sprintf('pop_vectors_%d_%c.pdf',ntrials,G),600);



%% population vector and correlation sensory air
while 1
    selected_property = 'untuned';
    an = 1;
    titles = {'Air','Light','NB-A','NB-AL','NB-A*'};
    si = [Ab_On Lb Abs_On Ar_On ArL_On Ars_On];
%     si = [Ab_Off Abs_Off Ar_Off ArL_Off Ars_Off];
%     si = [Ar_D ArL_D Ars_D Ar_T ArL_T];
    Rs = o.Rs(:,si);mR = o.mR(:,si);
    ntrials = 50;
    props1 = get_props_Rs(Rs,ntrials);
    respA = props1.vals;
%     respA = props1.good_zMI;
    resp1 = cell_list_op(respA(:,1:2),[],'or');
    resp2 = cell_list_op(respA(:,3:4),[],'or');
    resp = respA;%[resp1(:,1:2),resp2(:,1:2)];
%     eval(cmdTxt);
    ffM = makeFigureRowsCols(107,[1 0.5 4 1],'RowsCols',[2 5],...
        'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.08 0.11],'widthHeightAdjustment',...
        [10 -80]);    set(gcf,'color','w');    set(gcf,'Position',[10 3 5 2]);
    sInds = 1:5;
    [CRc,aCRc1,mRR] = find_population_vector_corr(Rs(:,sInds),mR(:,sInds),resp(:,sInds),0);
    ff = ffM; ff.axesPos = ffM.axesPos(:,sInds); ff.h_axes = ffM.h_axes(:,sInds);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,sInds),mRR(an,:),CRc(an,:),[],[],0);
    for ii = 1:length(ff.h_axes(1,:)) ht = get_obj(ff.h_axes(1,ii),'title'); set(ht,'String',titles{ii}); end
 
    save_pdf(ff.hf,mData.pdf_folder,sprintf('PV_air_%s.pdf',selected_property),600);

    % average correlation of all animals
    ffM = makeFigureRowsCols(108,[1 0.5 4 1],'RowsCols',[1 5],...
        'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.08 0.24],'widthHeightAdjustment',...
        [10 -300]);    set(gcf,'color','w');    set(gcf,'Position',[10 7 5 1]);
    sInds = 1:5;
    ff = ffM; ff.axesPos = ffM.axesPos(:,sInds); ff.h_axes = ffM.h_axes(:,sInds);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,sInds),[],aCRc1,[],[],0);
%     sInds = 3:4;
%     ff = ffM; ff.axesPos = ffM.axesPos(:,sInds); ff.h_axes = ffM.h_axes(:,sInds);
%     ff = show_population_vector_and_corr(mData,ff,Rs(an,sInds),[],aCRc1,[],[]);
    save_pdf(ff.hf,mData.pdf_folder,sprintf('aPV_air_%s.pdf',selected_property),600);
    %%
    break;
end

%% population vector and correlation sensory light
while 1
    selected_property = 'untuned';
    an = 1;
    titles = {'B-L','B-L*','NB-L'};
    si = [Lb Lbs ArL_L];
    Rs = o.Rs(:,si);mR = o.mR(:,si);
    ntrials = 50;
    props1 = get_props_Rs(Rs,ntrials);
    respA = props1.good_FR_and_untuned;
    resp = respA;%[resp1(:,1:2),resp2(:,1:2)];
%     eval(cmdTxt);
    ffM = makeFigureRowsCols(107,[1 0.5 4 1],'RowsCols',[2 3],...
        'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.08 0.11],'widthHeightAdjustment',...
        [10 -80]);    set(gcf,'color','w');    set(gcf,'Position',[10 3 2.75 2]);
    sInds = 1:3;
    [CRc,aCRc1,mRR] = find_population_vector_corr(Rs(:,sInds),mR(:,sInds),resp(:,sInds),0);
    ff = ffM; ff.axesPos = ffM.axesPos(:,sInds); ff.h_axes = ffM.h_axes(:,sInds);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,sInds),mRR(an,:),CRc(an,:),[],[],0);
    for ii = 1:length(ff.h_axes(1,:)) ht = get_obj(ff.h_axes(1,ii),'title'); set(ht,'String',titles{ii}); end
    
%     sInds = 3:4;
%     [CRc,aCRc2,mRR] = find_population_vector_corr(Rs(:,sInds),mR(:,sInds),resp(:,sInds),0);
%     ff = ffM; ff.axesPos = ffM.axesPos(:,sInds); ff.h_axes = ffM.h_axes(:,sInds);
%     ff = show_population_vector_and_corr(mData,ff,Rs(an,sInds),mRR(an,:),CRc(an,:),[],[]);
%     for ii = 1:length(ff.h_axes(1,:)) ht = get_obj(ff.h_axes(1,ii),'title'); set(ht,'String',titles{ii}); end
    
    
    save_pdf(ff.hf,mData.pdf_folder,sprintf('PV_air_%s.pdf',selected_property),600);

    % average correlation of all animals
    ffM = makeFigureRowsCols(108,[1 0.5 4 1],'RowsCols',[1 3],...
        'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.08 0.24],'widthHeightAdjustment',...
        [10 -300]);    set(gcf,'color','w');    set(gcf,'Position',[10 7 2.75 1]);
    sInds = 1:3;
    ff = ffM; ff.axesPos = ffM.axesPos(:,sInds); ff.h_axes = ffM.h_axes(:,sInds);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,sInds),[],aCRc1,[],[],0);

    save_pdf(ff.hf,mData.pdf_folder,sprintf('aPV_air_%s.pdf',selected_property),600);
    %%
    break;
end


%% population vector and correlation sensory light + Control light
while 1
    selected_property = 'untuned';
    an = 1;
    titles = {'B-L','B-L*','NB-L'};
    si = [Ab_On Ab_Offc Ar_On Ar_Offc];
    si = [Ab_Offc Ar_Offc];
    Rs = o.Rs(:,si);mR = o.mR(:,si);
    ntrials = 50;
    props1 = get_props_Rs(Rs,ntrials);
    respA = props1.vals;
    resp = respA;%[resp1(:,1:2),resp2(:,1:2)];
%     eval(cmdTxt);
    ffM = makeFigureRowsCols(107,[1 0.5 4 1],'RowsCols',[2 4],...
        'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.08 0.11],'widthHeightAdjustment',...
        [10 -80]);    set(gcf,'color','w');    set(gcf,'Position',[10 3 3.5 2]);
    sInds = 1:2;
    [CRc,aCRc1,mRR] = find_population_vector_corr(Rs(:,sInds),mR(:,sInds),resp(:,sInds),0);
    ff = ffM; ff.axesPos = ffM.axesPos(:,sInds); ff.h_axes = ffM.h_axes(:,sInds);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,sInds),mRR(an,:),CRc(an,:),[],[],0);
    for ii = 1:length(ff.h_axes(1,:)) ht = get_obj(ff.h_axes(1,ii),'title'); set(ht,'String',titles{ii}); end
    
%     sInds = 3:4;
%     [CRc,aCRc2,mRR] = find_population_vector_corr(Rs(:,sInds),mR(:,sInds),resp(:,sInds),0);
%     ff = ffM; ff.axesPos = ffM.axesPos(:,sInds); ff.h_axes = ffM.h_axes(:,sInds);
%     ff = show_population_vector_and_corr(mData,ff,Rs(an,sInds),mRR(an,:),CRc(an,:),[],[]);
%     for ii = 1:length(ff.h_axes(1,:)) ht = get_obj(ff.h_axes(1,ii),'title'); set(ht,'String',titles{ii}); end
    
    
    save_pdf(ff.hf,mData.pdf_folder,sprintf('PV_air_%s.pdf',selected_property),600);

    % average correlation of all animals
    ffM = makeFigureRowsCols(108,[1 0.5 4 1],'RowsCols',[1 3],...
        'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.08 0.24],'widthHeightAdjustment',...
        [10 -300]);    set(gcf,'color','w');    set(gcf,'Position',[10 7 2.75 1]);
    sInds = 1:3;
    ff = ffM; ff.axesPos = ffM.axesPos(:,sInds); ff.h_axes = ffM.h_axes(:,sInds);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,sInds),[],aCRc1,[],[],0);

    save_pdf(ff.hf,mData.pdf_folder,sprintf('aPV_air_%s.pdf',selected_property),600);
    %%
    break;
end


%% population vector and correlation distance and time
while 1
    an = 1;
    titles = {'B-AOn','B-AOff','B-Arb','NB-AOn','NB-AOff','NB-Arb'};
    si = [Ab_On Ab_Off Ab_Offc Ars_On Ars_Off Ars_Offc];
    Rs = o.Rs(:,si);mR = o.mR(:,si);
    ntrials = 50;
    props1 = get_props_Rs(Rs,ntrials);
    resp = props1.vals;
%     resp = respA;%[resp1(:,1:2),resp2(:,1:2)];
%     resp = [respDT.exc respDT.exc];
%     resp = [respDT.inh respDT.inh];
%     eval(cmdTxt);
    ffM = makeFigureRowsCols(107,[1 0.5 4 1],'RowsCols',[2 6],...
        'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.08 0.11],'widthHeightAdjustment',...
        [10 -80]);    set(gcf,'color','w');    set(gcf,'Position',[10 3 6 2]);
    sInds = 1:6;
    [CRc,aCRc1,mRR] = find_population_vector_corr(Rs(:,sInds),mR(:,sInds),resp(:,sInds),0);
    ff = ffM; ff.axesPos = ffM.axesPos(:,sInds); ff.h_axes = ffM.h_axes(:,sInds);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,sInds),mRR(an,:),CRc(an,:),[],[],0);
    for ii = 1:length(ff.h_axes(1,:)) ht = get_obj(ff.h_axes(1,ii),'title'); set(ht,'String',titles{ii}); end
    for ii = 4:6 set_obj(ff.h_axes(2,ii),{'xtick',[1 68 136],'xticklabels',[0 7.5 15]}); set_obj(ff.h_axes(2,ii),{'ytick',[1 68 136],'yticklabels',[0 7.5 15]}); end 
    save_pdf(ff.hf,mData.pdf_folder,sprintf('aPV_spatial_temp.pdf'),600);

    % average correlation of all animals
    ffM = makeFigureRowsCols(108,[1 0.5 4 1],'RowsCols',[1 6],...
        'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.08 0.24],'widthHeightAdjustment',...
        [10 -300]);    set(gcf,'color','w');    set(gcf,'Position',[10 7 6 1]);
    sInds = 1:6;
    ff = ffM; ff.axesPos = ffM.axesPos(:,sInds); ff.h_axes = ffM.h_axes(:,sInds);
    ff = show_population_vector_and_corr(mData,ff,Rs(an,sInds),[],aCRc1,[],[],0);
    for ii = 4:6 set_obj(ff.h_axes(1,ii),{'xtick',[1 68 136],'xticklabels',[0 7.5 15]}); set_obj(ff.h_axes(1,ii),{'ytick',[1 68 136],'yticklabels',[0 7.5 15]}); end 
    save_pdf(ff.hf,mData.pdf_folder,sprintf('aPV_spatial_temp.pdf'),600);
    %%
    break;
end
