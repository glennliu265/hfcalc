
%HF Calc Ens Average

% Calculate Ensemble Average 

%% ------------------------------
%  User Input
%  ------------------------------
%fluxes = {'RHFLX','THFLX','NHFLX','SHFLX','LHFLX','FSNS','FLNS'};
fluxes   = {'nhflx','thflx','rhflx'}  ; %Flux to calculate ensemble average
fluxtype = 'nhflx'

lonremap  = 1          ; %Set 1 to remap LON from 360 -> -180 (Change [bbox] acc.)
monwin    = 3          ; %Months considered (1 or 3)
mode      = 4          ; % (1) No Sig Testing, (2) SST testing (3) Flx testing (4) Both
deg5      = 0          ;% Set to 1 to use smoothed data
ensorem   = 1          ;% Indicate if ENSO was removed
dof_man   = 82
p         = 0.20      

% ---------------
% Plotting Options

% Bounding Boxes (lon1 lon2 lat1 lat2)
bbox         = [-100 20 -25 75]; % North Atlantic
caxischoose  = [-50 50];
plotpt       = 0
topleft      = 0 ;% For seasonal plots, plot the first subplot ont the left
seasonplot   = 0
use_contour  = 0
p2005_plot   = 1
%cint = [-48:8:-2,2:8:48];
cint = [-50:10:50];
% -----------------------
% Set Paths and load data
projpath = '/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/';
addpath('/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/00_Commons/01_Data'); % Path to common data
addpath('/Users/gliu/') % Path to scripts
outpath = [projpath,'02_Figures/20200629/'];

if deg5 == 1 
    datpath  = [projpath,'01_Data/5deg/'];
    load('CESM1_LATLON_5deg.mat');
    lonsize = 72;
    latsize = 36;
else
    datpath  = [projpath,'01_Data/'];
    load('CESM1_LATLON.mat')
    lonsize = 288;
    latsize = 192;
end


% Load Data %
outname = [datpath,'ensavg_',fluxtype,'damping'...
    '_monwin',num2str(monwin),...
    '_sig',num2str(p*100,'%03d'),...
    '_dof',num2str(dof_man,'%03d'),...
    '_mode',num2str(mode),...
    '.mat'];
load(outname) 


% -----------------------------
% Shouldn't change these, but I'll put these here for now...
monlist   = [1:12]     ; %List of months 
ensnum    = [1:42]     ; %List of ensemble members 

% Labeling Info
if monwin == 3
    monlab = {'DJF','JFM','FMA',...
        'MAM','AMJ','MJJ',...
        'JJA','JAS','ASO',...
        'SON','OND','NDJ'};
elseif monwin == 1
    monlab = {'Jan','Feb','Mar',...
    'Apr','May','Jun',...
    'Jul','Aug','Sep',...
    'Oct','Nov','Dec'};
end

% Run Startup Script
startup
%% ------------------------------------------------------------------------
%  Script Start
%% ------------------------------------------------------------------------



%% Type 1: Seasonal Plots (Separately For Each Month)
if seasonplot == 1
    close all
    
    % Set month order and counters
    monorder = [12,1:11]
    cnt1 = 1; cnt2 = 1; cnt3 = 1; cnt4 = 1;

    seasons = {'Winter','Spring','Summer','Fall'};
    
    % General subplots by month
    for imon = 1:12

        % Start with December to make things easier
        m = monorder(imon); %

        % Get data for that month
        vplot = squeeze(ensavg(:,:,m));

        if m < 3 || m > 11
            fign = 1;     
            cnt1 = cnt1 + 1;
            cnt = cnt1;
        elseif m < 6
            fign = 2;
            season = 'Spring';
            cnt2 = cnt2 + 1;
            cnt = cnt2;
        elseif m < 9
            fign = 3;
            season = 'Summer';
            cnt3 = cnt3 + 1;
            cnt = cnt3;
        elseif m < 12
            fign = 4;
            season = 'Fall';
            cnt4 = cnt4 + 1;
            cnt = cnt4;
        end

        % Open Corresponding Figure Number
        figure(fig`n)
        hold on

        if topleft == 1 && cnt == 2
            cnt = 1
        end

        subplot(2,2,cnt)

        % Make Projection and Plot
        m_proj('Miller','long',[bbox(1) bbox(2)],'lat',[bbox(3)  bbox(4)])  
        %set(gca,'Color','k')
        m_pcolor(LON1,LAT,vplot')
        %m_contourf(LON1,LAT,vplot',[-50:10:-10, -2,2,10:10:50])

        % Color Axis Options
        colormap(m_colmap('diverging',11));
        caxis(caxischoose)

        %cints=(caxischoose(2) - caxischoose(1))/10;
        %[tn,tl] = divergingticks(caxischoose,11,cints);
        %xcolorbar('Location','southoutside','Ticks',tn,'TickLabels',tl)

        m_coast('patch',[.7, .7, .7], 'edgecolor','black');
        m_grid('tickdir','out','linewi',2,'backgroundcolor','g');

        if plotpt == 1
            title(char(monlab(m)))
        else
            title(char(monlab(m)))
        end



        % pos = get(gcf,'Position')


        %set(gca,'Position',[-676 1402 300 300])
    end

    % Finalize plot for each season and save
    for s = 1:4
        season = char(seasons(s));
        figure(s)
        sgtitle([season,' ',upper(fluxtype),' Feedback'])

        set(gcf,'PaperUnits','inches','Position',[10 10 450 475])
        
        % (ex. SeasonalPlots_dampingnhflx_monwin1_sig005_dof082_mode4.png)
        figname = [outpath,season,...
            '_damping',fluxtype,...
            '_monwin',num2str(monwin),...
            '_sig',num2str(p*100,'%03d'),...
            '_dof',num2str(dof_man,'%03d'),...
            '_mode',num2str(mode),...
            '.png'];
        saveas(gcf,figname,'png')

    end
    
    % Make an annual avg plot to generate the colorbar
    figure(5)
    vplot = mean(ensavg,3);
    
    % Make Projection and Plot
    m_proj('Miller','long',[bbox(1) bbox(2)],'lat',[bbox(3)  bbox(4)])  
    m_pcolor(LON1,LAT,vplot')
    %m_contourf(LON1,LAT,vplot',[-50:10:-10, -2,2,10:10:50])

    % Color Axis Options
    colormap(m_colmap('diverging',10));
    caxis(caxischoose)
    %colorbar('Location','southoutside','Ticks',[-50:10:50])
    %colorbar('Location','southoutside
    m_coast('patch',[.7, .7, .7], 'edgecolor','black');
    m_grid('tickdir','out','linewi',2,'backgroundcolor','g');
    title(['Annual Average ',upper(fluxtype),' Feedback'],'FontSize',20)

    figname = [outpath,'AnnAvg',...
        '_damping',fluxtype,...
        '_monwin',num2str(monwin),...
        '_sig',num2str(p*100,'%03d'),...
        '_dof',num2str(dof_man,'%03d'),...
        '_mode',num2str(mode),...
        '.png'];
    saveas(gcf,figname,'png')
    
end


%% Type 2: Park et al style

if p2005_plot == 1
    
    figure(1)
    clf
    
    seasons = {'Winter','Spring','Summer','Fall','AnnAvg'};
    mlabel  = {'DJF','MAM','JJA','SON','Annual Average'};
    monorder = [12,1:11]
    mcnt = 1;
    for s = 1:5
        
        % Index in intervals of three
        if s < 5
            m0 = mcnt; m1   = mcnt+2;
            vplot  = squeeze(mean(ensavg(:,:,m0:m1),3));
            mcnt = m1 + 1; % Add to counter
        
            % Prepare plot
            subplot(2,2,s)
        elseif s == 5
            vplot = squeeze(mean(ensavg,3));
            
            % Prepare plot
            figure(2)
            clf
        end
        
        
        % Get appropriate labels
        season = char(seasons(s));
        mlab   = char(mlabel(s));
        
        % Make plot
        m_proj('Miller','long',[bbox(1) bbox(2)],'lat',[bbox(3)  bbox(4)])
        if use_contour == 1
            [C,h] = m_contourf(LON1,LAT,vplot',cint);
            clabel(C,h,'LabelSpacing',250,'FontSize',10);
            colormap(m_colmap('diverging',10));
        else
            m_pcolor(LON1,LAT,vplot')
            colormap(m_colmap('diverging',10));
        end
        
       
        % Color Axis Options
        
        if s == 5
            colorbar('Location','southoutside','Ticks',cint)
            caxis(caxischoose)
        else
            %colorbar('Location','southoutside','Ticks',cint)
            caxis(caxischoose)
            title(mlab)
            
        end
    
        % Coastline and Grid Options
        m_coast('patch',[.7, .7, .7], 'edgecolor','black');
        m_grid('tickdir','out','linewi',2,'backgroundcolor','y');
        
    end
    
    % Save Figure 1 (Seasonal Plots)
    figure(1)
    sgtitle([upper(fluxtype),' Feedback, Ensemble and Lag Average'],'FontSize',20)
    
    set(gcf,'PaperUnits','inches','Position',[10 10 450 475])

    % (ex. SeasonalPlots_dampingnhflx_monwin1_sig005_dof082_mode4.png)
    figname = [outpath,'Park_etal_',...
        '_damping',fluxtype,...
        '_monwin',num2str(monwin),...
        '_sig',num2str(p*100,'%03d'),...
        '_dof',num2str(dof_man,'%03d'),...
        '_mode',num2str(mode),...
        '.png'];
    saveas(gcf,figname,'png')
    
    %Save Figure 2 (Annual Avg)
    figure(2)
    title([upper(fluxtype),' Feedback, Annual Average'],'FontSize',20)
    
    set(gcf,'PaperUnits','inches','Position',[10 10 450 475])

    % (ex. SeasonalPlots_dampingnhflx_monwin1_sig005_dof082_mode4.png)
    figname = [outpath,'Park_etal_AnnAvg_',...
        '_damping',fluxtype,...
        '_monwin',num2str(monwin),...
        '_sig',num2str(p*100,'%03d'),...
        '_dof',num2str(dof_man,'%03d'),...
        '_mode',num2str(mode),...
        '.png'];
    saveas(gcf,figname,'png')
end
        
        
        
        

   

