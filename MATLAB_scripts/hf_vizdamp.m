% viz_hfdamploop_ens
% Visualize damping coefficient, and generate figure for all 40 ensemble
% members
% Can specify...
%   1) months to loop over
%   2) lags to average over
%   3) ensnums to average over
%   4) bounding box of map


% Dependencies
% lon360to180.m
% findcoords.m
% m_map


%% ------------------------------
%  User Input
%  ------------------------------
%fluxes = {'RHFLX','THFLX','NHFLX','SHFLX','LHFLX','FSNS','FLNS'};
fluxes   = {'flns'}  ; %Flux to Visualize

lonremap  = 1          ; %Set 1 to remap LON from 360 -> -180 (Change [bbox] acc.)
monlist   = [1]       ; %List of months to visualize
avgmon    = 0          ; %Set 1 to average across all months in monlist
lag       = [1:2]      ; %List of lags to average over/sum over
ensnum    = [1:42]      ; %List of ensemble members to average over
monwin    = 3          ; %Months considered (1mon vs 3mon)
ensorem   = 1          ; %Set to 1 if enso was removed
savedamp  = 0          ; % Set to 1 to save damping variables after applying sigtest
mode      = 4          ; % (1) No Sig Testing, (2) SST testing (3) Flx testing (4) Both
ensavgf   = 1          ; %Set to 1 if you want to take the ensemble average first
lag1_elim = 0          ;% Set to 1 to eliminate rest of lags if lag 1 fails
insigNaN  = 0          ;% Set to 1 to change insignificant feedbacks to NaN
deg5      = 0          ;% Set to 1 to use smoothed data
%(3) to msk cor, (4) for both, (5) for old threshold, (6) for custom
%threshold

% -------------------
% MODE6: SET SIGNIFICANCE TESTING OPTIONS
% -------------------
% T-Test parameters
p        = 0.20         ; % p-value that will be used
tails    = 2           ; % 1 or 2-tailed T-test
dof_man  = 82            ; % Manual DOF value

% -------------------
% Point Plot toggle
% -------------------
plotpt = 0;
lon_find = -72;
lat_find = 35;
caxischoose = [-50 50];% Good caxis threshold for THFLX, noisy pixel id

% -------------------
% Plot Opt
% -------------------
plot_allens = 0; % Set to 1 to create separate plot for all ens members
plot_ensavg = 1; % Set to 1 to create plot of ensemble average
plot_oneens = 0; 
plot_lagmap = 0; % Set to 1 to create plot of lags removed
plot_mask   = 0;

% Set Paths
projpath = '/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/01_hfdamping/';
outpath  = [projpath,'02_Figures/20200629/'];

% Add Paths

addpath('/Users/gliu/Downloads/02_Research/01_Projects/01_AMV/00_Commons/01_Data'); % Path to common data
addpath('/Users/gliu/') % Path to scripts

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

addpath(datpath) % Path to data
% Load SST autocorrelation coefficients
load(strcat(datpath,['SST_rho_ensorem',num2str(ensorem),'_monwin',num2str(monwin),'.mat']))

% Bounding Boxes (lon1 lon2 lat1 lat2)
bbox         = [-100 20 -25 75]; % North Atlantic

% Visualization settings
figtype = 0;
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

% Other Settings
mnum = [1:35,101:107] ;%All ensemble members

startup
%% ------------------------------------------------------------------------
%  Script Start

%% -------------------------------------------
%  Prepare Figure Labels
%  -------------------------------------------           
% Labeling Options (Detect if AVG or Single)
if length(ensnum) > 1
    stensnum = 'Avg';
else
    stensnum  = num2str(ensnum,'%03d');
end

if length(lag) > 1
    stlag = 'Avg';
else
    stlag = num2str(lag);
end


% Determine effective DOF
if dof_man == 0
    if monwin == 1
        n_eff = 84;
    elseif monwin == 3
        n_eff = 82;
    end
% Manually set dof
else
    n_eff = dof_man
end


% Calculate Correlation Threshold (fixed dof)
ptilde  = (1-p/tails);
critval = tinv(ptilde,n_eff);
corrthres = sqrt(1/((n_eff/critval^2)+1));
                
% Loop for each element of fluxes
for i = 1:length(fluxes)
    
    % Load variable and pull flux
    fluxtype = char(fluxes(i));
    fluxlab = upper(fluxtype);
    %eval(['v = d',fluxtype,';'])
    
    
    % Load damping variable
    load(strcat(datpath,[fluxlab,'_damping_ensorem',num2str(ensorem),'_monwin',num2str(monwin),'.mat']))
    v = damping;
    
    % Load correlation coefficients
    load(strcat(datpath,[fluxlab,'_rho_ensorem',num2str(ensorem),'_monwin',num2str(monwin),'.mat']))
    rflx = rho;
    
    % -----------------------------------
    % Masking
    % -----------------------------------
    % Make masks  
    if mode == 1
        % No mask (identity matrix)
        amask = ones(size(v));
        
    else
        
        % Create masks with threshold
        msst = rsst   > corrthres;
        mflx = rflx   > corrthres;

        
        % Lag 1-based elimination
        if lag1_elim == 1
            
            % Eliminate points that failed lag 1 for...
            % Autocorrelation
            msst(:,:,:,:,2) = msst(:,:,:,:,2) .* msst(:,:,:,:,1);
            msst(:,:,:,:,3) = msst(:,:,:,:,3) .* msst(:,:,:,:,1);
            
            % Cross Correlation
            mflx(:,:,:,:,2) = mflx(:,:,:,:,2) .* mflx(:,:,:,:,1);
            mflx(:,:,:,:,3) = mflx(:,:,:,:,3) .* mflx(:,:,:,:,1);            

        end
        
        % Selectively apply testing
        switch mode
            % Apply autocorrelation only
            case 2
                amask = msst;
            % Apply cross correlation only
            case 3
                amask = mflx;
            % Apply Both
            case 4
                amask = mflx .* msst;      
        end
        % End mask creation
    end
    
    % Turn 0 to NaNs if option is set
    if insigNaN == 1
        amask = double(amask);
        amask(amask==0) = NaN;
    end
    
    % Apply Mask
    v = v .* amask;
    
    % ----------------------------------------------
    % Further Adjustments
    % ----------------------------------------------
    
    % Check size and adjust variable for plotting
    chksize = size(v)
    refsize = [lonsize, latsize, 12, 42, 3];
    if isequal(chksize,refsize) ~= 1
        fprintf('Permuting variable to match script requirements')
        corrsize = zeros(1,length(refsize));
            for i = 1:length(refsize)
                corrsize(1,i) = find(chksize==refsize(i));
            end
        v = permute(v,corrsize);
    end
    
    % Remap LON to -180-180 if toggle is on
    if lonremap == 1
        if any(LON> 181)        
            [LON1,v] = lon360to180(LON,v);
            [~,amask] = lon360to180(LON,amask);
        else
            v = v;
            LON1 = LON;
        end
        
        % Find Corresponding Long/Lat
        [oidx,aidx] = findcoords(lon_find,lat_find,2,{LON1,LAT}) 
    else
        v = v;
        LON1 = LON;
    end
    
    
    % Find midpoint values
    midptx = LON1(oidx) + (LON1(oidx+1)-LON1(oidx))/2;
    midpty = LAT(aidx)  + (LAT(aidx+1)-LAT(aidx))/2;
    
    % Close any existing figures
    close all
    
    %% ---------------------------------
    %  Loop by month... 
    %  ---------------------------------
    for m = 1:length(monlist)
        
        % Average over all months if set.
        if avgmon == 1
            imon = monlist;
        else
            imon = monlist(m);
        end
        
        figure(m)
        clf
        
        % Save data if m = monlist
        if plot_oneens == 1
            
            ensavg = squeeze(v(:,:,imon,ensnum,lag));
        else
            ensavg = squeeze(v(:,:,imon,:,lag));
        end
        
        if ensavgf == 1
            fprintf("Ensavgf set to 1\n")
            while length(size(ensavg)) > 2 
                ensavg = squeeze(nanmean(ensavg,3));
            end        
        elseif ensavgf == 0
            fprintf("Ensavgf set to 0\n")
            ensavg = ensavg(:,:,:);
            ensavg = nanmean(ensavg,3);
        end
        
        %% ----------------------------------------
        %  Create Separate Plots for each ENS Mem.
        %  ----------------------------------------
        if plot_allens == 1
            pvals = ones(1,42);
            for e = 1:length(mnum)
                ensnum = mnum(e);

                %% -------------------------------------------
                %  1. Restrict/Avg Variable for Lag/Mon/Ensnum
                %  -------------------------------------------
                %Select by month, ensemble member, and lag
                if plot_mask == 1
                    movstd = squeeze(msst(:,:,e,imon,lag));
                else
                    movstd = squeeze(v(:,:,imon,e,lag));
                    %plotmask = squeeze(amask(:,:,imon,e,lag));
                end

                %Take mean along third dimension until only lat/lon remain
                while length(size(movstd)) > 2 
                    movstd = squeeze(nanmean(movstd,3));
                    %plotmask = squeeze(nanprod(plotmask,3));
                end 

                %% -------------------------------------------
                %  2. Create Figure (42 Member Plot)
                %  -------------------------------------------
                
                subplot(6,7,e)
                %subaxis(6,7,e,'Spacing',0.5)
                
                
                m_proj('Miller','long',[bbox(1) bbox(2)],'lat',[bbox(3)  bbox(4)]) 
                m_pcolor(LON1,LAT,movstd')
                colormap(m_colmap('diverging',10));
                caxis(caxischoose)
                m_coast('patch',[.7, .7, .7], 'edgecolor','black');

                % Plot stippling to show significance
                %stipple(LON1,LAT,plotmask'>0)
                
                if plotpt == 1           
                    pval = squeeze(movstd(oidx,aidx));
                    title(['e',num2str(ensnum,'%03d'),': ',num2str(pval,'%.02f')]);
                    pvals(:,e) = pval;
                else
                    title(['ens: ',num2str(ensnum,'%03d')]);
                end
                
                % Set ticks and figure size
                %m_grid('xticklabels',[],'yticklabels',[],'tickdir','out','linest','none','linewi',1,'backgroundcolor','g');
                set(gca,'Color','y')
                set(gca,'YTickLabel',[]);
                set(gca,'XTickLabel',[]);
                set(gcf,'Position',[388 -8 1094 804])
                set(gcf, 'InvertHardcopy', 'off')
                
                % Add Point Marker
                hold on

                if plotpt == 1
                    m_plot(midptx,midpty,'rx',...
                        'MarkerSize',25,'MarkerEdgeColor','y')

                    m_plot(midptx,midpty,'ro',...
                        'MarkerSize',25,'MarkerEdgeColor','y')
                end
                fprintf('\nCompleted month %s for lag %s, ENSNUM %s',num2str(imon,'%02d'),stlag,stensnum)
                % End loop for ann ensemble members
            end
            %% -------------------------------------------
            %  3. Save 42-member plot
            %  -------------------------------------------
            if plot_mask == 1
                figname = [outpath,fluxlab,'mask_mo',num2str(imon,'%02d'),...
                    '_lag',stlag,...
                    '_ensALL',...
                    '_monwin',num2str(monwin),...
                    '_ensorem',num2str(ensorem),...
                    '_mode',num2str(mode),...
                    '_sig',num2str(p*100,'%03d'),...
                    '.png'];       
            else
                figname = [outpath,fluxlab,'damping_mo',num2str(imon,'%02d'),...
                    '_lag',stlag,...
                    '_ensALL',...
                    '_monwin',num2str(monwin),...
                    '_ensorem',num2str(ensorem),...
                    '_mode',num2str(mode),...
                    '_sig',num2str(p*100,'%03d'),...
                    '.png'];
            end

            sgtitle(strcat(fluxlab,' Damping; Mon ',monlab(m),...
                '; Lag ',stlag,...
                '; All Ens'),...
                    'FontSize',20);
            set(gcf,'color','w')
            saveas(gcf,figname,'png')   
            fprintf('\nCompleted month %s',num2str(imon,'%02d'))
            % End conditional for plot allens
        end

        %% ----------------------------------------
        %  Create ENS Avg Plot
        %  ----------------------------------------
        if plot_ensavg == 1   
            figure(12+m)
            clf
            m_proj('Miller','long',[bbox(1) bbox(2)],'lat',[bbox(3)  bbox(4)]) 
            set(gca,'Color','k')
            m_pcolor(LON1,LAT,ensavg')
            colormap(m_colmap('diverging',10));
            caxis(caxischoose)
            
            %cints=(caxischoose(2) - caxischoose(1))/10
            %[tn,tl] = divergingticks(caxischoose,11,cints)
            
            
            %colorbar('Location','southoutside','Ticks',tn,'TickLabels',tl)
            colorbar('Location','southoutside','Ticks',[-50:10:50])
            m_coast('patch',[.7, .7, .7], 'edgecolor','black');
            m_grid('tickdir','out','linewi',2,'backgroundcolor','yellow');
            pval = ensavg(oidx,aidx);
            
            if plot_oneens == 1
                stensnum = num2str(ensnum(1),'%02d');
            end
            
            figname = strcat(outpath,fluxlab,'damping_mon',num2str(imon,'%02d'),...
                    '_lag',stlag,...
                    '_ens',stensnum,...
                    '_monwin',num2str(monwin),...
                    '_ensorem',num2str(ensorem),...
                    '_mode',num2str(mode),...
                    '_sig',num2str(p*100,'%03d'),...
                    '.png');
            % Add Point Marker
            hold on
            if plotpt == 1
                m_plot(midptx,midpty,'yx',...
                    'MarkerSize',25,'MarkerEdgeColor','y')

                m_plot(midptx,midpty,'yo',...
                    'MarkerSize',25,'MarkerEdgeColor','y')
                
                title({strcat(fluxlab,' Damping; Mon ',monlab(m),...
                    '; Lag ',stlag,...
                    '; Ens  ',stensnum);['Corr. Thres: ',num2str(corrthres,'%.04f'),'; Pt Value: ',num2str(pval,'%.02f')]},...
                    'FontSize',16)
            else
                title({strcat(fluxlab,' Damping; Mon ',num2str(imon,'%02d'),...
                    '; Lag ',stlag,...
                    '; Ens  ',stensnum)},...
                    'FontSize',16)
            end
            set(gcf,'PaperUnits','inches','Position',[10 10 450 475])
            saveas(gcf,figname,'png')
            fprintf('\nCompleted month %s',num2str(imon,'%02d'))


        end     
        
        % Exit loop if avgmon is selected
        if avgmon == 1
            break
        end
    end
end
