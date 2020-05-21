% -------------------------------------------------
% h1_enavgrm
% -------------------------------------------------
% First 2 steps of hfdamping calculation
% Load in variables for all 42 ensemble members
% Applies precalculated land/ice mask in maskloc, then...
% 
% 1) Remove Ensemble Mean
% 2) Remove Climatological Mean
% 3) Save Output (EX: outpath/LHFLX_ens035_proc.mat)
%
% Input:p
% 1) vars     - variables to process {'TS','LHFLX','SHFLX','FSNS','FLNS'}
% 2) expr     - experiment ('HTR', 'PIC',....)
%
% Output : [outpath,vname,'_ens',num2str(ensnum,'%03d'),'_proc.mat']
% 1) varenavg - 42-member Ensemble Average []
% 2) climmon  - Climatological Mean (mean over year dimension) 
% 3) varrm    - Processed data for the ensemble member (lon x lat x m x time)
% 4) step1    - Step 1 Flag (1 for enavg removal)
% 5) step2    - Step 2 Flag (1 for clim mean removal)
%
%% User Input ------------
vars = {'TS','LHFLX','SHFLX','FSNS','FLNS'};
expr = 'HTR';

% Path to data (monthly folder before the variable addition)
datpath = '/vortex/jetstream/climate/data1/yokwon/CESM1_LE/downloaded/atm/proc/tseries/monthly/';
outpath = '/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/hfdamping_matfiles/01_hf1out/';

% Other sets
yrs   = 1850:2005;
mnum  = [1:35,101:107]                   ;% Ensemble Members List

% Mask Location
maskloc  = '/home/glliu/01_Data/masks/'     ;% Location of Landice masks

% Toggles
step1 = 1; % Ens Avg Removal
step2 = 1; % Climatological Mean Removal
%% Script Start ------------
fprintf('Running hf1_enavgrm (%s)\n',datetime('now'))
allstart = datetime('now');

% Load startup
addpath('/home/glliu/')
startup

% Get Index of start year for first ensemble member
tstart = (find(yrs==1920)-1)*12+1;

% Get total size of time dimension
tsize = (yrs(end) - 1920 + 1) * 12;

% Loop by variable
for v = 1:length(vars)
    
    % Get Variable Name
    vname = char(vars(v));
    ncnames = build_e11name(expr,vname);
      
    % Preallocate
    varall = NaN(288,192,tsize,length(ncnames));
    
    for n = 1:length(ncnames)
        lstart = datetime('now');       
        
        % Get Ensemble Number
        ensnum = mnum(n);
        fprintf('%s: Loading ENS %s...',vname,num2str(ensnum,'%02d'))
        
        % Read in variable
        ncpath = strcat(datpath,vname,'/',ncnames(n));
        readvar = ncread(ncpath,vname);
               
        % Reduce dimensions if ensemble member 1
        if n == 1
            readvar = readvar(:,:,tstart:end);
        end
        
        % -------------------
        % Apply Land/Ice Mask
        % -------------------
        % Load in land/ice mask
        maskmat = [maskloc,'landicemask_ensnum',num2str(ensnum,'%03d'),'.mat'];
        load(maskmat)
        
        % Replicate mask to time period and apply max
        mask_landice = mask_landice(:,:,ones(1,tsize));
        readvar = readvar .* mask_landice;    
        
        % Store in master variable
        varall(:,:,:,n) = readvar;  
        
        % Timekeeping
        lend = datetime('now');
        elapsed = lend-lstart;
        fprintf('DONE! (%s | %s)\n',elapsed,lend-allstart)
    end
    
    %% ---------------------------
    %  Step 1. Ens Avg Removal
    %  ---------------------------
    if step1 == 1        
        % Calculate Ensemble Average
        varenavg = nanmean(varall,4);
    
        % Remove Ensemble Average
        varall = varall - varenavg;
    
    else
        % Create dummy ensavg variable
        varenavg = ones(size(varall(:,:,:,1)));
    end
    
    %% ---------------------------
    %  Step 2. Anomaly Calculation
    %  ---------------------------
    % Reshape to separate mon and year
    varanom = reshape(varall,288,192,12,size(varall,3)/12,length(ncnames));
    
    if step2 == 1
        % Calculate climatological mean
        climmean = nanmean(varanom,4);
    
        % Remove Climatological Mean
        varanom = varanom - climmean;
    else
        climmean = ones(size(varanom));
    end
    
    %% ---------------------------
    %  Save Variables
    %  ---------------------------
    % Loop by Ens number and save
    for n = 1:length(ncnames)
        sstart = datetime('now');
        
        
        % Get Ensemble Number
        ensnum = mnum(n);
        fprintf('%s: Saving ENS %s...',vname,num2str(ensnum))
        
        % Get data
        varrm   = squeeze(varanom(:,:,:,:,n));
        climmon = squeeze(climmean(:,:,:,:,n));
        
        % Save variables
        outname_var = [outpath,vname,'_ens',num2str(ensnum,'%03d'),'_proc.mat'];
        save(outname_var,'varenavg','varrm','climmon','step1','step2')
        
        sed = datetime('now');
        elapsed = sed-sstart;
        fprintf('DONE! (%s | %s)\n',elapsed,sed-allstart)
    end
    fprintf('Done with %s! (%s)\n',vname,sed-allstart)
end

allend = datetime('now');
elapsed = allend-allstart;
fprintf('hf1_enavgrm ran in %s (%s)\n',elapsed,datetime('now'))