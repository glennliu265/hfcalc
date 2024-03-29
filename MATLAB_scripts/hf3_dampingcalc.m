    % -------------------------------------------------
% hf3_dampingcalc (and correlation calculation)
% -------------------------------------------------
%
% This script is split into two parts
%
% Input
% 1. vars    - List of Fluxes/Variables {TS,'LHFLX','SHFLX','FSNS','FLNS'}
% 2. ensolag - Lag in months between Variable and PC before the regression
% 3. monwin  - Window of months to consider (1 or 3)
% 4. pcrem   - PCs to remove [1,2,3]
% 5. EOFcorr - Set to 1 to correct the sign of the EOF based on predef
% boxes
%
% Output
% 1. outfile1 - Damping Coefficients for each flux
% 2. outfile2 - covariance for each flux and SST autocovariance
%

% Note: All the files are automatically stacked into the same file.
%% User Input ------------
vars    = {'LHFLX','SHFLX','FLNS','FSNS','RHFLX','THFLX','NHFLX'}; % Note, remove SST calculation
lags    = [1:3]         ; % Lag between variable and SST [1,2,3]
monwin  = 3             ; % Months to include (ex: 1 = D, 3 = "NDJ")
skipto  = 1             ; % Skip to specific member
stopat  = 43            ; % Stop Iteration at given number
timerng = 1921:2004     ; % Range of Years to Include
startyr = 1921          ; %
deg5    = 1             ; % Toggle for the smoothed data

% T-Test parameters
pval     = 0.05          ; % p-value that will be used
dof_type = 4             ; % Type of edof used
ensorem  = 1             ; % 1 = use variables with enso rem
tails    = 2             ; % 1 or 2-tailed T-test
dof_man  = 82            ; % Manual DOF value

% Path to data (monthly folder before the variable addition)
if deg5 == 1
     
    datpath = '/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/hfdamping_matfiles/02_hf2out/5deg/';
    latlon  = '/home/glliu/01_Data/CESM1_LATLON_5deg.mat';
    
    outpath = ['/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/hfdamping_matfiles/03_hf3out/5deg/',num2str(monwin),'mon/'];
    lonsize = 72;
    latsize = 36;
else
    if ensorem == 1
        datpath = '/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/hfdamping_matfiles/02_hf2out/ensorem/';
    else
        datpath = '/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/hfdamping_matfiles/01_hf1out/';
    end
    latlon  = '/home/glliu/01_Data/CESM1_LATLON.mat';
    
    % OutPath
    outpath = ['/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/hfdamping_matfiles/03_hf3out/',num2str(monwin),'mon/'];
    
    lonsize = 288;
    latsize = 192;
end



% Other sets
mnum  = [1:35,101:107];% Ensemble Members List

% Toggles
savedamping = 1;% Set to 1 to save damping values
savecov     = 0;% Set to 1 to save covariance and autocov
savetest    = 0;% Set to 1 to save significance test results (rsst)
saverho     = 1;% Set to 1 to save map of cross-correlation coefficients
savesst     = 1;% Set to 1 to save SST autcorrelation
%% -----------------------------------------------------------------------
%  Script Start/Setup
%  -----------------------------------------------------------------------

allstart = datetime('now');
fprintf('Now running hf3_dampingcalc (%s)',allstart)
fprintf('\n Save Options -- Damping:%i | Cov/AutoCov:%i | SigTest:%i | Corr:%i | SST:%i',savedamping,savecov,savetest,saverho,savesst)
fprintf('\n\t Save Loc: %s',outpath)

% Load startup
addpath('/home/glliu/')
startup

% Load Lat/lon
load(latlon)

% Set windowsize
if monwin == 3
    winsize = 1;
elseif monwin == 1
    winsize = 0;
end

% Preallocate Arrays for SST autocorrelation 
if savesst == 1
    % Preallocate for SST Autocorrelation Testing Results
    if savetest == 1
        sigval    = NaN(1,  length(vars)+1)       ; % Critical Correlation Value
        st_aut    = NaN(lonsize,latsize,12,length(lags))  ; % Significance Test Result
    end
    
    % Preallocate for SST Autocorrelation Coefficients
    if saverho == 1
          rho_aut   = NaN(lonsize,latsize,12,length(lags))  ; % Correlation Coefficient
    end
end

% Loop by ensemble member
for n = 1:length(mnum)
    
    % Skip to specific ensemble member if option is set
    if n < skipto
        continue
    end
    
    % Stop at specific ensemble member if option is set
    if n >= stopat
        fprintf('\nStopping now at %i (%s)',stopat,datetime('now'))
        break
    end
    
    
    nstart = datetime('now');
    ensnum = mnum(n);
    
    % Load file depending on if enso was removed
    if ensorem == 1
        % Load in SST variable [12 84 lonsize latsize]
        matname = [datpath,'TS_noenso_ENS',num2str(ensnum,'%03d'),'.mat'];
        load(matname)
    
        % Reassign variable to SST 
        sst = var_noenso;
    else
        % Load from hf1 output
        datpath = '/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/hfdamping_matfiles/01_hf1out/';
        matname = [datpath,'TS_ens',num2str(ensnum,'%03d'),'_proc.mat'];
        load(matname)
        
        % Reassign variable to SST
        sst = varrm;
        sst = permute(sst,[3,4,1,2]);
    end
    
    % -------------------------------
    % Restrict to desired time period
    % -------------------------------
    % Make Year Vector
    numyrs = size(sst,2);
    yrfull = startyr-1+[1:numyrs];
    
    % Find Years
    yrsel  = find(yrfull >= timerng(1) & yrfull <= timerng(end));
    
    % Select accordingly
    sst    = sst(:,yrsel,:,:);
    
    for v = 1:length(vars)
        vname = char(vars(v));
        
        % Preallocate arrays for cross correlation
        if savetest == 1       
            st_cov     = NaN(lonsize,latsize,12,length(lags))  ; % Significance Test Result 
        end
        
        % Preallocate arrays for cross-corr coefficients
        if saverho == 1
            rho_cov   = NaN(lonsize,latsize,12,length(lags))  ; % Correlation Coefficient
        end
        
        % Note: load in and sum fluxes (if NHFLX RHFLX or THFLX)
        if strcmp(vname,'NHFLX') == 1
            flx = NaN(12,84,lonsize,latsize,4);
            for vi = 1:4
                mult = 1;
                switch vi
                    case 1
                        vn = 'LHFLX';
                    case 2
                        vn = 'SHFLX';
                    case 3
                        vn = 'FLNS';
                    case 4
                        vn = 'FSNS';
                        mult = 1*-1;
                end         
                matname = [datpath,vn,'_noenso_ENS',num2str(ensnum,'%03d'),'.mat'];
                load(matname)
                if mult < 1
                    fprintf('\t Flipping sign for %s',vn)
                end
                flx(:,:,:,:,vi) = var_noenso*mult;
                clear var_noenso
            end
            flx = nansum(flx,5);
            
        elseif strcmp(vname,'RHFLX') == 1
            flx = NaN(12,84,lonsize,latsize,4);
            for vi = 3:4
                mult = 1;
                switch vi
                    case 3
                        vn = 'FLNS';
                    case 4
                        vn = 'FSNS';
                        mult = 1*-1;
                end         
                matname = [datpath,vn,'_noenso_ENS',num2str(ensnum,'%03d'),'.mat'];
                load(matname)
                if mult < 1
                    fprintf('\t Flipping sign for %s',vn)
                end
                flx(:,:,:,:,vi) = var_noenso*mult;
                clear var_noenso
            end
            flx = nansum(flx,5);
            
        elseif strcmp(vname,'THFLX') == 1
            flx = NaN(12,84,lonsize,latsize,4);
            for vi = 1:2
                switch vi
                    case 1
                        vn = 'LHFLX';
                    case 2
                        vn = 'SHFLX';
                end         
                matname = [datpath,vn,'_noenso_ENS',num2str(ensnum,'%03d'),'.mat'];
                load(matname)
                flx(:,:,:,:,vi) = var_noenso*mult;
                clear var_noenso   
            end
            flx = nansum(flx,5);
            
        else
            % Load in flux variable + reassign
            if ensorem == 1
                matname = [datpath,vname,'_noenso_ENS',num2str(ensnum,'%03d'),'.mat'];
                load(matname)
                flx = var_noenso;
            else
                matname = [datpath,vname,'_ens',num2str(ensnum,'%03d'),'_proc.mat'];
                load(matname)
                flx = varrm;
                flx = permute(flx,[3,4,1,2]);
            end
        end

        % Restrict to time period
        flx = flx(:,yrsel,:,:);
            
        % Preallocate Arrays (Oof this might take a lot of mem?)
        if savedamping == 1
            hfdamping = NaN(lonsize,latsize,12,3);
        end
        if savecov == 1 
            covall    = NaN(lonsize,latsize,12,3);
            autall    = NaN(lonsize,latsize,12,3);
        end       

        % Loop by lag
        for lag = 1:length(lags)
            %% ----------------------------------------------------------------
            %  Indexing
            %  ----------------------------------------------------------------
            l = lags(lag);
            
            % Crossing Flags           
            tflag   = 0;         
            tshift  = 0;
            
            % If windowsize > 0, reduce leading end to allow for year
            % crossing acommodation
            if winsize > 0
                lshift  = -1;
            else
                lshift  = 0;
            end
            
             % Total time series length
            totyr  = size(sst,2);
            
            % Print Message
            fprintf('\nComputing Damping for ENS:%s Var:%s Lag:%i...',...
                    num2str(ensnum,'%03d'),vname,l)
                
            % Loop by month
            for m = 1:12
                mstart = datetime('now');
                
                
                % Get Lag Month
                lm = mod(m - l,12);
                if lm == 0
                    lm = 12;
                end
                
                % When leading edge of the base variable's window crosses 
                % into next year, shift the end of the period forward
                % (When varidx -> ND|J)
                if m + winsize > 12
                    lshift = lshift + 1;
                end
                
                % When the trailing end of the lagged variable crosses into
                % the central period, shift the start of the period forward
                
                if tflag == 1
                    tshift = tshift + 1;
                    tflag = 0;
                end
                
                % Flag for tshift based on lag month and window size
                % (When pcidx -> D|JF)
                if winsize > 0
                    if lm - winsize == 0 && tflag == 0
                        tflag = 1;                 
                    end
                else
                    if lm == 12
                        tflag = 1;
                    end
                end
                
                % Produce time period
                base_ts = 1+tshift:totyr+lshift;
                
                % Index according to lag
                [flxidx,basflag] = laglead_indexer(flx,m,winsize,base_ts);
                [sstidx,~] = laglead_indexer(sst,m,winsize,base_ts);
                [sstlagidx,lagflag] = laglead_indexer(sst,lm,winsize,base_ts);
                
                % Reduce yr dim (correct time periods) based on flag
                % If both months are NOT in the interior
                if lagflag+basflag ~= 0 || tshift == 0
                    % Take the interior period if month has crossed the initial
                    % threshold (trailing end of window no longer in the 
                    % previous year) (sstlagidx -> D|JF ; sst/flxidx -> |JFM)
                    if basflag == 0
                        flxidx = flxidx(:,2:end,:,:);
                        sstidx = sstidx(:,2:end,:,:);
                    end
                
                    % Take the interior period once the leading edge of the
                    % window of the base variable has crossed into the next
                    % year (sstlagidx -> OND| ; sst/flxidx -> ND|J)
                    if lagflag == 0
                        sstlagidx = sstlagidx(:,1:end-1,:,:);
                    end
                end
                
                % Combine time dimensions to create 2D matrix
                vsize = size(sstidx);
                sst2d = reshape(sstidx,prod(vsize(1:2)),prod(vsize(3:4)));
                flx2d = reshape(flxidx,prod(vsize(1:2)),prod(vsize(3:4)));
                sstlag2d = reshape(sstlagidx,prod(vsize(1:2)),prod(vsize(3:4)));
                
                                
                %% ---------------------------------------------
                %  Compute Correlation Coefficients and P Values
                %  ---------------------------------------------
                
                % Cross Correlation Testing
                if savetest == 1             
                    % Compute cross correlation coefficients and test
                    % results
                    [~,~,cc_test,~,~,cc_thres] = pcorr(flx2d,sstlag2d,1,pval,tails,dof_type,dof_man);
                    
                    % Store values        
                    st_cov(:,:,m,l)      = reshape(cc_test,lonsize,latsize);
                    sigval(1,v)              = cc_thres;
                end
                
                % Cross Correlation Coefficients
                if saverho == 1
                    [cc_rho,~,~,~,~,~] = pcorr(flx2d,sstlag2d,1,pval,tails,dof_type,dof_man);
                    rho_cov(:,:,m,l)     = reshape(cc_rho,lonsize,latsize);
                end
                
                % ---------------------------------------------------------
                % SST Autocorrelation: Do on first iteration if opt is set
                % ---------------------------------------------------------
                if savesst == 1 && v == 1
                    if savetest == 1
                        [~,~,ac_test,~,~,ac_thres]    = pcorr(sst2d,sstlag2d,1,pval,tails,dof_type,dof_man);
                        st_aut(:,:,m,l)        = reshape(ac_test,lonsize,latsize);
                        sigval(1,length(vars)+1) = ac_thres;
                    end
                    
                    if saverho == 1
                        [ac_rho,~,~,~,~,~]    = pcorr(sst2d,sstlag2d,1,pval,tails,dof_type,dof_man);
                        rho_aut(:,:,m,l)     = reshape(ac_rho,lonsize,latsize);
                    end
                end
                
                %% ---------------------------------------------
                %  Compute Damping
                %  ---------------------------------------------
                if savedamping ==1
                    % -------------------------------------------
                    % Remove time mean
                    % -------------------------------------------
                    ssta = sst2d - nanmean(sst2d,1);
                    flxa = flx2d - nanmean(flx2d,1);
                    sstlaga = sstlag2d - nanmean(sstlag2d,1);

    
                    % ---------------------------------------
                    % Calculate covariance and autocovariance
                    % ---------------------------------------
                    cov = nansum(flxa.*sstlaga)./size(flxa,1);
                    aut = nansum(ssta.*sstlaga)./size(ssta,1);
                
                    % ---------------------------------------
                    % Calculate damping coefficient
                    % ---------------------------------------
                    damping = cov ./ aut;
                    
                    % ---------------------------------------
                    % Store variables
                    % ---------------------------------------
                    hfdamping(:,:,m,l) = reshape(damping,lonsize,latsize)     ; % Damping Coeffs
                    
                    % Also save covariance if option is set
                    if savecov == 1
                        covall(:,:,m,l)    = reshape(cov,lonsize,latsize)     ; % Covariance (FLX, SSTlag)
                        
                        % Save SST autocovariance on first iteration
                        if v == 1
                           autall(:,:,m,l)    = reshape(aut,lonsize,latsize)  ; % Autocovariance (SST, SSTlag)
                        end
                    end
                end

                
            end
            elapsed = datetime('now') - mstart;
            fprintf('Complete! (%s | %s)',elapsed,datetime('now')-allstart)
        end
        
        % Store in named variables
        if savedamping == 1
            eval([vname,'_damping =hfdamping;'])
        end
        
        if savecov == 1
            eval([vname,'_cov =covall;'])
        end
        
        if savetest == 1
            eval([vname,'_test =st_cov;'])  
        end
        
        if saverho == 1
            eval([vname,'_rho =rho_cov;']) 
        end
        
        % Only Store SST values on the first loop    
        if savesst == 1 && v == 1
            if savecov == 1
                eval(['SST_aut =autall;'])
            end
            
            if savetest == 1
                eval(['SST_test = st_aut;'])  
            end
            
            if saverho == 1
                eval(['SST_rho = rho_aut;'])  
            end     
        end
            
        % End loop for variable
    end
    
    % Save damping variables
    if savedamping == 1
        
        % Make cell array of names
        dampnames = strcat(vars,'_damping');
        
        
        % Save (ex: hfdamping_ENS001_ensorem1_monwin1.mat)
        outname1 = [outpath,'hfdamping_ENS',num2str(ensnum,'%03d'),'_ensorem',num2str(ensorem),'_monwin',num2str(monwin),'.mat'];
        save(outname1,dampnames{:})
        
    end
    
    % Save intermediary covariance/autocovariance variables
    if savecov == 1
        
        % Get variable names
        covnames = strcat(vars,'_cov');
        
        % Save (ex: covars_ENS002_ensorem1.mat)
    	outname2 = [outpath,'covars_ENS',num2str(ensnum,'%03d'),'_ensorem',num2str(ensorem),'.mat'];  
        save(outname2,'LHFLX_cov','SHFLX_cov','FSNS_cov','FLNS_cov','SST_aut')
        
    end
    
    % Save correlation and pvalues
    if savetest == 1
        
        % Get variable names and add significance testing value
        testnames = strcat(vars,'_test');
        testnames{end+1} = 'sigval';
        
        % ex: (test_ENS001_ensorem1_monwin1.mat)
        outname3 = [outpath,'test_ENS',num2str(ensnum,'%03d'),'_ensorem',num2str(ensorem),'_monwin',num2str(monwin),'.mat'];
        save(outname3,testnames{:})
        
       
    end
    
    
    % (4) Save Correlation Coefficients
    if saverho == 1
        
        % Get variable names
        rhonames = strcat(vars,'_rho');
        
        % Save (ex: corrENS_ensorem1_monwin1.mat)
        outname4= [outpath,'corr_ENS',num2str(ensnum,'%03d'),'_ensorem',num2str(ensorem),'_monwin',num2str(monwin),'.mat'];  
        save(outname4,rhonames{:})
        
    end
    
    
    % (5) Save SST Autocorrelation Variables
    if savesst == 1
        
        % Allocate cell array
        sstsave = {};
        
        % SST Covariance
        if savecov == 1
            sstsave{end+1} = 'SST_aut';
        end
        
        % SST Significance Testing Results
        if savetest == 1
            sstsave{end+1} = 'SST_test';
        end
        
        % SST Autocorrelation Coefficients
        if saverho == 1
            sstsave{end+1} = 'SST_rho';
        end
        
        % Save (ex: sstAuto_ENS001_ensorem1_monwin1.mat)
        outname5 = [outpath,'sstAuto_ENS',num2str(ensnum,'%03d'),'_ensorem',num2str(ensorem),'_monwin',num2str(monwin),'.mat'];  
        save(outname5,sstsave{:})
    end
            
            
    
    nend = datetime('now');
    elapsed = nend-nstart;
    fprintf('\nENS%s Completed in %s (%s)',num2str(ensnum,'%03d'),elapsed,nend-allstart)
    
end