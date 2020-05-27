%% ---------------------
%  hf4_sumfluxes
%  ---------------------
% Sum fluxes together
% Restructured to output file for different members

%% User Inputs
vars   = {'LHFLX','SHFLX','FLNS','FSNS','RHFLX','THFLX','NHFLX'};
mnum   = [1:35,101:107] ;% Ensemble Members List
lags   = [1:3]          ;% Number of lags to include
deg5   = 1              ;% Set to 1 to use smoothed data


monwin = 3              ;% Month Window
% Note: Net outputs heat flux damping, test results, and correlation coefficients
% automatically, regardless of the toggle options below!

% Paths to data and output
if deg5 == 1
    datpath = ['/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/hfdamping_matfiles/03_hf3out/5deg/',num2str(monwin),'mon/'];
    outpath = ['/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/hfdamping_matfiles/04_hf4out/5deg/',num2str(monwin),'mon/'];
    lonsize = 72;
    latsize = 36;
else
    datpath = ['/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/hfdamping_matfiles/03_hf3out/',num2str(monwin),'mon/'];
    outpath = ['/stormtrack/data3/glliu/01_Data/02_AMV_Project/01_hfdamping/hfdamping_matfiles/04_hf4out/',num2str(monwin),'mon/'];
    lonsize = 288;
    latsize = 192;
end

% Enso Removal?
ensorem     = 1; % Set to 1 if ENSO was removed

% Saving Options
savedamping = 0; % Set to 1 to save damping
savetest    = 0; % Set to 1 to save significance testing
saverho     = 1; % Set to 1 to save correlation coefficients
savesst     = 1; % Set to 1 to save SST


%% Script Start
allstart = datetime('now');
fprintf('Now running hf4_sumfluxes for monwin %i (%s)',monwin,allstart)
fprintf('\n Save Options -- Damping:%i | SigTest:%i | Corr:%i | SST:%i | deg5:%i',savedamping,savetest,saverho,savesst,deg5)
fprintf('\n\t Save Loc: %s',outpath)

% -----------------------------------------------
% Preallocate (LON x LAT x ENS x MON x LAG X VAR)
% -----------------------------------------------
numvars = length(vars);
pasize = [lonsize,latsize,42,12,length(lags),numvars];

if savedamping == 1
    dampall = NaN(pasize);
end

if savetest == 1
    testall = NaN(pasize);
    %sigvall = NaN(42,numvars);
    
    if savesst == 1
        tsst = NaN(lonsize,latsize,42,12,length(lags));
    end
end

if saverho == 1
    corrall = NaN(pasize);
    
    
    if savesst == 1
        rsst = NaN(lonsize,latsize,42,12,length(lags));
    end
end

% Loop by ensemble member to combine heat fluxes
for n = 1:length(mnum)
    ensnum = mnum(n);
    
     % Load and sum damping variables
    if savedamping == 1
        
        % Construct matfile name (ex: hfdamping_ENS001_ensorem1_monwin1.mat)
        matname = [datpath,'hfdamping_ENS',num2str(ensnum,'%03d'),'_ensorem',num2str(ensorem),'_monwin',num2str(monwin),'.mat'];
         
        % Load in variables
        load(matname)
        dampnames = strcat(vars,'_damping');
        
        for v = 1:length(dampnames)
            eval(['dampall(:,:,n,:,:,v)=',dampnames{v},'(:,:,:,1:3);']);
        end
    end
    
    % Load and sum significance test results
    if savetest == 1
        
        % Construct matfile name (ex: test_ENS001_ensorem1_monwin1.mat)
        matname = [datpath,'test_ENS',num2str(ensnum,'%03d'),'_ensorem',num2str(ensorem),'_monwin',num2str(monwin),'.mat'];
        
        % Load variables
        load(matname)
        testnames = strcat(vars,'_test');
        
        for v = 1:length(testnames)
            eval(['testall(:,:,n,:,:,v)=',testnames{v},'(:,:,:,1:3);']);
        end
        
    end
    
    
    % Load and sum correlation maps
    if saverho == 1
        
        % Construct matfile name (ex: corr_ENS001_ensorem1_monwin1.mat)
        matname = [datpath,'corr_ENS',num2str(ensnum,'%03d'),'_ensorem',num2str(ensorem),'_monwin',num2str(monwin),'.mat'];
        
        % Load variables
        load(matname)
        corrnames = strcat(vars,'_rho');
        
        for v = 1:length(corrnames)
            eval(['corrall(:,:,n,:,:,v)=',corrnames{v},'(:,:,:,1:3);']);
        end
        
    end
    
    
    % SST Autocorrelation Options
    if savesst == 1
        
        matname = [datpath,'sstAuto_ENS',num2str(ensnum,'%03d'),'_ensorem',num2str(ensorem),'_monwin',num2str(monwin),'.mat'];  
        load(matname)
        if savetest == 1
            tsst(:,:,n,:,:) = SST_test;
        end
        
        if saverho == 1
            rsst(:,:,n,:,:) = SST_rho;
        end
        
    end
    
    fprintf('\n ENS%s loaded! (%s)',num2str(ensnum,'%03d'),datetime('now')-allstart)
end


% Assign variables to their separate summed files
for v = 1:length(vars)
    
    
    % Get lowercase form of the name
    vname = char(vars(v));
    
    
    % Save Damping Variables Separately
    if savedamping ==1     
        % Get corresponding variable
        damping = dampall(:,:,:,:,:,v);
        
        % Output ex: (LHFLX_damping_ensorem1_monwin1.mat)
        outname = [outpath,vname,'_damping_ensorem',num2str(ensorem),'_monwin',num2str(monwin),'.mat'];
        save(outname,'vname','damping');     
        
        fprintf('\n Saved damping for %s! (%s)',vname,datetime('now')-allstart)
    end
    
    % Save Test Variables Separately
    if savetest == 1     
        % Get corresponding variable (crosscorrelation test)
        test = testall(:,:,:,:,:,v);
        
        
        % Output ex: (LHFLX_test_ensorem1_monwin1.mat)
        outname = [outpath,vname,'_test_ensorem',num2str(ensorem),'_monwin',num2str(monwin),'.mat'];
        save(outname,'vname','test');
        fprintf('\n Saved significance test results for %s! (%s)',vname,datetime('now')-allstart)
        
        % Save SST testing if option is set
        if savesst == 1 && v == 1
            outname = [outpath,'SST_test_ensorem',num2str(ensorem),'_monwin',num2str(monwin),'.mat']
            save(outname,'tsst');
            
            fprintf('\n Saved significance test results for SST! (%s)',datetime('now')-allstart)
        end
        
        
    end
        
        

    % Save Cross Correlation Coefficients
    if saverho == 1      
        % Get corresponding variable (crosscorrelation coefficients)
        rho = corrall(:,:,:,:,:,v);
        
        % Output ex: (LHFLX_rho_ensorem1_monwin1.mat)
        outname = [outpath,vname,'_rho_ensorem',num2str(ensorem),'_monwin',num2str(monwin),'.mat'];
        save(outname,'vname','rho');
        
        fprintf('\n Saved correlation coefficients for %s! (%s)',vname,datetime('now')-allstart)
        
        % Save SST Autocorrelation Coefficients if option is set
        if savesst == 1 && v == 1
            
            outname = [outpath,'SST_rho_ensorem',num2str(ensorem),'_monwin',num2str(monwin),'.mat'];
            save(outname,'rsst');
            fprintf('\n Saved correlation coefficients for SST! (%s)',datetime('now')-allstart)
        end  
    end
 
end
