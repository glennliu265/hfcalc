function [rho,T,sigtest,critval,n_eff,corrthres] = pcorr(A,B,dim,p,tails,dof_type,dof_man)

% Calculates correlation coefficient between two matrices A and B, of equal
% size and also returns the significance test result (0 for fail, 1 for
% pass)
%   
% dof_type choices (1,2,3)
%     1 - Assume n-2 degrees of freedom
%     2 - N* = N* (1-rho^2)/(1+rho^2) [Recommended by Yu-Chiao]
%     3 - N* = N* (1-rho)/(1+rho)     [From PSU website]
%     4 - Input (manually a dof)
% dim indicates the dimension along which the correlation value will be
% computed
%
% tail indicates whether it is one or two sided
% Based on von Storch and Zwiers, Statistical Analysis in Climate Research 
% (2003)

%% Script Start

% Student's T Distribution for 2 Degrees of freedom
% Upper Tail Critical Values
% (from Appendix F p.423)
%pct  = [0.750 0.900 0.950 0.975 0.990 0.995  0.999];
%df2  = [0.816 1.886 2.920 4.303 6.965 9.925 22.327];

% Find Anomaly
Aanom = A - nanmean(A,dim);
Banom = B - nanmean(B,dim);

% Find elementwise product of A and B
AB    = Aanom .* Banom;

% Find the squared variance
A2    = Aanom.^2;
B2    = Banom.^2;

% Compute Pearson's Correlation Coefficient
rho   = nansum(AB,dim) ./ sqrt(nansum(A2,dim) .* nansum(B2,dim));

% Compute effective degrees of freedom
n     = size(A,dim);
switch dof_type
    case 1
        n_eff = n-2;
    case 2
        n_eff = n .* (1-rho.^2)./(1+rho.^2);
    case 3
        n_eff = n .* (1-rho)./(1+rho);
    case 4        
        % Manually input dof
        n_eff = dof_man;
end

% Compute p-value based on tails
ptilde  = (1-p/tails);

if dof_type == 1 || dof_type == 4
    % Compute T at each point
    T     = rho .* sqrt(n_eff ./ (1-rho.^2));
    
    
    % Compute threshold critical value
    critval = tinv(ptilde,n_eff); % Note that this does not consider degrees of freedom

    % H0: Correlation is not significantly different from 0
    % Ha: Correlation IS significantly different from 0 (keep)
    % So if the critical value exceeds critval, then we take 
    sigtest = abs(T) > critval;
    
    % Compute correlation threshold (only holds for single case atm)
    corrthres = sqrt(1/((n_eff/critval^2)+1));
       
else
    % Compute T @ each point
    T       = rho .* sqrt(n_eff ./ (1-rho.^2));
    
    % Compute threshold critical value @ each point
    critval = tinv(ptilde,n_eff); % Note that this does not consider degrees of freedom
    
    % Perform significance testing
    sigtest = T > critval;
    
    % Compute correlation threshold (only holds for single case atm)
    corrthres = sqrt(1/((n_eff./critval.^2)+1));
    
end



