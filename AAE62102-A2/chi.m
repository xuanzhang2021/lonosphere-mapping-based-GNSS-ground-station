function [T_threshold,PL] = chi()

    P_fa = 1e-2;       
    sigma = 3;         

    n_sat = 5;         

    dof = n_sat - 4;  

    T_threshold = chi2inv(1 - P_fa, dof) * sigma^2; % ≈ 82.89
    fprintf('T_threshold = %.2f meters\n', T_threshold);

    P_md = 1e-7;
    
    fun = @(PL) ncx2cdf(T_threshold, dof, PL.^2 / sigma^2) - (1 - P_md);
    PL_guess = 50; 
    PL = fzero(fun, PL_guess); 
    fprintf('3D Protection Level = %.2f meters\n', PL);
end