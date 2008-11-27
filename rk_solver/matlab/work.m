function work()

n_master = 10;
[t1, rate] = rk(@rate_equations, 0, 50000, 0.5, 0,params);
[t2, moment] = rk(@moment_equations, 0, 50000, 0.1, [0,0],params);
[t3, master] = rk(@master_equations, 0, 10000, 0.1, ...    
    [1, zeros(1,n_master-1)],params);
master_exp = master * [0:(n_master - 1)]';
master_moment2 = master * ([0:(n_master - 1)].^2)';