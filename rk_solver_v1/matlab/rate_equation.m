function res = rate_equation(t, n, params)
res = params.F - (params.W * n) - (2 * params.K * (n^2));
return;