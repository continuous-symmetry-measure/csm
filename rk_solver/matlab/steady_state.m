function res = steady_state(params)
res = (-params.W + sqrt((params.W)^2 + 8*params.K * params.F))/ (4*params.K);
return;