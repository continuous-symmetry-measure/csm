function res = moment_equations(t, moments, params)
res = [];
res(1) = params.F + (2*params.K - params.W) * moments(1) - ...
    2*params.K * moments(2);
res(2) = params.F + (4*params.K + 2*params.F + params.W) * moments(1) + ...
    (-4*params.K - 2*params.W) * moments(2);
return