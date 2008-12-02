function res = steady_state_moments(params)
A = [(2*params.K - params.W) ,(-2*params.K);...
    (4*params.K+2*params.F+params.W), (-4*params.K-2*params.W)];
res = inv(A)*[-params.F; -params.F];
return;