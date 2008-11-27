function [times, results] = rk(func, a, b, dt, y0, params)

intervals = floor((b-a)/dt) + 1;
results = zeros(intervals, length(y0));
times = a:dt:b;
results(1,:) = y0;
t = a;
for ii = 2:intervals
    prev = results(ii - 1,:);
    k1 = func(t, prev, params);
    k2 = func(t + dt / 2, prev + dt/2*k1, params);
    k3 = func(t + dt / 2, prev + dt/2*k2, params);
    k4 = func(t + dt, prev + dt * k3, params);
    deriv = (k1 + 2 * k2 + 2 * k3 + k4) / 6;
    results(ii, :) = prev + dt*deriv;
    t = t + dt;   
    if (norm(deriv) < 1e-12)                          
        times = times(1:ii);
        results = results(1:ii,:);        
        return;
    end
end
disp('did not converge');
return;