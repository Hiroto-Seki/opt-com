% STMを求める

function STM = getSTM(x0,mu,t0, tf)
    tspan = [t0 tf];
    y0     = [x0;reshape(eye(6),36,1)];
    opts   = odeset('Reltol',1e-12,'AbsTol',1e-12);
    odefun = @(t,x) Spacecraft.twobody_stateAndSTM(t,x,mu);
    % nominal trajectory
    [t,ynom] = ode113(odefun,tspan,y0,opts);
    STM = reshape(ynom(size(ynom,1),7:end),6,6);

end