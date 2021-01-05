% Dt時間後の宇宙機の位置・速度を求める
% RK4の時間間隔はdtとする
function xvsc = timeUpdate_sc(xvsc,mu, Dt, dt)
   % 初期化
   timeStep = dt;
   Dt_rest = Dt;
   %% RK4で伝搬する
   while Dt * Dt_rest > 1e-10
       if abs(Dt_rest) < timeStep
           timeStep = Dt_rest;
           if Dt_rest > 0
               Dt_rest = Dt_rest - dt;
           else
               Dt_rest = Dt_rest + dt;
           end
       else
           if Dt_rest > 0
               timeStep = dt;
               Dt_rest = Dt_rest - dt;
           else
               timeStep = -dt;
               Dt_rest = Dt_rest + dt;
           end

       end
   k1sc = CelestialBody.twobody(xvsc,mu,0);
   k2sc = CelestialBody.twobody(xvsc+0.5*timeStep*k1sc,mu,0);
   k3sc = CelestialBody.twobody(xvsc+0.5*timeStep*k2sc,mu,0);
   k4sc = CelestialBody.twobody(xvsc+timeStep*k3sc,mu,0);
   xvsc = xvsc + timeStep/6*(k1sc+2*k2sc+2*k3sc+k4sc);
   end
   
   
end