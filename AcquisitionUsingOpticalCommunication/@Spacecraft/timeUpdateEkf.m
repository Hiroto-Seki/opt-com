% 探査機のEKFの際，伝搬時間Dtを与えてリファレンスの状態量，誤差共分散行列の更新をする
% RK4の時間間隔はdtとする
function [X, P] = timeUpdateEkf(X, P, constant, Dt, dt,error)
   % 初期化
   timeStep = dt;
   Dt_rest = Dt;
   % system noise
   % system function Error 
   w = [  1e-7;  1e-7 * ones(3,1) ;  error.dynamics * ones(3,1) ];
   Q_t = w * w.';
%    Q_t = [(error.randomClock*1e0)^2, zeros(1,6); zeros(3,7) ; zeros(3,4),eye(3)*((error.dynamics*1e0)^2)]; 
%    Q_t = [zeros(4,7);zeros(3,4),eye(3)*(error.dynamics^2)]; 
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
   STM  = eye(7);
   xvsc = X(2:7);
   k1sc = CelestialBody.twobody(xvsc,constant.sunMu,0);
   k2sc = CelestialBody.twobody(xvsc+0.5*timeStep*k1sc,constant.sunMu,0);
   k3sc = CelestialBody.twobody(xvsc+0.5*timeStep*k2sc,constant.sunMu,0);
   k4sc = CelestialBody.twobody(xvsc+timeStep*k3sc,constant.sunMu,0);
   % STMの伝搬
   k1stm = Spacecraft.delFdelX(xvsc,constant.sunMu)*STM;
   k2stm = Spacecraft.delFdelX(xvsc+0.5*timeStep*k1sc,constant.sunMu)*(STM+0.5*timeStep*k1stm);
   k3stm = Spacecraft.delFdelX(xvsc+0.5*timeStep*k2sc,constant.sunMu)*(STM+0.5*timeStep*k2stm);
   k4stm = Spacecraft.delFdelX(xvsc+timeStep*k3sc,constant.sunMu)*(STM+timeStep*k3stm);
   STM = STM + timeStep/6*(k1stm+2*k2stm+2*k3stm+k4stm);
   % STMを用いてリファレンスの状態量と誤差共分散行列の更新
   X = [X(1) ;xvsc + timeStep/6*(k1sc+2*k2sc+2*k3sc+k4sc)];
   P = STM * P * STM.' +  STM *  Q_t *  STM.' * dt;
   end
   
  
   
end