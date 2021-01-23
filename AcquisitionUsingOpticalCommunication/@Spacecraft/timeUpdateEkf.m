% 探査機のEKFの際，伝搬時間Dtを与えてリファレンスの状態量，誤差共分散行列の更新をする
% RK4の時間間隔はdtとする
function [X, P] = timeUpdateEkf(X, P, constant, Dt, dt,error)
   % 初期化
   timeStep = dt;
   Dt_rest = Dt;
   % system noise
   % system function Error 
%    w = [  0;  [0;0;0] ;  error.dynamics * ones(3,1) ];
%    Q_t = w * w.';
   Q = zeros(7,7);
   Q(5:7,5:7) = error.dynamics^2 * eye(3);
   Q(1,1)     = 1e-18; % 1e-16
   Q(2:4,2:4) = 1e-18 * eye(3); % 1e-16
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
%    delFdelX = Spacecraft.delFdelX(X(2:7),constant.sunMu);
%    STM     = expm(delFdelX * dt);
   
   
   % STMを用いてリファレンスの状態量と誤差共分散行列の更新
   X = [X(1) ;xvsc + timeStep/6*(k1sc+2*k2sc+2*k3sc+k4sc)];
   P = STM * P * STM.' +  STM *  Q *  STM.' * dt;
   end
   

end