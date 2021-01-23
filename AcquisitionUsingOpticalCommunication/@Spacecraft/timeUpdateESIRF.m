% 探査機のEKFの際，伝搬時間Dtを与えてリファレンスの状態量，誤差共分散行列の更新をする
% RK4の時間間隔はdtとする
function [X, P] = timeUpdateESIRF(X, P, constant, Dt, dt,error)
   % 初期化
   timeStep = dt;
   Dt_rest = Dt;
   % system noise
   % system function Error % コレスキー分解するために，clockと位置にも微笑の値を入れている
   Q = zeros(7,7);
   Q(5:7,5:7) = error.dynamics^2 * eye(3);
   % 分解できるように，逆行列を計算できるように数値計算誤差程度の値を加える. 若干適当．．
   Q(1,1) = 1e-10;
   Q(2:4,2:4) = 1e-10 * eye(3);
   % 特別な形なので分解できる
   W          = chol(Q);
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
       
       %% Pの伝搬
       delFdelX = Spacecraft.delFdelX(X(2:7),constant.sunMu);
       STM      = expm(delFdelX * dt);
       odefun = @(tau,gamma) Spacecraft.calcDGamma(tau,delFdelX,dt);
       tspan  = [0 dt];
       gamma0 = zeros(7^2,1);
       [~, gamma]    = ode113(odefun,tspan,gamma0);
       gamma = gamma(end,:);
       gamma = reshape(gamma,[7,7]);
       S = chol(P);
       matrixBeforePropagation = [inv(W), zeros(7,7), zeros(7,1);
                           -inv(S) * inv(STM)*gamma, inv(S)*inv(STM),zeros(7,1)];
       [~,matrixAfterPropagation] = qr(matrixBeforePropagation);
       S = inv(matrixAfterPropagation(8:14,8:14));
       P = S.' * S;
       % Xの更新
       xvsc = X(2:7);
       k1sc = CelestialBody.twobody(xvsc,constant.sunMu,0);
       k2sc = CelestialBody.twobody(xvsc+0.5*timeStep*k1sc,constant.sunMu,0);
       k3sc = CelestialBody.twobody(xvsc+0.5*timeStep*k2sc,constant.sunMu,0);
       k4sc = CelestialBody.twobody(xvsc+timeStep*k3sc,constant.sunMu,0);
       X = [X(1) ;xvsc + timeStep/6*(k1sc+2*k2sc+2*k3sc+k4sc)];
   end
   

end


% function dGamma = calcDGamma(tau,delFdelX,dt)
%     dGamma = expm(delFdelX * (dt-tau));
%     
% end