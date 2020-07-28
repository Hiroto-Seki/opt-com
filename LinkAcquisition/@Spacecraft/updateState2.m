% 地上局のEKFの際，伝搬時間を与えてリファレンスの状態量，誤差共分散行列の更新をする
function [X_bar, P_bar] = updateState2(X_hat,P, dt,mu)
   STM = eye(6);
   xvsc = X_hat;
   k1sc = Spacecraft.twobody(xvsc,mu,0);
   k2sc = Spacecraft.twobody(xvsc+0.5*dt*k1sc,mu,0);
   k3sc = Spacecraft.twobody(xvsc+0.5*dt*k2sc,mu,0);
   k4sc = Spacecraft.twobody(xvsc+dt*k3sc,mu,0);
   % STMの伝搬
   k1stm = Spacecraft.delFdelX2(xvsc,mu)*STM;
   k2stm = Spacecraft.delFdelX2(xvsc+0.5*dt*k1sc,mu)*(STM+0.5*dt*k1stm);
   k3stm = Spacecraft.delFdelX2(xvsc+0.5*dt*k2sc,mu)*(STM+0.5*dt*k2stm);
   k4stm = Spacecraft.delFdelX2(xvsc+dt*k3sc,mu)*(STM+dt*k3stm);
   STM = STM + dt/6*(k1stm+2*k2stm+2*k3stm+k4stm);
   % STMを用いてリファレンスの状態量と誤差共分散行列の更新
   X_bar = STM * X_hat;
   P_bar = STM * P * STM.';
end