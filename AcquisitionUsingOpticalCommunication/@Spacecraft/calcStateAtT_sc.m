% ある時刻tのscの状態量を計算する
function  xvAtT = calcStateAtT_sc(obj,t,time,constant)
mu = constant.sunMu;
% 最終時刻までのリスト内に収まっている時
if t < obj.t(length(obj.t))
    closeTimeIndex = floor( (t - obj.t(1) )/ time.simDt ) + 1;
    closeTimeOffset = t - obj.t(closeTimeIndex);      % 一番近い時刻から伝搬しなければいけない時間
    % spacecraftについて伝搬時刻後(temp)時間後の状態量を得る
    xvsc = obj.state(:,closeTimeIndex);
    k1sc = CelestialBody.twobody(xvsc,mu,0);
    k2sc = CelestialBody.twobody(xvsc+0.5*closeTimeOffset*k1sc,mu,0);
    k3sc = CelestialBody.twobody(xvsc+0.5*closeTimeOffset*k2sc,mu,0);
    k4sc = CelestialBody.twobody(xvsc+closeTimeOffset*k3sc,mu,0);
    xvsc = xvsc + closeTimeOffset/6*(k1sc+2*k2sc+2*k3sc+k4sc); 
 % 最終時刻までのリストに収まっていない時
else
    restTime = t -  obj.t(length(obj.t)); 
    % 最終時刻から伝搬していく
    xvsc =  obj.state(:,length(obj.t));
    % time.simDtずつ伝搬する．
    while restTime > 0
        if restTime > time.simDt
            timeStep = time.simDt;
        else
            timeStep = restTime;
        end
        k1sc = CelestialBody.twobody(xvsc,mu,0);
        k2sc = CelestialBody.twobody(xvsc+0.5*timeStep*k1sc,mu,0);
        k3sc = CelestialBody.twobody(xvsc+0.5*timeStep*k2sc,mu,0);
        k4sc = CelestialBody.twobody(xvsc+timeStep*k3sc,mu,0);
        xvsc = xvsc + timeStep/6*(k1sc+2*k2sc+2*k3sc+k4sc);
        restTime = restTime - time.simDt;
    end           
end
 xvAtT = xvsc;
     
end