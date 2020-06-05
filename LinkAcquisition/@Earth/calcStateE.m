% ある時刻の状態量を計算する
function  xvAtT = calcStateE(earth,t,time)
% 最終時刻までのリスト内に収まっている時
if t < earth.t(length(earth.t))
    closeTimeIndex = floor( (t - earth.t(1) )/ time.simDt ) + 1;
    closeTimeOffset = t - earth.t(closeTimeIndex);      % 一番近い時刻から伝搬しなければいけない時間
    % spacecraftについて伝搬時刻後(temp)時間後の状態量を得る
    xvsc = earth.state(:,closeTimeIndex);
    k1sc = earth.twobody(xvsc,earth.mu,0);
    k2sc = earth.twobody(xvsc+0.5*closeTimeOffset*k1sc,earth.mu,0);
    k3sc = earth.twobody(xvsc+0.5*closeTimeOffset*k2sc,earth.mu,0);
    k4sc = earth.twobody(xvsc+closeTimeOffset*k3sc,earth.mu,0);
    xvsc = xvsc + closeTimeOffset/6*(k1sc+2*k2sc+2*k3sc+k4sc); 
 % 最終時刻までのリストに収まっていない時
else
    restTime = t -  earth.t(length(earth.t)); 
    % 最終時刻から伝搬していく
    xvsc =  earth.state(:,length(earth.t));
    % time.simDtずつ伝搬する．
    while restTime > 0
        if restTime > time.simDt
            timeStep = time.simDt;
        else
            timeStep = restTime;
        end
        k1sc = earth.twobody(xvsc,earth.mu,0);
        k2sc = earth.twobody(xvsc+0.5*timeStep*k1sc,earth.mu,0);
        k3sc = earth.twobody(xvsc+0.5*timeStep*k2sc,earth.mu,0);
        k4sc = earth.twobody(xvsc+timeStep*k3sc,earth.mu,0);
        xvsc = xvsc + timeStep/6*(k1sc+2*k2sc+2*k3sc+k4sc);
        restTime = restTime - time.simDt;
    end           
end
 xvAtT = xvsc;
    
    
end