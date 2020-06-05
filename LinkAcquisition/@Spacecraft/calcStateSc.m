% ある時刻tのscの状態量を計算する
function  xvAtT = calcStateSc(sc,t,time)
% 最終時刻までのリスト内に収まっている時
if t < sc.t(length(sc.t))
    closeTimeIndex = floor( (t - sc.t(1) )/ time.simDt ) + 1;
    closeTimeOffset = t - sc.t(closeTimeIndex);      % 一番近い時刻から伝搬しなければいけない時間
    % spacecraftについて伝搬時刻後(temp)時間後の状態量を得る
    xvsc = sc.state(:,closeTimeIndex);
    k1sc = sc.twobody(xvsc,sc.mu,0);
    k2sc = sc.twobody(xvsc+0.5*closeTimeOffset*k1sc,sc.mu,0);
    k3sc = sc.twobody(xvsc+0.5*closeTimeOffset*k2sc,sc.mu,0);
    k4sc = sc.twobody(xvsc+closeTimeOffset*k3sc,sc.mu,0);
    xvsc = xvsc + closeTimeOffset/6*(k1sc+2*k2sc+2*k3sc+k4sc); 
 % 最終時刻までのリストに収まっていない時
else
    restTime = t -  sc.t(length(sc.t)); 
    % 最終時刻から伝搬していく
    xvsc =  sc.state(:,length(sc.t));
    % time.simDtずつ伝搬する．
    while restTime > 0
        if restTime > time.simDt
            timeStep = time.simDt;
        else
            timeStep = restTime;
        end
        k1sc = sc.twobody(xvsc,sc.mu,0);
        k2sc = sc.twobody(xvsc+0.5*timeStep*k1sc,sc.mu,0);
        k3sc = sc.twobody(xvsc+0.5*timeStep*k2sc,sc.mu,0);
        k4sc = sc.twobody(xvsc+timeStep*k3sc,sc.mu,0);
        xvsc = xvsc + timeStep/6*(k1sc+2*k2sc+2*k3sc+k4sc);
        restTime = restTime - time.simDt;
    end           
end
 xvAtT = xvsc;
    
    
end