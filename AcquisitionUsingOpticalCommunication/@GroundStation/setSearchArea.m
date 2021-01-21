% 1回の観測(探索)にかかる時間
function [gs,time,directionAccuracy_ut] = setSearchArea(time,gs,SSD,scEstByGsSeqP,error)
    % 推定値の分散から，探索が必要な範囲を決定
    posError = (scEstByGsSeqP(2,2) + scEstByGsSeqP(3,3) + scEstByGsSeqP(4,4))^0.5;
    gs.searchArea     = max([(posError/(SSD.AU*10)) * 3, 10*1e-6, gs.searchArea]) ; %10AUくらいを想定．3sigmaをカバーする. 
    gs.searchStep     = min(0.66e-6, ceil(gs.searchArea/50*1e8)*1e-8 ); %探索時の1stepあたりの間隔(rad). 0.1urad~0.8uradくらいで探索する
    %% 送る情報量から，1回のuplinkにかかる時間を求める
    % 今は仮で情報量を置いておく
    dataVolume = 13 * 64; % bit
    gs.searchTimeStep = gs.switchTime + dataVolume/gs.bps_up;
    % 探索1回にかかる時間
    time.obs = (2 * gs.searchArea/gs.searchStep)^2 * gs.searchTimeStep;    
    % 探索1回にかかる時間がシミュレーションの何stepに相当するか
    time.obsStep = ceil(time.obs/time.simDt);
    % 観測誤差共分散の計算
    directionAccuracy_ut = ( gs.searchStep ^2 + error.gsPoint^2)^0.5; 
    
end