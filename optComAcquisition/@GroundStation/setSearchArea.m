% 1回の観測(探索)にかかる時間
function [time,R1wSc,R2wSc] = setSearchArea(time,gs,SSD,scEstByGsSeqP,R1wSc,R2wSc,error)
    % 推定値の分散から，探索が必要な範囲を決定
    posError = (scEstByGsSeqP(2,2) + scEstByGsSeqP(3,3) + scEstByGsSeqP(4,4))^0.5;
    gs.searchArea     = (posError/(SSD.AU*10)) * 3 ; %10AUくらいを想定．3sigmaをカバーする
    gs.searchStep     = min(2e-6, ceil(gs.searchArea/20*1e8)*1e-8 ); %探索時の1stepあたりの間隔(rad)
    % 探索1回にかかる時間
    time.obs = (2 * gs.searchArea/gs.searchStep)^2 * gs.searchTimeStep;    
    % 探索1回にかかる時間がシミュレーションの何stepに相当するか
    time.obsStep = ceil(time.obs/time.simDt);
    % Rの書き換え
    R1wSc(6:7,6:7) = (gs.searchStep^2+error.gsPoint^2)*eye(2);
    R2wSc(6:7,6:7) = (gs.searchStep^2+error.gsPoint^2)*eye(2);
    
end