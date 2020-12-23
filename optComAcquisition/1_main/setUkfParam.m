function ukf= setUkfParam()
    ukf.n = 7; % 推定する状態量の数(クロックオフセット，位置3・速度3)
    ukf.alpha = 1; % 0~1の間
    % 棄却シグマ値レベル．どのくらいが妥当なのか？？→とりあえず全て3シグマにしておく
    ukf.sigmaN = 3;
    
    ukf.kappa  = 3 - ukf.n;
    ukf.lambda = ukf.alpha^2 * (ukf.n + ukf.kappa) - ukf.n;
    ukf.beta   = 2;
    ukf.w0_m   = ukf.lambda/(ukf.n + ukf.lambda);
    ukf.wi_m   = 1/(2*(ukf.n + ukf.lambda));
    ukf.w0_c   = ukf.lambda/(ukf.n + ukf.lambda) + (1- ukf.alpha^2 + ukf.beta);
    ukf.wi_c   = 1/(2*(ukf.n + ukf.lambda));  
end