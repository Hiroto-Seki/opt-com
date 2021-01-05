% シグマポイントの計算
% 出力: obj.x_sc M x N (M:状態変数の数，N:シグマ点の数)
function x_sp = calcSigmaPoint(X,P,ukf)
n = ukf.n;
x_sp0 = X;

% 誤差共分散行列の平方根行列
sqrtP = sqrtm(P);

% 虚数部分は無視する．(これがどこまでやっていいのか？)
sqrtP = real(sqrtP);


x_spi = zeros(n, 2*n);
for i = 1:n
    x_spi(:,i) = X + (n+ukf.lambda)^0.5 * sqrtP(:,i);
    x_spi(:,i+n) = X - (n+ukf.lambda)^0.5 * sqrtP(:,i);
end


x_sp = [x_sp0, x_spi];


end
