% compare true orbit and estimated orbit in sun-earth fix coordinate
% especially compare when covariance is small. 



% plot the last part of orbit
scTrue_RotFrame = RotFrame2D(earth_afterLos.state, scTrue_afterLos.state);
scEst_RotFrame = RotFrame2D(earth_afterLos.state, scEstByScEkf_afterLos.state);



figure('Name','comparison between estimated orbit and true orbit')
plot(scTrue_RotFrame(1,end-1:end), scTrue_RotFrame(2,end-1:end))
hold on 
plot(scEst_RotFrame(1,end-1:end), scEst_RotFrame(2,end-1:end))
legend("true orbit","estimated orbit")




function yRotFrame2D = RotFrame2D(yRef,y)  % yRefは地球の座標. 地球を中心にする
    yRotFrame2D= [(y(1,:).*yRef(1,:) + y(2,:).*yRef(2,:))./(yRef(1,:).^2 + yRef(2,:).^2).^0.5;(y(2,:).*yRef(1,:) - y(1,:).*yRef(2,:))./(yRef(1,:).^2 + yRef(2,:).^2).^0.5];
    yRotFrame2D = yRotFrame2D - yRef(1:2,:);
end
