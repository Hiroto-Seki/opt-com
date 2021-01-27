% result.mを読み込んで，グラフを作る
% グラフを保存する場所
resultFile = '2021012614'; %書き換える
resultPath = ['~/Documents/lab/master/research/result/',resultFile];
addpath(resultPath);

% Downlinkの成功率の分布グラフ
f1 = figure('visible', 'on');
f1.OuterPosition = [100 100 800 400];
% グラフの範囲を決める
minAvailability = min(reshape([result(1).downAvail(:);result(3).downAvail(:)],1,[] ));
minXLabel_f1 = floor(minAvailability/10)*10;  
subplot(1,2,1)
histogram(result(1).downAvail,(100-minXLabel_f1)*5)
xlim([minXLabel_f1 100])
xlabel("downlink availability [%]")
ylabel("number of result in 100 simulations")
title({'downlink availability during link acquisition'; 'before long communication disconnection'})
subplot(1,2,2)
histogram(result(3).downAvail,(100-minXLabel_f1)*5)
xlim([minXLabel_f1 100])
xlabel("downlink availability [%]")
ylabel("number of result in 100 simulations")
title({'downlink availability during link acquisition'; 'after long communication disconnection'})
file1 = [resultPath,'/downlinkAvailability.png'];
saveas(f1, file1)


% 最終時刻の軌道決定精度の分布のグラフ
f2 = figure('visible', 'on');
f2.OuterPosition = [100 100 800 400];
minPosError = min([result(1).posErrorSc(4,:),result(3).posErrorSc(4,:)]);
maxPosError = max([result(1).posErrorSc(4,:),result(3).posErrorSc(4,:)]);
minXLabel_f2 = floor(log10(minPosError));
maxXLabel_f2 = ceil(log10(maxPosError));
edge_f2      = logspace(minXLabel_f2,maxXLabel_f2,(maxXLabel_f2-minXLabel_f2)*5 + 1);
subplot(1,2,1)
histogram(result(1).posErrorSc(4,:), edge_f2)
set(gca,'xscale','log')
xlim([10^minXLabel_f2 10^maxXLabel_f2])
xlabel("position error [km]")
ylabel("number of result in 100 simulations")
title({'position error after orbit determination'; 'before long communication disconnection'})
subplot(1,2,2)
histogram(result(3).posErrorSc(4,:), edge_f2)
set(gca,'xscale','log')
xlim([10^minXLabel_f2 10^maxXLabel_f2])
xlabel("position error [km]")
ylabel("number of result in 100 simulations")
title({'position error after orbit determination'; 'after long communication disconnection'})
file2 = [resultPath,'/positionErrorSammary.png'];
saveas(f2, file2)


f3 = figure('visible', 'on');
f3.OuterPosition = [100 100 800 400];
minVelError = min([result(1).velErrorSc(4,:),result(3).velErrorSc(4,:)]);
maxVelError = max([result(1).velErrorSc(4,:),result(3).velErrorSc(4,:)]);
minXLabel_f3 = floor(log10(minVelError));
maxXLabel_f3 = ceil(log10(maxVelError));
edge_f3      = logspace(minXLabel_f3,maxXLabel_f3,(maxXLabel_f3-minXLabel_f3)*5 + 1);
subplot(1,2,1)
histogram(result(1).velErrorSc(4,:), edge_f3)
set(gca,'xscale','log')
xlim([10^minXLabel_f3 10^maxXLabel_f3])
xlabel("velocity error [km]")
ylabel("number of result in 100 simulations")
title({'velocity error after orbit determination'; 'before long communication disconnection'})
subplot(1,2,2)
histogram(result(3).velErrorSc(4,:), edge_f3)
set(gca,'xscale','log')
xlim([10^minXLabel_f3 10^maxXLabel_f3])
xlabel("velocity error [km]")
ylabel("number of result in 100 simulations")
title({'velocity error after orbit determination'; 'after long communication disconnection'})
file3 = [resultPath,'/velocityErrorSammary.png'];
saveas(f3, file3)


% クロックのオフセット
f4 = figure('visible', 'on');
f4.OuterPosition = [100 100 800 400];
minClockError = min([result(1).clockErrorSc(1,:),result(3).clockErrorSc(1,:)]);
maxClockError = max([result(1).clockErrorSc(1,:),result(3).clockErrorSc(1,:)]);
minXLabel_f4 = floor(log10(minClockError));
maxXLabel_f4 = ceil(log10(maxClockError));
edge_f4      = logspace(minXLabel_f4,maxXLabel_f4,(maxXLabel_f4-minXLabel_f4)*5 + 1);
subplot(1,2,1)
histogram(result(1).clockErrorSc(1,:), edge_f4)
set(gca,'xscale','log')
xlim([10^minXLabel_f4 10^maxXLabel_f4])
xlabel("clock bias error [km]")
ylabel("number of result in 100 simulations")
title({'clock bias error after orbit determination'; 'before long communication disconnection'})
subplot(1,2,2)
histogram(result(3).clockErrorSc(1,:), edge_f4)
set(gca,'xscale','log')
xlim([10^minXLabel_f4 10^maxXLabel_f4])
xlabel("clock bias error")
ylabel("number of result in 100 simulations")
title({'clock bias error after orbit determination'; 'after long communication disconnection'})
file4 = [resultPath,'/clcokBiasErrorSammary.png'];
saveas(f4, file4)

% 誤差共分散と誤差の比較
% clockのoffset
f5 = figure('visible', 'on');
f5.OuterPosition = [100 100 800 400];
errorBySigma_clockOffest1 = result(1).clockErrorSc(1,:)./result(1).clockErrorSc(2,:);
errorBySigma_clockOffest3 = result(3).clockErrorSc(1,:)./result(3).clockErrorSc(2,:);
minErrorBySigma_ClockOffset = min([errorBySigma_clockOffest1,errorBySigma_clockOffest3]);
maxErrorBySigma_ClockOffset = max([errorBySigma_clockOffest1,errorBySigma_clockOffest3]);
minXLabel_f5 = floor(log10(minErrorBySigma_ClockOffset));
maxXLabel_f5 = ceil(log10(maxErrorBySigma_ClockOffset));
edge_f5      = logspace(minXLabel_f5,maxXLabel_f5,(maxXLabel_f5-minXLabel_f5)*5 + 1);
subplot(1,2,1)
histogram(errorBySigma_clockOffest1, edge_f5)
set(gca,'xscale','log')
xlim([10^minXLabel_f5 10^maxXLabel_f5])
xlabel("clock bias error/(clock bias error covariance)^{0.5}")
ylabel("number of result in 100 simulations")
title({'clock bias error to sigma ratio after orbit determination'; 'before long communication disconnection'})
subplot(1,2,2)
histogram(errorBySigma_clockOffest3, edge_f5)
set(gca,'xscale','log')
xlim([10^minXLabel_f5 10^maxXLabel_f5])
xlabel("clock bias error/(clock bias error covariance)^{0.5}")
ylabel("number of result in 100 simulations")
title({'clock bias error to sigma ratio after orbit determination'; 'after long communication disconnection'})
file5 = [resultPath,'/clcokBiasError2SigmaRatioSammary.png'];
saveas(f5, file5)

f6 = figure('visible', 'on');
f6.OuterPosition = [100 100 800 400];
errorBySigma_position1 = result(1).posErrorSc(4,:)./result(1).posErrorSc(2,:);
errorBySigma_position3 = result(3).posErrorSc(4,:)./result(3).posErrorSc(2,:);
minErrorBySigma_position = min([errorBySigma_position1,errorBySigma_position3]);
maxErrorBySigma_position = max([errorBySigma_position1,errorBySigma_position3]);
minXLabel_f6 = floor(log10(minErrorBySigma_position));
maxXLabel_f6 = ceil(log10(maxErrorBySigma_position));
edge_f6      = logspace(minXLabel_f6,maxXLabel_f6,(maxXLabel_f6-minXLabel_f6)*5 + 1);
subplot(1,2,1)
histogram(errorBySigma_position1, edge_f6)
set(gca,'xscale','log')
xlim([10^minXLabel_f6 10^maxXLabel_f6])
xlabel("position error/(position error covariance)^{0.5}")
ylabel("number of result in 100 simulations")
title({'position error to sigma ratio after orbit determination'; 'before long communication disconnection'})
subplot(1,2,2)
histogram(errorBySigma_position3, edge_f6)
set(gca,'xscale','log')
xlim([10^minXLabel_f6 10^maxXLabel_f6])
xlabel("position error/(position error covariance)^{0.5}")
ylabel("number of result in 100 simulations")
title({'position error to sigma ratio after orbit determination'; 'after long communication disconnection'})
file6 = [resultPath,'/positionError2SigmaRatioSammary.png'];
saveas(f6, file6)


f7 = figure('visible', 'on');
f7.OuterPosition = [100 100 800 400];
errorBySigma_velocity1 = result(1).velErrorSc(4,:)./result(1).velErrorSc(2,:);
errorBySigma_velocity3 = result(3).velErrorSc(4,:)./result(3).velErrorSc(2,:);
minErrorBySigma_velocity = min([errorBySigma_velocity1,errorBySigma_velocity3]);
maxErrorBySigma_velocity = max([errorBySigma_velocity1,errorBySigma_velocity3]);
minXLabel_f7 = floor(log10(minErrorBySigma_velocity));
maxXLabel_f7 = ceil(log10(maxErrorBySigma_velocity));
edge_f7      = logspace(minXLabel_f7,maxXLabel_f7,(maxXLabel_f7-minXLabel_f7)*5 + 1);
subplot(1,2,2)
histogram(errorBySigma_velocity1, edge_f7)
set(gca,'xscale','log')
xlim([10^minXLabel_f7 10^maxXLabel_f7])
xlabel("velocity error/(velocity error covariance)^{0.5}")
ylabel("number of result in 100 simulations")
title({'velocity error to sigma ratio after orbit determination'; 'before long communication disconnection'})
subplot(1,2,1)
histogram(errorBySigma_velocity3, edge_f7)
set(gca,'xscale','log')
xlim([10^minXLabel_f7 10^maxXLabel_f7])
xlabel("velocity error/(velocity error covariance)^{0.5}")
ylabel("number of result in 100 simulations")
title({'velocity error to sigma ratio after orbit determination'; 'after long communication disconnection'})
file7 = [resultPath,'/velocityError2SigmaRatioSammary.png'];
saveas(f7, file7)


