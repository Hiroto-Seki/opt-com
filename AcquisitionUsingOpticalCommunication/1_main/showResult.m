function showResult(scTrue,scEstByScEkf,scEstByGsEkf,error)

% 各要素ごとに出力
figure(1)
title("clock error")
semilogy((scEstByScEkf.t-scEstByScEkf.t(1))/60/60, abs(scEstByScEkf.clockError),...
         (scEstByGsEkf.t-scEstByScEkf.t(1))/60/60, abs(scEstByGsEkf.clockError),...
         (scEstByScEkf.t-scEstByScEkf.t(1))/60/60, 1 * reshape(scEstByScEkf.P_list(1,1,:), 1, length(scEstByScEkf.t)).^0.5, ...
         (scEstByGsEkf.t-scEstByScEkf.t(1))/60/60, 1 * reshape(scEstByGsEkf.P_list(1,1,:), 1, length(scEstByScEkf.t)).^0.5)
%          (scEstByScEkf.t-scEstByScEkf.t(1))/60/60, 3 * reshape(scEstByScEkf.P_list(1,1,:), 1, length(scEstByScEkf.t)).^0.5,...
%          (scEstByGsEkf.t-scEstByScEkf.t(1))/60/60, 3 * reshape(scEstByGsEkf.P_list(1,1,:), 1, length(scEstByScEkf.t)).^0.5)
xlabel('time[h]')
legend('estimated by sc', 'estimated by gs', '1\sigma estimated by sc', '1\sigma estimated by gs')

figure(2)
tiledlayout(1,4)
nexttile
hold on
title("position X error")
plot((scEstByScEkf.t-scEstByScEkf.t(1))/60/60, (scEstByScEkf.state(1,:) - scTrue.state(1,:)),'-r')
plot((scEstByGsEkf.t-scEstByScEkf.t(1))/60/60, (scEstByGsEkf.state(1,:) - scTrue.state(1,:)),'-g')
plot((scEstByScEkf.t-scEstByScEkf.t(1))/60/60,  1 * reshape(scEstByScEkf.P_list(2,2,:), 1, length(scEstByScEkf.t)).^0.5, '--r')
plot((scEstByScEkf.t-scEstByScEkf.t(1))/60/60, -1 * reshape(scEstByScEkf.P_list(2,2,:), 1, length(scEstByScEkf.t)).^0.5, '--r')
plot((scEstByScEkf.t-scEstByScEkf.t(1))/60/60,  3 * reshape(scEstByScEkf.P_list(2,2,:), 1, length(scEstByScEkf.t)).^0.5, ':r')
plot((scEstByScEkf.t-scEstByScEkf.t(1))/60/60, -3 * reshape(scEstByScEkf.P_list(2,2,:), 1, length(scEstByScEkf.t)).^0.5, ':r')
plot((scEstByGsEkf.t-scEstByScEkf.t(1))/60/60,  1 * reshape(scEstByGsEkf.P_list(2,2,:), 1, length(scEstByScEkf.t)).^0.5, '--g')
plot((scEstByGsEkf.t-scEstByScEkf.t(1))/60/60, -1 * reshape(scEstByGsEkf.P_list(2,2,:), 1, length(scEstByScEkf.t)).^0.5, '--g')
plot((scEstByGsEkf.t-scEstByScEkf.t(1))/60/60,  3 * reshape(scEstByGsEkf.P_list(2,2,:), 1, length(scEstByScEkf.t)).^0.5, ':g')
plot((scEstByGsEkf.t-scEstByScEkf.t(1))/60/60, -3 * reshape(scEstByGsEkf.P_list(2,2,:), 1, length(scEstByScEkf.t)).^0.5, ':g')
legend('error estimated by sc', 'error estimated by gs')
ylim(error.scPosSigma*[-1, 1])
xlabel('time[h]')
ylabel('position error [km]')
nexttile
hold on
title("position Y error")
plot((scEstByScEkf.t-scEstByScEkf.t(1))/60/60, (scEstByScEkf.state(2,:) - scTrue.state(2,:)),'-r')
plot((scEstByGsEkf.t-scEstByScEkf.t(1))/60/60, (scEstByGsEkf.state(2,:) - scTrue.state(2,:)),'-g')
plot((scEstByScEkf.t-scEstByScEkf.t(1))/60/60,  1 * reshape(scEstByScEkf.P_list(3,3,:), 1, length(scEstByScEkf.t)).^0.5, '--r')
plot((scEstByScEkf.t-scEstByScEkf.t(1))/60/60, -1 * reshape(scEstByScEkf.P_list(3,3,:), 1, length(scEstByScEkf.t)).^0.5, '--r')
plot((scEstByScEkf.t-scEstByScEkf.t(1))/60/60,  3 * reshape(scEstByScEkf.P_list(3,3,:), 1, length(scEstByScEkf.t)).^0.5, ':r')
plot((scEstByScEkf.t-scEstByScEkf.t(1))/60/60, -3 * reshape(scEstByScEkf.P_list(3,3,:), 1, length(scEstByScEkf.t)).^0.5, ':r')
plot((scEstByGsEkf.t-scEstByScEkf.t(1))/60/60,  1 * reshape(scEstByGsEkf.P_list(3,3,:), 1, length(scEstByScEkf.t)).^0.5, '--g')
plot((scEstByGsEkf.t-scEstByScEkf.t(1))/60/60, -1 * reshape(scEstByGsEkf.P_list(3,3,:), 1, length(scEstByScEkf.t)).^0.5, '--g')
plot((scEstByGsEkf.t-scEstByScEkf.t(1))/60/60,  3 * reshape(scEstByGsEkf.P_list(3,3,:), 1, length(scEstByScEkf.t)).^0.5, ':g')
plot((scEstByGsEkf.t-scEstByScEkf.t(1))/60/60, -3 * reshape(scEstByGsEkf.P_list(3,3,:), 1, length(scEstByScEkf.t)).^0.5, ':g')
% legend('error estimated by sc', 'error estimated by gs')
ylim(error.scPosSigma*[-1, 1])
xlabel('time[h]')
ylabel('position error [km]')
nexttile
hold on
title("position Z error")
plot((scEstByScEkf.t-scEstByScEkf.t(1))/60/60, (scEstByScEkf.state(3,:) - scTrue.state(3,:)),'-r')
plot((scEstByGsEkf.t-scEstByScEkf.t(1))/60/60, (scEstByGsEkf.state(3,:) - scTrue.state(3,:)),'-g')
plot((scEstByScEkf.t-scEstByScEkf.t(1))/60/60,  1 * reshape(scEstByScEkf.P_list(4,4,:), 1, length(scEstByScEkf.t)).^0.5, '--r')
plot((scEstByScEkf.t-scEstByScEkf.t(1))/60/60, -1 * reshape(scEstByScEkf.P_list(4,4,:), 1, length(scEstByScEkf.t)).^0.5, '--r')
plot((scEstByScEkf.t-scEstByScEkf.t(1))/60/60,  3 * reshape(scEstByScEkf.P_list(4,4,:), 1, length(scEstByScEkf.t)).^0.5, ':r')
plot((scEstByScEkf.t-scEstByScEkf.t(1))/60/60, -3 * reshape(scEstByScEkf.P_list(4,4,:), 1, length(scEstByScEkf.t)).^0.5, ':r')
plot((scEstByGsEkf.t-scEstByScEkf.t(1))/60/60,  1 * reshape(scEstByGsEkf.P_list(4,4,:), 1, length(scEstByScEkf.t)).^0.5, '--g')
plot((scEstByGsEkf.t-scEstByScEkf.t(1))/60/60, -1 * reshape(scEstByGsEkf.P_list(4,4,:), 1, length(scEstByScEkf.t)).^0.5, '--g')
plot((scEstByGsEkf.t-scEstByScEkf.t(1))/60/60,  3 * reshape(scEstByGsEkf.P_list(4,4,:), 1, length(scEstByScEkf.t)).^0.5, ':g')
plot((scEstByGsEkf.t-scEstByScEkf.t(1))/60/60, -3 * reshape(scEstByGsEkf.P_list(4,4,:), 1, length(scEstByScEkf.t)).^0.5, ':g')
% legend('error estimated by sc', 'error estimated by gs')
ylim(error.scPosSigma*[-1, 1])
xlabel('time[h]')
ylabel('position error [km]')
nexttile
hold on
title("position error")
plot((scEstByScEkf.t-scEstByScEkf.t(1))/60/60, ((scEstByScEkf.state(1,:) - scTrue.state(1,:)).^2 + (scEstByScEkf.state(2,:) - scTrue.state(2,:)).^2 + (scEstByScEkf.state(3,:) - scTrue.state(3,:)).^2).^0.5 ,'-r')
plot((scEstByGsEkf.t-scEstByScEkf.t(1))/60/60, ((scEstByGsEkf.state(1,:) - scTrue.state(1,:)).^2 + (scEstByGsEkf.state(2,:) - scTrue.state(2,:)).^2 + (scEstByGsEkf.state(3,:) - scTrue.state(3,:)).^2).^0.5,'-g')
plot((scEstByScEkf.t-scEstByScEkf.t(1))/60/60,  1 * (reshape(scEstByScEkf.P_list(2,2,:), 1, length(scEstByScEkf.t)) + reshape(scEstByScEkf.P_list(3,3,:), 1, length(scEstByScEkf.t))+ reshape(scEstByScEkf.P_list(4,4,:), 1, length(scEstByScEkf.t))).^0.5, '--r')
% plot((scEstByScEkf.t-scEstByScEkf.t(1))/60/60, -1 * (reshape(scEstByScEkf.P_list(2,2,:), 1, length(scEstByScEkf.t)) + reshape(scEstByScEkf.P_list(3,3,:), 1, length(scEstByScEkf.t))+ reshape(scEstByScEkf.P_list(4,4,:), 1, length(scEstByScEkf.t))).^0.5, '--r')
plot((scEstByScEkf.t-scEstByScEkf.t(1))/60/60,  3 * (reshape(scEstByScEkf.P_list(2,2,:), 1, length(scEstByScEkf.t)) + reshape(scEstByScEkf.P_list(3,3,:), 1, length(scEstByScEkf.t))+ reshape(scEstByScEkf.P_list(4,4,:), 1, length(scEstByScEkf.t))).^0.5, ':r')
% plot((scEstByScEkf.t-scEstByScEkf.t(1))/60/60, -3 * (reshape(scEstByScEkf.P_list(2,2,:), 1, length(scEstByScEkf.t)) + reshape(scEstByScEkf.P_list(3,3,:), 1, length(scEstByScEkf.t))+ reshape(scEstByScEkf.P_list(4,4,:), 1, length(scEstByScEkf.t))).^0.5, ':r')
plot((scEstByGsEkf.t-scEstByScEkf.t(1))/60/60,  1 * (reshape(scEstByGsEkf.P_list(2,2,:), 1, length(scEstByScEkf.t)) + reshape(scEstByGsEkf.P_list(3,3,:), 1, length(scEstByScEkf.t))+ reshape(scEstByGsEkf.P_list(4,4,:), 1, length(scEstByScEkf.t))).^0.5, '--g')
% plot((scEstByGsEkf.t-scEstByScEkf.t(1))/60/60, -1 * (reshape(scEstByGsEkf.P_list(2,2,:), 1, length(scEstByScEkf.t)) + reshape(scEstByGsEkf.P_list(3,3,:), 1, length(scEstByScEkf.t))+ reshape(scEstByGsEkf.P_list(4,4,:), 1, length(scEstByScEkf.t))).^0.5, '--g')
plot((scEstByGsEkf.t-scEstByScEkf.t(1))/60/60,  3 * (reshape(scEstByGsEkf.P_list(2,2,:), 1, length(scEstByScEkf.t)) + reshape(scEstByGsEkf.P_list(3,3,:), 1, length(scEstByScEkf.t))+ reshape(scEstByGsEkf.P_list(4,4,:), 1, length(scEstByScEkf.t))).^0.5, ':g')
% plot((scEstByGsEkf.t-scEstByScEkf.t(1))/60/60, -3 * (reshape(scEstByGsEkf.P_list(2,2,:), 1, length(scEstByScEkf.t)) + reshape(scEstByGsEkf.P_list(3,3,:), 1, length(scEstByScEkf.t))+ reshape(scEstByGsEkf.P_list(4,4,:), 1, length(scEstByScEkf.t))).^0.5, ':g')
% legend('error estimated by sc', 'error estimated by gs')
ylim(error.scPosSigma*[0, 1])
xlabel('time[h]')
ylabel('position error [km]')




figure(3)
tiledlayout(1,4)
nexttile
hold on
title("velocity X error")
plot((scEstByScEkf.t-scEstByScEkf.t(1))/60/60, (scEstByScEkf.state(4,:) - scTrue.state(4,:)),'-r')
plot((scEstByGsEkf.t-scEstByScEkf.t(1))/60/60, (scEstByGsEkf.state(4,:) - scTrue.state(4,:)),'-g')
plot((scEstByScEkf.t-scEstByScEkf.t(1))/60/60,  1 * reshape(scEstByScEkf.P_list(5,5,:), 1, length(scEstByScEkf.t)).^0.5, '--r')
plot((scEstByScEkf.t-scEstByScEkf.t(1))/60/60, -1 * reshape(scEstByScEkf.P_list(5,5,:), 1, length(scEstByScEkf.t)).^0.5, '--r')
plot((scEstByScEkf.t-scEstByScEkf.t(1))/60/60,  3 * reshape(scEstByScEkf.P_list(5,5,:), 1, length(scEstByScEkf.t)).^0.5, ':r')
plot((scEstByScEkf.t-scEstByScEkf.t(1))/60/60, -3 * reshape(scEstByScEkf.P_list(5,5,:), 1, length(scEstByScEkf.t)).^0.5, ':r')
plot((scEstByGsEkf.t-scEstByScEkf.t(1))/60/60,  1 * reshape(scEstByGsEkf.P_list(5,5,:), 1, length(scEstByScEkf.t)).^0.5, '--g')
plot((scEstByGsEkf.t-scEstByScEkf.t(1))/60/60, -1 * reshape(scEstByGsEkf.P_list(5,5,:), 1, length(scEstByScEkf.t)).^0.5, '--g')
plot((scEstByGsEkf.t-scEstByScEkf.t(1))/60/60,  3 * reshape(scEstByGsEkf.P_list(5,5,:), 1, length(scEstByScEkf.t)).^0.5, ':g')
plot((scEstByGsEkf.t-scEstByScEkf.t(1))/60/60, -3 * reshape(scEstByGsEkf.P_list(5,5,:), 1, length(scEstByScEkf.t)).^0.5, ':g')
legend('error estimated by sc', 'error estimated by gs')
ylim(error.scVelSigma* [-1, 1])
xlabel('time[h]')
ylabel('velocity error [km/s]')
nexttile
hold on
title("velocity Y error")
plot((scEstByScEkf.t-scEstByScEkf.t(1))/60/60, (scEstByScEkf.state(5,:) - scTrue.state(5,:)),'-r')
plot((scEstByGsEkf.t-scEstByScEkf.t(1))/60/60, (scEstByGsEkf.state(5,:) - scTrue.state(5,:)),'-g')
plot((scEstByScEkf.t-scEstByScEkf.t(1))/60/60,  1 * reshape(scEstByScEkf.P_list(6,6,:), 1, length(scEstByScEkf.t)).^0.5, '--r')
plot((scEstByScEkf.t-scEstByScEkf.t(1))/60/60, -1 * reshape(scEstByScEkf.P_list(6,6,:), 1, length(scEstByScEkf.t)).^0.5, '--r')
plot((scEstByScEkf.t-scEstByScEkf.t(1))/60/60,  3 * reshape(scEstByScEkf.P_list(6,6,:), 1, length(scEstByScEkf.t)).^0.5, ':r')
plot((scEstByScEkf.t-scEstByScEkf.t(1))/60/60, -3 * reshape(scEstByScEkf.P_list(6,6,:), 1, length(scEstByScEkf.t)).^0.5, ':r')
plot((scEstByGsEkf.t-scEstByScEkf.t(1))/60/60,  1 * reshape(scEstByGsEkf.P_list(6,6,:), 1, length(scEstByScEkf.t)).^0.5, '--g')
plot((scEstByGsEkf.t-scEstByScEkf.t(1))/60/60, -1 * reshape(scEstByGsEkf.P_list(6,6,:), 1, length(scEstByScEkf.t)).^0.5, '--g')
plot((scEstByGsEkf.t-scEstByScEkf.t(1))/60/60,  3 * reshape(scEstByGsEkf.P_list(6,6,:), 1, length(scEstByScEkf.t)).^0.5, ':g')
plot((scEstByGsEkf.t-scEstByScEkf.t(1))/60/60, -3 * reshape(scEstByGsEkf.P_list(6,6,:), 1, length(scEstByScEkf.t)).^0.5, ':g')
% legend('error estimated by sc', 'error estimated by gs')
ylim(error.scVelSigma* [-1, 1])
xlabel('time[h]')
ylabel('velocity error [km/s]')
nexttile
hold on
title("velocity Z error")
plot((scEstByScEkf.t-scEstByScEkf.t(1))/60/60, (scEstByScEkf.state(6,:) - scTrue.state(6,:)),'-r')
plot((scEstByGsEkf.t-scEstByScEkf.t(1))/60/60, (scEstByGsEkf.state(6,:) - scTrue.state(6,:)),'-g')
plot((scEstByScEkf.t-scEstByScEkf.t(1))/60/60,  1 * reshape(scEstByScEkf.P_list(7,7,:), 1, length(scEstByScEkf.t)).^0.5, '--r')
plot((scEstByScEkf.t-scEstByScEkf.t(1))/60/60, -1 * reshape(scEstByScEkf.P_list(7,7,:), 1, length(scEstByScEkf.t)).^0.5, '--r')
plot((scEstByScEkf.t-scEstByScEkf.t(1))/60/60,  3 * reshape(scEstByScEkf.P_list(7,7,:), 1, length(scEstByScEkf.t)).^0.5, ':r')
plot((scEstByScEkf.t-scEstByScEkf.t(1))/60/60, -3 * reshape(scEstByScEkf.P_list(7,7,:), 1, length(scEstByScEkf.t)).^0.5, ':r')
plot((scEstByGsEkf.t-scEstByScEkf.t(1))/60/60,  1 * reshape(scEstByGsEkf.P_list(7,7,:), 1, length(scEstByScEkf.t)).^0.5, '--g')
plot((scEstByGsEkf.t-scEstByScEkf.t(1))/60/60, -1 * reshape(scEstByGsEkf.P_list(7,7,:), 1, length(scEstByScEkf.t)).^0.5, '--g')
plot((scEstByGsEkf.t-scEstByScEkf.t(1))/60/60,  3 * reshape(scEstByGsEkf.P_list(7,7,:), 1, length(scEstByScEkf.t)).^0.5, ':g')
plot((scEstByGsEkf.t-scEstByScEkf.t(1))/60/60, -3 * reshape(scEstByGsEkf.P_list(7,7,:), 1, length(scEstByScEkf.t)).^0.5, ':g')
ylim(error.scVelSigma* [-1, 1])
xlabel('time[h]')
ylabel('velocity error [km/s]')
nexttile
hold on
title("velocity error")
plot((scEstByScEkf.t-scEstByScEkf.t(1))/60/60, ((scEstByScEkf.state(4,:) - scTrue.state(4,:)).^2 + (scEstByScEkf.state(5,:) - scTrue.state(5,:)).^2 + (scEstByScEkf.state(6,:) - scTrue.state(6,:)).^2).^0.5 ,'-r')
plot((scEstByGsEkf.t-scEstByScEkf.t(1))/60/60, ((scEstByGsEkf.state(4,:) - scTrue.state(4,:)).^2 + (scEstByGsEkf.state(5,:) - scTrue.state(5,:)).^2 + (scEstByGsEkf.state(6,:) - scTrue.state(6,:)).^2).^0.5,'-g')
plot((scEstByScEkf.t-scEstByScEkf.t(1))/60/60,  1 * (reshape(scEstByScEkf.P_list(5,5,:), 1, length(scEstByScEkf.t)) + reshape(scEstByScEkf.P_list(6,6,:), 1, length(scEstByScEkf.t))+ reshape(scEstByScEkf.P_list(7,7,:), 1, length(scEstByScEkf.t))).^0.5, '--r')
% plot((scEstByScEkf.t-scEstByScEkf.t(1))/60/60, -1 * (reshape(scEstByScEkf.P_list(5,5,:), 1, length(scEstByScEkf.t)) + reshape(scEstByScEkf.P_list(6,6,:), 1, length(scEstByScEkf.t))+ reshape(scEstByScEkf.P_list(7,7,:), 1, length(scEstByScEkf.t))).^0.5, '--r')
plot((scEstByScEkf.t-scEstByScEkf.t(1))/60/60,  3 * (reshape(scEstByScEkf.P_list(5,5,:), 1, length(scEstByScEkf.t)) + reshape(scEstByScEkf.P_list(6,6,:), 1, length(scEstByScEkf.t))+ reshape(scEstByScEkf.P_list(7,7,:), 1, length(scEstByScEkf.t))).^0.5, ':r')
% plot((scEstByScEkf.t-scEstByScEkf.t(1))/60/60, -3 * (reshape(scEstByScEkf.P_list(5,5,:), 1, length(scEstByScEkf.t)) + reshape(scEstByScEkf.P_list(6,6,:), 1, length(scEstByScEkf.t))+ reshape(scEstByScEkf.P_list(7,7,:), 1, length(scEstByScEkf.t))).^0.5, ':r')
plot((scEstByGsEkf.t-scEstByScEkf.t(1))/60/60,  1 * (reshape(scEstByGsEkf.P_list(5,5,:), 1, length(scEstByScEkf.t)) + reshape(scEstByGsEkf.P_list(6,6,:), 1, length(scEstByScEkf.t))+ reshape(scEstByGsEkf.P_list(7,7,:), 1, length(scEstByScEkf.t))).^0.5, '--g')
% plot((scEstByGsEkf.t-scEstByScEkf.t(1))/60/60, -1 * (reshape(scEstByGsEkf.P_list(5,5,:), 1, length(scEstByScEkf.t)) + reshape(scEstByGsEkf.P_list(6,6,:), 1, length(scEstByScEkf.t))+ reshape(scEstByGsEkf.P_list(7,7,:), 1, length(scEstByScEkf.t))).^0.5, '--g')
plot((scEstByGsEkf.t-scEstByScEkf.t(1))/60/60,  3 * (reshape(scEstByGsEkf.P_list(5,5,:), 1, length(scEstByScEkf.t)) + reshape(scEstByGsEkf.P_list(6,6,:), 1, length(scEstByScEkf.t))+ reshape(scEstByGsEkf.P_list(7,7,:), 1, length(scEstByScEkf.t))).^0.5, ':g')
% plot((scEstByGsEkf.t-scEstByScEkf.t(1))/60/60, -3 * (reshape(scEstByGsEkf.P_list(5,5,:), 1, length(scEstByScEkf.t)) + reshape(scEstByGsEkf.P_list(6,6,:), 1, length(scEstByScEkf.t))+ reshape(scEstByGsEkf.P_list(7,7,:), 1, length(scEstByScEkf.t))).^0.5, ':g')
ylim(error.scVelSigma* [0, 1])
xlabel('time[h]')
ylabel('velocity error [km/s]')




end