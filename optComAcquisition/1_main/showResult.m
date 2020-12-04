function showResult()

% % 時計誤差の時間履歴
% figure(1)
% title('clock error of spacecraft')
% hold on
% plot((time.list-time.list(1))/60/60, scEst.resClockError)
% xlabel('time [h]')
% ylabel('clock error [s]')
% hold off
% 
% % 軌道決定精度をplotする
% % 探査機位置誤差の時間履歴
% figure(2)
% tiledlayout(3,1)
% nexttile
% title('position error of estimated value')
% hold on
% plot((time.list-time.list(1))/60/60, scEst.state(1,:) - scTrue.state(1,:))
% plot((time.list-time.list(1))/60/60, scEst.state(2,:) - scTrue.state(2,:))
% plot((time.list-time.list(1))/60/60, scEst.state(3,:) - scTrue.state(3,:))
% plot((time.list-time.list(1))/60/60,...
%     ( (scEst.state(1,:) - scTrue.state(1,:)).^2 + (scEst.state(2,:) - scTrue.state(2,:)).^2 + (scEst.state(3,:) - scTrue.state(3,:)).^2).^0.5)
% xlabel('time [h]')
% ylabel('position error [km]')
% legend('x', 'y', 'z','length')
% % ylim([-1000 1000])
% hold off
% nexttile
% title('velocity error of estimated value')
% hold on
% plot((time.list-time.list(1))/60/60, scEst.state(4,:) - scTrue.state(4,:))
% plot((time.list-time.list(1))/60/60, scEst.state(5,:) - scTrue.state(5,:))
% plot((time.list-time.list(1))/60/60, scEst.state(6,:) - scTrue.state(6,:))
% plot((time.list-time.list(1))/60/60,...
%     ( (scEst.state(4,:) - scTrue.state(4,:)).^2 + (scEst.state(5,:) - scTrue.state(5,:)).^2 + (scEst.state(6,:) - scTrue.state(6,:)).^2).^0.5)
% xlabel('time [h]')
% ylabel('position error [km]')
% legend('x', 'y', 'z','speed')
% % ylim([-1000 1000])
% hold off
% nexttile
% title('downlink direction error')
% hold on
% plot((time.list(1:length(time.list)) - time.list(1))/60/60, scEst.azmDown - scTrue.azmDown)
% plot((time.list(1:length(time.list)) - time.list(1))/60/60, scEst.elvDown - scTrue.elvDown)
% plot((time.list(1:length(time.list)) - time.list(1))/60/60, ((scEst.azmDown - scTrue.azmDown).^2 + (scEst.elvDown - scTrue.elvDown).^2).^0.5);
% hold off
% xlabel('time [h]')
% ylabel('angle [rad]')
% legend('azimuth', 'elevation', 'angle')
% 
% % 地上局の推定値
% figure(3)
% tiledlayout(2,1)
% nexttile
% title('position error of estimated value')
% hold on
% plot((time.list-time.list(1))/60/60, scEstGs.state(1,:) - scTrue.state(1,:))
% plot((time.list-time.list(1))/60/60, scEstGs.state(2,:) - scTrue.state(2,:))
% plot((time.list-time.list(1))/60/60, scEstGs.state(3,:) - scTrue.state(3,:))
% plot((time.list-time.list(1))/60/60,...
%     ( (scEstGs.state(1,:) - scTrue.state(1,:)).^2 + (scEstGs.state(2,:) - scTrue.state(2,:)).^2 + (scEstGs.state(3,:) - scTrue.state(3,:)).^2).^0.5)
% xlabel('time [h]')
% ylabel('position error [km]')
% legend('x', 'y', 'z','length')
% nexttile
% title('velocity error of estimated value')
% hold on
% plot((time.list-time.list(1))/60/60, scEstGs.state(4,:) - scTrue.state(4,:))
% plot((time.list-time.list(1))/60/60, scEstGs.state(5,:) - scTrue.state(5,:))
% plot((time.list-time.list(1))/60/60, scEstGs.state(6,:) - scTrue.state(6,:))
% plot((time.list-time.list(1))/60/60,...
%     ( (scEstGs.state(4,:) - scTrue.state(4,:)).^2 + (scEstGs.state(5,:) - scTrue.state(5,:)).^2 + (scEstGs.state(6,:) - scTrue.state(6,:)).^2).^0.5)
% xlabel('time [h]')
% ylabel('velocity error [km]')
% legend('x', 'y', 'z','speed')
% 
% figure(4)
% title('estemated/true position of spacecraft')
% hold on
% plot3(scTrue.state(1,:),scTrue.state(2,:),scTrue.state(3,:))
% plot3(scEst.state(1,:),scEst.state(2,:),scEst.state(3,:))
% plot3(scEstGs.state(1,:),scEstGs.state(2,:),scEstGs.state(3,:))
% xlabel('x [km]')
% ylabel('y [km]')
% zlabel('z [km]')
% legend('true', 'estimated(sc)', 'estimated(gs)')
% hold off
end