function showResult(scTrue,scEstByScEkf,scEstByGsEkf,error,a,gsTrue,gs,resultPath)

% 各要素ごとに出力
f1 = figure('visible', 'off');
f1.OuterPosition = [100 100 500 500];
semilogy((scEstByScEkf.t-scEstByScEkf.t(1))/60/60/24, abs(scEstByScEkf.clockError),...
         (scEstByGsEkf.t-scEstByScEkf.t(1))/60/60/24, abs(scEstByGsEkf.clockError),...
         (scEstByScEkf.t-scEstByScEkf.t(1))/60/60/24, 1 * reshape(scEstByScEkf.P_list(1,1,:), [], size(scEstByScEkf.P_list,3)).^0.5, ...
         (scEstByGsEkf.t-scEstByScEkf.t(1))/60/60/24, 1 * reshape(scEstByGsEkf.P_list(1,1,:), [], size(scEstByScEkf.P_list,3)).^0.5)
xlabel('day')
ylim([1e-8 1e0])
legend('estimated by sc', 'estimated by gs', '1\sigma estimated by sc', '1\sigma estimated by gs')
title("clock error befor long communication disconnection ")
file1 = [resultPath,'/clockError',num2str(a+1),'.png'];
saveas(f1, file1)

f2 = figure('visible', 'off');
f2.OuterPosition = [100 100 1200 600];
subplot(1,4,1)
% nexttile

semilogy((scEstByScEkf.t-scEstByScEkf.t(1))/60/60/24, abs(scEstByScEkf.state(1,:) - scTrue.state(1,:)),'-r',...
         (scEstByGsEkf.t-scEstByScEkf.t(1))/60/60/24, abs(scEstByGsEkf.state(1,:) - scTrue.state(1,:)),'-g',...
         (scEstByScEkf.t-scEstByScEkf.t(1))/60/60/24,  1 * reshape(scEstByScEkf.P_list(2,2,:),  [], size(scEstByScEkf.P_list,3)).^0.5, '--r',...
         (scEstByScEkf.t-scEstByScEkf.t(1))/60/60/24,  3 * reshape(scEstByScEkf.P_list(2,2,:),  [], size(scEstByScEkf.P_list,3)).^0.5, ':r',...
         (scEstByGsEkf.t-scEstByScEkf.t(1))/60/60/24,  1 * reshape(scEstByGsEkf.P_list(2,2,:),  [], size(scEstByScEkf.P_list,3)).^0.5, '--g',...
         (scEstByGsEkf.t-scEstByScEkf.t(1))/60/60/24,  3 * reshape(scEstByGsEkf.P_list(2,2,:),  [], size(scEstByScEkf.P_list,3)).^0.5, ':g')

legend('error estimated by sc', 'error estimated by gs', '1\sigma of P estimated by sc', '3\sigma of P estimated by sc', '1\sigma of P estimated by gs', '3\sigma of P estimated by gs')
ylim([1e-1 1e6])
xlabel('day')
ylabel('position error [km]')
title("position X error")
% nexttile
% hold on
subplot(1,4,2)
semilogy((scEstByScEkf.t-scEstByScEkf.t(1))/60/60/24, abs(scEstByScEkf.state(2,:) - scTrue.state(2,:)),'-r',...
         (scEstByGsEkf.t-scEstByScEkf.t(1))/60/60/24, abs(scEstByGsEkf.state(2,:) - scTrue.state(2,:)),'-g',...
         (scEstByScEkf.t-scEstByScEkf.t(1))/60/60/24,  1 * reshape(scEstByScEkf.P_list(3,3,:),  [], size(scEstByScEkf.P_list,3)).^0.5, '--r',...
         (scEstByScEkf.t-scEstByScEkf.t(1))/60/60/24,  3 * reshape(scEstByScEkf.P_list(3,3,:),  [], size(scEstByScEkf.P_list,3)).^0.5, ':r',...
         (scEstByGsEkf.t-scEstByScEkf.t(1))/60/60/24,  1 * reshape(scEstByGsEkf.P_list(3,3,:),  [], size(scEstByScEkf.P_list,3)).^0.5, '--g',...
         (scEstByGsEkf.t-scEstByScEkf.t(1))/60/60/24,  3 * reshape(scEstByGsEkf.P_list(3,3,:),  [], size(scEstByScEkf.P_list,3)).^0.5, ':g')
% legend('error estimated by sc', 'error estimated by gs')
ylim([1e-1 1e6])
xlabel('day')
ylabel('position error [km]')
title("position Y error")
% nexttile
subplot(1,4,3)
semilogy((scEstByScEkf.t-scEstByScEkf.t(1))/60/60/24, abs(scEstByScEkf.state(3,:) - scTrue.state(3,:)),'-r',...
         (scEstByGsEkf.t-scEstByScEkf.t(1))/60/60/24, abs(scEstByGsEkf.state(3,:) - scTrue.state(3,:)),'-g',...
         (scEstByScEkf.t-scEstByScEkf.t(1))/60/60/24,  1 * reshape(scEstByScEkf.P_list(4,4,:),  [], size(scEstByScEkf.P_list,3)).^0.5, '--r',...
         (scEstByScEkf.t-scEstByScEkf.t(1))/60/60/24,  3 * reshape(scEstByScEkf.P_list(4,4,:),  [], size(scEstByScEkf.P_list,3)).^0.5, ':r',...
         (scEstByGsEkf.t-scEstByScEkf.t(1))/60/60/24,  1 * reshape(scEstByGsEkf.P_list(4,4,:),  [], size(scEstByScEkf.P_list,3)).^0.5, '--g',...
         (scEstByGsEkf.t-scEstByScEkf.t(1))/60/60/24,  3 * reshape(scEstByGsEkf.P_list(4,4,:),  [], size(scEstByScEkf.P_list,3)).^0.5, ':g')
xlabel('day')
ylabel('position error [km]')
ylim([1e-1 1e6])
title("position Z error")
% nexttile
subplot(1,4,4)
semilogy((scEstByScEkf.t-scEstByScEkf.t(1))/60/60/24, ((scEstByScEkf.state(1,:) - scTrue.state(1,:)).^2 + (scEstByScEkf.state(2,:) - scTrue.state(2,:)).^2 + (scEstByScEkf.state(3,:) - scTrue.state(3,:)).^2).^0.5, '-r',...
         (scEstByGsEkf.t-scEstByScEkf.t(1))/60/60/24, ((scEstByGsEkf.state(1,:) - scTrue.state(1,:)).^2 + (scEstByGsEkf.state(2,:) - scTrue.state(2,:)).^2 + (scEstByGsEkf.state(3,:) - scTrue.state(3,:)).^2).^0.5, '-g',...
         (scEstByScEkf.t-scEstByScEkf.t(1))/60/60/24,  1 * (reshape(scEstByScEkf.P_list(2,2,:),  [], size(scEstByScEkf.P_list,3)) + reshape(scEstByScEkf.P_list(3,3,:), [], size(scEstByScEkf.P_list,3))+ reshape(scEstByScEkf.P_list(4,4,:),  [], size(scEstByScEkf.P_list,3))).^0.5, '--r',...
         (scEstByScEkf.t-scEstByScEkf.t(1))/60/60/24,  3 * (reshape(scEstByScEkf.P_list(2,2,:),  [], size(scEstByScEkf.P_list,3)) + reshape(scEstByScEkf.P_list(3,3,:), [], size(scEstByScEkf.P_list,3))+ reshape(scEstByScEkf.P_list(4,4,:),  [], size(scEstByScEkf.P_list,3))).^0.5, ':r',...
         (scEstByGsEkf.t-scEstByScEkf.t(1))/60/60/24,  1 * (reshape(scEstByGsEkf.P_list(2,2,:),  [], size(scEstByScEkf.P_list,3)) + reshape(scEstByGsEkf.P_list(3,3,:), [], size(scEstByScEkf.P_list,3))+ reshape(scEstByGsEkf.P_list(4,4,:),  [], size(scEstByScEkf.P_list,3))).^0.5, '--g',...
         (scEstByGsEkf.t-scEstByScEkf.t(1))/60/60/24,  3 * (reshape(scEstByGsEkf.P_list(2,2,:),  [], size(scEstByScEkf.P_list,3)) + reshape(scEstByGsEkf.P_list(3,3,:), [], size(scEstByScEkf.P_list,3))+ reshape(scEstByGsEkf.P_list(4,4,:),  [], size(scEstByScEkf.P_list,3))).^0.5, ':g')

xlabel('day')
ylabel('position error [km]')
ylim([1e-1 1e6])
title("position error")
file2 = [resultPath,'/positionError',num2str(a+1),'.png'];
saveas(f2, file2)



f3 = figure('visible', 'off');
f3.OuterPosition = [100 100 1200 600];
subplot(1,4,1)
% nexttile

semilogy((scEstByScEkf.t-scEstByScEkf.t(1))/60/60/24, abs(scEstByScEkf.state(4,:) - scTrue.state(4,:)),'-r',...
         (scEstByGsEkf.t-scEstByScEkf.t(1))/60/60/24, abs(scEstByGsEkf.state(4,:) - scTrue.state(4,:)),'-g',...
         (scEstByScEkf.t-scEstByScEkf.t(1))/60/60/24,  1 * reshape(scEstByScEkf.P_list(5,5,:),  [], size(scEstByScEkf.P_list,3)).^0.5, '--r',...
         (scEstByScEkf.t-scEstByScEkf.t(1))/60/60/24,  3 * reshape(scEstByScEkf.P_list(5,5,:),  [], size(scEstByScEkf.P_list,3)).^0.5, ':r',...
         (scEstByGsEkf.t-scEstByScEkf.t(1))/60/60/24,  1 * reshape(scEstByGsEkf.P_list(5,5,:),  [], size(scEstByScEkf.P_list,3)).^0.5, '--g',...
         (scEstByGsEkf.t-scEstByScEkf.t(1))/60/60/24,  3 * reshape(scEstByGsEkf.P_list(5,5,:),  [], size(scEstByScEkf.P_list,3)).^0.5, ':g')
legend('error estimated by sc', 'error estimated by gs', '1\sigma of P estimated by sc', '3\sigma of P estimated by sc', '1\sigma of P estimated by gs', '3\sigma of P estimated by gs')
xlabel('day')
ylabel('velocity error [km/s]')
ylim([1e-6 1e0])
title("velocity X error")
% nexttile
subplot(1,4,2)
semilogy((scEstByScEkf.t-scEstByScEkf.t(1))/60/60/24, abs(scEstByScEkf.state(5,:) - scTrue.state(5,:)),'-r',...
         (scEstByGsEkf.t-scEstByScEkf.t(1))/60/60/24, abs(scEstByGsEkf.state(5,:) - scTrue.state(5,:)),'-g',...
         (scEstByScEkf.t-scEstByScEkf.t(1))/60/60/24,  1 * reshape(scEstByScEkf.P_list(6,6,:), [], size(scEstByScEkf.P_list,3)).^0.5, '--r',...
         (scEstByScEkf.t-scEstByScEkf.t(1))/60/60/24,  3 * reshape(scEstByScEkf.P_list(6,6,:), [], size(scEstByScEkf.P_list,3)).^0.5, ':r',...
         (scEstByGsEkf.t-scEstByScEkf.t(1))/60/60/24,  1 * reshape(scEstByGsEkf.P_list(6,6,:), [], size(scEstByScEkf.P_list,3)).^0.5, '--g',...
         (scEstByGsEkf.t-scEstByScEkf.t(1))/60/60/24,  3 * reshape(scEstByGsEkf.P_list(6,6,:), [], size(scEstByScEkf.P_list,3)).^0.5, ':g')
xlabel('day')
ylabel('velocity error [km/s]')
ylim([1e-6 1e0])
title("velocity Y error")
% nexttile
subplot(1,4,3)
semilogy((scEstByScEkf.t-scEstByScEkf.t(1))/60/60/24, abs(scEstByScEkf.state(6,:) - scTrue.state(6,:)),'-r',...
         (scEstByGsEkf.t-scEstByScEkf.t(1))/60/60/24, abs(scEstByGsEkf.state(6,:) - scTrue.state(6,:)),'-g',...
         (scEstByScEkf.t-scEstByScEkf.t(1))/60/60/24,  1 * reshape(scEstByScEkf.P_list(7,7,:),  [], size(scEstByScEkf.P_list,3)).^0.5, '--r',...
         (scEstByScEkf.t-scEstByScEkf.t(1))/60/60/24,  3 * reshape(scEstByScEkf.P_list(7,7,:),  [], size(scEstByScEkf.P_list,3)).^0.5, ':r',...
         (scEstByGsEkf.t-scEstByScEkf.t(1))/60/60/24,  1 * reshape(scEstByGsEkf.P_list(7,7,:),  [], size(scEstByScEkf.P_list,3)).^0.5, '--g',...
         (scEstByGsEkf.t-scEstByScEkf.t(1))/60/60/24,  3 * reshape(scEstByGsEkf.P_list(7,7,:),  [], size(scEstByScEkf.P_list,3)).^0.5, ':g')
xlabel('day')
ylabel('velocity error [km/s]')
ylim([1e-6 1e0])
title("velocity Z error")
% nexttile
subplot(1,4,4)
semilogy((scEstByScEkf.t-scEstByScEkf.t(1))/60/60/24, ((scEstByScEkf.state(4,:) - scTrue.state(4,:)).^2 + (scEstByScEkf.state(5,:) - scTrue.state(5,:)).^2 + (scEstByScEkf.state(6,:) - scTrue.state(6,:)).^2).^0.5, '-r',...
         (scEstByGsEkf.t-scEstByScEkf.t(1))/60/60/24, ((scEstByGsEkf.state(4,:) - scTrue.state(4,:)).^2 + (scEstByGsEkf.state(5,:) - scTrue.state(5,:)).^2 + (scEstByGsEkf.state(6,:) - scTrue.state(6,:)).^2).^0.5, '-g',...
         (scEstByScEkf.t-scEstByScEkf.t(1))/60/60/24,  1 * (reshape(scEstByScEkf.P_list(5,5,:),  [], size(scEstByScEkf.P_list,3)) + reshape(scEstByScEkf.P_list(6,6,:),  [], size(scEstByScEkf.P_list,3))+ reshape(scEstByScEkf.P_list(7,7,:),  [], size(scEstByScEkf.P_list,3))).^0.5, '--r',...
         (scEstByScEkf.t-scEstByScEkf.t(1))/60/60/24,  3 * (reshape(scEstByScEkf.P_list(5,5,:),  [], size(scEstByScEkf.P_list,3)) + reshape(scEstByScEkf.P_list(6,6,:),  [], size(scEstByScEkf.P_list,3))+ reshape(scEstByScEkf.P_list(7,7,:),  [], size(scEstByScEkf.P_list,3))).^0.5, ':r',...
         (scEstByGsEkf.t-scEstByScEkf.t(1))/60/60/24,  1 * (reshape(scEstByGsEkf.P_list(5,5,:),  [], size(scEstByScEkf.P_list,3)) + reshape(scEstByGsEkf.P_list(6,6,:),  [], size(scEstByScEkf.P_list,3))+ reshape(scEstByGsEkf.P_list(7,7,:),  [], size(scEstByScEkf.P_list,3))).^0.5, '--g',...
         (scEstByGsEkf.t-scEstByScEkf.t(1))/60/60/24,  3 * (reshape(scEstByGsEkf.P_list(5,5,:),  [], size(scEstByScEkf.P_list,3)) + reshape(scEstByGsEkf.P_list(6,6,:),  [], size(scEstByScEkf.P_list,3))+ reshape(scEstByGsEkf.P_list(7,7,:),  [], size(scEstByScEkf.P_list,3))).^0.5, ':g')
xlabel('day')
ylabel('velocity error [km/s]')
ylim([1e-6 1e0])
title("velocity error")
file3 = [resultPath,'/velocityError',num2str(a+1),'.png'];
saveas(f3, file3)

% 通信成立性のplot

if mod(a,3) == 0 || mod(a,3) == 2
    f4 = figure('visible', 'off');
    f4.OuterPosition = [100 100 500 500];
    semilogy((gsTrue.t_drList-scTrue.t(1))/24/60/60, scTrue.targetError_dtList,'p',...
          (gsTrue.t_drList-scTrue.t(1))/24/60/60,scTrue.pointError_dtList, 'o',...
          (gsTrue.t_drList-scTrue.t(1))/24/60/60,5.3e-6 * ones(1,length(gsTrue.t_drList)), '-')
    xlabel("time downlink is received")
    ylabel("downlink Pointing error")
    legend("target error", "pointing error", "required pointing accuracy")

    %通信ができた割合を算出
    sucsessDown = gsTrue.snr_drList > gs.reqSnr_down;
    sucsessDownRate = sum(sucsessDown)/length(sucsessDown) * 100;
    title(['Downlink Availability is ',num2str(sucsessDownRate,4), '%'])
    ylim([1e-8 5e-4])
    file4 = [resultPath,'/pointingError',num2str(a+1),'.png'];
    saveas(f4, file4)
end
end