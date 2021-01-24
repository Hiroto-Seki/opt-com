% 結果を格納していく
% 長期通信断絶前
result(1).posErrorSc(1,n) = abs(scEstByScEkf.state(1,end) - scTrue.state(1,end));
result(1).posErrorSc(2,n) = abs(scEstByScEkf.state(2,end) - scTrue.state(2,end));
result(1).posErrorSc(3,n) = abs(scEstByScEkf.state(3,end) - scTrue.state(3,end));
result(1).posErrorSc(4,n) = (result(1).posErrorSc(1,n)^2 + result(1).posErrorSc(2,n)^2 + result(1).posErrorSc(3,n)^2)^0.5;
result(1).posErrorSc(5,n) = scEstByScEkf.P_list(2,2,end)^0.5;
result(1).posErrorSc(6,n) = scEstByScEkf.P_list(3,3,end)^0.5;
result(1).posErrorSc(7,n) = scEstByScEkf.P_list(4,4,end)^0.5;
result(1).posErrorSc(8,n) = (result(1).posErrorSc(5,n)^2 + result(1).posErrorSc(6,n)^2 + result(1).posErrorSc(7,n)^2)^0.5;
result(1).velErrorSc(1,n) = abs(scEstByScEkf.state(4,end) - scTrue.state(4,end));
result(1).velErrorSc(2,n) = abs(scEstByScEkf.state(5,end) - scTrue.state(5,end));
result(1).velErrorSc(3,n) = abs(scEstByScEkf.state(6,end) - scTrue.state(6,end));
result(1).velErrorSc(4,n) = (result(1).velErrorSc(1,n)^2 + result(1).velErrorSc(2,n)^2 + result(1).velErrorSc(3,n)^2)^0.5;
result(1).velErrorSc(5,n) = scEstByScEkf.P_list(5,5,end)^0.5;
result(1).velErrorSc(6,n) = scEstByScEkf.P_list(6,6,end)^0.5;
result(1).velErrorSc(7,n) = scEstByScEkf.P_list(7,7,end)^0.5;
result(1).velErrorSc(8,n) = (result(1).velErrorSc(5,n)^2 + result(1).velErrorSc(6,n)^2 + result(1).velErrorSc(7,n)^2)^0.5;
result(1).clockErrorSc(1,n) = abs(scEstByScEkf.clockError(end));
result(1).clockErrorSc(2,n) = scEstByScEkf.P_list(1,1,end)^0.5;

result(1).posErrorGs(1,n) = abs(scEstByGsEkf.state(1,end) - scTrue.state(1,end));
result(1).posErrorGs(2,n) = abs(scEstByGsEkf.state(2,end) - scTrue.state(2,end));
result(1).posErrorGs(3,n) = abs(scEstByGsEkf.state(3,end) - scTrue.state(3,end));
result(1).posErrorGs(4,n) = (result(1).posErrorGs(1,n)^2 + result(1).posErrorGs(2,n)^2 + result(1).posErrorGs(3,n)^2)^0.5;
result(1).posErrorGs(5,n) = scEstByGsEkf.P_list(2,2,end)^0.5;
result(1).posErrorGs(6,n) = scEstByGsEkf.P_list(3,3,end)^0.5;
result(1).posErrorGs(7,n) = scEstByGsEkf.P_list(4,4,end)^0.5;
result(1).posErrorGs(8,n) = (result(1).posErrorGs(5,n)^2 + result(1).posErrorGs(6,n)^2 + result(1).posErrorGs(7,n)^2)^0.5;
result(1).velErrorGs(1,n) = abs(scEstByGsEkf.state(4,end) - scTrue.state(4,end));
result(1).velErrorGs(2,n) = abs(scEstByGsEkf.state(5,end) - scTrue.state(5,end));
result(1).velErrorGs(3,n) = abs(scEstByGsEkf.state(6,end) - scTrue.state(6,end));
result(1).velErrorGs(4,n) = (result(1).velErrorGs(1,n)^2 + result(1).velErrorGs(2,n)^2 + result(1).velErrorGs(3,n)^2)^0.5;
result(1).velErrorGs(5,n) = scEstByGsEkf.P_list(5,5,end)^0.5;
result(1).velErrorGs(6,n) = scEstByGsEkf.P_list(6,6,end)^0.5;
result(1).velErrorGs(7,n) = scEstByGsEkf.P_list(7,7,end)^0.5;
result(1).velErrorGs(8,n) = (result(1).velErrorGs(5,n)^2 + result(1).velErrorGs(6,n)^2 + result(1).velErrorGs(7,n)^2)^0.5;
result(1).clockErrorGs(1,n) = abs(scEstByGsEkf.clockError(end));
result(1).clockErrorGs(2,n) = scEstByGsEkf.P_list(1,1,end)^0.5;

result(1).downAvail(1,n) = sum(gsTrue.snr_drList > gs.reqSnr_down)/length(gsTrue.snr_drList) * 100;

% 長期通信断絶中
result(2).posErrorSc(1,n) = abs(scEstBySc_los.state(1,end) - scTrue_los.state(1,end));
result(2).posErrorSc(2,n) = abs(scEstBySc_los.state(2,end) - scTrue_los.state(2,end));
result(2).posErrorSc(3,n) = abs(scEstBySc_los.state(3,end) - scTrue_los.state(3,end));
result(2).posErrorSc(4,n) = (result(2).posErrorSc(1,n)^2 + result(2).posErrorSc(2,n)^2 + result(2).posErrorSc(3,n)^2)^0.5;
result(2).posErrorSc(5,n) = scEstBySc_los.P_list(2,2,end)^0.5;
result(2).posErrorSc(6,n) = scEstBySc_los.P_list(3,3,end)^0.5;
result(2).posErrorSc(7,n) = scEstBySc_los.P_list(4,4,end)^0.5;
result(2).posErrorSc(8,n) = (result(2).posErrorSc(5,n)^2 + result(2).posErrorSc(6,n)^2 + result(2).posErrorSc(7,n)^2)^0.5;
result(2).velErrorSc(1,n) = abs(scEstBySc_los.state(4,end) - scTrue_los.state(4,end));
result(2).velErrorSc(2,n) = abs(scEstBySc_los.state(5,end) - scTrue_los.state(5,end));
result(2).velErrorSc(3,n) = abs(scEstBySc_los.state(6,end) - scTrue_los.state(6,end));
result(2).velErrorSc(4,n) = (result(2).velErrorSc(1,n)^2 + result(2).velErrorSc(2,n)^2 + result(2).velErrorSc(3,n)^2)^0.5;
result(2).velErrorSc(5,n) = scEstBySc_los.P_list(5,5,end)^0.5;
result(2).velErrorSc(6,n) = scEstBySc_los.P_list(6,6,end)^0.5;
result(2).velErrorSc(7,n) = scEstBySc_los.P_list(7,7,end)^0.5;
result(2).velErrorSc(8,n) = (result(2).velErrorSc(5,n)^2 + result(2).velErrorSc(6,n)^2 + result(2).velErrorSc(7,n)^2)^0.5;
result(2).clockErrorSc(1,n) = abs(scEstBySc_los.clockError(end));
result(2).clockErrorSc(2,n) = scEstBySc_los.P_list(1,1,end)^0.5;

result(2).posErrorGs(1,n) = abs(scEstByGs_los.state(1,end) - scTrue_los.state(1,end));
result(2).posErrorGs(2,n) = abs(scEstByGs_los.state(2,end) - scTrue_los.state(2,end));
result(2).posErrorGs(3,n) = abs(scEstByGs_los.state(3,end) - scTrue_los.state(3,end));
result(2).posErrorGs(4,n) = (result(2).posErrorGs(1,n)^2 + result(2).posErrorGs(2,n)^2 + result(2).posErrorGs(3,n)^2)^0.5;
result(2).posErrorGs(5,n) = scEstByGs_los.P_list(2,2,end)^0.5;
result(2).posErrorGs(6,n) = scEstByGs_los.P_list(3,3,end)^0.5;
result(2).posErrorGs(7,n) = scEstByGs_los.P_list(4,4,end)^0.5;
result(2).posErrorGs(8,n) = (result(2).posErrorGs(5,n)^2 + result(2).posErrorGs(6,n)^2 + result(2).posErrorGs(7,n)^2)^0.5;
result(2).velErrorGs(1,n) = abs(scEstByGs_los.state(4,end) - scTrue_los.state(4,end));
result(2).velErrorGs(2,n) = abs(scEstByGs_los.state(5,end) - scTrue_los.state(5,end));
result(2).velErrorGs(3,n) = abs(scEstByGs_los.state(6,end) - scTrue_los.state(6,end));
result(2).velErrorGs(4,n) = (result(2).velErrorGs(1,n)^2 + result(2).velErrorGs(2,n)^2 + result(2).velErrorGs(3,n)^2)^0.5;
result(2).velErrorGs(5,n) = scEstByGs_los.P_list(5,5,end)^0.5;
result(2).velErrorGs(6,n) = scEstByGs_los.P_list(6,6,end)^0.5;
result(2).velErrorGs(7,n) = scEstByGs_los.P_list(7,7,end)^0.5;
result(2).velErrorGs(8,n) = (result(2).velErrorGs(5,n)^2 + result(2).velErrorGs(6,n)^2 + result(2).velErrorGs(7,n)^2)^0.5;
result(2).clockErrorGs(1,n) = abs(scEstByGs_los.clockError(end));
result(2).clockErrorGs(2,n) = scEstByGs_los.P_list(1,1,end)^0.5;

result(2).downAvail(1,n) = 0;
% 長期通信断絶後
result(3).posErrorSc(1,n) = abs(scEstByScEkf_afterLos.state(1,end) - scTrue_afterLos.state(1,end));
result(3).posErrorSc(2,n) = abs(scEstByScEkf_afterLos.state(2,end) - scTrue_afterLos.state(2,end));
result(3).posErrorSc(3,n) = abs(scEstByScEkf_afterLos.state(3,end) - scTrue_afterLos.state(3,end));
result(3).posErrorSc(4,n) = (result(3).posErrorSc(1,n)^2 + result(3).posErrorSc(2,n)^2 + result(3).posErrorSc(3,n)^2)^0.5;
result(3).posErrorSc(5,n) = scEstByScEkf_afterLos.P_list(2,2,end)^0.5;
result(3).posErrorSc(6,n) = scEstByScEkf_afterLos.P_list(3,3,end)^0.5;
result(3).posErrorSc(7,n) = scEstByScEkf_afterLos.P_list(4,4,end)^0.5;
result(3).posErrorSc(8,n) = (result(3).posErrorSc(5,n)^2 + result(3).posErrorSc(6,n)^2 + result(3).posErrorSc(7,n)^2)^0.5;
result(3).velErrorSc(1,n) = abs(scEstByScEkf_afterLos.state(4,end) - scTrue_afterLos.state(4,end));
result(3).velErrorSc(2,n) = abs(scEstByScEkf_afterLos.state(5,end) - scTrue_afterLos.state(5,end));
result(3).velErrorSc(3,n) = abs(scEstByScEkf_afterLos.state(6,end) - scTrue_afterLos.state(6,end));
result(3).velErrorSc(4,n) = (result(3).velErrorSc(1,n)^2 + result(3).velErrorSc(2,n)^2 + result(3).velErrorSc(3,n)^2)^0.5;
result(3).velErrorSc(5,n) = scEstByScEkf_afterLos.P_list(5,5,end)^0.5;
result(3).velErrorSc(6,n) = scEstByScEkf_afterLos.P_list(6,6,end)^0.5;
result(3).velErrorSc(7,n) = scEstByScEkf_afterLos.P_list(7,7,end)^0.5;
result(3).velErrorSc(8,n) = (result(3).velErrorSc(5,n)^2 + result(3).velErrorSc(6,n)^2 + result(3).velErrorSc(7,n)^2)^0.5;
result(3).clockErrorSc(1,n) = abs(scEstByScEkf_afterLos.clockError(end));
result(3).clockErrorSc(2,n) = scEstByScEkf_afterLos.P_list(1,1,end)^0.5;

result(3).posErrorGs(1,n) = abs(scEstByGsEkf_afterLos.state(1,end) - scTrue_afterLos.state(1,end));
result(3).posErrorGs(2,n) = abs(scEstByGsEkf_afterLos.state(2,end) - scTrue_afterLos.state(2,end));
result(3).posErrorGs(3,n) = abs(scEstByGsEkf_afterLos.state(3,end) - scTrue_afterLos.state(3,end));
result(3).posErrorGs(4,n) = (result(3).posErrorGs(1,n)^2 + result(3).posErrorGs(2,n)^2 + result(3).posErrorGs(3,n)^2)^0.5;
result(3).posErrorGs(5,n) = scEstByGsEkf_afterLos.P_list(2,2,end)^0.5;
result(3).posErrorGs(6,n) = scEstByGsEkf_afterLos.P_list(3,3,end)^0.5;
result(3).posErrorGs(7,n) = scEstByGsEkf_afterLos.P_list(4,4,end)^0.5;
result(3).posErrorGs(8,n) = (result(3).posErrorGs(5,n)^2 + result(3).posErrorGs(6,n)^2 + result(3).posErrorGs(7,n)^2)^0.5;
result(3).velErrorGs(1,n) = abs(scEstByGsEkf_afterLos.state(4,end) - scTrue_afterLos.state(4,end));
result(3).velErrorGs(2,n) = abs(scEstByGsEkf_afterLos.state(5,end) - scTrue_afterLos.state(5,end));
result(3).velErrorGs(3,n) = abs(scEstByGsEkf_afterLos.state(6,end) - scTrue_afterLos.state(6,end));
result(3).velErrorGs(4,n) = (result(3).velErrorGs(1,n)^2 + result(3).velErrorGs(2,n)^2 + result(3).velErrorGs(3,n)^2)^0.5;
result(3).velErrorGs(5,n) = scEstByGsEkf_afterLos.P_list(5,5,end)^0.5;
result(3).velErrorGs(6,n) = scEstByGsEkf_afterLos.P_list(6,6,end)^0.5;
result(3).velErrorGs(7,n) = scEstByGsEkf_afterLos.P_list(7,7,end)^0.5;
result(3).velErrorGs(8,n) = (result(3).velErrorGs(5,n)^2 + result(3).velErrorGs(6,n)^2 + result(3).velErrorGs(7,n)^2)^0.5;
result(3).clockErrorGs(1,n) = abs(scEstByGsEkf_afterLos.clockError(end));
result(3).clockErrorGs(2,n) = scEstByGsEkf_afterLos.P_list(1,1,end)^0.5;

result(3).downAvail(1,n) = sum(gsTrue_afterLos.snr_drList > gs.reqSnr_down)/length(gsTrue_afterLos.snr_drList) * 100;


