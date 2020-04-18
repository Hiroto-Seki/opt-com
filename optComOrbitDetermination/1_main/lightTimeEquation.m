% solve light time eqaution estimating the orbit of spacecraft, ground
% station

clear all; close all; clc

%% setting parameter
rng('default');
rng(1)
% �萔
constant.sunMu         = 1.32712438e11;   % km^3/s^2 ���z�̏d�͒萔
constant.earthMu       = 3.98600441e5;    % km^3/s^2 �n���̏d�͒萔
constant.saturnMu      = 3.7931187e7;     %[km/s^2]
constant.au            = 149597870.7;     % km ���z-�n���ԋ���
constant.lightSpeed    = 299792.458;      % [km/s] ���� 
constant.earthRadius   = 6.371e3;         % [km] �n�����a
constant.earthRotation = 2*pi/24/60/60;   % [rad/s] �n�����]���x
constant.eathAxis      = 23.4/180*pi;     % �n���̌X�� 

%% ���ԂɊւ���p�����[�^�[
% simulation timeStep[s]
time.simDt = 100;
% ������́Ctime.list�̒��ŉ��Ԗڂ�
time.refId = 600; 
% �����0�Ƃ��Ăǂ͈̔͂̎������v�Z���邩(-6000�b~6000�b���f�t�H���g�Ƃ���)
time.list = linspace((-time.refId+1)*time.simDt, (time.refId)*time.simDt,time.refId*2);

% ������̐ݒ�
time.t0 = juliandate(2030,1,1,0,0,0);

%% �덷�Ɋւ���p�����[�^�[
% �������v�덷
error.clock = 10   * randn; %���v�덷(�b)
error.time0 = error.clock/ (60*60*24); %���v�덷(��)
% �����T���@�O���덷[km]. �n������݂�1��rad���炢�̌덷�ɂ���
% error.scPos0 = randn(3,1) *  150;
error.scPos0 = randn(3,1) *  0;
% �K����0.1km/s���x�̌덷�Ƃ���
% error.scVel0 = randn(3,1) *  0.1;
error.scVel0 = randn(3,1) *  0;
%%  ��������l�̐ݒ�
% �n���̐��莞����ephemeris
earthEst = orbitalState('Earth','estimate','SolarSystem');
[earthEst.pos0,earthEst.vel0]=planetEphemeris(time.t0,'SolarSystem','Earth');
earthEst.pos0 =  reshape(earthEst.pos0,[3,1]);
earthEst.vel0 =  reshape(earthEst.vel0,[3,1]);
% �Ƃ肠�����C�T���@�O���̐���l = �y���O���̐���l�Ƃ���
scEst = orbitalState('spacecraft','estimate','SolarSystem');
[scEst.pos0,scEst.vel0]=planetEphemeris(time.t0,'SolarSystem','Saturn');
scEst.pos0 =  reshape(scEst.pos0,[3,1]);
scEst.vel0 =  reshape(scEst.vel0,[3,1]);

% �n���
gsEst = groundState('estimate');
% �n��ǂ̐ݒ�(ECI���W�n��)(�K��)
gsEst.pos0 =  [constant.earthRadius;0;0];

%% �^�l�̐ݒ�
%�������ephemeris�̐^�l
earthTrue = orbitalState('Earth','true','SolarSystem');
[earthTrue.pos0,earthTrue.vel0]=planetEphemeris(time.t0,'SolarSystem','Earth');
earthTrue.pos0 =  reshape(earthTrue.pos0,[3,1]);
earthTrue.vel0 =  reshape(earthTrue.vel0,[3,1]);
% [saturnPosTrue0,saturnVelTrue0]=planetEphemeris(time.t0+error.time0,'SolarSystem','Saturn');
% �T���@�̓y�����S�ł̈ʒu���x
scTrue = orbitalState('spacecraft','true','SolarSystem');
scTrue.pos0 = scEst.pos0 + error.scPos0;
scTrue.vel0 = scEst.vel0 + error.scVel0;

% �n��ǂ̐^�l
gsTrue = groundState('true');
[gsTrue.pos0,~] = groundState.earthRotation(gsEst.pos0, 0, constant);

%% �T���@�͋O���`���ɂ����time.list���̈ʒu�C���x�𓾂�
scEst.getOrbitTwoBody(time, constant);
scTrue.getOrbitTwoBody(time,constant);

%% �n���̈ʒu�C���x��`�����ċ��߂�
earthEst.getOrbitTwoBody(time, constant);
earthTrue.getOrbitTwoBody(time,constant);

%% �n��ǂ̈ʒu�C���x�����]�̉^��(t)�̊֐��ɂ���ċ��߂�
gsEst.getTrajectoryEarthRotation(time,constant)
gsTrue.getTrajectoryEarthRotation(time,constant);

%% �O���������̊֐��ɂ���(�V�~�����[�V������Ԃ�ʂ��Ă̊֐��Ƃ���)
% �񎚋Ȑ��ߎ� or r(t+dt)=r(t)+v(t)dt+0.5a(t)*dt^2�C Chebyshev coefficients�@�����Ƃ܂킵

%% �ϑ�����(�ϑ��͐^�l��`�����ċ��߂Ă���)
observe = observeState(time);
% �����\�z���v�덷
clockErrorCorrection = 0;
for i = 1:length(time.list)
    % �܂��͊ϑ��ʂ𓾂�
    observe.fromSc2Gs(time,earthTrue,gsTrue,scTrue,constant,error.clock,i);
    % �ϑ��ʂ��玞�v�덷�𐄑�����
    observe.calcClockOffset(time,earthEst,gsEst,scEst,constant,i,clockErrorCorrection)
    clockErrorCorrection = observe.clockError(i);
end

hold on
plot(time.list/60/60/24, observe.clockError)
plot(time.list/60/60/24, error.clock * ones(1, length(time.list)))
xlabel('year')
ylabel('light time delay [sec]')
legend('estimated','true')
hold off


% �n���ʒu��

%% �ϑ�����`���x�������߂�(�����ł͋O���͎��Ԃ̊֐��Ƃ������̂�p����)�����v�덷�����ς��肽��
% observe.calcClockOffset(time,earthEst,gsEst,scEst,constant)


% �����܂ł��������ɂ��

%% ���ς������`���x������C�w��������(�����O��������Ă���΂����Ǝw���ł��Ȃ����C�^�l�Ɛ���l�������v�Z���āC������������)


%%% �����܂ł�RG�܂łɂ�肽��
%% �����I�ɂ́C
% 1. �ϑ��𑝂₵�ċO���𐄒肷��

% 2. ���肵���O�������Ƃɒn��ǂɎw������
% �����ŁC�ǂ̂悤�ɒT�����邩����(�������ϒl�̕����Ɍ������C�T����)

% 3. �n��ǂŎ�M�ł���
% �ϑ����d�˂Ēn��ǂ��C�T���@�̋O���𐄒肷��D(�T���@�̎��v�̂���C�O��)
% �����ł́C���S��two way���m���ł����Ɖ��肵�C�T���@�̎��v�̂���͊��S�ɂ킩�����Ƃ���

% 4. �n��ǂ̃A�b�v�����N�ɂ���āC�T���@�͎��v(�n��ǂƈ�v)�ƋO����␳����D

% 5 �ŏI�I�Ȏw�����x�����ς���D



