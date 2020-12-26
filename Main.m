%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%此文件用于仿真单个主用户的情况；仿真图的横纵坐标分别为：最大可用发射功率和加权
%和速率；在仿真时，令各个次级用户权重均为1；
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clear;
clc;

%% 初始化参数
func = normalFuncSet; %导入函数集
n_monte = 100; %蒙特卡洛仿真次数
n_pow = 10; %功率迭代的次数     
[scene,dist] = func.init();%初始化场景参数
rate_vec = zeros(1,n_pow);
rate_vec1 = zeros(1,n_pow);

%% 进行蒙特卡洛仿真
tic
for cnt_monte = 1:n_monte
    disp(['Iteration = ', num2str(cnt_monte)]);
    channel = func.setChannel(scene, dist);
    for cnt_pow = 1:n_pow
        scene.max_pow = cnt_pow;
        disp(['    The max Power is ', num2str(scene.max_pow),'W.']);
        %% 初始化过程
        %初始化反射系数矩阵和预编码矩阵
        [precode_mat,reflect_mat] = func.initPrecodeAndReflectMat(scene);       
        %计算联合信道
        [g_AP_PU,g_AP_SUs] = func.getUnionChannel(channel,reflect_mat);
        %计算每个次级用户对应的功率矩阵和干扰协方差矩阵
        [sig_mat,jam_mat] = func.getSigAndJamMat(g_AP_SUs,precode_mat,scene.noise_SU);
        %计算所有次级用户的速率和
        sum_rate = func.getWeightSumRate(sig_mat,jam_mat);
        disp(['The front sum of rates is ',num2str(sum_rate),' bps.']);
        
        %% 基于交替优化求解问题
        for cnt_iter = 1:1000
            sum_rate_tmp = sum_rate;
            %计算解码矩阵和辅助矩阵
            [decode_mat,weight_mat] = getDecodeAndWeightMat(sig_mat,jam_mat,g_AP_SUs,precode_mat);
            %给定其他参数的情况下计算预编码矩阵
            precode_mat = getPrecodeMat(scene,g_AP_PU,g_AP_SUs,decode_mat,weight_mat,precode_mat);
            if(cnt_iter == 1)
                %计算每个次级用户对应的功率矩阵和干扰协方差矩阵
                [sig_mat,jam_mat] = func.getSigAndJamMat(g_AP_SUs,precode_mat,scene.noise_SU);
                %计算所有次级用户的速率和
                sum_rate = func.getWeightSumRate(sig_mat,jam_mat);
                rate_vec1(cnt_pow) = rate_vec1(cnt_pow) + sum_rate;
            end
            %给定其他参数的情况下计算反射系数矩阵
            reflect_mat = getReflectMat(scene,channel,precode_mat,reflect_mat,decode_mat,weight_mat);          
            %计算联合信道
            [g_AP_PU,g_AP_SUs] = func.getUnionChannel(channel,reflect_mat);
            %计算每个次级用户对应的功率矩阵和干扰协方差矩阵
            [sig_mat,jam_mat] = func.getSigAndJamMat(g_AP_SUs,precode_mat,scene.noise_SU);
            %计算所有次级用户的速率和
            sum_rate = func.getWeightSumRate(sig_mat,jam_mat);
            %disp(['The back sum of rates is ',num2str(sum_rate),' bps.']);
            if(abs(sum_rate - sum_rate_tmp) < 0.0001)
                break;
            end
        end
        rate_vec(cnt_pow) = rate_vec(cnt_pow) + sum_rate;
        disp(['The back sum of rates is ',num2str(sum_rate),' bps.']);
    end
end
rate_vec = rate_vec/n_monte;
rate_vec1 = rate_vec1/n_monte;
toc