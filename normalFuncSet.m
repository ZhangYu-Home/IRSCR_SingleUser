%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 相当于C语言中的函数声明，在此处进行声明的函数，均是需在外部使用的函数；
%% 若函数只在本文件其他函数内调用，可以不在此处声明，外部也就无法调用，相当
%% 于私有函数；
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Func = normalFuncSet
    %需要对外展示的接口
    Func.init = @init;
    Func.setChannel = @setChannel;
    Func.getUnionChannel = @getUnionChannel;
    Func.initPrecodeAndReflectMat = @initPrecodeAndReflectMat;
    Func.getWeightSumRate = @getWeightSumRate;
    Func.getSigAndJamMat = @getSigAndJamMat;
    Func.getDecodeAndWeightMat = @getDecodeAndWeightMat;
    Func.getPrecodeMat = @getPrecodeMat;
    Func.getReflectMat = @getReflectMat;
end

%% 初始化函数
function [scene,dist] = init()
    scene = initScene();
    pos = initNodePos(scene.n_SU);%pos作为中间变量无需返回
    dist = initNodeDist(pos);
end

%% 初始化只包含单个主用户时的场景参数
function scene = initScene()
    %scene中用于保存需要在后续函数中传递的参数
    scene.n_SU = 5;%次级用户数
    scene.n_ante_AP = 6;%次级接入点天线数
    scene.n_ante_PU = 7;%主用户天线数
    scene.n_ante_SU = 8;%次级用户天线数
    scene.m_IRS = 5;%IRS反射单元个数
    scene.noise_SU = 1e-11;%次级用户处噪声平均功率
    scene.max_pow = 0.2;%定义最大功率
    scene.leak_pow = 1e-6;%定义主用户的干扰泄漏约束阈值
    scene.n_data = min(scene.n_ante_AP,scene.n_ante_SU);%发射的数据流数量
end

%% 初始化各个节点位置参数
function pos = initNodePos(n_SU)
    %各个节点的位置
    pos.AP.x = 0;%次级接入点位置
    pos.AP.y = 0;
    pos.PU.x = 50;%单个主用户的位置
    pos.PU.y = 0;
    pos.SUs.x = -5 + 10*rand(1,n_SU);%各个次级用户的位置
    pos.SUs.y = 20 + 10*rand(1,n_SU);
    pos.IRS.x = 40;%IRS的位置
    pos.IRS.y = 10;
end

%% 初始化各个节点间距离参数
function dist = initNodeDist(pos)
    dist.AP_PU = sqrt((pos.AP.x-pos.PU.x)^2+(pos.AP.y-pos.PU.y)^2);
    dist.AP_SUs = sqrt((pos.AP.x-pos.SUs.x).^2+(pos.AP.y-pos.SUs.y).^2);
    dist.AP_IRS = sqrt((pos.AP.x-pos.IRS.x)^2+(pos.AP.y-pos.IRS.y)^2);
    dist.IRS_PU = sqrt((pos.IRS.x-pos.PU.x)^2+(pos.IRS.y-pos.PU.y)^2);
    dist.IRS_SUs = sqrt((pos.IRS.x-pos.SUs.x).^2+(pos.IRS.y-pos.SUs.y).^2);
end

%% 初始化各个节点间的路损
function pathloss = initPathLoss(dist)
    %路损指数：
    pl_exp = 2;
    %路损求解
    pathloss.AP_PU = 10.^((-30-10*pl_exp*log10(dist.AP_PU))/10).^0.5;
    pathloss.AP_SUs = 10.^((-30-10*pl_exp*log10(dist.AP_SUs))/10).^0.5;
    pathloss.AP_IRS = 10.^((-30-10*pl_exp*log10(dist.AP_IRS))/10).^0.5;
    pathloss.IRS_PU = 10.^((-30-10*pl_exp*log10(dist.IRS_PU))/10).^0.5;
    pathloss.IRS_SUs = 10.^((-30-10*pl_exp*log10(dist.IRS_SUs))/10).^0.5;
    
end

%% 设置各个节点间的信道
%% 当前所使用信道均为瑞丽衰落信道，后续可以进行修改
function channel = setChannel(scene, dist)
    pathloss = initPathLoss(dist);
    h_AP_PU = randn(scene.n_ante_PU,scene.n_ante_AP)+1j*randn(scene.n_ante_PU,scene.n_ante_AP); 
    h_AP_SUs = randn(scene.n_ante_SU,scene.n_ante_AP,scene.n_SU)+1j*randn(scene.n_ante_SU,scene.n_ante_AP,scene.n_SU); 
    h_AP_IRS = randn(scene.m_IRS,scene.n_ante_AP)+1j*randn(scene.m_IRS,scene.n_ante_AP); 
    h_IRS_PU = randn(scene.n_ante_PU,scene.m_IRS)+1j*randn(scene.n_ante_PU,scene.m_IRS); 
    h_IRS_SUs = randn(scene.n_ante_SU,scene.m_IRS,scene.n_SU)+1j*randn(scene.n_ante_SU,scene.m_IRS,scene.n_SU); 

    %对与次级用户相关的信道单独处理
    for i = 1:scene.n_SU
        h_AP_SUs(:,:,i) = h_AP_SUs(:,:,i)*pathloss.AP_SUs(i);
        h_IRS_SUs(:,:,i) = h_IRS_SUs(:,:,i)*pathloss.IRS_SUs(i);
    end

    channel.h_AP_PU = h_AP_PU*pathloss.AP_PU;
    channel.h_AP_SUs = h_AP_SUs;
    channel.h_AP_IRS = h_AP_IRS*pathloss.AP_IRS;
    channel.h_IRS_PU = h_IRS_PU*pathloss.IRS_PU;
    channel.h_IRS_SUs = h_IRS_SUs;
 
end

%% 计算联合信道
function [g_AP_PU,g_AP_SUs] = getUnionChannel(channel,reflect_mat)
    g_AP_PU = channel.h_AP_PU+channel.h_IRS_PU*reflect_mat*channel.h_AP_IRS;
    g_AP_SUs = zeros(size(channel.h_AP_SUs));
    n_SU = size(channel.h_AP_SUs,3);
    for i = 1:n_SU
        g_AP_SUs(:,:,i) = channel.h_AP_SUs(:,:,i)+channel.h_IRS_SUs(:,:,i)*reflect_mat*channel.h_AP_IRS;
    end
end

%% 初始化预编码矩阵和反射系数矩阵
function [precode_mat,reflect_mat] = initPrecodeAndReflectMat(scene)
    %初始化预编码矩阵，按照等功率分配至每个数据流
    precode_mat = zeros(scene.n_ante_AP,scene.n_data,scene.n_SU);
    tmp_coeff = sqrt(scene.max_pow/scene.n_SU/scene.n_data);
    for i = 1:scene.n_SU
        precode_mat(:,:,i) = tmp_coeff*eye(scene.n_ante_AP,scene.n_data);
    end
    %初始化反射系数矩阵，令所有反射单元相移为0；
    reflect_mat = eye(scene.n_SU);
end

%% 计算信号功率矩阵与干扰协方差矩阵
function [sig_mat,jam_mat] = getSigAndJamMat(g_AP_SUs,precode_mat,noise_SU)
    [n_ante_SU,n_ante_AP,n_SU] = size(g_AP_SUs);
    I_SU = eye(n_ante_SU);%单位矩阵
    
    %临时变量：用于存储所有预编码矩阵与其共轭转置乘积的和
    Q = zeros(n_ante_AP,n_ante_AP);
    for i = 1:n_SU
        Q = Q + precode_mat(:,:,i)*precode_mat(:,:,i)';
    end
    
    sig_mat = zeros(n_ante_SU,n_ante_SU,n_SU);
    jam_mat = zeros(n_ante_SU,n_ante_SU,n_SU);
    for i = 1:n_SU
        tmp_mat = g_AP_SUs(:,:,i)*precode_mat(:,:,i);
        sig_mat(:,:,i) = tmp_mat*tmp_mat';
        jam_mat(:,:,i) = g_AP_SUs(:,:,i)*Q*g_AP_SUs(:,:,i)'+noise_SU*I_SU-sig_mat(:,:,i);
    end
end

%% 计算所有次级用户加权和速率
function sum_rate = getWeightSumRate(sig_mat,jam_mat)
    [n_ante_SU,~,n_SU] = size(sig_mat);
    I_SU = eye(n_ante_SU);%单位矩阵，用于求解速率
    
    sum_rate = 0;
    for i = 1:n_SU
        %此处是否要加入与最大功率max_pow相乘，需要后续考虑
        sum_rate = sum_rate + real(log(det(I_SU+sig_mat(:,:,i)*inv(jam_mat(:,:,i)))));
    end
end

