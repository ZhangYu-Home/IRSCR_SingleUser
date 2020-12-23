%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% �൱��C�����еĺ����������ڴ˴����������ĺ��������������ⲿʹ�õĺ�����
%% ������ֻ�ڱ��ļ����������ڵ��ã����Բ��ڴ˴��������ⲿҲ���޷����ã��൱
%% ��˽�к�����
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Func = normalFuncSet
    %��Ҫ����չʾ�Ľӿ�
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

%% ��ʼ������
function [scene,dist] = init()
    scene = initScene();
    pos = initNodePos(scene.n_SU);%pos��Ϊ�м�������践��
    dist = initNodeDist(pos);
end

%% ��ʼ��ֻ�����������û�ʱ�ĳ�������
function scene = initScene()
    %scene�����ڱ�����Ҫ�ں��������д��ݵĲ���
    scene.n_SU = 5;%�μ��û���
    scene.n_ante_AP = 6;%�μ������������
    scene.n_ante_PU = 7;%���û�������
    scene.n_ante_SU = 8;%�μ��û�������
    scene.m_IRS = 5;%IRS���䵥Ԫ����
    scene.noise_SU = 1e-11;%�μ��û�������ƽ������
    scene.max_pow = 0.2;%���������
    scene.leak_pow = 1e-6;%�������û��ĸ���й©Լ����ֵ
    scene.n_data = min(scene.n_ante_AP,scene.n_ante_SU);%���������������
end

%% ��ʼ�������ڵ�λ�ò���
function pos = initNodePos(n_SU)
    %�����ڵ��λ��
    pos.AP.x = 0;%�μ������λ��
    pos.AP.y = 0;
    pos.PU.x = 50;%�������û���λ��
    pos.PU.y = 0;
    pos.SUs.x = -5 + 10*rand(1,n_SU);%�����μ��û���λ��
    pos.SUs.y = 20 + 10*rand(1,n_SU);
    pos.IRS.x = 40;%IRS��λ��
    pos.IRS.y = 10;
end

%% ��ʼ�������ڵ��������
function dist = initNodeDist(pos)
    dist.AP_PU = sqrt((pos.AP.x-pos.PU.x)^2+(pos.AP.y-pos.PU.y)^2);
    dist.AP_SUs = sqrt((pos.AP.x-pos.SUs.x).^2+(pos.AP.y-pos.SUs.y).^2);
    dist.AP_IRS = sqrt((pos.AP.x-pos.IRS.x)^2+(pos.AP.y-pos.IRS.y)^2);
    dist.IRS_PU = sqrt((pos.IRS.x-pos.PU.x)^2+(pos.IRS.y-pos.PU.y)^2);
    dist.IRS_SUs = sqrt((pos.IRS.x-pos.SUs.x).^2+(pos.IRS.y-pos.SUs.y).^2);
end

%% ��ʼ�������ڵ���·��
function pathloss = initPathLoss(dist)
    %·��ָ����
    pl_exp = 2;
    %·�����
    pathloss.AP_PU = 10.^((-30-10*pl_exp*log10(dist.AP_PU))/10).^0.5;
    pathloss.AP_SUs = 10.^((-30-10*pl_exp*log10(dist.AP_SUs))/10).^0.5;
    pathloss.AP_IRS = 10.^((-30-10*pl_exp*log10(dist.AP_IRS))/10).^0.5;
    pathloss.IRS_PU = 10.^((-30-10*pl_exp*log10(dist.IRS_PU))/10).^0.5;
    pathloss.IRS_SUs = 10.^((-30-10*pl_exp*log10(dist.IRS_SUs))/10).^0.5;
    
end

%% ���ø����ڵ����ŵ�
%% ��ǰ��ʹ���ŵ���Ϊ����˥���ŵ����������Խ����޸�
function channel = setChannel(scene, dist)
    pathloss = initPathLoss(dist);
    h_AP_PU = randn(scene.n_ante_PU,scene.n_ante_AP)+1j*randn(scene.n_ante_PU,scene.n_ante_AP); 
    h_AP_SUs = randn(scene.n_ante_SU,scene.n_ante_AP,scene.n_SU)+1j*randn(scene.n_ante_SU,scene.n_ante_AP,scene.n_SU); 
    h_AP_IRS = randn(scene.m_IRS,scene.n_ante_AP)+1j*randn(scene.m_IRS,scene.n_ante_AP); 
    h_IRS_PU = randn(scene.n_ante_PU,scene.m_IRS)+1j*randn(scene.n_ante_PU,scene.m_IRS); 
    h_IRS_SUs = randn(scene.n_ante_SU,scene.m_IRS,scene.n_SU)+1j*randn(scene.n_ante_SU,scene.m_IRS,scene.n_SU); 

    %����μ��û���ص��ŵ���������
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

%% ���������ŵ�
function [g_AP_PU,g_AP_SUs] = getUnionChannel(channel,reflect_mat)
    g_AP_PU = channel.h_AP_PU+channel.h_IRS_PU*reflect_mat*channel.h_AP_IRS;
    g_AP_SUs = zeros(size(channel.h_AP_SUs));
    n_SU = size(channel.h_AP_SUs,3);
    for i = 1:n_SU
        g_AP_SUs(:,:,i) = channel.h_AP_SUs(:,:,i)+channel.h_IRS_SUs(:,:,i)*reflect_mat*channel.h_AP_IRS;
    end
end

%% ��ʼ��Ԥ�������ͷ���ϵ������
function [precode_mat,reflect_mat] = initPrecodeAndReflectMat(scene)
    %��ʼ��Ԥ������󣬰��յȹ��ʷ�����ÿ��������
    precode_mat = zeros(scene.n_ante_AP,scene.n_data,scene.n_SU);
    tmp_coeff = sqrt(scene.max_pow/scene.n_SU/scene.n_data);
    for i = 1:scene.n_SU
        precode_mat(:,:,i) = tmp_coeff*eye(scene.n_ante_AP,scene.n_data);
    end
    %��ʼ������ϵ�����������з��䵥Ԫ����Ϊ0��
    reflect_mat = eye(scene.n_SU);
end

%% �����źŹ��ʾ��������Э�������
function [sig_mat,jam_mat] = getSigAndJamMat(g_AP_SUs,precode_mat,noise_SU)
    [n_ante_SU,n_ante_AP,n_SU] = size(g_AP_SUs);
    I_SU = eye(n_ante_SU);%��λ����
    
    %��ʱ���������ڴ洢����Ԥ����������乲��ת�ó˻��ĺ�
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

%% �������дμ��û���Ȩ������
function sum_rate = getWeightSumRate(sig_mat,jam_mat)
    [n_ante_SU,~,n_SU] = size(sig_mat);
    I_SU = eye(n_ante_SU);%��λ���������������
    
    sum_rate = 0;
    for i = 1:n_SU
        %�˴��Ƿ�Ҫ�����������max_pow��ˣ���Ҫ��������
        sum_rate = sum_rate + real(log(det(I_SU+sig_mat(:,:,i)*inv(jam_mat(:,:,i)))));
    end
end

