%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%���ļ����ڷ��浥�����û������������ͼ�ĺ�������ֱ�Ϊ�������÷��书�ʺͼ�Ȩ
%�����ʣ��ڷ���ʱ��������μ��û�Ȩ�ؾ�Ϊ1��
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clear;
clc;

%% ��ʼ������
func = normalFuncSet; %���뺯����
n_monte = 100; %���ؿ���������
n_pow = 10; %���ʵ����Ĵ���     
[scene,dist] = func.init();%��ʼ����������
rate_vec = zeros(1,n_pow);
rate_vec1 = zeros(1,n_pow);

%% �������ؿ������
tic
for cnt_monte = 1:n_monte
    disp(['Iteration = ', num2str(cnt_monte)]);
    channel = func.setChannel(scene, dist);
    for cnt_pow = 1:n_pow
        scene.max_pow = cnt_pow;
        disp(['    The max Power is ', num2str(scene.max_pow),'W.']);
        %% ��ʼ������
        %��ʼ������ϵ�������Ԥ�������
        [precode_mat,reflect_mat] = func.initPrecodeAndReflectMat(scene);       
        %���������ŵ�
        [g_AP_PU,g_AP_SUs] = func.getUnionChannel(channel,reflect_mat);
        %����ÿ���μ��û���Ӧ�Ĺ��ʾ���͸���Э�������
        [sig_mat,jam_mat] = func.getSigAndJamMat(g_AP_SUs,precode_mat,scene.noise_SU);
        %�������дμ��û������ʺ�
        sum_rate = func.getWeightSumRate(sig_mat,jam_mat);
        disp(['The front sum of rates is ',num2str(sum_rate),' bps.']);
        
        %% ���ڽ����Ż��������
        for cnt_iter = 1:1000
            sum_rate_tmp = sum_rate;
            %����������͸�������
            [decode_mat,weight_mat] = getDecodeAndWeightMat(sig_mat,jam_mat,g_AP_SUs,precode_mat);
            %������������������¼���Ԥ�������
            precode_mat = getPrecodeMat(scene,g_AP_PU,g_AP_SUs,decode_mat,weight_mat,precode_mat);
            if(cnt_iter == 1)
                %����ÿ���μ��û���Ӧ�Ĺ��ʾ���͸���Э�������
                [sig_mat,jam_mat] = func.getSigAndJamMat(g_AP_SUs,precode_mat,scene.noise_SU);
                %�������дμ��û������ʺ�
                sum_rate = func.getWeightSumRate(sig_mat,jam_mat);
                rate_vec1(cnt_pow) = rate_vec1(cnt_pow) + sum_rate;
            end
            %������������������¼��㷴��ϵ������
            reflect_mat = getReflectMat(scene,channel,precode_mat,reflect_mat,decode_mat,weight_mat);          
            %���������ŵ�
            [g_AP_PU,g_AP_SUs] = func.getUnionChannel(channel,reflect_mat);
            %����ÿ���μ��û���Ӧ�Ĺ��ʾ���͸���Э�������
            [sig_mat,jam_mat] = func.getSigAndJamMat(g_AP_SUs,precode_mat,scene.noise_SU);
            %�������дμ��û������ʺ�
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