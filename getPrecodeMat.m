%% ����Ԥ�������
function precode_mat = getPrecodeMat(scene,g_AP_PU,g_AP_SUs,decode_mat,weight_mat,precode_mat)   
    %% ���뺯����
    func = normalFuncSet;
    %% ������ز���
    X_0 = zeros(scene.n_ante_AP,scene.n_ante_AP);
    for i = 1:scene.n_SU
        X_0 = X_0 + g_AP_SUs(:,:,i)'*decode_mat(:,:,i)*weight_mat(:,:,i)*decode_mat(:,:,i)'*g_AP_SUs(:,:,i);
    end

    X_p = g_AP_PU'*g_AP_PU;
    Z_p = max(eig(X_p))*eye(scene.n_ante_AP);

    Y = zeros(scene.n_data,scene.n_ante_AP,scene.n_SU);
    for i = 1:scene.n_SU
        Y(:,:,i) = weight_mat(:,:,i)*decode_mat(:,:,i)'*g_AP_SUs(:,:,i);
    end

    while(1)
        %precode_mat_tmpΪ����͹����ÿһ�ֵ����ĳ�ʼ��,leak_pow_tmpΪ��Ӧ�Ľ��Ƹ���й©��ֵ
        precode_mat_tmp = precode_mat;
        leak_pow_tmp = scene.leak_pow;
        for i = 1:scene.n_SU
            leak_pow_tmp = leak_pow_tmp - real(trace(precode_mat_tmp(:,:,i)'*(Z_p-X_p)*precode_mat_tmp(:,:,i)));
        end
        
       %% ����һ������Լ��Ϊ��ЧԼ��
        mu = 0;
        precode_mat = calPrecodeMat1(X_0,Z_p,mu,Y,X_p,precode_mat_tmp);
        val_J = calInterferenceLeak(Z_p,X_p,precode_mat,precode_mat_tmp,leak_pow_tmp);
        if(val_J > 0)         
            %����ָ������ȷ��mu���Ͻ���½�
            mu = 0.1;
            while(1)
                precode_mat = calPrecodeMat1(X_0,Z_p,mu,Y,X_p,precode_mat_tmp);
                val_J = calInterferenceLeak(Z_p,X_p,precode_mat,precode_mat_tmp,scene.leak_pow);
                if(val_J > 0) mu = mu * 2; else break; end
            end
            mu_low = mu/2;
            mu_up = mu;
            
            %������õ��Ͻ���½磬���ö�����������Ͻ���½�
            while(1)
                mu = (mu_low + mu_up)/2;
                precode_mat = calPrecodeMat1(X_0,Z_p,mu,Y,X_p,precode_mat_tmp);
                val_J = calInterferenceLeak(Z_p,X_p,precode_mat,precode_mat_tmp,scene.leak_pow);
                if(val_J > 0) mu_low = mu; else mu_up = mu; end
                if(mu_up - mu_low < 0.001)
                    %ȡmu���½磬��֤Լ����������
                    precode_mat = calPrecodeMat1(X_0,Z_p,mu_low,Y,X_p,precode_mat_tmp);
                    break;
                end
            end
        end
        
        %% �ж��Ƿ����㹦��Լ����������㣬��һ����������������һ�֣����������ζ�
        if(isSatifyPowerCon(scene,precode_mat))
            continue;
        end

        %% ���ζ�������Լ��Ϊ��ЧԼ��
        leak_pow_tmp = max(eig(X_p))*scene.max_pow - leak_pow_tmp;
        lambda = 0;
        mu = 0;
        precode_mat = calPrecodeMat2(X_0,Z_p,lambda,mu,Y,X_p,precode_mat_tmp);
        if(isSatifyPowerCon(scene,precode_mat))
            
        end
        
        
        while(1)
            
        end 
        
        %% �ж��Ƿ�����
        [sig_mat,jam_mat] = func.getSigAndJamMat(g_AP_SUs,precode_mat_tmp,scene.noise_SU);
        sum_rate_tmp = func.getWeightSumRate(sig_mat,jam_mat);
        [sig_mat,jam_mat] = func.getSigAndJamMat(g_AP_SUs,precode_mat,scene.noise_SU);
        sum_rate = func.getWeightSumRate(sig_mat,jam_mat);
        
        if(abs(sum_rate-sum_rate_tmp) < 1e-6)
            disp(['mu_up = ',num2str(mu_up)]);
            disp(['mu_low = ',num2str(mu_low)]);
            break;
        end
    end
end

%% ������һ�У���mu�仯ʱ������Ԥ�������ʱ�õ��ĺ���
function precode_mat = calPrecodeMat1(X_0,Z_p,mu,Y,X_p,precode_mat)
    n_SU = size(precode_mat,3);
    tmp_coeff1 = inv(X_0+mu*Z_p);
    tmp_coeff2 = mu*(Z_p - X_p);
    for i = 1:n_SU
        precode_mat(:,:,i) = tmp_coeff1*(Y(:,:,i)-tmp_coeff2*precode_mat(:,:,i));
    end
end

%% ������һ�У���mu�仯ʱ���������й©Լ��
function val_J = calInterferenceLeak1(Z_p,X_p,precode_mat,precode_mat_tmp,leak_pow_tmp)
    n_SU = size(precode_mat,3);
    val_J = 0;
    for i = 1:n_SU
        val_J = val_J + real(trace(precode_mat(:,:,i)'*Z_p*precode_mat(:,:,i)));
        val_J = val_J - 2*real(trace(precode_mat_tmp(:,:,i)'*(Z_p-X_p)*precode_mat(:,:,i)));
    end
    val_J = val_J - leak_pow_tmp;
end

%% �ж��Ƿ����㹦��Լ��
function Flag = isSatifyPowerCon(scene,precode_mat)
    sum_pow = 0;
    for i = 1:scene.n_SU
        sum_pow = sum_pow + trace(precode_mat(:,:,i)'*precode_mat(:,:,i));
    end
    Flag = scene.max_pow - sum_pow >= 0;
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%���ζ���غ���%%%%%%%%%%%%%%%%%%%%%% 
%% �����ζ��У���mu��lambda�仯ʱ������Ԥ�������ʱ�õ��ĺ���
function precode_mat = calPrecodeMat2(X_0,Z_p,lambda,mu,Y,X_p,precode_mat)
    n_SU = size(precode_mat,3);
    tmp_coeff1 = inv(X_0+lambda*eye(n_SU));
    tmp_coeff2 = mu*(Z_p - X_p);
    for i = 1:n_SU
        precode_mat(:,:,i) = tmp_coeff1*(Y(:,:,i)-tmp_coeff2*precode_mat(:,:,i));
    end
end

%% ������һ�У���mu�仯ʱ���������й©Լ��
function val_J = calInterferenceLeak2(Z_p,X_p,precode_mat,precode_mat_tmp,leak_pow_tmp)
    n_SU = size(precode_mat,3);
    val_J = 0;
    for i = 1:n_SU
        val_J = 2*real(trace(precode_mat_tmp(:,:,i)'*(Z_p-X_p)*precode_mat(:,:,i)));
    end
    val_J = val_J - leak_pow_tmp;
end
%% ����muֵ
function mu = calMu2()
    mu = 0;
end
