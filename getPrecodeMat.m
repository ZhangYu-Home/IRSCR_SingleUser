%% 计算预编码矩阵
function precode_mat = getPrecodeMat(scene,g_AP_PU,g_AP_SUs,decode_mat,weight_mat,precode_mat)   
    %% 引入函数集
    func = normalFuncSet;
    %% 计算相关参数
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
        %precode_mat_tmp为连续凸近似每一轮迭代的初始点,leak_pow_tmp为相应的近似干扰泄漏阈值
        precode_mat_tmp = precode_mat;
        leak_pow_tmp = scene.leak_pow;
        for i = 1:scene.n_SU
            leak_pow_tmp = leak_pow_tmp - real(trace(precode_mat_tmp(:,:,i)'*(Z_p-X_p)*precode_mat_tmp(:,:,i)));
        end
        
       %% 情形一：功率约束为无效约束
        mu = 0;
        precode_mat = calPrecodeMat1(X_0,Z_p,mu,Y,X_p,precode_mat_tmp);
        val_J = calInterferenceLeak1(Z_p,X_p,precode_mat,precode_mat_tmp);
        if(val_J > leak_pow_tmp)         
            %利用指数步进确定mu的上界和下界
            mu = 0.1;
            while(1)
                precode_mat = calPrecodeMat1(X_0,Z_p,mu,Y,X_p,precode_mat_tmp);
                val_J = calInterferenceLeak1(Z_p,X_p,precode_mat,precode_mat_tmp);
                if(val_J > leak_pow_tmp) mu = mu * 2; else break; end
            end
            mu_low = 0;
            mu_up = mu;
            
            %基于求得的上界和下界，利用二分搜索求解上界和下界
            while(1)
                mu = (mu_low + mu_up)/2;
                precode_mat = calPrecodeMat1(X_0,Z_p,mu,Y,X_p,precode_mat_tmp);
                val_J = calInterferenceLeak1(Z_p,X_p,precode_mat,precode_mat_tmp);
                if(val_J > leak_pow_tmp) mu_low = mu; else mu_up = mu; end
                if(mu_up - mu_low < 0.001)
                    %取mu的下界，保证约束条件满足
                    precode_mat = calPrecodeMat1(X_0,Z_p,mu_up,Y,X_p,precode_mat_tmp);
                    break;
                end
            end
        end
        
        %% 判断是否满足功率约束，如果不满足，考虑情形二
        if(calTotalPower(precode_mat) > scene.max_pow)
            %% 情形二：功率约束为有效约束
            leak_pow_tmp = max(eig(X_p))*scene.max_pow - leak_pow_tmp;
            lambda = 0;
            mu = 0;
            precode_mat = calPrecodeMat2(X_0,Z_p,lambda,mu,Y,X_p,precode_mat_tmp);
            %如果两个约束条件有一个不满足，则继续计算
            if(calTotalPower(precode_mat) > scene.max_pow || calInterferenceLeak2(Z_p,X_p,precode_mat,precode_mat_tmp) < leak_pow_tmp)
                mu = calMu2(leak_pow_tmp,X_0,Z_p,X_p,Y,precode_mat,precode_mat_tmp,lambda);
                precode_mat = calPrecodeMat2(X_0,Z_p,lambda,mu,Y,X_p,precode_mat_tmp);
                %如果功率约束条件不满足，则继续计算
                if(calTotalPower(precode_mat) > scene.max_pow)
                    %利用指数步进寻找lambda的上界
                    lambda = 0.1;
                    while(1)
                        mu = calMu2(leak_pow_tmp,X_0,Z_p,X_p,Y,precode_mat,precode_mat_tmp,lambda);
                        precode_mat = calPrecodeMat2(X_0,Z_p,lambda,mu,Y,X_p,precode_mat_tmp);
                        if(calTotalPower(precode_mat) > scene.max_pow)
                            lambda = lambda*2;
                        else
                            break;
                        end
                    end
                    lambda_low = 0;
                    lambda_up = lambda;
                    
                    %利用二分搜索求解最优lambda
                    while(1)
                        lambda = (lambda_up + lambda_low)/2;
                        mu = calMu2(leak_pow_tmp,X_0,Z_p,X_p,Y,precode_mat,precode_mat_tmp,lambda);
                        precode_mat = calPrecodeMat2(X_0,Z_p,lambda,mu,Y,X_p,precode_mat_tmp);
                        val_P = calTotalPower(precode_mat);
                        if(val_P > scene.max_pow) 
                            lambda_up = lambda; 
                        else
                            lambda_low = lambda; 
                        end
                        if(lambda_up-lambda_low < 0.001)
                            lambda = lambda_up;
                            mu = calMu2(leak_pow_tmp,X_0,Z_p,X_p,Y,precode_mat,precode_mat_tmp,lambda);
                            precode_mat = calPrecodeMat2(X_0,Z_p,lambda,mu,Y,X_p,precode_mat_tmp);
                            break;
                        end
                    end
                end 
            end
        end

        
        
        %% 判断是否跳出
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

%% 在情形一中，当mu变化时，计算预编码矩阵时用到的函数
function precode_mat = calPrecodeMat1(X_0,Z_p,mu,Y,X_p,precode_mat)
    n_SU = size(precode_mat,3);
    mat_coeff1 = inv(X_0+mu*Z_p);
    mat_coeff2 = mu*(Z_p - X_p);
    for i = 1:n_SU
        precode_mat(:,:,i) = mat_coeff1*(Y(:,:,i)-mat_coeff2*precode_mat(:,:,i));
    end
end

%% 在情形一中，当mu变化时，计算干扰泄漏
function val_J = calInterferenceLeak1(Z_p,X_p,precode_mat,precode_mat_tmp)
    n_SU = size(precode_mat,3);
    val_J = 0;
    for i = 1:n_SU
        val_J = val_J + real(trace(precode_mat(:,:,i)'*Z_p*precode_mat(:,:,i)));
        val_J = val_J - 2*real(trace(precode_mat_tmp(:,:,i)'*(Z_p-X_p)*precode_mat(:,:,i)));
    end
end

%% 判断是否满足功率约束
function val_P = calTotalPower(precode_mat)
    val_P = 0;
    for i = 1:size(precode_mat,3)
        val_P = val_P + trace(precode_mat(:,:,i)'*precode_mat(:,:,i));
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%情形二相关函数%%%%%%%%%%%%%%%%%%%%%% 
%% 在情形二中，当mu和lambda变化时，计算预编码矩阵时用到的函数
function precode_mat = calPrecodeMat2(X_0,Z_p,lambda,mu,Y,X_p,precode_mat)
    tmp_coeff1 = inv(X_0+lambda*eye(size(X_0,1)));
    tmp_coeff2 = mu*(Z_p - X_p);
    for i = 1:size(precode_mat,3)
        precode_mat(:,:,i) = tmp_coeff1*(Y(:,:,i)-tmp_coeff2*precode_mat(:,:,i));
    end
end

%% 在情形一中，当mu变化时，计算干扰泄漏约束
function val_J = calInterferenceLeak2(Z_p,X_p,precode_mat,precode_mat_tmp)
    val_J = 0;
    for i = 1:size(precode_mat,3)
        val_J = 2*real(trace(precode_mat_tmp(:,:,i)'*(Z_p-X_p)*precode_mat(:,:,i)));
    end
end
%% 计算mu值
function mu = calMu2(leak_pow_tmp,X_0,Z_p,X_p,Y,precode_mat,precode_mat_tmp,lambda)
    n_SU = size(precode_mat,3);
    mat_coe1 = (Z_p-X_p)*inv(X_0+lambda*eye(size(X_0,1)));
    mat_coe2 = mat_coe1*(Z_p-X_p);
    tmp_coeff1 = leak_pow_tmp;
    for i = 1:n_SU
        tmp_coeff1 = tmp_coeff1 - 2*real(trace(precode_mat_tmp(:,:,i)'*mat_coe1*Y(:,:,i)'));
    end
    tmp_coeff2 = 0;
    for i = 1:n_SU
        tmp_coeff2 = tmp_coeff2 + 2*real(trace(precode_mat_tmp(:,:,i)'*mat_coe2*precode_mat_tmp(:,:,i)));
    end

    mu = tmp_coeff1/tmp_coeff2;
end
