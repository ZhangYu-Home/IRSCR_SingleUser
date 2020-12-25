%% 计算预编码矩阵
function precode_mat = getPrecodeMat(scene,g_AP_PU,g_AP_SUs,decode_mat,weight_mat,precode_mat)   
    %% 计算相关参数
    X_0 = zeros(scene.n_ante_AP,scene.n_ante_AP);
    for i = 1:scene.n_SU
        X_0 = X_0 + g_AP_SUs(:,:,i)'*decode_mat(:,:,i)*weight_mat(:,:,i)*decode_mat(:,:,i)'*g_AP_SUs(:,:,i);
    end
    X_p = g_AP_PU'*g_AP_PU;
    Z_p = max(real(eig(X_p)))*eye(scene.n_ante_AP);

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
            %利用指数步进确定mu的上界
            mu_low = 0;
            mu_up = 0.1;
            while(val_J > leak_pow_tmp)
                mu_up = mu_up * 2;
%                 val_J_tmp = val_J;
                precode_mat = calPrecodeMat1(X_0,Z_p,mu_up,Y,X_p,precode_mat_tmp);
                val_J = calInterferenceLeak1(Z_p,X_p,precode_mat,precode_mat_tmp);
%                 if(val_J > val_J_tmp)
%                      disp('出错');
%                      pause;
%                 end
            end
            
            %基于求得的上界和下界，利用二分搜索求解上界和下界
            while(mu_up - mu_low > 0.001)
                mu = (mu_low + mu_up)/2;
                precode_mat = calPrecodeMat1(X_0,Z_p,mu,Y,X_p,precode_mat_tmp);
                val_J = calInterferenceLeak1(Z_p,X_p,precode_mat,precode_mat_tmp);
                if(val_J > leak_pow_tmp) mu_low = mu; else mu_up = mu; end
            end
            %取mu的上界，保证约束条件满足
            precode_mat = calPrecodeMat1(X_0,Z_p,mu_up,Y,X_p,precode_mat_tmp);
        end
        %% 判断是否满足功率约束，如果不满足，考虑情形二
        if(calTotalPower(precode_mat) > scene.max_pow)
            %% 情形二：功率约束为有效约束
            leak_pow_tmp = max(eig(X_p))*scene.max_pow - leak_pow_tmp;
            lambda = 0;
            mu = 0;
            precode_mat = calPrecodeMat2(X_0,Z_p,lambda,mu,Y,X_p,precode_mat_tmp);
            %如果两个约束条件有一个不满足，则继续计算
            if(calInterferenceLeak2(Z_p,X_p,precode_mat,precode_mat_tmp) < leak_pow_tmp)
                mu = calMu2(leak_pow_tmp,X_0,Z_p,X_p,Y,precode_mat,precode_mat_tmp,lambda);
                precode_mat = calPrecodeMat2(X_0,Z_p,lambda,mu,Y,X_p,precode_mat_tmp);
            end
            %如果功率约束条件不满足，则继续计算
            %val_P = calTotalPower(precode_mat);
            if(calTotalPower(precode_mat) > scene.max_pow)
                %利用指数步进寻找lambda的上界
                lambda_low = 0;
                lambda_up = 0.1;
                while(calTotalPower(precode_mat) > scene.max_pow)
                    lambda_up = lambda_up*2;
                    mu = 0;
                    precode_mat = calPrecodeMat2(X_0,Z_p,lambda_up,mu,Y,X_p,precode_mat_tmp);
                    if(calInterferenceLeak2(Z_p,X_p,precode_mat,precode_mat_tmp) < leak_pow_tmp)
                        mu = calMu2(leak_pow_tmp,X_0,Z_p,X_p,Y,precode_mat,precode_mat_tmp,lambda_up);
                        precode_mat = calPrecodeMat2(X_0,Z_p,lambda_up,mu,Y,X_p,precode_mat_tmp);
                    end
                end               

                %利用二分搜索求解最优lambda
                while(lambda_up-lambda_low > 0.001)
                    lambda = (lambda_up + lambda_low)/2;
                    mu = 0;
                    precode_mat = calPrecodeMat2(X_0,Z_p,lambda,mu,Y,X_p,precode_mat_tmp);
                    %如果两个约束条件有一个不满足，则继续计算
                    if(calInterferenceLeak2(Z_p,X_p,precode_mat,precode_mat_tmp) < leak_pow_tmp)
                        mu = calMu2(leak_pow_tmp,X_0,Z_p,X_p,Y,precode_mat,precode_mat_tmp,lambda);
                        precode_mat = calPrecodeMat2(X_0,Z_p,lambda,mu,Y,X_p,precode_mat_tmp);
                    end

                    val_P = calTotalPower(precode_mat);
                    if(val_P > scene.max_pow) 
                        lambda_low = lambda; 
                    else
                        lambda_up = lambda;
                    end
                end
                %取上界，保证约束条件满足
                lambda = lambda_up;
                mu = 0;
                precode_mat = calPrecodeMat2(X_0,Z_p,lambda,mu,Y,X_p,precode_mat_tmp);
                %如果两个约束条件有一个不满足，则继续计算
                if(calInterferenceLeak2(Z_p,X_p,precode_mat,precode_mat_tmp) < leak_pow_tmp)
                    mu = calMu2(leak_pow_tmp,X_0,Z_p,X_p,Y,precode_mat,precode_mat_tmp,lambda);
                    precode_mat = calPrecodeMat2(X_0,Z_p,lambda,mu,Y,X_p,precode_mat_tmp);
                end
%                 val_P = calTotalPower(precode_mat);
%                 disp(['val_P = ',num2str(val_P)]);
%                 mu = calMu2(leak_pow_tmp,X_0,Z_p,X_p,Y,precode_mat,precode_mat_tmp,lambda_up);
%                 precode_mat = calPrecodeMat2(X_0,Z_p,lambda_up,mu,Y,X_p,precode_mat_tmp);
%                 val_P = calTotalPower(precode_mat);
%                 disp(['val_P = ',num2str(val_P)]);
%                 disp(['lambda = ',num2str(lambda)]);
            end 
        end
       
        %% 判断是否跳出
        val_obj_tmp = calValObjFunc(X_0,Y,precode_mat_tmp);
        val_obj = calValObjFunc(X_0,Y,precode_mat);
        if(abs(val_obj-val_obj_tmp) < 1e-6)
            break;
        end
    end
end

%% 在情形一中，当mu变化时，计算预编码矩阵时用到的函数
function precode_mat = calPrecodeMat1(X_0,Z_p,mu,Y,X_p,precode_mat_tmp)
     
%     [U,V] = eig(X_0);
%     V = diag(1./(diag(V)+mu*Z_p(1,1)));
%     mat_coeff1 = U*V*inv(U);
    mat_coeff1 = inv(X_0+mu*Z_p);
    mat_coeff2 = mu*(Z_p - X_p);
    precode_mat = zeros(size(precode_mat_tmp));
    for i = 1:size(precode_mat_tmp,3)
        precode_mat(:,:,i) = mat_coeff1*(Y(:,:,i)'+mat_coeff2*precode_mat_tmp(:,:,i));
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
        val_P = val_P + real(trace(precode_mat(:,:,i)'*precode_mat(:,:,i)));
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%情形二相关函数%%%%%%%%%%%%%%%%%%%%%% 
%% 在情形二中，当mu和lambda变化时，计算预编码矩阵时用到的函数
function precode_mat = calPrecodeMat2(X_0,Z_p,lambda,mu,Y,X_p,precode_mat)
    
%     [U,V] = eig(X_0);
%     for i = 1:size(X_0,1)
%         V(i,i) = 1/(V(i,i)+lambda);
%     end
%     tmp_coeff1 = U*V*inv(U);
    tmp_coeff1 = inv(X_0+lambda*eye(size(X_0,1)));
    tmp_coeff2 = mu*(Z_p - X_p);
    for i = 1:size(precode_mat,3)
        precode_mat(:,:,i) = tmp_coeff1*(Y(:,:,i)'+tmp_coeff2*precode_mat(:,:,i));
    end
end

%% 在情形一中，当mu变化时，计算干扰泄漏约束
function val_J = calInterferenceLeak2(Z_p,X_p,precode_mat,precode_mat_tmp)
    val_J = 0;
    for i = 1:size(precode_mat,3)
        val_J = 2*real(trace(precode_mat_tmp(:,:,i)'*(Z_p-X_p)*precode_mat(:,:,i)));
    end
end
%% 在情形二中，给定lamdba下，计算mu值
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

%% 计算目标函数值
function val_obj = calValObjFunc(X_0,Y,precode_mat)
    val_obj = 0;
    for i = 1:size(precode_mat,3)
        val_obj = val_obj + real(trace(precode_mat(:,:,i)'*X_0*precode_mat(:,:,i)));
        val_obj = val_obj - real(trace(Y(:,:,i)*precode_mat(:,:,i)));
        val_obj = val_obj - real(trace(Y(:,:,i)'*precode_mat(:,:,i)'));
    end
end