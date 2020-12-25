%% ����Ԥ�������
function precode_mat = getPrecodeMat(scene,g_AP_PU,g_AP_SUs,decode_mat,weight_mat,precode_mat)   
    %% ������ز���
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
        %precode_mat_tmpΪ����͹����ÿһ�ֵ����ĳ�ʼ��,leak_pow_tmpΪ��Ӧ�Ľ��Ƹ���й©��ֵ
        precode_mat_tmp = precode_mat;
        leak_pow_tmp = scene.leak_pow;
        for i = 1:scene.n_SU
            leak_pow_tmp = leak_pow_tmp - real(trace(precode_mat_tmp(:,:,i)'*(Z_p-X_p)*precode_mat_tmp(:,:,i)));
        end
        
       %% ����һ������Լ��Ϊ��ЧԼ��
        mu = 0;
        precode_mat = calPrecodeMat1(X_0,Z_p,mu,Y,X_p,precode_mat_tmp);
        val_J = calInterferenceLeak1(Z_p,X_p,precode_mat,precode_mat_tmp);
        if(val_J > leak_pow_tmp)         
            %����ָ������ȷ��mu���Ͻ�
            mu_low = 0;
            mu_up = 0.1;
            while(val_J > leak_pow_tmp)
                mu_up = mu_up * 2;
%                 val_J_tmp = val_J;
                precode_mat = calPrecodeMat1(X_0,Z_p,mu_up,Y,X_p,precode_mat_tmp);
                val_J = calInterferenceLeak1(Z_p,X_p,precode_mat,precode_mat_tmp);
%                 if(val_J > val_J_tmp)
%                      disp('����');
%                      pause;
%                 end
            end
            
            %������õ��Ͻ���½磬���ö�����������Ͻ���½�
            while(mu_up - mu_low > 0.001)
                mu = (mu_low + mu_up)/2;
                precode_mat = calPrecodeMat1(X_0,Z_p,mu,Y,X_p,precode_mat_tmp);
                val_J = calInterferenceLeak1(Z_p,X_p,precode_mat,precode_mat_tmp);
                if(val_J > leak_pow_tmp) mu_low = mu; else mu_up = mu; end
            end
            %ȡmu���Ͻ磬��֤Լ����������
            precode_mat = calPrecodeMat1(X_0,Z_p,mu_up,Y,X_p,precode_mat_tmp);
        end
        %% �ж��Ƿ����㹦��Լ������������㣬�������ζ�
        if(calTotalPower(precode_mat) > scene.max_pow)
            %% ���ζ�������Լ��Ϊ��ЧԼ��
            leak_pow_tmp = max(eig(X_p))*scene.max_pow - leak_pow_tmp;
            lambda = 0;
            mu = 0;
            precode_mat = calPrecodeMat2(X_0,Z_p,lambda,mu,Y,X_p,precode_mat_tmp);
            %�������Լ��������һ�������㣬���������
            if(calInterferenceLeak2(Z_p,X_p,precode_mat,precode_mat_tmp) < leak_pow_tmp)
                mu = calMu2(leak_pow_tmp,X_0,Z_p,X_p,Y,precode_mat,precode_mat_tmp,lambda);
                precode_mat = calPrecodeMat2(X_0,Z_p,lambda,mu,Y,X_p,precode_mat_tmp);
            end
            %�������Լ�����������㣬���������
            %val_P = calTotalPower(precode_mat);
            if(calTotalPower(precode_mat) > scene.max_pow)
                %����ָ������Ѱ��lambda���Ͻ�
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

                %���ö��������������lambda
                while(lambda_up-lambda_low > 0.001)
                    lambda = (lambda_up + lambda_low)/2;
                    mu = 0;
                    precode_mat = calPrecodeMat2(X_0,Z_p,lambda,mu,Y,X_p,precode_mat_tmp);
                    %�������Լ��������һ�������㣬���������
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
                %ȡ�Ͻ磬��֤Լ����������
                lambda = lambda_up;
                mu = 0;
                precode_mat = calPrecodeMat2(X_0,Z_p,lambda,mu,Y,X_p,precode_mat_tmp);
                %�������Լ��������һ�������㣬���������
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
       
        %% �ж��Ƿ�����
        val_obj_tmp = calValObjFunc(X_0,Y,precode_mat_tmp);
        val_obj = calValObjFunc(X_0,Y,precode_mat);
        if(abs(val_obj-val_obj_tmp) < 1e-6)
            break;
        end
    end
end

%% ������һ�У���mu�仯ʱ������Ԥ�������ʱ�õ��ĺ���
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

%% ������һ�У���mu�仯ʱ���������й©
function val_J = calInterferenceLeak1(Z_p,X_p,precode_mat,precode_mat_tmp)
    n_SU = size(precode_mat,3);
    val_J = 0;
    for i = 1:n_SU
        val_J = val_J + real(trace(precode_mat(:,:,i)'*Z_p*precode_mat(:,:,i)));
        val_J = val_J - 2*real(trace(precode_mat_tmp(:,:,i)'*(Z_p-X_p)*precode_mat(:,:,i)));
    end
end

%% �ж��Ƿ����㹦��Լ��
function val_P = calTotalPower(precode_mat)
    val_P = 0;
    for i = 1:size(precode_mat,3)
        val_P = val_P + real(trace(precode_mat(:,:,i)'*precode_mat(:,:,i)));
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%���ζ���غ���%%%%%%%%%%%%%%%%%%%%%% 
%% �����ζ��У���mu��lambda�仯ʱ������Ԥ�������ʱ�õ��ĺ���
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

%% ������һ�У���mu�仯ʱ���������й©Լ��
function val_J = calInterferenceLeak2(Z_p,X_p,precode_mat,precode_mat_tmp)
    val_J = 0;
    for i = 1:size(precode_mat,3)
        val_J = 2*real(trace(precode_mat_tmp(:,:,i)'*(Z_p-X_p)*precode_mat(:,:,i)));
    end
end
%% �����ζ��У�����lamdba�£�����muֵ
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

%% ����Ŀ�꺯��ֵ
function val_obj = calValObjFunc(X_0,Y,precode_mat)
    val_obj = 0;
    for i = 1:size(precode_mat,3)
        val_obj = val_obj + real(trace(precode_mat(:,:,i)'*X_0*precode_mat(:,:,i)));
        val_obj = val_obj - real(trace(Y(:,:,i)*precode_mat(:,:,i)));
        val_obj = val_obj - real(trace(Y(:,:,i)'*precode_mat(:,:,i)'));
    end
end