%% 计算IRS反射系数
function reflect_mat = getReflectMat(scene,channel,precode_mat,reflect_mat)
    % 计算所有预编码矩阵与其共轭转职乘积的和
    Q_s = zeros(scene.n_ante_AP,scene.n_ante_AP);
    for i = 1:scene.n_SU
        Q_s = Q_s + precode_mat(:,:,i)*precode_mat(:,:,i)';
    end
    
    % 计算涉及到的参数
    B_0 = zeros(scene.m_IRS,scene.m_IRS);
    B_p = channel.h_IRS_PU'*channel.h_IRS_PU;
    C = channel.h_AP_IRS*Q_s*channel.h_AP_IRS';
    D_0 = zeros(scene.m_IRS,scene.m_IRS);
    D_p = channel.h_AP_IRS*Q_s*channel.h_AP_PU'*channel.h_IRS_PU;
    for i = 1:scene.n_SU
        B_0 = B_0 + channel.h_IRS_SUs(:,:,i)'*decode_mat(:,:,i)*weight_mat(:,:,i)*decode_mat(:,:,i)'*channel.h_IRS_SUs(:,:,i);
        D_0 = D_0 + channel.h_AP_IRS*Q_s*channel.h_AP_SUs(:,:,i)'*decode_mat(:,:,i)*weight_mat(:,:,i)*decode_mat(:,:,i)'*channel.h_IRS_SUs(:,:,i);
        D_0 = D_0 - channel.h_AP_IRS*precode_mat(:,:,i)*weight_mat(:,:,i)*decode_mat(:,:,i)'*channel.h_IRS_SUs(:,:,i);
    end
    
    leak_pow_tmp = scene.leak_pow-real(trace(channel.h_AP_PU*Q_s*channel.h_AP_PU'));
    gamma_0 = B_0.*C;
    gamma_p = B_p.*C;
    lambda_0 = max(eig(gamma_0));
    lambda_p = max(eig(gamma_p));
    d_0 = diag(D_0);
    d_p = diag(D_p);
    d_0_conj = conj(d_0);
    d_p_conj = conj(d_p);
    
    % 利用连续凸近似和二分搜索求解
    reflect_vec = diag(reflect_mat);
    while(1)
        reflect_vec_tmp = reflect_vec;
        q_0_n = (lambda_0*eye(scene.m_IRS)-gamma_0)*reflect_vec_tmp - d_0_conj;
        q_p_n = (lambda_p*eye(scene.m_IRS)-gamma_p)*reflect_vec_tmp - d_p_conj;
        omga_p_n = (scene.m_IRS*lambda_p+reflect_vec_tmp'*(X_p-gamma_p)*reflect_vec_tmp-leak_pow_tmp)/2;

        reflect_vec = exp(1j*angle(q_0_n));
        if(real(reflect_vec'*q_p_n) < omga_p_n)
            %确定alpha的上下界
            alpha = 0.1;
            while(1)
                reflect_vec = exp(1j*angle(q_0_n+alpha*q_p_n));
                if(real(reflect_vec'*q_p_n) < omga_p_n)
                    alpha = alpha *2;
                else
                    break;
                end
            end
            alpha_low = 0;
            alpha_up = alpha;

            %二分搜索确定alpha取值
            while(1)
                alpha = (alpha_up - alpha_low)/2;
                reflect_vec = exp(1j*angle(q_0_n+alpha*q_p_n));
                if(real(reflect_vec'*q_p_n) < omga_p_n)
                    alpha_low = alpha;
                else
                    alpha_up = alpha;
                end
                           
                if(alpha_up - alpha_low < 0.001)
                    break;
                end
            end
            
        end
        
        % 判断是否跳出循环
        tmp1 = reflect_vec_tmp'*gamma_0*reflect_vec_tmp + 2*real(reflect_vec_tmp'*d_0_conj);
        tmp2 = reflect_vec'*gamma_0*reflect_vec + 2*real(reflect_vec'*d_0_conj);
        if(tmp1 - tmp2 < 0.001)
            break;
        end
    end
    reflect_mat = diag(reflect_vec);
end