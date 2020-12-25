%% ����IRS����ϵ��
function reflect_mat = getReflectMat(scene,channel,precode_mat,reflect_mat,decode_mat,weight_mat)
    disp('���뷴��ϵ�����');
    % ��������Ԥ����������乲��תְ�˻��ĺ�
    Q_s = zeros(scene.n_ante_AP,scene.n_ante_AP);
    for i = 1:scene.n_SU
        Q_s = Q_s + precode_mat(:,:,i)*precode_mat(:,:,i)';
    end
    
    % �����漰���Ĳ���
    B_p = channel.h_IRS_PU'*channel.h_IRS_PU;
    C = channel.h_AP_IRS*Q_s*channel.h_AP_IRS';
    D_p = channel.h_AP_IRS*Q_s*channel.h_AP_PU'*channel.h_IRS_PU;
    B_0 = zeros(scene.m_IRS,scene.m_IRS);
    D_0 = zeros(scene.m_IRS,scene.m_IRS);
    for i = 1:scene.n_SU
        B_0 = B_0 + channel.h_IRS_SUs(:,:,i)'*decode_mat(:,:,i)*weight_mat(:,:,i)*decode_mat(:,:,i)'*channel.h_IRS_SUs(:,:,i);
        D_0 = D_0 + channel.h_AP_IRS*Q_s*channel.h_AP_SUs(:,:,i)'*decode_mat(:,:,i)*weight_mat(:,:,i)*decode_mat(:,:,i)'*channel.h_IRS_SUs(:,:,i);
        D_0 = D_0 - channel.h_AP_IRS*precode_mat(:,:,i)*weight_mat(:,:,i)*decode_mat(:,:,i)'*channel.h_IRS_SUs(:,:,i);
    end
    
    leak_pow_tmp = scene.leak_pow-real(trace(channel.h_AP_PU*Q_s*channel.h_AP_PU'));
    gamma_0 = B_0.*(C.');
    gamma_p = B_p.*(C.');
    lambda_0 = max(real(eig(gamma_0)));
    lambda_p = max(real(eig(gamma_p)));
    d_0 = diag(D_0);
    d_p = diag(D_p);
    d_0_conj = conj(d_0);
    d_p_conj = conj(d_p);
    
    % ��������͹���ƺͶ����������
    reflect_vec = diag(reflect_mat);
    %while(1)
    %���Ƶ�������������100��
    for cnt_iter = 1:100
        reflect_vec_tmp = reflect_vec;
        q_0_n = (lambda_0*eye(scene.m_IRS)-gamma_0)*reflect_vec_tmp - d_0_conj;
        q_p_n = (lambda_p*eye(scene.m_IRS)-gamma_p)*reflect_vec_tmp - d_p_conj;
        leak_pow_tmp = real((scene.m_IRS*lambda_p+reflect_vec_tmp'*(lambda_p*eye(scene.m_IRS)-gamma_p)*reflect_vec_tmp-leak_pow_tmp)/2);
        disp(['leak_pow_tmp = ',num2str(leak_pow_tmp)]);
        reflect_vec = exp(1j*angle(q_0_n));
        if(real(reflect_vec'*q_p_n) < leak_pow_tmp)
            %ȷ��alpha�����½�
            alpha_low = 0;
            alpha_up = 0.1;
            while(real(reflect_vec'*q_p_n) < leak_pow_tmp)
                alpha_up = alpha_up *2;
                reflect_vec = exp(1j*angle(q_0_n+alpha_up*q_p_n));
                %disp(['alpha_up = ',num2str(alpha_up)]);
            end
            disp('ѭ��1');
            
            %��������ȷ��alphaȡֵ
            while(alpha_up - alpha_low > 0.001)
                alpha = (alpha_up + alpha_low)/2;
                reflect_vec = exp(1j*angle(q_0_n+alpha*q_p_n));
                if(real(reflect_vec'*q_p_n) < leak_pow_tmp)
                    alpha_low = alpha;
                else
                    alpha_up = alpha;
                end
            end
            disp('ѭ��2');
        end
        
        % �ж��Ƿ�����ѭ��
        tmp1 = real(reflect_vec_tmp'*gamma_0*reflect_vec_tmp) + 2*real(reflect_vec_tmp'*d_0_conj);
        tmp2 = real(reflect_vec'*gamma_0*reflect_vec) + 2*real(reflect_vec'*d_0_conj);
        %disp(['tmp1 = ',num2str(tmp1),';tmp2 = ',num2str(tmp2)]);
        if(abs(tmp2 - tmp1) < 0.001)
            break;
        end
    end
    reflect_mat = diag(reflect_vec);
    disp('��������ϵ�����');
end