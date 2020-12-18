%% 计算解码矩阵和权重矩阵
function [decode_mat,weight_mat] = getDecodeAndWeightMat(sig_mat,jam_mat,g_AP_SUs,precode_mat)
    [n_ante_SU,~,n_SU] = size(g_AP_SUs);
    n_data = size(precode_mat,2);
    I_d_d = eye(n_data);
    decode_mat = zeros(n_ante_SU,n_data,n_SU);
    weight_mat = zeros(n_data,n_data,n_SU);
    for i = 1:n_SU
        decode_mat(:,:,i) = inv(jam_mat(:,:,i)+sig_mat(:,:,i))*g_AP_SUs(:,:,i)*precode_mat(:,:,i);
        weight_mat(:,:,i) = inv(I_d_d-precode_mat(:,:,i)'*g_AP_SUs(:,:,i)'*inv(jam_mat(:,:,i)+sig_mat(:,:,i))*g_AP_SUs(:,:,i)*precode_mat(:,:,i));
    end
end