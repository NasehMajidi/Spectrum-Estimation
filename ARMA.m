function [f,P_ARMA] = ARMA(Xn, p_algo, q_algo)
    f = 0:0.001:0.5;
    Rx = R_estimation(Xn, p_algo+q_algo);
    cntr = (p_algo+q_algo)+1;
    R = fliplr(Rx(cntr + q_algo - p_algo + 1 : cntr + q_algo));
    for i = 1 : p_algo-1
       R = [R; fliplr(Rx(cntr + q_algo - p_algo + 1 + i : cntr + q_algo + i))];
    end
    P = -Rx(cntr + q_algo + 1 : cntr + q_algo + p_algo)';
    a = inv(R)*P;
    Vn = filter([1,a'],[1],Xn);
    Rv = R_estimation(Vn,q_algo);
    matrix_v = exp(1i*2*pi*f*q_algo);
    for i=1 : 2*q_algo
        matrix_v = [matrix_v; matrix_v(i,:).*exp(-1i*2*pi*f)];
    end
    Pv = Rv*matrix_v;
    matrix_a = ones(1,length(f));
    for i = 1:p_algo
        matrix_a = [matrix_a;  matrix_a(i,:).*exp(-1i*2*pi*f)];
    end
    Pa = abs([1,a']*matrix_a).^2;
    P_ARMA = real(Pv)./real(Pa);
end