function [R] = R_estimation(x,num)
    cntr = num+1;
    for n = 0 : num
        temp = 0;
        for i = 1 : length(x)-n
            temp = temp + x(i) * x(n+i);
        end
        R(cntr - n) =temp;
    end
    R = (1/length(x))*[R(1:cntr) flip(R(1:cntr-1))];
end