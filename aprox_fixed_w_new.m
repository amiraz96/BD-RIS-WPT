function [coeff, f0] = aprox_fixed_w_new(w0, h, K2, K4, sub_number, N, MODE_P)

%h should be a row vector
coeff = zeros(sub_number, N);

w0 = abs(w0);
h = abs(h);

for n = 1:sub_number
    first_term = (h(n)*w0(n))*h(n);
    sum1 = 0;
    sum3 = 0;
    for n1 = 1:sub_number
        if n1 ~= n
            sum1 = sum1 + 8*((norm(w0(n1)))^2)*((norm(h(n)))^2)*((norm(h(n1)))^2)*w0(n).';
        end
%         sum1 = sum1 + 8*((norm(w0(n1)))^2)*((norm(h(n)))^2)*((norm(h(n1)))^2)*w0(n).';
        sum2 = 0;
        for n2 = 1:sub_number
            for n3 = 1:sub_number
                if (n2 ~= n3) && (n2 + n3 == 2*n)
                    sum2 = sum2 + 2*(h(n2)*w0(n2))'*(h(n)*w0(n))*(h(n3)*w0(n3))'.*h(n) + ...
                        2*(h(n2)*w0(n2))*(h(n)*w0(n))'*(h(n3)*w0(n3)).*conj(h(n));
                end
                if (-n1 +n2 + n3 == n) && (n ~= n1) && (n ~= n2) && (n ~= n3) ...
                        && (n1 ~= n2) && (n1 ~= n3) && (n2 ~= n3)
                    sum3 = sum3 + 2*(h(n1)*w0(n1))*(h(n2)*w0(n2))'*(h(n3)*w0(n3))'.*h(n) + ...
                        2*(h(n1)*w0(n1))*(h(n2)*w0(n2))*(h(n3)*w0(n3))'.*conj(h(n));
                end
            end
        end
    end
    term1 = 4.*((norm(w0(n)))^2)*((norm(h(n)))^4)*w0(n).';
    if MODE_P == 1
        coeff(n) = K2*(first_term) + (3/8)*K4*(term1 + sum1 + sum2 + sum3);
    end
    if MODE_P == 2
        coeff(n) = (first_term);
    end
end

expw = 0;
expw2 = 0;

for nnn = 1:sub_number
    expw = expw + 0.5*K2*norm(h(nnn)*w0(nnn))^2;
    for nnn1 = 1:sub_number
        for nnn2 = 1:sub_number
            for nnn3 = 1:sub_number
                if nnn + nnn1 == nnn2 + nnn3
                    expw2 = expw2 +  ((3/8)*K4*(h(nnn)*w0(nnn))*(h(nnn1)*w0(nnn1))*...
                            (h(nnn2)*w0(nnn2))'*(h(nnn3)*w0(nnn3))');
                end
            end
        end
    end
end
% expw
% expw2
if MODE_P == 1
    f0 = expw + expw2;
end
if MODE_P == 2
    f0 = expw/K2;
end
% f0 = expw/K2;

end
