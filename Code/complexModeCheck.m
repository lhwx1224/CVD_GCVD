function modes_hat_checked =  complexModeCheck(modes_hat)

[m, n] = size(modes_hat);
modes_hat_checked = zeros(size(modes_hat));
for k = 1:n/2
    if isreal(modes_hat(:,2*k-1))
        warning(['Mode #',num2str(k),' is real-valued, random perturbation is imposed!'])
        modes_hat_checked(:, 2*k-1) = modes_hat(:, 2*k-1) + 1i*1e-6*randn(size(modes_hat(:,2*k-1)));
        modes_hat_checked(:, 2*k) = real(modes_hat_checked(:, 2*k-1)) - 1i*imag(modes_hat_checked(:, 2*k-1));
    else
        modes_hat_checked(:, 2*k-1) = modes_hat(:, 2*k-1);
        modes_hat_checked(:, 2*k) = modes_hat(:, 2*k);
    end
end