function [coeffs, adjusted_points] = minimos_quadrados(x, y, n)

    % deslocando a curva para passar na origem
    y_zero = y(1);
    for i = 1:size(y, 2)
        y(i) = y(i) - y_zero;
    end

    % Definição da Matriz A
    A = zeros(length(x), n);

    for i = 1:n
        A(:, i) = x.^(n-i+1); % Preenche cada coluna com x^(n-i+1)
    end

    newX = (A' * A) \ (A' * y');

    coeffs = [newX; 0];
    adjusted_points = polyval(coeffs, x);
    
    % deslocando a curva de volta
    for i = 1:size(adjusted_points, 2)
        adjusted_points(i) = adjusted_points(i) + y_zero;
    end
end