function value = gaussian( x, mu, sigma )
% Calculates the value of the unnormalized Gaussian function with variable
% x, center at mu and width sigma

value = exp(-(x-mu).^2./(2*sigma.^2));

end

