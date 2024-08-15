function spectrum = generate_spectrum(soln_hat)



[Nx, Ny] = size(soln_hat);
spectrum = zeros(1, Nx/2);


for j = 1:Ny/2
    for i = 1:Nx/2
        modes = floor(sqrt(i^2 + j^2));
        if modes <= Nx/2
            spectrum(modes) = spectrum(modes) + abs(soln_hat(i,j))^2;
        end
    end
end
modes = 1:Nx/2;
spectrum = sqrt(spectrum)/Nx^2;
end