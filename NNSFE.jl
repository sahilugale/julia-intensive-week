using FFTW, CairoMakie;  # Libraries for FFT and plotting

# Variables
N = 1000;  # Number of spatial points
tfinal = 5.0;  # Final time
dt = tfinal / 10000;  # Time step
M = Int(tfinal / dt);  # Total number of time steps
L = 50.0;  # Spatial domain length
x_grid = 20 * LinRange(-1, 1, N);  # Spatial grid

# Fourier Transform Setup
Lx = N;  # Grid length
Lt = M;  # Total time steps
k = 2 * π * (fftshift(-Lx/2:Lx/2 - 1)) / (x_grid[end] - x_grid[1]);  # Wavenumbers

# Time evolution function
function time_evolution(psi, alpha, dt, k, g, s, V)
    """
    Evolves the wave function `psi` in time using the specified parameters.

    # Paramters
    - `psi`: The current wave function (complex array).
    - `dt`: Time step size (float).
    - `k`: Wavenumbers in Fourier space (vector).
    - `g`: Nonlinearity parameter (float).
    - `s`: Saturation parameter (float).
    - `V`: Potential term (float).

    # Returns
    - The evolved wave function after one time step (complex array).
    """
    psi_nonlinear = fft(exp.(-dt * 1im * ((g .* psi .* conj(psi)) ./ (1 .+ s .* psi .* conj(psi)) .+ V)) .* psi)
    psi_linear = ifft(exp.(-0.5 * dt * 1im * abs.(k).^alpha) .* psi_nonlinear)
    return psi_linear
end

# Potential functions

# Step potential
function step_potential(x_grid, x0, V0)
    V = zeros(length(x_grid))  # Initialize the potential array with zeros
    for i in 1:length(x_grid)
        if x_grid[i] >= x0
            V[i] = V0  # Assign V0 for x >= x0
        end
    end
    return V
end

# Harmonic potential
function harmonic_potential(x_grid)
    return x_grid .^ 2  # Harmonic potential V(x) = x^2
end

# quartic potential
function quartic_potential(x_grid)
    return x_grid .^ 4  # Quartic potential V(x) = x^4
end