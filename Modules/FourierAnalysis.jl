function FS_analysis(fx_vec, x_vec; N=20)
    p = x_vec[end] - x_vec[1]
    n_vec = 0:N-1
    an_vec = zeros(N)
    bn_vec = zeros(N)
    for (i_n, n) in enumerate(n_vec)
        problem = SampledIntegralProblem((@. fx_vec * cos((2*pi*n/p)*x_vec)), x_vec)
        sol = solve(problem, TrapezoidalRule()) 
        an_vec[i_n] = 2/p * sol.u
        problem = SampledIntegralProblem((@. fx_vec * sin((2*pi*n/p)*x_vec)), x_vec)
        sol = solve(problem, TrapezoidalRule())
        bn_vec[i_n] = 2/p * sol.u
    end
    return an_vec, bn_vec, n_vec
end

function FS_synthesis(an_vec, bn_vec, p)
    f(x) = an_vec[1]/2 + sum( an_vec[n+1]*cos((2*pi*n/p)*x) + bn_vec[n+1]*sin((2*pi*n/p)*x) for n in range(1, length(bn_vec)-1) )
end