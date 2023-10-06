using Distributions
f = function (μ,ρ,σ,m,N)
    nd = Normal()
    y = zeros(N)
    y[1] = μ - m*sqrt(σ^2/(1-ρ^2))
    y[N] = μ + m*sqrt(σ^2/(1-ρ^2))
    m = zeros(N)
    Δ = (y[N] - y[1])/(N-1)
    π = zeros(N,N)
    for i in 1:N
        y[i] = y[1] + (i-1)*Δ
        m[i] = y[i] + Δ/2
    end
    for i in 1:N
        for j in 1:N
            if j == 1
                π[i,j] = cdf(nd, (m[j] - ρ*y[i] + -μ*(1-ρ))/σ)
            elseif j == N
                π[i,j] = 1 - cdf(nd, (m[j-1] - ρ*y[i] + -μ*(1-ρ))/σ)
            else
                π[i,j] = cdf(nd, (m[j] - ρ*y[i] + -μ*(1-ρ))/σ) - cdf(nd, (m[j-1] - ρ*y[i] + -μ*(1-ρ))/σ)
            end
        end
    end
    print("yᵢ ∈ {")
    for i in 1:N print(y[i] ,", ") end
    print("}")
    print("\n")
    print("π = ",π)         
end
f(1.0,0.9,0.5,3,3)