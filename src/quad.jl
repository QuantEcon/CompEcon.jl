function qnwmonomial1(vcv::AbstractMatrix)
    n = size(vcv, 1)
    @assert n == size(vcv, 2) "Variance covariance matrix must be square"
    n_nodes = 2n

    z1 = zeros(n_nodes, n)

    # In each node, random variable i takes value either 1 or -1, and
    # all other variables take value 0. For example, for N = 2,
    # z1 = [1 0; -1 0; 0 1; 0 -1]
    for i=1:n
        z1[2*(i-1)+1:2*i, i] = [1, -1]
    end

    sqrt_vcv = chol(vcv)
    R = sqrt(n)*sqrt_vcv
    ϵj = z1*R
    ωj = ones(n_nodes) ./ n_nodes
    ϵj, ωj
end


function qnwmonomial2(vcv::AbstractMatrix)
    n = size(vcv, 1)
    @assert n == size(vcv, 2) "Variance covariance matrix must be square"
    n_nodes = 2n^2 + 1
    z0 = zeros(1, n)

    z1 = zeros(2n, n)
    # In each node, random variable i takes value either 1 or -1, and
    # all other variables take value 0. For example, for N = 2,
    # z1 = [1 0; -1 0; 0 1; 0 -1]
    for i=1:n
        z1[2*(i-1)+1:2*i, i] = [1, -1]
    end

    z2 = zeros(2n*(n-1), n)
    i = 0

    # In each node, a pair of random variables (p,q) takes either values
    # (1,1) or (1,-1) or (-1,1) or (-1,-1), and all other variables take
    # value 0. For example, for N = 2, `z2 = [1 1; 1 -1; -1 1; -1 1]`
    for p=1:n-1
        for q=p+1:n
            i += 1
            z2[4*(i-1)+1:4*i, p] = [1, -1, 1, -1]
            z2[4*(i-1)+1:4*i, q] = [1, 1, -1, -1]
        end
    end

    sqrt_vcv = chol(vcv)
    R = sqrt(n+2)*sqrt_vcv
    S = sqrt((n+2)/2)*sqrt_vcv
    ϵj = [z0; z1*R; z2*S]
    ωj = vcat(2/(n+2) * ones(size(z0, 1)),
              (4-n)/(2*(n+2)^2) * ones(size(z1, 1)),
               1/(n+2)^2 * ones(size(z2, 1)))
    return ϵj, ωj
end
