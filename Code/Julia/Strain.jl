
#using Gadfly
using Graphs
@everywhere using DataFrames
@everywhere using DistributedArrays
@everywhere using PTools

@everywhere function scaled_ones(n)
    eye(n) - (1/n) * ones(n,n)
end

@everywhere function scale_matrix(A)
    return scaled_ones(size(A,1)) * A
end

@everywhere function rsq(A,B)
    tbatab = transpose(B) * A * transpose(A) * B
    return sum(sqrt(abs(eigvals(tbatab))))
end

@everywhere function proc_reflect(A,B)
    Abar = scale_matrix(A)
    Bbar = scale_matrix(B)
    return rsq(Abar,Bbar) / sum(diag(transpose(Bbar) * Bbar))
end

@everywhere function proc_rotate_reflect(A,B)
    X = transpose(scale_matrix(A)) * scale_matrix(B)
    Xsvd = svdfact!(X)
    return transpose(Xsvd.Vt) * transpose(Xsvd.U)
end

@everywhere function proc_OSS(A,B)
    R = proc_rotate_reflect(A,B)
    s = proc_reflect(A,B)
    Ahat = scale_matrix(A)
    Bhat = scale_matrix(B) * R * s
    resids = Ahat - Bhat
    return sum(diag(transpose(resids) * resids))
end

@everywhere function get_M_matrix(X)
    size_M = size(X,1)
    Xsvd = svdfact(X)
    Xhat1 = Xsvd.U[:,1:size_M] * diagm(sqrt(Xsvd.S[1:size_M]))
    Xhat2 = transpose(Xsvd.Vt[:,1:size_M]) * diagm(sqrt(Xsvd.S[1:size_M]))
    return hcat(Xhat1,Xhat2)
end

function strain(X,ind,max_rank)
    size_M = size(X,1)
    arr = deleteat!([1:size_M],ind)
    M_obs = get_M_matrix(X)
    M_ind = get_M_matrix(X[arr,arr])
    strain_ind = dzeros(max_rank)
    for rank in 1:1:max_rank
        strain_ind[rank] = @spawn proc_OSS(M_ind[:,1:rank], M_obs[arr,1:rank])
    end
    fetch(strain_ind)
    return strain_ind
end

function strain_ranks_p(X,max_rank)
    size_M = size(X,1)
    strains = SharedArray(Float64, (max_rank, size_M))
    @parallel for species in 1:size_M
        strains[:,species] = strain(X,species,max_rank)
    end
    return convert(DataFrame,strains)
end

function strain_ranks(X,max_rank)
    size_M = size(X,1)
    strains = fill(0.0,max_rank,size_M)
    for species in 1:size_M
        strains[:,species] = strain(X,species,max_rank)
    end
    return convert(DataFrame,strains)
end

M = rand(0:1,100,100);
@time strains = strain_ranks(M,10);
#plot(melt(strains),x="variable",y="value",Geom.boxplot)

M = rand(0:1,500,500);
strains = strain_ranks(M,10);
plot(melt(strains),x="variable",y="value",Geom.boxplot)

S = SharedArray(Int, (5, 8))

DataFrames.mean(strains)

test_graph = erdos_renyi_graph(10,0.4,is_directed=false)
 u_er = erdos_renyi_graph(n, m, is_directed=false)


test_niche = niche_model_graph(100,0.5)

ar[1.0]


