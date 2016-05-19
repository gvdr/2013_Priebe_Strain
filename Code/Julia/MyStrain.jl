module MyStrain

using DataFrames

export scaled_ones, scale_matrix, rsq, proc_reflect, proct_rotate_reflect, proc_OSS, get_M_matrix, strain, strain_ranks 


function scaled_ones(n)
	eye(n) - (1/n) * ones(n,n)
end

function scale_matrix(A)
	return scaled_ones(size(A,1)) * A
end

function rsq(A,B)
	tbatab = transpose(B) * A * transpose(A) * B
	return sum(sqrt(abs(eigvals(tbatab))))
end

function proc_reflect(A,B)
	Abar = scale_matrix(A)
	Bbar = scale_matrix(B)
	return rsq(Abar,Bbar) / sum(diag(transpose(Bbar) * Bbar))
end

function proc_rotate_reflect(A,B)
	X = transpose(scale_matrix(A)) * scale_matrix(B)
	Xsvd = svdfact!(X)
	return transpose(Xsvd.Vt) * transpose(Xsvd.U)
end

function proc_OSS(A,B)
	R = proc_rotate_reflect(A,B)
	s = proc_reflect(A,B)
	Ahat = scale_matrix(A)
	Bhat = scale_matrix(B) * R * s
	resids = Ahat - Bhat
	return sum(diag(transpose(resids) * resids))
end

function get_M_matrix(X)
	size_M = size(X,1)
	Xsvd = svdfact(X)
	Xhat1 = Xsvd.U[:,1:size_M] * diagm(sqrt(Xsvd.S[1:size_M]))
	Xhat2 = transpose(Xsvd.Vt[:,1:size_M]) * diagm(sqrt(Xsvd.S[1:size_M]))
	return hcat(Xhat1,Xhat2)
end


function strain(X,ind,max_rank)
	size_M = size(X,1)
	arr = deleteat!([1:size_M;],ind)
	M_obs = get_M_matrix(X)
	M_ind = get_M_matrix(X[arr,arr])
	strain_ind = zeros(max_rank)
	for rank in 1:1:max_rank
		strain_ind[rank] = proc_OSS(M_ind[:,1:rank], M_obs[arr,1:rank])
	end
	return strain_ind
end


function strain_ranks(X,max_rank)
	size_M = size(X,1)
	strains = Array(Float64, (max_rank, size_M))
	strains = @parallel (hcat) for species in 1:size_M
		strain(X,species,max_rank)
	end
	return convert(DataFrame,strains)
end

end
