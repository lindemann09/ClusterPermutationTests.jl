###
### CPCollection
###

struct CPCollection{M}
	iv::Symbol # name of the to be shuffled independent variable
	mass_fnc::Function # cluster mass function
	cc::TClusterCritODef # cluster definition

	m::Vector{M} # fitted models of initial fit
	ts::TParameterVector # time series statistics of the initial fit
	S::Vector{TParameterVector} # time series stats for the permutations
end;

CPCollection{M}(iv::SymbolOString, mass_fnc::Function, cluster_criterium::TClusterCritODef) where {M} =
	CPCollection(Symbol(iv), mass_fnc, cluster_criterium, M[], TParameterVector(), TParameterVector[])

npermutations(x::CPCollection) = length(x.S)


function perm_table(list_para_vector::Vector{TParameterVector})
	time_points = Int32[]
	perm = Int32[]
	parameters = Vector{Float64}[]
	for (cnt, para_vec) in enumerate(list_para_vector)
		for tp in para_vec
			push!(time_points, tp.t)
			push!(perm, cnt)
			push!(parameters, [tp.z..., 67, 34])
		end
	end
	matrix = reduce(hcat, parameters)
    table = (; perm=perm, time = time_points, (Symbol("v", i) => matrix[i, :] for i in 1:size(matrix, 1))...)
	return (time_points, table)
end


function _permutation_stats(x::CPCollection, cl_ranges::Vector{TClusterRange}, parameter_id::Integer)::Matrix{Float64}
	if length(x.S) == 0
		return zeros(Float64, 0, 0)
	else
		rtn = Vector{Float64}[]
		for cl in cl_ranges
			l = length(cl)
			cl_mass_para = Float64[]
			for para_vector in x.S
				# get all parameter values of this cluster of this permutation
				cl_para = get_parameter(para_vector, cl, parameter_id) # FIXME optimize maybe, create first a structure that easy to search
				mp = length(cl_para) == l ? x.mass_fnc(cl_para) : NaN # check if parameter for all samples are found
				push!(cl_mass_para, mp)
			end
			push!(rtn, cl_mass_para)
		end
		return reduce(hcat, rtn)
	end
end

###
### ClusterPermutationTest
###

abstract type ClusterPermutationTest end
#requires
#	cpc::CPCollection
#	dat::CPData

nepochs(x::ClusterPermutationTest) = nepochs(x.dat)
epoch_length(x::ClusterPermutationTest) = epoch_length(x.dat)
design_table(x::ClusterPermutationTest) = design_table(x.dat)
StudyDesigns.unit_observation(x::ClusterPermutationTest) = unit_observation(x.dat.design.uo)

cluster_criterium(x::ClusterPermutationTest) = x.cpc.cc
initial_fits(x::ClusterPermutationTest) = x.cpc.m
npermutations(x::ClusterPermutationTest) = npermutations(x.cpc)


function StatsAPI.summary(x::ClusterPermutationTest)
	println("$(test_info(x))")
	println("  data: $(nepochs(x)) x $(epoch_length(x))")
	pt = pretty_table(String, cluster_table(x);
		show_subheader = false,
		formatters = (ft_printf("%0.2f", [5, 6]),
			ft_printf("%0.3f", [7])),
		vlines = :none,
		tf = tf_unicode_rounded)
	print(pt)
	return println("n permutations: $(npermutations(x))")
end;

function Base.show(io::IO, mime::MIME"text/plain", x::ClusterPermutationTest)
	println(io, "$(test_info(x))")
	println(io, "  data: $(nepochs(x)) x $(epoch_length(x))")
	clr = TClusterRange[]
	try
		clr = cluster_ranges(x)
	catch e
		println(io, "  WARNING: Could not compute cluster ranges. Effect probably not found.")
	end
	cc = cluster_criterium(x)
	if cc isa ClusterDefinition
		println(io, "  $(length(clr)) cluster (ranges=$(cc.ranges)):")
	else
		println(io,
			"  $(length(clr)) cluster (threshold=$(cc.threshold), min_size=$(cc.min_size)):")
	end
	return println(io, "  $(npermutations(x)) permutations")
end;
