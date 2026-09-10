module CPTestMakieExt

using Makie
using ClusterPermutationTests

export 	plot_time_series_stats!,
		plot_cluster_nhd!

function ClusterPermutationTests.plot_time_series_stats!(ax::Axis,
	cpt::ClusterPermutationTest,
	effect::Union{Integer, Symbol, String}=1)
	para = time_series_stats(cpt, effect)
	xs = 1:(length(para))
	lines!(ax, xs, para, color = :red)
	return ax
end

function ClusterPermutationTests.plot_time_series_stats!(fig::Figure, cpt::ClusterPermutationTest;
	effect::Union{Integer, Symbol, String}=1)
	plot_time_series_stats!(Axis(fig[1, 1]), cpt, effect)
	return fig
end;

function ClusterPermutationTests.plot_cluster_nhd!(ax::Axis,
	cpt::ClusterPermutationTest,
	effect::Union{Integer, Symbol, String},
	cluster_id::Integer;
	bins=100, color = (:blue, 0.3))

	cs = cluster_mass_stats(cpt, effect)
	dist = cluster_nhd(cpt, effect)
	hist!(ax, dist[:, cluster_id]; bins, color)
	#density!(ax, dist[:, i], color=(:blue, 0.3))
	vlines!(ax, [cs[cluster_id]]; color = :black,
			linewidth = 3, linestyle = :dash)
	return ax
end;

function ClusterPermutationTests.plot_cluster_nhd!(fig::Figure,
	cpt::ClusterPermutationTest,
	effect::Union{Integer, Symbol, String};
	xlabel = "test statistics", bins=100,
	color = (:blue, 0.3))

	dist = cluster_nhd(cpt, effect)
	for i in 1:size(dist, 2)
		ax = Axis(fig[i, 1]; xlabel, ylabel = "count")
		plot_cluster_nhd!(ax, cpt, effect, i; bins, color)
	end
	return fig
end;

function ClusterPermutationTests.plot_cluster_nhd!(fig::Figure,
	cpt::ClusterPermutationTest;
	bins=100,
	color = (:blue, 0.3))

	for e in 1:length(coefnames(cpt))
		for c in 1:length(cluster(cpt, e))
			ax = Axis(fig[e, c], title = "Effect $e, Cluster $c")
			plot_cluster_nhd!(ax, cpt, e, c; bins, color)
		end
	end
	return fig
end;


end # module
