module CPTestMakieExt

using Makie
using ClusterPermutationTests

export 	plot_time_series_stats!,
		plot_cluster_nhd!

function ClusterPermutationTests.plot_time_series_stats!(ax::Axis,
	cpt::ClusterPermutationTest)
	para = time_series_stats(cpt)
	xs = 1:(length(para))
	lines!(ax, xs, para, color = :red)
	return ax
end

function ClusterPermutationTests.plot_time_series_stats!(fig::Figure, cpt::ClusterPermutationTest)
	plot_time_series_stats!(Axis(fig[1, 1]), cpt)
	return fig
end;

function ClusterPermutationTests.plot_cluster_nhd!(fig::Figure, cpt::ClusterPermutationTest;
	xlabel = "test statistics", bins=100)
	cs = cluster_mass_stats(cpt)
	dist = cluster_nhd(cpt)
	for i in 1:size(dist, 2)
		ax = Axis(fig[i, 1]; xlabel, ylabel = "count")
		hist!(ax, dist[:, i]; bins)
		#density!(ax, dist[:, i], color=(:blue, 0.9))
		vlines!(ax, [cs[i]]; color = :black,
			linewidth = 3, linestyle = :dash)
	end
	return fig
end;


end # module
