module CPTestMakieExt

using Makie
using ClusterPermutationTests

export 	plot_sample_stats!,
		plot_sample_distribution!

function ClusterPermutationTests.plot_sample_stats!(ax::Axis,
	cpt::ClusterPermutationTest)
	para = sample_stats(cpt)
	xs = 1:(length(para))
	lines!(ax, xs, para, color = :red)
	return ax
end

function ClusterPermutationTests.plot_sample_stats!(fig::Figure, cpt::ClusterPermutationTest)
	plot_sample_stats!(Axis(fig[1, 1]), cpt)
	return fig
end;

function ClusterPermutationTests.plot_sample_distribution!(fig::Figure, cpt::ClusterPermutationTest;
	xlabel = "test statistics", bins=100)
	cs = cluster_stats(cpt)
	dist = permutation_stats(cpt)
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
