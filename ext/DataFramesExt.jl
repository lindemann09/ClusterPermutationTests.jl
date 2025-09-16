module DataFramesExt

# convenience functions for using DataFrames in ClusterPermutationTests

using DataFrames
using ClusterPermutationTests

DataFrames.DataFrame(perm_design::PermutationDesign) =
    DataFrame(design_table(perm_design))

end;