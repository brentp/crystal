import crystal
import crystal.utils as cu
covs, cluster = cu.real_cluster()
formula = "methylation ~ age + gender"
c = crystal.wrapper(crystal.zscore_cluster, formula, cluster, covs, "gender")
crystal.plot.factorplot_cluster(c, covs)
