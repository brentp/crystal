import crystal
import crystal.utils as cu
covs, cluster = cu.real_count_cluster()
covs.head()
formula = "methylation ~ ko"
c = crystal.wrapper(crystal.zscore_cluster, formula, cluster, covs, "ko")
crystal.plot.spaghetti_plot(c, covs)
