
import toolshed as ts
import pandas as pd

covariates_file = "notebooks/covariates.csv"

covs_full = pd.read_csv(covariates_file, index_col=0)

covs_full['id'] = covs_full.index
covs_full['plate'] = covs_full.Basename.map(lambda s: s.split("_")[0])

# find plates with between 5 and 7 females
gplates = covs_full.groupby('plate')['gender'].agg(lambda p : 5 <= (p == 'F').sum() <= 7)
good_plates = gplates[gplates].index

dplates = ['6042324097', '6042324119']
vplates = set(good_plates) - set(dplates)


meth = pd.read_table('/drive/450k/norm.beta.txt.gz', compression='gzip', index_col=0)

for l, iplates in (('v', vplates), ('d', dplates)):

    for i, plate in enumerate(iplates):

        covs = covs_full[covs_full.plate == plate]
        covs.to_csv('notebooks/%s%i-450K-12.covariates.csv' % (l, i))

        smeth = meth.ix[:, covs['id']]
        smeth.to_csv(ts.nopen('notebooks/%s%i-450K-12.meth.txt.gz' % (l, i), 'w'),
        index_label='probe', index=True, sep="\t", float_format='%.3f')
