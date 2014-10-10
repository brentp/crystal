import sys
import toolshed as ts
from crystal.evaluate import evaluate
from matplotlib import pyplot as plt
import seaborn as sns

rbed = sys.argv[1]


fig, ax = plt.subplots(1)

for pbed in sys.argv[2:]:
    bediter = [(d.get('chrom', d.get('chr', d.get('seqnames'))), int(float(d['start'])), int(float(d['end'])),
                   float(d.get('p', d.get('median.p', d.get('z_p', d.get('slk_p'))))), 1) for d in ts.reader(pbed)]

    evaluate(bediter, rbed, ax=ax, label=pbed.rsplit(".", 1)[0].split("/", 1)[-1])


ax.plot([0, 1], [0, 1], 'k--')
ax.legend()
plt.show()

