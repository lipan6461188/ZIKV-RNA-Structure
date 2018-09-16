

import icSHAPE



ids = """
'KU740184|China|2016|Zika'
'KU820898|China|2016|Zika'
'KU501215.1|Puerto_Rico|2015|Zika'
'KU365778|Brazil|2015|Zika'
'KU312312|Suriname|2015|Zika'
'KU365777|Brazil|2015|Zika'
'KU365780|Brazil|2015|Zika'
'KU707826|Brazil|2015|Zika'
'KU365779|Brazil|2015|Zika'
'KU321639|Brazil|2015|Zika'
'KU509998|Haiti|2014|Zika'
'KU870645|USA|2016|Zika'
'KU501216|Guatemala|2015|Zika'
'KU501217|Guatemala|2015|Zika'
'KU497555|Brazil|2015|Zika'
'KU922960|Mexico|2016|Zika'
'KU922923|Mexico|2016|Zika'
'KU647676|Martinique|2015|Zika'
'KU820897|Colombia|2015|Zika'
'KU527068|Brazil|2015|Zika'
'KU866423|China|2016|Zika'
'KJ776791|French|2013|Zika'
'KX813683|Singapore|2016|Zika'
'KX827309|Singapore|2016|Zika'
'KU681081|Thailand|2014|Zika'
'KU955593|Cambodia|2010|Zika'
'KU681082|Philippines|2012|Zika'
'EU545988|Micronesia|2007|Zika'
'HQ234499|Malaysia|1966|Zika'
'KF383119|Senegal|2001|Zika'
'AY632535.2|Uganda|1947|Zika'
'KU955594|Uganda|1947|Zika'
'KU955591|Senegal|1984|Zika'
'KU955592|Senegal|1984|Zika'
'KU955595|Senegal|1984|Zika'
"""

ids = ids.strip().split('\n')
ids = [item[1:-1] for item in ids]

def init_matrix(dim):
    import pandas as pd
    import numpy as np
    return pd.DataFrame(np.zeros((dim, dim)))


sequence = icSHAPE.readSeq("/Users/lee/Desktop/Projects/Virus/Virus_Genome/phylogenetic_tree/Zike_Z01.afa", removeVersion=False)
Seq = icSHAPE.readSeq("/Users/lee/Desktop/Projects/Virus/Virus_Genome/final_virus_genome.fa", removeVersion=False)

def sequence_identity(seq_1, seq_2):
    assert(len(seq_1) == len(seq_2))
    total = 0
    identity = 0
    for base_1, base_2 in zip(list(seq_1), list(seq_2)):
        if base_1 == '-' and base_2 == '-':
            pass
        elif base_1 == '-' or base_2 == '-':
            pass
        elif base_1 == base_2:
            identity += 1
            total += 1
        else:
            total += 1
    if total > 200:
        return 1.0 * identity / total
    else:
        return 0.98



cds_start_59 = sequence['KU501215.1|Puerto_Rico|2015|Zika'].find(Seq['KU501215.1'][106:113])
cds_end_59 = sequence['KU501215.1|Puerto_Rico|2015|Zika'].find(Seq['KU501215.1'][10375:10386])



######### Full Length

similarity_matrix = init_matrix(len(ids))
for idx_1 in range(len(ids)):
    for idx_2 in range(idx_1, len(ids)):
        similarity_matrix[idx_1][idx_2] = sequence_identity(sequence[ids[idx_1]], sequence[ids[idx_2]])
        similarity_matrix[idx_2][idx_1] = similarity_matrix[idx_1][idx_2]

sns.heatmap(similarity_matrix, cmap=sns.color_palette("BuGn", 100), vmin=0.85, vmax=1.0)
plt.savefig("/tmp/similarity_full.pdf")
plt.show()

######### UTR

similarity_matrix = init_matrix(len(ids))
for idx_1 in range(len(ids)):
    for idx_2 in range(idx_1, len(ids)):
        seq_1 = sequence[ids[idx_1]][:cds_start_59]+sequence[ids[idx_1]][cds_end_59:]
        seq_2 = sequence[ids[idx_2]][:cds_start_59]+sequence[ids[idx_1]][cds_end_59:]
        similarity_matrix[idx_1][idx_2] = sequence_identity(seq_1, seq_2)
        similarity_matrix[idx_2][idx_1] = similarity_matrix[idx_1][idx_2]


sns.heatmap(similarity_matrix, cmap=sns.color_palette("BuGn", 100), vmin=0.85, vmax=1.0)
plt.savefig("/tmp/similarity_utr.pdf")
plt.show()


######### CDS

similarity_matrix = init_matrix(len(ids))
for idx_1 in range(len(ids)):
    for idx_2 in range(idx_1, len(ids)):
        seq_1 = sequence[ids[idx_1]][cds_start_59:cds_end_59]
        seq_2 = sequence[ids[idx_2]][cds_start_59:cds_end_59]
        similarity_matrix[idx_1][idx_2] = sequence_identity(seq_1, seq_2)
        similarity_matrix[idx_2][idx_1] = similarity_matrix[idx_1][idx_2]




sns.heatmap(similarity_matrix, cmap=sns.cubehelix_palette(30), vmin=0.88, vmax=1.0)
plt.savefig("/tmp/similarity_cds.pdf")
plt.show()




