
"""

    Sequence identity and icSHAPE correlation

"""


def read_align(align_file):
    align_dict = {}
    IN = open(align_file)
    line = IN.readline()
    while line:
        data = line.strip().split()
        align_dict[data[0]] = data[1]
        line = IN.readline()
    return align_dict

def sequence_identity(seq_1, seq_2):
    assert(len(seq_1) == len(seq_2))
    same_count = 0
    for idx in range(len(seq_1)):
        if seq_1[idx] == seq_2[idx]:
            same_count += 1
    return 1.0 * same_count / len(seq_1)

def shape_correlation(shape_1, shape_2):
    assert(len(shape_1) == len(shape_2))
    vs_1 = []; vs_2 = []
    for s1, s2 in zip(shape_1, shape_2):
        if s1 != 'NULL' and s2 != 'NULL':
            vs_1.append( float(s1) )
            vs_2.append( float(s2) )
    return scipy.stats.pearsonr(vs_1, vs_2)[0]


def Seq_Shape_Coorrelation(shape, align, window_size=200, step=50):
    key_1 = "KU501215.1"; key_2 = "AY632535.2"
    seq_1 = align[key_1]; seq_2 = align[key_2]; shape_1 = []; shape_2 = [];
    idx_1 = -1; idx_2 = -1;
    for idx in range(len(align[key_1])):
        if align[key_1][idx] != '-':
            idx_1 += 1
        if align[key_2][idx] != '-':
            idx_2 += 1
        if align[key_1][idx] == '-':
            shape_1.append( 'NULL' )
        else:
            shape_1.append( shape[key_1][idx_1] )
        if align[key_2][idx] == '-':
            shape_2.append( 'NULL' )
        else:
            shape_2.append( shape[key_2][idx_2] )
    print len(seq_1), len(seq_2), len(shape_1), len(shape_2)
    start = 0; Cor = []; idx = 0
    while len(seq_1) - start > window_size / 2:
        seq_cor = sequence_identity(seq_1[start:start+window_size], seq_2[start:start+window_size])
        shape_cor = shape_correlation(shape_1[start:start+window_size], shape_2[start:start+window_size])
        Cor.append([idx, seq_cor, shape_cor])
        start += step
        idx += 1
    return Cor



from icSHAPE import *

seq = readSeq("/Users/lee/Desktop/Projects/Virus/Virus_Genome/final_virus_genome.fa", removeVersion=False)

v_59 = seq['KU501215.1']
v_766 = seq['AY632535.2']

shape = loadicSHAPE("/Users/lee/Desktop/Projects/Virus/Virus_Genome/final_virus_icSHAPE.out", removeVersion=False)
align = read_align("/Users/lee/Desktop/Projects/Virus/Virus_Genome/final_align.stoch")

window_size = 100
window_step = 20
Cor = Seq_Shape_Coorrelation(shape, align, window_size=window_size, step=window_step)

Cor = pd.DataFrame(Cor, columns=['index', 'seq_identity', 'shape_Cor'])

plt.figure(figsize=(15,3))

plt.subplot(2,1,1)
plt.bar(Cor['index'], Cor['seq_identity'], width=1, linewidth=0, color='black')
plt.ylim(0, 1.0)

plt.subplot(2,1,2)
plt.bar(Cor['index'], Cor['shape_Cor'], width=1, linewidth=0, color='black')
plt.ylim(0, 1.0)

plt.show()


