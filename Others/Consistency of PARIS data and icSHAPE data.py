

"""

Consistency of PARIS data and icSHAPE data

"""


import icSHAPE, structure, visual, random, virus


def read_tab_dg(inFile):
    dg_list = []
    IN = open(inFile)
    line = IN.readline()
    while line:
        if line[0] == '>':
            data = line.strip().split()
            l_start, l_end, r_start, r_end, support = int(data[4]), int(data[5]), int(data[8]), int(data[9]), int(data[10])
            dg_list.append( (l_start, l_end, r_start, r_end, support) )
        line = IN.readline()
    IN.close()
    return dg_list


def find_pair_out_base(dot):
    stack = []
    out_base_list = []
    for idx, symbol in enumerate(list(dot)):
        if symbol == ')':
            if len(stack) == 0:
                out_base_list.append(idx+1)
            else:
                stack.pop()
        elif symbol == '(':
            stack.append(idx+1)
    out_base_list += stack
    return out_base_list


def get_dg_bp_shape(l_start, l_end, r_start, r_end, sequence, shape, visualization=False):
    if r_start - l_end < 20:
        return None, None
    assert(len(sequence) == len(shape))
    
    # expand
    l_start = max(l_start-10, 1)
    l_end = min(len(sequence), l_end+10)
    r_start = max(r_start-10, 1)
    r_end = min(len(sequence), r_end+10)
    
    # get sequence and predict structure
    seq_1 = sequence[l_start-1:l_end]
    seq_2 = sequence[r_start-1:r_end]
    try:
        pred_ss = structure.bi_fold(seq_1, seq_2, local_pairing=True, mfe=True)
    except ValueError, e:
        return None, None
    
    # get shape and local structure
    cur_seq = seq_1 + "III" + seq_2
    cur_shape = shape[l_start-1:l_end] + ['NULL']*3 + shape[r_start-1:r_end]
    left_ss = pred_ss[:len(seq_1)]
    right_ss = pred_ss[len(seq_1)+3:]
    
    # find stem from 
    stems = virus.find_stem(structure.dot2ct(pred_ss), max_stem_gap=1, min_stem_len=7)
    stems_list = []
    for index in range(stems.size()):
        if stems.at(index).l_end <= len(seq_1) and stems.at(index).r_start >= len(seq_1)+3:
            stems_list.append( (stems.at(index).l_start, stems.at(index).l_end, stems.at(index).r_start, stems.at(index).r_end) )
    if len(stems_list) == 0:
        return None, None
    if visualization:
        stems.show()
    
    # find interaction base pairs
    left_list = find_pair_out_base(left_ss)
    right_list = find_pair_out_base(right_ss)
    bias = len(seq_1) + 3
    right_list = [ index+bias for index in right_list ]
    
    # collect shape scores
    try:
        interaction_bs_shape_list = []
        for base_idx in left_list+right_list:
            if cur_shape[base_idx-1] != 'NULL':
                for cur_stem in stems_list:
                    if cur_stem[0] < base_idx < cur_stem[1] or cur_stem[2] < base_idx < cur_stem[3]:
                        interaction_bs_shape_list.append( float(cur_shape[base_idx-1]) )
                        break
    except IndexError:
        print cur_shape, left_list, right_list
    
    if visualization:
        print interaction_bs_shape_list
        visual.Plot_RNAStructure_Shape(cur_seq, pred_ss, cur_shape, mode='fill')
    
    # return
    if len(interaction_bs_shape_list) > 5:
        rand_shape = [ float(item) for item in random.sample(shape, len(interaction_bs_shape_list)) if item != "NULL" ]
        return 1.0*sum(interaction_bs_shape_list)/len(interaction_bs_shape_list), sum(rand_shape)/len(rand_shape)
    else:
        return None, None

def collect_dg_shape(dg_list, trans_id, visualization=True):
    short_shape_list = []; rand_list = []; long_shape_list = []
    random.shuffle(dg_list)
    for dg in dg_list:
        cur_shape, cur_rand = get_dg_bp_shape(dg[0], dg[1], dg[2], dg[3], sequence[trans_id], shape[trans_id], visualization=visualization)
        #cur_shape, cur_rand = get_dg_bp_shape(dg[0], dg[1], dg[2], dg[3], sequence['KU501215.1'], shape['KU501215.1'], visualization=visualization)
        if cur_shape:
            if dg[2] - dg[1] > 1000:
                long_shape_list.append(cur_shape)
            else:
                short_shape_list.append(cur_shape)
            rand_list.append(cur_rand)
    return long_shape_list, short_shape_list, rand_list


sequence = icSHAPE.readSeq("/Users/lee/Desktop/Projects/Virus/Virus_Genome/final_virus_genome.fa", removeVersion=False)
shape = icSHAPE.loadicSHAPE("/Users/lee/Desktop/Projects/Virus/Virus_Genome/final_virus_icSHAPE.out", removeVersion=False)

dg_list = read_tab_dg("/tmp/virus_paris/766/766.tab")

long_shape_scores, short_shape_scores, rand_scores = collect_dg_shape(dg_list, visualization=False, trans_id="KU501215.1")
sum(long_shape_scores)/len(long_shape_scores)
sum(short_shape_scores)/len(short_shape_scores)
sum(rand_scores)/len(rand_scores)


################  Density plot

sns.distplot(rand_scores, rug=True, bins=50)
sns.distplot(short_shape_scores, rug=True, bins=50)
sns.distplot(long_shape_scores, rug=True, bins=50)
plt.savefig("/tmp/paris_shape.pdf")

plt.show()

################  significance test

# mannwhitney U test

scipy.stats.mannwhitneyu(rand_scores, short_shape_scores, use_continuity=True, alternative=None)
scipy.stats.mannwhitneyu(rand_scores, long_shape_scores, use_continuity=True, alternative=None)

# KS test

from scipy.stats import ks_2samp
import statsmodels, statsmodels.sandbox.stats
p1 = ks_2samp(long_shape_scores, rand_scores)[1]
p2 = ks_2samp(short_shape_scores, rand_scores)[1]
statsmodels.sandbox.stats.multicomp.multipletests((p1, p2), method="bonferroni")


################  CDF plot

visual.cdf(short_shape_scores, bins=100, color="red")
visual.cdf(long_shape_scores, bins=100, color="green")
visual.cdf(rand_scores, bins=100, color="blue")
plt.xlim(0, 1.0)
plt.savefig("/tmp/766.pdf")
plt.show()












