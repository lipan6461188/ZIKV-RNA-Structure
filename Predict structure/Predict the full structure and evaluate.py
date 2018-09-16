

"""

Predict the virus full structure model

"""

import structure,visual,virus,tools,covariation

def read_dg(inFile):
    dg_list = []
    IN = open(inFile)
    line = IN.readline()
    while line:
        if line[0] == '>':
            data = line.strip().split()
            if float(data[13]) > 0.01:
                dg_list.append( (int(data[4]), int(data[5]), int(data[8]), int(data[9])) )
        line = IN.readline()
    IN.close()
    return dg_list

def filter_intra_stem(stem, break_point):
    new_dg_list = []
    for idx in range(stem.size()):
        ls,le,rs,re = stem.at(idx).l_start, stem.at(idx).l_end, stem.at(idx).r_start, stem.at(idx).r_end
        if ls<le<break_point<rs<re:
            new_dg_list.append( (ls,le,rs,re) )
    return new_dg_list

def accept_single_bp_list(bp_list, min_single_shape):
    if not bp_list:
        return True
    bp_list = sorted(bp_list, reverse=True)
    if len(bp_list) < 3:
        accept_top = True
        #accept_top = all([ it > 0.1 for it in bp_list ])
    else:
        min_accept_num = (len(bp_list)-3) / 3 + 3
        accept_top = all([ it > 0.6 for it in bp_list[:min_accept_num] ])
    accept_mean = sum(bp_list)/len(bp_list)
    return accept_top and accept_mean > min_single_shape

def accept_double_bp_list(bp_list, max_pair_shape):
    if not bp_list:
        return False
    bp_list = sorted(bp_list, reverse=True)
    step = len(bp_list) / 5
    accept_top = all([ it < 0.6 for it in bp_list[step-1:2*step] ])
    accept_mean = sum(bp_list)/len(bp_list)
    return accept_top and accept_mean < max_pair_shape


def dg_predict_stem(dg, sequence, shape, extend, min_single_shape=0.3, max_pair_shape=0.6, test=False):
    import virus
    assert( len(sequence) == len(shape) )
    length = len(sequence)
    ls,le,rs,re = dg[:4]
    if rs-le<1:
        return []
    
    ls = max(1, ls-extend)
    re = min(length, re+extend)
    if rs-le<extend*2:
        le = le + extend/2
        re = re - extend/2
    
    seq_1 = sequence[ls-1:le]
    seq_2 = sequence[rs-1:re]
    shape_1 = shape[ls-1:le]
    shape_2 = shape[rs-1:re]
    
    try:
        ss = structure.bi_fold(seq_1,seq_2,local_pairing=False)
        #ss = structure.predictStructure(seq_1+"III"+seq_2, shape_1+['NULL']*3+shape_2)
    except:
        return []
    stem = virus.find_stem(structure.dot2ct(ss), max_stem_gap=2, min_stem_len=7)
    intra_stem = filter_intra_stem(stem, len(seq_1)+1)
    
    preserve_stem = []
    pair_list = []
    single_list = []
    for cur_stem in intra_stem:
        if test: print cur_stem
        sls,sle,srs,sre = cur_stem[:4]
        if test: print (sls+ls,sle+ls,srs+rs-len(seq_1)-3,sre+rs-len(seq_1)-3)
        for i in range(sls-1+1, sle-1):
            global_coor = ls + i
            cur_shape = shape[global_coor-1]
            if ss[i] == '.':
                if cur_shape != 'NULL':
                    single_list.append( float(cur_shape) )
            else:
                #print sequence[global_coor-1]
                if cur_shape != 'NULL':
                    pair_list.append( float(cur_shape) )
        #print "<==>"
        for i in range(srs-1+1, sre-1):
            global_coor = rs + (i - len(seq_1) - 3)
            cur_shape = shape[global_coor-1]
            if ss[i] == '.':
                if cur_shape != 'NULL':
                    single_list.append( float(cur_shape) )
            else:
                #print global_coor, sequence[global_coor-1]
                if cur_shape != 'NULL':
                    pair_list.append( float(cur_shape) )
        if test: print single_list, pair_list
        #visual.Plot_RNAStructure_Shape(seq_1+"III"+seq_2, ss, shape_1+['NULL']*3+shape_2, mode='fill')
        
        if accept_single_bp_list(single_list, min_single_shape) and accept_double_bp_list(pair_list, max_pair_shape):
            preserve_stem.append( [sls+ls-1,sle+ls-1,srs+rs-len(seq_1)-4,sre+rs-len(seq_1)-4] )
            """
            if (990,996,1409,1415) == (sls+ls-1,sle+ls-1,srs+rs-len(seq_1)-4,sre+rs-len(seq_1)-4):
                print sls,sle,srs,sre
                print Sequence['KU501215.1'][2035-1:2042], Sequence['KU501215.1'][2396-1:2403]
                visual.Plot_RNAStructure_Shape(seq_1+"III"+seq_2, ss, shape_1+['NULL']*3+shape_2, mode='fill')
                return None
            """
        if test: 
            if accept_single_bp_list(single_list, min_single_shape):
                print "accept_single_bp_list"
            if accept_double_bp_list(pair_list, max_pair_shape):
                print "accept_double_bp_list"
        if test:
            visual.Plot_RNAStructure_Shape(seq_1+"III"+seq_2, ss, shape_1+['NULL']*3+shape_2, mode='fill')
        """
        if single_list:
            if sum(single_list)/len(single_list) < min_single_shape:
                continue
        if pair_list and sum(pair_list)/len(pair_list) < max_pair_shape:
            preserve_stem.append( (sls+ls,sle+ls,srs+rs-len(seq_1)-3,sre+rs-len(seq_1)-3) )
        #visual.Plot_RNAStructure_Shape(seq_1+"III"+seq_2, ss, shape_1+['NULL']*3+shape_2, mode='fill')
        """
    
    return preserve_stem

def read_paris_reads(inFile):
    Reads = []
    IN = open(inFile)
    line = IN.readline()
    while line:
        data = line.strip().split()
        Reads.append( (int(data[5]), int(data[6]), int(data[12]), int(data[13])) )
        line = IN.readline()
    IN.close()
    return Reads

def stem_reads_num(stem, Reads):
    num = 0
    for read in Reads:
        if read[0] < stem[1] and read[1] > stem[0] and read[2] < stem[3] and read[3] > stem[2]:
            num += 1
    return num

def is_dg_in_domain(dg, domain, tolerant=0):
    for idx in range(domain.size()):
        start, end = domain.at(idx)
        if start-tolerant <= dg[0] < dg[3] <= end+tolerant:
            return idx+1
    return False

def predict_duplex_dg(dg_list, domain, key, Reads):
    dg_in_domain = []
    dg_out_domain = []
    for i in range(domain.size()): dg_in_domain.append([])
    dg_index = 0
    for dg in dg_list:
        dg_index += 1
        cur_stem = dg_predict_stem(dg, Sequence[key], Shape[key], extend=10, min_single_shape=0.05, max_pair_shape=0.4)
        if cur_stem:
            for stem in cur_stem:
                domain_index = is_dg_in_domain(stem, domain, tolerant=3)
                support_reads = stem_reads_num(stem, Reads)
                if domain_index:
                    dg_in_domain[domain_index-1].append(tuple(stem+[dg_index, support_reads]))
                else:
                    dg_out_domain.append(tuple(stem+[dg_index, support_reads]))
    return dg_in_domain, dg_out_domain


def prune_structure(ss, shape, Reads, base_coor=1):
    assert( len(ss)==len(shape) )
    
    ctList = structure.dot2ct(ss)
    stems = virus.find_stem(ctList, max_stem_gap=2, min_stem_len=1)
    refused_stem_list = []; force_refuse = []
    new_ss = list(ss)
    for idx in range(stems.size()):
        stem = stems.at(idx)
        ls,le,rs,re = stem.l_start, stem.l_end, stem.r_start, stem.r_end
        if le-ls+1 <= 2 or re-rs+1 <= 2:
            #print "1", ls,le,rs,re
            refused_stem_list.append( (ls,le,rs,re) )
            continue
        pair_shape_list = []
        for idx in range(ls-1, le):
            if ss[idx] != '.' and shape[idx] != 'NULL':
                pair_shape_list.append( float(shape[idx]) )
        for idx in range(rs-1, re):
            if ss[idx] != '.' and shape[idx] != 'NULL':
                pair_shape_list.append( float(shape[idx]) )
        pair_shape_list.sort(reverse=True)
        large_shape_num = sum([it>=0.5 for it in pair_shape_list])
        great_large_shape_num = sum([it>=0.7 for it in pair_shape_list])
        #print "0", ls,le,rs,re, pair_shape_list
        if large_shape_num > 1.0*len(pair_shape_list)/3-0.01:
            #print "2", ls,le,rs,re
            refused_stem_list.append( (ls,le,rs,re) )
            if great_large_shape_num >= 1.0*len(pair_shape_list)/2-0.01:
                force_refuse.append((ls,le,rs,re))
    
    #print refused_stem_list
    for ls,le,rs,re in refused_stem_list:
        support = stem_reads_num((ls+base_coor-1,le+base_coor-1,rs+base_coor-1,re+base_coor-1), Reads)
        if support > 20 and (ls,le,rs,re) not in force_refuse:
            continue
        for i in range(ls-1, le):
            new_ss[i] = '.'
        for i in range(rs-1, re):
            new_ss[i] = '.'
    
    return "".join(new_ss)


def evaluate_with_shape(ss, shape_list):
    assert( len(ss) == len(shape_list) )
    single_strand = []
    double_strand = []
    for idx in range(len(ss)):
        if ss[idx] == '.':
            if shape_list[idx] != 'NULL':
                single_strand.append(float(shape_list[idx]))
        else:
            if shape_list[idx] != 'NULL':
                double_strand.append(float(shape_list[idx]))
    plt.subplot(2,1,1)
    tools.plt.violinplot((single_strand, double_strand), (1,2), showmeans=True)
    plt.subplot(2,1,2)
    tools.plt.boxplot((single_strand, double_strand), (1,2), showmeans=True)
    #tools.plt.show()



def evaluate_with_paris(ss, Reads, base_coor, OUT=sys.stdout):
    ctList = structure.dot2ct(ss)
    ctList = [ (it[0]+base_coor-1, it[1]+base_coor-1) for it in ctList ]
    stems = virus.find_stem(ctList, max_stem_gap=3, min_stem_len=6)
    support_list = []
    for idx in range(stems.size()):
        stem = stems.at(idx)
        ls,le,rs,re = stem.l_start, stem.l_end, stem.r_start, stem.r_end
        support = stem_reads_num((ls,le,rs,re), Reads)
        support_list.append( (ls,le,rs,re,support) )
        print >>OUT, "%s-%s <==> %s-%s  %s-%s <==> %s-%s  %s" % (ls,le,rs,re, ls-base_coor+1,le-base_coor+1,rs-base_coor+1,re-base_coor+1, support)
    plot_support_list = []
    tools.sns.swarmplot(data=np.log10([1+it[4] for it in support_list]))
    return support_list

"""

seq = Sequence[key766][5057-1:5266]
shape = Shape[key766][5057-1:5266]
ss = structure.predictStructure(seq, shape)
visual.Plot_RNAStructure_Shape(seq, ss, shape, mode='fill')

"""

############
#  读取数据 Mac
############

key59, key766 = ['KU501215.1', 'AY632535.2']

Sequence = tools.readSeq("/Users/lee/Desktop/Projects/Virus/Virus_Genome/final_virus_genome.fa")
Shape = tools.loadicSHAPE("/Users/lee/Desktop/Projects/Virus/Virus_Genome/final_virus_icSHAPE.out")

dg_list_59 = read_dg("/Users/lee/Desktop/Projects/Virus/Virus_Genome/PARIS/virus_paris/59/59.tab")
dg_list_766 = read_dg("/Users/lee/Desktop/Projects/Virus/Virus_Genome/PARIS/virus_paris/766/766.tab")

reads_59 = read_paris_reads("/Users/lee/Desktop/Projects/Virus/Virus_Genome/PARIS/virus_paris/59/59.dg")
reads_766 = read_paris_reads("/Users/lee/Desktop/Projects/Virus/Virus_Genome/PARIS/virus_paris/766/766.dg")

domain_59 = virus.read_domain("/Users/lee/pCloud Drive/Projects/Virus/re-domain/59_domain.txt")
domain_59.load_seq(Sequence[key59])
domain_59.load_shape(Shape[key59])

domain_766 = virus.read_domain("/Users/lee/pCloud Drive/Projects/Virus/re-domain/766_domain.txt")
domain_766.load_seq(Sequence[key766])
domain_766.load_shape(Shape[key766])


############
#  读取数据 Linux
############

key59, key766 = ['KU501215.1', 'AY632535.2']

Sequence = tools.readSeq("/Share/home/zhangqf8/lipan/virus/final_virus_genome.fa")
Shape = tools.loadicSHAPE("/Share/home/zhangqf8/lipan/virus/final_virus_icSHAPE.out")

#dg_list_59 = read_dg("/tmp/virus_paris/59/59.tab")
#dg_list_766 = read_dg("/tmp/virus_paris/766/766.tab")

reads_59 = read_paris_reads("/tmp/virus_paris/59/59.dg")
reads_766 = read_paris_reads("/tmp/virus_paris/766/766.dg")

domain_59 = virus.read_domain("/Share/home/zhangqf8/lipan/virus/Re-Domain/59_domain.txt")
domain_59.load_seq(Sequence[key59])
domain_59.load_shape(Shape[key59])

domain_766 = virus.read_domain("/Share/home/zhangqf8/lipan/virus/Re-Domain/766_domain.txt")
domain_766.load_seq(Sequence[key766])
domain_766.load_shape(Shape[key766])

seqdb_fa = "/Share/home/zhangqf8/lipan/virus/select_virus_genome/Flaviviridae-4256.fasta"



############
#  寻找PARIS/icSHAPE共支持的pair
############

dg_in_domain_59, dg_out_domain_59 = predict_duplex_dg(dg_list_59, domain_59, key59, reads_59)
dg_in_domain_766, dg_out_domain_766 = predict_duplex_dg(dg_list_766, domain_766, key766)


############
#  分别对每一个domain预测结构并对每一次stem进行评估
############

def make_single(ss, region_list):
    ss = list(ss)
    for i in range(len(region_list)):
        start, end = region_list[i]
        for j in range(start-1, end):
            ss[j] = '.'
    return "".join(ss)


def print_warning(domain, domain_index):
    RED = "\033[31m"
    DEF = "\033[39m"
    print RED
    if domain.at(domain_index-1)[0] % 10 != 1:
        print "Warning: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        print "Warning: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        print "Warning: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        print "Warning: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    print DEF

def correct_ss_info(seq, ss, shape, base_coor):
    base_coor = base_coor % 10
    if base_coor == 0:
        app_num = 9
    else:
        app_num = base_coor - 1
    seq = "I"*app_num + seq
    ss = "."*app_num + ss
    shape = ['NULL']*app_num + shape
    return seq, ss, shape


#####################################
#######   PRVABC59
#####################################

###################
# 1
###################


domain_index = 1

print_warning(domain_59, domain_index)
seq, shape = domain_59.sequence(domain_index-1), domain_59.shape(domain_index-1)
ss = structure.predictStructure(seq, shape, si=-0.4, sm=2.0)
prune_ss = prune_structure(ss, shape, reads_59, domain_59.at(domain_index-1)[0])

seq, prune_ss, shape = correct_ss_info(seq, prune_ss, shape, domain_59.at(domain_index-1)[0])
visual.Plot_RNAStructure_Shape(seq, prune_ss, shape, wait=False, mode='fill', title='59_domain_%s_after'%(domain_index, ))

handle_1 = covariation.submit_get_covariate_bp(seq, ss, sequence_db_fasta=seqdb_fa, cluster='Z-ZQF', cpu=26, clean=False, randID='59_domain_%s'%(domain_index, ), memory=100000, tmp_outdir="/Share/home/zhangqf8/.covariation")


###################
# 2
###################


domain_index = 2

print_warning(domain_59, domain_index)
seq, shape = domain_59.sequence(domain_index-1), domain_59.shape(domain_index-1)
ss = structure.predictStructure(seq, shape, si=-0.4, sm=2.0)
prune_ss = prune_structure(ss, shape, reads_59, domain_59.at(domain_index-1)[0])

handle_2 = covariation.submit_get_covariate_bp(seq, ss, sequence_db_fasta=seqdb_fa, cluster='Z-ZQF', cpu=26, clean=False, randID='59_domain_%s'%(domain_index, ), memory=100000, tmp_outdir="/Share/home/zhangqf8/.covariation")


seq, prune_ss, shape = correct_ss_info(seq, prune_ss, shape, domain_59.at(domain_index-1)[0])
visual.Plot_RNAStructure_Shape(seq, prune_ss, shape, wait=False, mode='fill', title='59_domain_%s_after'%(domain_index, ))


###################
# 3
###################


domain_index = 3

print_warning(domain_59, domain_index)
seq, shape = domain_59.sequence(domain_index-1), domain_59.shape(domain_index-1)
ss = structure.predictStructure(seq, shape, si=-0.4, sm=2.0)
prune_ss = prune_structure(ss, shape, reads_59, domain_59.at(domain_index-1)[0])

seq, prune_ss, shape = correct_ss_info(seq, prune_ss, shape, domain_59.at(domain_index-1)[0])
visual.Plot_RNAStructure_Shape(seq, prune_ss, shape, wait=False, mode='fill', title='59_domain_%s_after'%(domain_index, ))



###################
# 4
###################


domain_index = 4

print_warning(domain_59, domain_index)
seq, shape = domain_59.sequence(domain_index-1), domain_59.shape(domain_index-1)
ss = structure.predictStructure(seq, shape, si=-0.4, sm=2.0)
prune_ss = prune_structure(ss, shape, reads_59, domain_59.at(domain_index-1)[0])
prune_ss = make_single(prune_ss, [(1,10),(474,480)])

seq, prune_ss, shape = correct_ss_info(seq, prune_ss, shape, domain_59.at(domain_index-1)[0])
visual.Plot_RNAStructure_Shape(seq, prune_ss, shape, wait=False, mode='fill', title='59_domain_%s_after'%(domain_index, ))


###################
# 5
###################


domain_index = 5

print_warning(domain_59, domain_index)
seq, shape = domain_59.sequence(domain_index-1), domain_59.shape(domain_index-1)
ss = structure.predictStructure(seq, shape, si=-0.4, sm=2.0)
prune_ss = prune_structure(ss, shape, reads_59, domain_59.at(domain_index-1)[0])
#prune_ss = make_single(prune_ss, [(1,10),(474,480)])

seq, prune_ss, shape = correct_ss_info(seq, prune_ss, shape, domain_59.at(domain_index-1)[0])
visual.Plot_RNAStructure_Shape(seq, prune_ss, shape, wait=False, mode='fill', title='59_domain_%s_after'%(domain_index, ))




###################
# 6
###################


domain_index = 6

print_warning(domain_59, domain_index)
seq, shape = domain_59.sequence(domain_index-1), domain_59.shape(domain_index-1)
ss = structure.predictStructure(seq, shape, si=-0.4, sm=2.0)
prune_ss = prune_structure(ss, shape, reads_59, domain_59.at(domain_index-1)[0])
prune_ss = make_single(prune_ss, [(1,13),(542,581)])

seq, prune_ss, shape = correct_ss_info(seq, prune_ss, shape, domain_59.at(domain_index-1)[0])
visual.Plot_RNAStructure_Shape(seq, prune_ss, shape, wait=False, mode='fill', title='59_domain_%s_after'%(domain_index, ))


###################
# 7
###################


domain_index = 7

print_warning(domain_59, domain_index)
seq, shape = domain_59.sequence(domain_index-1), domain_59.shape(domain_index-1)
ss = structure.predictStructure(seq, shape, si=-0.4, sm=2.0)
prune_ss = prune_structure(ss, shape, reads_59, domain_59.at(domain_index-1)[0])
#prune_ss = make_single(prune_ss, [(1,13),(542,581)])

seq, prune_ss, shape = correct_ss_info(seq, prune_ss, shape, domain_59.at(domain_index-1)[0])
visual.Plot_RNAStructure_Shape(seq, prune_ss, shape, wait=False, mode='fill', title='59_domain_%s_after'%(domain_index, ))


###################
# 8
###################


domain_index = 8

print_warning(domain_59, domain_index)
seq, shape = domain_59.sequence(domain_index-1), domain_59.shape(domain_index-1)
ss = structure.predictStructure(seq, shape, si=-0.4, sm=2.0)
prune_ss = prune_structure(ss, shape, reads_59, domain_59.at(domain_index-1)[0])
#prune_ss = make_single(prune_ss, [(1,13),(542,581)])

seq, prune_ss, shape = correct_ss_info(seq, prune_ss, shape, domain_59.at(domain_index-1)[0])
visual.Plot_RNAStructure_Shape(seq, prune_ss, shape, wait=False, mode='fill', title='59_domain_%s_after'%(domain_index, ))


###################
# 9
###################


domain_index = 9

print_warning(domain_59, domain_index)
seq, shape = domain_59.sequence(domain_index-1), domain_59.shape(domain_index-1)
ss = structure.predictStructure(seq, shape, si=-0.4, sm=2.0)
prune_ss = prune_structure(ss, shape, reads_59, domain_59.at(domain_index-1)[0])
#prune_ss = make_single(prune_ss, [(1,13),(542,581)])

seq, prune_ss, shape = correct_ss_info(seq, prune_ss, shape, domain_59.at(domain_index-1)[0])
visual.Plot_RNAStructure_Shape(seq, prune_ss, shape, wait=False, mode='fill', title='59_domain_%s_after'%(domain_index, ))




###################
# 10
###################


domain_index = 10

print_warning(domain_59, domain_index)
seq, shape = domain_59.sequence(domain_index-1), domain_59.shape(domain_index-1)
ss = structure.predictStructure(seq, shape, si=-0.4, sm=2.0)
prune_ss = prune_structure(ss, shape, reads_59, domain_59.at(domain_index-1)[0])
#prune_ss = make_single(prune_ss, [(1,13),(542,581)])

seq, prune_ss, shape = correct_ss_info(seq, prune_ss, shape, domain_59.at(domain_index-1)[0])
visual.Plot_RNAStructure_Shape(seq, prune_ss, shape, wait=False, mode='fill', title='59_domain_%s_after'%(domain_index, ))




###################
# 11
###################


domain_index = 11

print_warning(domain_59, domain_index)
seq, shape = domain_59.sequence(domain_index-1), domain_59.shape(domain_index-1)
ss = structure.predictStructure(seq, shape, si=-0.4, sm=2.0)
prune_ss = prune_structure(ss, shape, reads_59, domain_59.at(domain_index-1)[0])
#prune_ss = make_single(prune_ss, [(1,13),(542,581)])

seq, prune_ss, shape = correct_ss_info(seq, prune_ss, shape, domain_59.at(domain_index-1)[0])
visual.Plot_RNAStructure_Shape(seq, prune_ss, shape, wait=False, mode='fill', title='59_domain_%s_after'%(domain_index, ))




###################
# 12
###################


domain_index = 12

print_warning(domain_59, domain_index)
seq, shape = domain_59.sequence(domain_index-1), domain_59.shape(domain_index-1)
ss = structure.predictStructure(seq, shape, si=-0.4, sm=2.0)
prune_ss = prune_structure(ss, shape, reads_59, domain_59.at(domain_index-1)[0])
#prune_ss = make_single(prune_ss, [(1,13),(542,581)])

seq, prune_ss, shape = correct_ss_info(seq, prune_ss, shape, domain_59.at(domain_index-1)[0])
visual.Plot_RNAStructure_Shape(seq, prune_ss, shape, wait=False, mode='fill', title='59_domain_%s_after'%(domain_index, ))




###################
# 13
###################


domain_index = 13

print_warning(domain_59, domain_index)
seq, shape = domain_59.sequence(domain_index-1), domain_59.shape(domain_index-1)
ss = structure.predictStructure(seq, shape, si=-0.4, sm=2.0)
prune_ss = prune_structure(ss, shape, reads_59, domain_59.at(domain_index-1)[0])
#prune_ss = make_single(prune_ss, [(1,13),(542,581)])

seq, prune_ss, shape = correct_ss_info(seq, prune_ss, shape, domain_59.at(domain_index-1)[0])
visual.Plot_RNAStructure_Shape(seq, prune_ss, shape, wait=False, mode='fill', title='59_domain_%s_after'%(domain_index, ))



###################
# 14
###################


domain_index = 14

print_warning(domain_59, domain_index)
seq, shape = domain_59.sequence(domain_index-1), domain_59.shape(domain_index-1)
ss = structure.predictStructure(seq, shape, si=-0.4, sm=2.0)
prune_ss = prune_structure(ss, shape, reads_59, domain_59.at(domain_index-1)[0])
prune_ss = make_single(prune_ss, [(20-5,31-5),(258-5,270-5)])

seq, prune_ss, shape = correct_ss_info(seq, prune_ss, shape, domain_59.at(domain_index-1)[0])
visual.Plot_RNAStructure_Shape(seq, prune_ss, shape, wait=False, mode='fill', title='59_domain_%s_after'%(domain_index, ))




###################
# 15
###################


domain_index = 15

print_warning(domain_59, domain_index)
seq, shape = domain_59.sequence(domain_index-1), domain_59.shape(domain_index-1)
ss = structure.predictStructure(seq, shape, si=-0.4, sm=2.0)
prune_ss = prune_structure(ss, shape, reads_59, domain_59.at(domain_index-1)[0])
prune_ss = make_single(prune_ss, [(643,669)])

seq, prune_ss, shape = correct_ss_info(seq, prune_ss, shape, domain_59.at(domain_index-1)[0])
visual.Plot_RNAStructure_Shape(seq, prune_ss, shape, wait=False, mode='fill', title='59_domain_%s_after'%(domain_index, ))



###################
# 16
###################


domain_index = 16

print_warning(domain_59, domain_index)
seq, shape = domain_59.sequence(domain_index-1), domain_59.shape(domain_index-1)
ss = structure.predictStructure(seq, shape, si=-0.4, sm=2.0)
prune_ss = prune_structure(ss, shape, reads_59, domain_59.at(domain_index-1)[0])
prune_ss = make_single(prune_ss, [(10-5,16-5), (398-5, 401-5)])

seq, prune_ss, shape = correct_ss_info(seq, prune_ss, shape, domain_59.at(domain_index-1)[0])
visual.Plot_RNAStructure_Shape(seq, prune_ss, shape, wait=False, mode='fill', title='59_domain_%s_after'%(domain_index, ))


###################
# 17
###################


domain_index = 17

print_warning(domain_59, domain_index)
seq, shape = domain_59.sequence(domain_index-1), domain_59.shape(domain_index-1)
ss = structure.predictStructure(seq, shape, si=-0.4, sm=2.0)
prune_ss = prune_structure(ss, shape, reads_59, domain_59.at(domain_index-1)[0])
#prune_ss = make_single(prune_ss, [(10-5,16-5), (398-5, 401-5)])

seq, prune_ss, shape = correct_ss_info(seq, prune_ss, shape, domain_59.at(domain_index-1)[0])
visual.Plot_RNAStructure_Shape(seq, prune_ss, shape, wait=False, mode='fill', title='59_domain_%s_after'%(domain_index, ))




###################
# 18
###################


domain_index = 18

print_warning(domain_59, domain_index)
seq, shape = domain_59.sequence(domain_index-1), domain_59.shape(domain_index-1)
ss = structure.predictStructure(seq, shape, si=-0.4, sm=2.0)
prune_ss = prune_structure(ss, shape, reads_59, domain_59.at(domain_index-1)[0])
#prune_ss = make_single(prune_ss, [(10-5,16-5), (398-5, 401-5)])

seq, prune_ss, shape = correct_ss_info(seq, prune_ss, shape, domain_59.at(domain_index-1)[0])
visual.Plot_RNAStructure_Shape(seq, prune_ss, shape, wait=False, mode='fill', title='59_domain_%s_after'%(domain_index, ))




###################
# 19
###################


domain_index = 19

print_warning(domain_59, domain_index)
seq, shape = domain_59.sequence(domain_index-1), domain_59.shape(domain_index-1)
ss = structure.predictStructure(seq, shape, si=-0.4, sm=2.0)
prune_ss = prune_structure(ss, shape, reads_59, domain_59.at(domain_index-1)[0])
#prune_ss = make_single(prune_ss, [(10-5,16-5), (398-5, 401-5)])

seq, prune_ss, shape = correct_ss_info(seq, prune_ss, shape, domain_59.at(domain_index-1)[0])
visual.Plot_RNAStructure_Shape(seq, prune_ss, shape, wait=False, mode='fill', title='59_domain_%s_after'%(domain_index, ))



###################
# 20
###################


domain_index = 20

print_warning(domain_59, domain_index)
seq, shape = domain_59.sequence(domain_index-1), domain_59.shape(domain_index-1)
ss = structure.predictStructure(seq, shape, si=-0.4, sm=2.0)
prune_ss = prune_structure(ss, shape, reads_59, domain_59.at(domain_index-1)[0])
#prune_ss = make_single(prune_ss, [(10-5,16-5), (398-5, 401-5)])

seq, prune_ss, shape = correct_ss_info(seq, prune_ss, shape, domain_59.at(domain_index-1)[0])
visual.Plot_RNAStructure_Shape(seq, prune_ss, shape, wait=False, mode='fill', title='59_domain_%s_after'%(domain_index, ))




###################
# 21
###################


domain_index = 21

print_warning(domain_59, domain_index)
seq, shape = domain_59.sequence(domain_index-1), domain_59.shape(domain_index-1)
ss = structure.predictStructure(seq, shape, si=-0.4, sm=2.0)
prune_ss = prune_structure(ss, shape, reads_59, domain_59.at(domain_index-1)[0])
#prune_ss = make_single(prune_ss, [(10-5,16-5), (398-5, 401-5)])

seq, prune_ss, shape = correct_ss_info(seq, prune_ss, shape, domain_59.at(domain_index-1)[0])
visual.Plot_RNAStructure_Shape(seq, prune_ss, shape, wait=False, mode='fill', title='59_domain_%s_after'%(domain_index, ))




###################
# 22
###################


domain_index = 22

print_warning(domain_59, domain_index)
seq, shape = domain_59.sequence(domain_index-1), domain_59.shape(domain_index-1)
ss = structure.predictStructure(seq, shape, si=-0.4, sm=2.0)
prune_ss = prune_structure(ss, shape, reads_59, domain_59.at(domain_index-1)[0])
prune_ss = make_single(prune_ss, [(60-5,70-5), (200-5, 202-5)])

seq, prune_ss, shape = correct_ss_info(seq, prune_ss, shape, domain_59.at(domain_index-1)[0])
visual.Plot_RNAStructure_Shape(seq, prune_ss, shape, wait=False, mode='fill', title='59_domain_%s_after'%(domain_index, ))



###################
# 23
###################


domain_index = 23

print_warning(domain_59, domain_index)
seq, shape = domain_59.sequence(domain_index-1), domain_59.shape(domain_index-1)
ss = structure.predictStructure(seq, shape, si=-0.4, sm=2.0)
prune_ss = prune_structure(ss, shape, reads_59, domain_59.at(domain_index-1)[0])
#prune_ss = make_single(prune_ss, [(60-5,70-5), (200-5, 202-5)])

seq, prune_ss, shape = correct_ss_info(seq, prune_ss, shape, domain_59.at(domain_index-1)[0])
visual.Plot_RNAStructure_Shape(seq, prune_ss, shape, wait=False, mode='fill', title='59_domain_%s_after'%(domain_index, ))



###################
# 24
###################


domain_index = 24

print_warning(domain_59, domain_index)
seq, shape = domain_59.sequence(domain_index-1), domain_59.shape(domain_index-1)
ss = structure.predictStructure(seq, shape, si=-0.4, sm=2.0)
prune_ss = prune_structure(ss, shape, reads_59, domain_59.at(domain_index-1)[0])
prune_ss = make_single(prune_ss, [(10,12), (137,139)])

seq, prune_ss, shape = correct_ss_info(seq, prune_ss, shape, domain_59.at(domain_index-1)[0])
visual.Plot_RNAStructure_Shape(seq, prune_ss, shape, wait=False, mode='fill', title='59_domain_%s_after'%(domain_index, ))




###################
# Predict full structure
###################

## 59

full_ss = ""
for domain_index in range(1,25):
    seq, shape = domain_59.sequence(domain_index-1), domain_59.shape(domain_index-1)
    ss = structure.predictStructure(seq, shape, si=-0.4, sm=2.0)
    full_ss += ss

### evaluate with SHAPE reactivities
evaluate_with_shape(full_ss, Shape[key59])
plt.savefig("/tmp/59_ss.pdf")
plt.close()

### evaluate with PARIS data
OUT = open("/tmp/59_paris.txt", 'w')
paris_support = evaluate_with_paris(full_ss, reads_59, base_coor=1, OUT=OUT)
OUT.close()
plt.savefig("/tmp/59_paris.pdf")
plt.show()

###################
# Call covariation base pairs
###################

################## 59

handle_list = []
for domain_index in range(1,25):
    seq, shape = domain_59.sequence(domain_index-1), domain_59.shape(domain_index-1)
    ss = structure.predictStructure(seq, shape, si=-0.4, sm=2.0)
    handle = covariation.submit_get_covariate_bp(seq, ss, sequence_db_fasta=seqdb_fa, cluster='Z-ZQF', cpu=20, clean=False, randID='59_domain_%s'%(domain_index, ), memory=100000, tmp_outdir="/Share/home/zhangqf8/.covariation")
    handle_list.append(handle)

for handle in handle_list:
    handle.wait()

covariation_list = []
for domain_index in range(1, len(handle_list)+1):
    seq, shape = domain_59.sequence(domain_index-1), domain_59.shape(domain_index-1)
    ss = structure.predictStructure(seq, shape, si=-0.4, sm=2.0)
    domain_start = domain_59.at(domain_index-1)[0]
    cv_list = handle_list[domain_index-1].get_output()['R-scape']
    prune_ss = prune_structure(ss, shape, reads_59, domain_59.at(domain_index-1)[0])
    visual.Plot_RNAStructure_highlight(seq, prune_ss, covariation.bp2blist(cv_list))
    raw_input("next")
    for bp in cv_list:
        covariation_list.append((bp[0]+domain_start-1, bp[1]+domain_start-1))

covariation_list = [(4, 68), (7, 65), (8, 64), (12, 59), (13, 58), (20, 42), 
                    (21, 41), (22, 40), (23, 39), (24, 38), (25, 37), (27, 36), 
                    (28, 35), (29, 34), (45, 52), (46, 51), (88, 101), (103, 116), 
                    (127, 144), (128, 143), (129, 142), (130, 141), (133, 138), (167, 206), 
                    (168, 205), (170, 203), (171, 202), (172, 201), (177, 196), (178, 195), 
                    (179, 194), (180, 193), (181, 192), (182, 191), (184, 189), (214, 246), 
                    (258, 338), (259, 337), (265, 331), (340, 345), (526, 998), (527, 997), 
                    (528, 996), (529, 995), (646, 765), (697, 711), (714, 735), (715, 734), 
                    (767, 991), (772, 814), (837, 980), (1024, 1221), (1237, 1247), (1381, 1388), 
                    (1399, 1441), (1445, 1455), (1530, 1800), (1531, 1799), (1532, 1798), (1807, 1831), 
                    (1811, 1827), (1812, 1826), (1813, 1825), (1997, 2108), (2091, 2098), (2155, 2509), 
                    (2525, 2859), (2602, 2707), (2657, 2696), (2866, 3067), (2885, 3046), (3222, 3306), 
                    (3225, 3302), (3227, 3300), (3348, 3358), (3349, 3357), (3516, 3717), (3517, 3716), 
                    (3518, 3715), (3519, 3714), (3520, 3712), (3521, 3711), (3522, 3710), (3523, 3709),
                    (3524, 3708), (3525, 3707), (3527, 3705), (3530, 3702), (3531, 3701), (3532, 3700), 
                    (3537, 3694), (3539, 3692), (3540, 3686), (3541, 3685), (3543, 3682), (3544, 3681), 
                    (3866, 3948), (3871, 3943), (4452, 4619), (4453, 4618), (4658, 5443), (4782, 4791), 
                    (4873, 4930), (5073, 5279), (5080, 5272), (5081, 5271), (5089, 5263), (5090, 5262), 
                    (5103, 5246), (5106, 5244), (5109, 5241), (5116, 5239), (5124, 5230), (5128, 5226), 
                    (5133, 5220), (5138, 5215), (5160, 5196), (5513, 5532), (5694, 5703), (5783, 5928), 
                    (5874, 5910), (5960, 6350), (6365, 6387), (6452, 6468), (6491, 6506), (6512, 6678), 
                    (6514, 6676), (6515, 6675), (6516, 6674), (6528, 6660), (6590, 6597), (6685, 6802), 
                    (6694, 6785), (6708, 6771), (6710, 6769), (6850, 6899), (6871, 6880), (6925, 7074), 
                    (6987, 6993), (6988, 6992), (7091, 7371), (7118, 7125), (7180, 7193), (7182, 7191), 
                    (7301, 7314), (7384, 7390), (7385, 7389), (7401, 7487), (7402, 7486), (7403, 7485), 
                    (7489, 7567), (7492, 7564), (7498, 7519), (7499, 7518), (7502, 7515), (7580, 7634), 
                    (7583, 7631), (7650, 7789), (7671, 7741), (7861, 7945), (7864, 7942), (8349, 8364), 
                    (8350, 8363), (8369, 8719), (8459, 8470), (8460, 8469), (8484, 8502), (8486, 8500), 
                    (8782, 8999), (8784, 8997), (8868, 8910), (9110, 9363), (9111, 9362), (9565, 9702), 
                    (9566, 9701), (9573, 9620), (10218, 10356), (10219, 10355), (10236, 10287), (10252, 10270), 
                    (10335, 10341), (10370, 10382), (10371, 10381), (10372, 10380), (10373, 10379), 
                    (10374, 10378), (10396, 10727), (10397, 10726), (10401, 10528), (10402, 10527), 
                    (10415, 10435), (10416, 10434), (10417, 10433), (10418, 10432), (10419, 10431), 
                    (10450, 10465), (10451, 10464), (10452, 10463), (10453, 10462), (10454, 10461), 
                    (10482, 10523), (10485, 10520), (10498, 10514), (10499, 10513), (10500, 10512), 
                    (10531, 10546), (10532, 10545), (10533, 10544), (10534, 10543), (10535, 10542), 
                    (10576, 10608), (10577, 10607), (10616, 10685), (10617, 10684), (10618, 10683), 
                    (10619, 10682), (10620, 10681), (10624, 10656), (10628, 10652), (10629, 10651), 
                    (10630, 10650), (10631, 10649), (10632, 10648), (10633, 10647), (10634, 10646), 
                    (10711, 10722), (10712, 10721), (10729, 10802), (10730, 10801), (10731, 10800), 
                    (10733, 10798), (10737, 10794), (10738, 10793), (10742, 10790), (10743, 10789), 
                    (10745, 10787), (10746, 10786), (10747, 10785), (10748, 10784), (10749, 10783), 
                    (10750, 10782), (10751, 10781), (10752, 10780), (10753, 10779), (10754, 10778), (10757, 10769)]


#####################################
#######   766
#####################################

###################
# 766 1
###################

domain_index = 1

print_warning(domain_766, domain_index)
seq, shape = domain_766.sequence(domain_index-1), domain_766.shape(domain_index-1)
ss = structure.predictStructure(seq, shape, si=-0.2, sm=0.6)
prune_ss = prune_structure(ss, shape, reads_766, domain_766.at(domain_index-1)[0])

seq, prune_ss, shape = correct_ss_info(seq, prune_ss, shape, domain_766.at(domain_index-1)[0])
visual.Plot_RNAStructure_Shape(seq, prune_ss, shape, wait=False, mode='fill', title='766_domain_%s_after'%(domain_index, ))


###################
# 766 2
###################

domain_index = 2

print_warning(domain_766, domain_index)
seq, shape = domain_766.sequence(domain_index-1), domain_766.shape(domain_index-1)
ss = structure.predictStructure(seq, shape, si=-0.2, sm=0.6)
prune_ss = prune_structure(ss, shape, reads_766, domain_766.at(domain_index-1)[0])

seq, prune_ss, shape = correct_ss_info(seq, prune_ss, shape, domain_766.at(domain_index-1)[0])
visual.Plot_RNAStructure_Shape(seq, prune_ss, shape, wait=False, mode='fill', title='766_domain_%s_after'%(domain_index, ))


###################
# 766 3
###################

domain_index = 3

print_warning(domain_766, domain_index)
seq, shape = domain_766.sequence(domain_index-1), domain_766.shape(domain_index-1)
ss = structure.predictStructure(seq, shape, si=-0.2, sm=0.6)
prune_ss = prune_structure(ss, shape, reads_766, domain_766.at(domain_index-1)[0])

seq, prune_ss, shape = correct_ss_info(seq, prune_ss, shape, domain_766.at(domain_index-1)[0])
visual.Plot_RNAStructure_Shape(seq, prune_ss, shape, wait=False, mode='fill', title='766_domain_%s_after'%(domain_index, ))


###################
# 766 4
###################

domain_index = 4

print_warning(domain_766, domain_index)
seq, shape = domain_766.sequence(domain_index-1), domain_766.shape(domain_index-1)
ss = structure.predictStructure(seq, shape, si=-0.2, sm=0.6)
prune_ss = prune_structure(ss, shape, reads_766, domain_766.at(domain_index-1)[0])

seq, prune_ss, shape = correct_ss_info(seq, prune_ss, shape, domain_766.at(domain_index-1)[0])
visual.Plot_RNAStructure_Shape(seq, prune_ss, shape, wait=False, mode='fill', title='766_domain_%s_after'%(domain_index, ))



###################
# 766 5
###################

domain_index = 5

print_warning(domain_766, domain_index)
seq, shape = domain_766.sequence(domain_index-1), domain_766.shape(domain_index-1)
ss = structure.predictStructure(seq, shape, si=-0.2, sm=0.6)
prune_ss = prune_structure(ss, shape, reads_766, domain_766.at(domain_index-1)[0])

seq, prune_ss, shape = correct_ss_info(seq, prune_ss, shape, domain_766.at(domain_index-1)[0])
visual.Plot_RNAStructure_Shape(seq, prune_ss, shape, wait=False, mode='fill', title='766_domain_%s_after'%(domain_index, ))



###################
# 766 6
###################

domain_index = 6

print_warning(domain_766, domain_index)
seq, shape = domain_766.sequence(domain_index-1), domain_766.shape(domain_index-1)
ss = structure.predictStructure(seq, shape, si=-0.2, sm=0.6)
prune_ss = prune_structure(ss, shape, reads_766, domain_766.at(domain_index-1)[0])

seq, prune_ss, shape = correct_ss_info(seq, prune_ss, shape, domain_766.at(domain_index-1)[0])
visual.Plot_RNAStructure_Shape(seq, prune_ss, shape, wait=False, mode='fill', title='766_domain_%s_after'%(domain_index, ))




###################
# 766 7
###################

domain_index = 7

print_warning(domain_766, domain_index)
seq, shape = domain_766.sequence(domain_index-1), domain_766.shape(domain_index-1)
ss = structure.predictStructure(seq, shape, si=-0.2, sm=0.6)
prune_ss = prune_structure(ss, shape, reads_766, domain_766.at(domain_index-1)[0])
prune_ss = make_single(prune_ss, [(10,20), (420,430)])

seq, prune_ss, shape = correct_ss_info(seq, prune_ss, shape, domain_766.at(domain_index-1)[0])
visual.Plot_RNAStructure_Shape(seq, prune_ss, shape, wait=False, mode='fill', title='766_domain_%s_after'%(domain_index, ))


###################
# 766 8
###################

domain_index = 8

print_warning(domain_766, domain_index)
seq, shape = domain_766.sequence(domain_index-1), domain_766.shape(domain_index-1)
ss = structure.predictStructure(seq, shape, si=-0.2, sm=0.6)
prune_ss = prune_structure(ss, shape, reads_766, domain_766.at(domain_index-1)[0])
#prune_ss = make_single(prune_ss, [(10,20), (420,430)])

seq, prune_ss, shape = correct_ss_info(seq, prune_ss, shape, domain_766.at(domain_index-1)[0])
visual.Plot_RNAStructure_Shape(seq, prune_ss, shape, wait=False, mode='fill', title='766_domain_%s_after'%(domain_index, ))


###################
# 766 9
###################

domain_index = 9

print_warning(domain_766, domain_index)
seq, shape = domain_766.sequence(domain_index-1), domain_766.shape(domain_index-1)
ss = structure.predictStructure(seq, shape, si=-0.2, sm=0.6)
prune_ss = prune_structure(ss, shape, reads_766, domain_766.at(domain_index-1)[0])
#prune_ss = make_single(prune_ss, [(10,20), (420,430)])

seq, prune_ss, shape = correct_ss_info(seq, prune_ss, shape, domain_766.at(domain_index-1)[0])
visual.Plot_RNAStructure_Shape(seq, prune_ss, shape, wait=False, mode='fill', title='766_domain_%s_after'%(domain_index, ))




###################
# 766 10
###################

domain_index = 10

print_warning(domain_766, domain_index)
seq, shape = domain_766.sequence(domain_index-1), domain_766.shape(domain_index-1)
ss = structure.predictStructure(seq, shape, si=-0.2, sm=0.6)
prune_ss = prune_structure(ss, shape, reads_766, domain_766.at(domain_index-1)[0])
#prune_ss = make_single(prune_ss, [(10,20), (420,430)])

seq, prune_ss, shape = correct_ss_info(seq, prune_ss, shape, domain_766.at(domain_index-1)[0])
visual.Plot_RNAStructure_Shape(seq, prune_ss, shape, wait=False, mode='fill', title='766_domain_%s_after'%(domain_index, ))




###################
# 766 11
###################

domain_index = 11

print_warning(domain_766, domain_index)
seq, shape = domain_766.sequence(domain_index-1), domain_766.shape(domain_index-1)
ss = structure.predictStructure(seq, shape, si=-0.2, sm=0.6)
prune_ss = prune_structure(ss, shape, reads_766, domain_766.at(domain_index-1)[0])
prune_ss = make_single(prune_ss, [(1, 38)])

seq, prune_ss, shape = correct_ss_info(seq, prune_ss, shape, domain_766.at(domain_index-1)[0])
visual.Plot_RNAStructure_Shape(seq, prune_ss, shape, wait=False, mode='fill', title='766_domain_%s_after'%(domain_index, ))




###################
# 766 12
###################

domain_index = 12

print_warning(domain_766, domain_index)
seq, shape = domain_766.sequence(domain_index-1), domain_766.shape(domain_index-1)
ss = structure.predictStructure(seq, shape, si=-0.2, sm=0.6)
prune_ss = prune_structure(ss, shape, reads_766, domain_766.at(domain_index-1)[0])
#prune_ss = make_single(prune_ss, [(1, 38)])

seq, prune_ss, shape = correct_ss_info(seq, prune_ss, shape, domain_766.at(domain_index-1)[0])
visual.Plot_RNAStructure_Shape(seq, prune_ss, shape, wait=False, mode='fill', title='766_domain_%s_after'%(domain_index, ))




###################
# 766 13
###################

domain_index = 13

print_warning(domain_766, domain_index)
seq, shape = domain_766.sequence(domain_index-1), domain_766.shape(domain_index-1)
ss = structure.predictStructure(seq, shape, si=-0.2, sm=0.6)
prune_ss = prune_structure(ss, shape, reads_766, domain_766.at(domain_index-1)[0])
#prune_ss = make_single(prune_ss, [(1, 38)])

seq, prune_ss, shape = correct_ss_info(seq, prune_ss, shape, domain_766.at(domain_index-1)[0])
visual.Plot_RNAStructure_Shape(seq, prune_ss, shape, wait=False, mode='fill', title='766_domain_%s_after'%(domain_index, ))



###################
# 766 14
###################

domain_index = 14

print_warning(domain_766, domain_index)
seq, shape = domain_766.sequence(domain_index-1), domain_766.shape(domain_index-1)
ss = structure.predictStructure(seq, shape, si=-0.2, sm=0.6)
prune_ss = prune_structure(ss, shape, reads_766, domain_766.at(domain_index-1)[0])
#prune_ss = make_single(prune_ss, [(1, 38)])

seq, prune_ss, shape = correct_ss_info(seq, prune_ss, shape, domain_766.at(domain_index-1)[0])
visual.Plot_RNAStructure_Shape(seq, prune_ss, shape, wait=False, mode='fill', title='766_domain_%s_after'%(domain_index, ))



###################
# 766 15
###################

domain_index = 15

print_warning(domain_766, domain_index)
seq, shape = domain_766.sequence(domain_index-1), domain_766.shape(domain_index-1)
ss = structure.predictStructure(seq, shape, si=-0.2, sm=0.6)
prune_ss = prune_structure(ss, shape, reads_766, domain_766.at(domain_index-1)[0])
#prune_ss = make_single(prune_ss, [(1, 38)])

seq, prune_ss, shape = correct_ss_info(seq, prune_ss, shape, domain_766.at(domain_index-1)[0])
visual.Plot_RNAStructure_Shape(seq, prune_ss, shape, wait=False, mode='fill', title='766_domain_%s_after'%(domain_index, ))


###################
# 766 16
###################

domain_index = 16

print_warning(domain_766, domain_index)
seq, shape = domain_766.sequence(domain_index-1), domain_766.shape(domain_index-1)
ss = structure.predictStructure(seq, shape, si=-0.2, sm=0.6)
prune_ss = prune_structure(ss, shape, reads_766, domain_766.at(domain_index-1)[0])
prune_ss = make_single(prune_ss, [(166-2, 168-2), (105-2, 107-2)])

seq, prune_ss, shape = correct_ss_info(seq, prune_ss, shape, domain_766.at(domain_index-1)[0])
visual.Plot_RNAStructure_Shape(seq, prune_ss, shape, wait=False, mode='fill', title='766_domain_%s_after'%(domain_index, ))



###################
# 766 17
###################

domain_index = 17

print_warning(domain_766, domain_index)
seq, shape = domain_766.sequence(domain_index-1), domain_766.shape(domain_index-1)
ss = structure.predictStructure(seq, shape, si=-0.2, sm=0.6)
prune_ss = prune_structure(ss, shape, reads_766, domain_766.at(domain_index-1)[0])
#prune_ss = make_single(prune_ss, [(166-2, 168-2), (105-2, 107-2)])

seq, prune_ss, shape = correct_ss_info(seq, prune_ss, shape, domain_766.at(domain_index-1)[0])
visual.Plot_RNAStructure_Shape(seq, prune_ss, shape, wait=False, mode='fill', title='766_domain_%s_after'%(domain_index, ))




###################
# 766 18
###################

domain_index = 18

print_warning(domain_766, domain_index)
seq, shape = domain_766.sequence(domain_index-1), domain_766.shape(domain_index-1)
ss = structure.predictStructure(seq, shape, si=-0.2, sm=0.6)
prune_ss = prune_structure(ss, shape, reads_766, domain_766.at(domain_index-1)[0])
prune_ss = make_single(prune_ss, [(160, 169), (585,587)])

seq, prune_ss, shape = correct_ss_info(seq, prune_ss, shape, domain_766.at(domain_index-1)[0])
visual.Plot_RNAStructure_Shape(seq, prune_ss, shape, wait=False, mode='fill', title='766_domain_%s_after'%(domain_index, ))



###################
# 766 19
###################

domain_index = 19

print_warning(domain_766, domain_index)
seq, shape = domain_766.sequence(domain_index-1), domain_766.shape(domain_index-1)
ss = structure.predictStructure(seq, shape, si=-0.2, sm=0.6)
prune_ss = prune_structure(ss, shape, reads_766, domain_766.at(domain_index-1)[0])
#prune_ss = make_single(prune_ss, [(160, 169), (585,587)])

seq, prune_ss, shape = correct_ss_info(seq, prune_ss, shape, domain_766.at(domain_index-1)[0])
visual.Plot_RNAStructure_Shape(seq, prune_ss, shape, wait=False, mode='fill', title='766_domain_%s_after'%(domain_index, ))

###################
# 766 20
###################

domain_index = 20

print_warning(domain_766, domain_index)
seq, shape = domain_766.sequence(domain_index-1), domain_766.shape(domain_index-1)
ss = structure.predictStructure(seq, shape, si=-0.2, sm=0.6)
prune_ss = prune_structure(ss, shape, reads_766, domain_766.at(domain_index-1)[0])
#prune_ss = make_single(prune_ss, [(160, 169), (585,587)])

seq, prune_ss, shape = correct_ss_info(seq, prune_ss, shape, domain_766.at(domain_index-1)[0])
visual.Plot_RNAStructure_Shape(seq, prune_ss, shape, wait=False, mode='fill', title='766_domain_%s_after'%(domain_index, ))


###################
# 766 21
###################

domain_index = 21

print_warning(domain_766, domain_index)
seq, shape = domain_766.sequence(domain_index-1), domain_766.shape(domain_index-1)
ss = structure.predictStructure(seq, shape, si=-0.2, sm=0.6)
prune_ss = prune_structure(ss, shape, reads_766, domain_766.at(domain_index-1)[0])
#prune_ss = make_single(prune_ss, [(160, 169), (585,587)])

seq, prune_ss, shape = correct_ss_info(seq, prune_ss, shape, domain_766.at(domain_index-1)[0])
visual.Plot_RNAStructure_Shape(seq, prune_ss, shape, wait=False, mode='fill', title='766_domain_%s_after'%(domain_index, ))



###################
# 766 22
###################

domain_index = 22

print_warning(domain_766, domain_index)
seq, shape = domain_766.sequence(domain_index-1), domain_766.shape(domain_index-1)
ss = structure.predictStructure(seq, shape, si=-0.2, sm=0.6)
prune_ss = prune_structure(ss, shape, reads_766, domain_766.at(domain_index-1)[0])
#prune_ss = make_single(prune_ss, [(160, 169), (585,587)])

seq, prune_ss, shape = correct_ss_info(seq, prune_ss, shape, domain_766.at(domain_index-1)[0])
visual.Plot_RNAStructure_Shape(seq, prune_ss, shape, wait=False, mode='fill', title='766_domain_%s_after'%(domain_index, ))


###################
# 766 23
###################

domain_index = 23

print_warning(domain_766, domain_index)
seq, shape = domain_766.sequence(domain_index-1), domain_766.shape(domain_index-1)
ss = structure.predictStructure(seq, shape, si=-0.2, sm=0.6)
prune_ss = prune_structure(ss, shape, reads_766, domain_766.at(domain_index-1)[0])
prune_ss = make_single(prune_ss, [(9-8, 16-8), (106-8, 113-8), (136-8, 142-8), (314-8, 320-8)])

seq, prune_ss, shape = correct_ss_info(seq, prune_ss, shape, domain_766.at(domain_index-1)[0])
visual.Plot_RNAStructure_Shape(seq, prune_ss, shape, wait=False, mode='fill', title='766_domain_%s_after'%(domain_index, ))


#### Predict the full structure

full_ss = ""
for domain_index in range(1,24):
    seq, shape = domain_766.sequence(domain_index-1), domain_766.shape(domain_index-1)
    ss = structure.predictStructure(seq, shape, si=-0.2, sm=0.6)
    full_ss += ss

### evaluate with SHAPE reactivities
evaluate_with_shape(full_ss, Shape[key766])
plt.savefig("/tmp/766.pdf")
plt.show()
plt.close()

### evaluate with PARIS data
OUT = open("/tmp/766_paris.txt", 'w')
paris_support = evaluate_with_paris(full_ss, reads_766, base_coor=1, OUT=OUT)
OUT.close()
plt.savefig("/tmp/766_paris.pdf")
plt.show()


### Show those duplex with minimun support 20
for domain_index in range(1,24):
    seq, shape = domain_766.sequence(domain_index-1), domain_766.shape(domain_index-1)
    ss = structure.predictStructure(seq, shape, si=-0.2, sm=0.6)
    prune_ss = prune_structure(ss, shape, reads_766, domain_766.at(domain_index-1)[0])
    start, end = domain_766.at(domain_index-1)
    highlight_list = []
    for stem in paris_support:
        if stem[4] >= 20:
            if start <= stem[0] < stem[3] <= end:
                for idx in range(stem[0], stem[1]+1):
                    highlight_list.append( idx-start+1 )
                for idx in range(stem[2], stem[3]+1):
                    highlight_list.append( idx-start+1 )
                print stem
    visual.Plot_RNAStructure_highlight(seq, prune_ss, highlight_list)
    raw_input("Next")


###################
# Call covariation base pairs
###################

################## 766

handle_list = []
for domain_index in range(1,24):
    seq, shape = domain_766.sequence(domain_index-1), domain_766.shape(domain_index-1)
    ss = structure.predictStructure(seq, shape, si=-0.2, sm=0.6)
    handle = covariation.submit_get_covariate_bp(seq, ss, sequence_db_fasta=seqdb_fa, cluster='Z-ZQF', cpu=20, clean=False, randID='766_domain_%s'%(domain_index, ), memory=100000, tmp_outdir="/Share/home/zhangqf8/.covariation")
    handle_list.append(handle)

for handle in handle_list:
    handle.wait()

covariation_list = []
for domain_index in range(1, len(handle_list)+1):
    seq, shape = domain_766.sequence(domain_index-1), domain_766.shape(domain_index-1)
    ss = structure.predictStructure(seq, shape, si=-0.2, sm=0.6)
    domain_start = domain_766.at(domain_index-1)[0]
    cv_list = handle_list[domain_index-1].get_output()['R-scape']
    prune_ss = prune_structure(ss, shape, reads_766, domain_766.at(domain_index-1)[0])
    visual.Plot_RNAStructure_highlight(seq, prune_ss, covariation.bp2blist(cv_list))
    raw_input("next")
    for bp in cv_list:
        covariation_list.append((bp[0]+domain_start-1, bp[1]+domain_start-1))



covariation_list = [(5, 69), (6, 68), (8, 66), (9, 65), (11, 62), (13, 60), (14, 59), (20, 44), 
                    (21, 43), (26, 38), (28, 37), (29, 36), (46, 53), (47, 52), (126, 145), 
                    (127, 144), (128, 143), (129, 142), (130, 141), (132, 139), (133, 138), 
                    (259, 337), (593, 943), (645, 766), (652, 759), (696, 713), (719, 728), 
                    (735, 757), (738, 754), (772, 784), (774, 781), (802, 876), (816, 838), 
                    (1021, 1225), (1024, 1221), (1027, 1218), (1029, 1216), (1033, 1208), (1108, 1203), 
                    (1111, 1200), (1233, 1251), (1234, 1250), (1237, 1247), (1238, 1246), (1261, 1692), 
                    (1262, 1691), (1263, 1690), (1276, 1677), (1335, 1351), (1336, 1350), (1341, 1345), 
                    (1381, 1665), (1382, 1664), (1387, 1659), (1390, 1656), (1394, 1408), (1395, 1407), 
                    (1397, 1405), (1410, 1650), (1411, 1649), (1413, 1648), (1418, 1643), (1420, 1641), 
                    (1422, 1640), (1423, 1639), (1425, 1628), (1426, 1627), (1428, 1625), (1442, 1616), 
                    (1443, 1615), (1452, 1568), (1453, 1567), (1454, 1566), (1705, 1804), (1706, 1803), 
                    (1817, 2007), (1906, 1954), (2052, 2465), (2053, 2464), (2054, 2463), (2055, 2462), 
                    (2056, 2461), (2083, 2128), (2093, 2118), (2138, 2156), (2198, 2365), (2292, 2348), 
                    (2404, 2458), (2408, 2453), (2409, 2452), (2410, 2451), (2490, 2527), (2531, 2556), 
                    (2562, 2587), (2846, 3061), (2854, 3055), (2862, 3045), (2873, 3034), (2874, 3032), 
                    (2879, 3027), (2880, 3026), (2883, 3024), (2884, 3023), (2885, 3022), (2891, 3013), 
                    (2899, 2907), (2900, 2906), (2901, 2905), (2910, 2985), (2913, 2982), (2914, 2981), 
                    (2915, 2980), (3063, 3324), (3110, 3231), (3235, 3274), (3335, 3437), (3347, 3425), 
                    (3365, 3406), (3439, 3677), (3441, 3675), (3446, 3670), (3499, 3564), (3526, 3545), 
                    (3532, 3541), (3533, 3540), (3576, 3643), (3590, 3630), (4050, 4131), (4422, 4442), 
                    (4423, 4441), (4424, 4439), (4425, 4438), (4427, 4436), (4428, 4435), (4446, 4798), 
                    (4450, 4507), (4452, 4505), (4509, 4518), (4510, 4517), (4526, 4660), (4527, 4659), 
                    (4532, 4654), (4533, 4653), (4602, 4626), (4805, 5479), (4822, 5461), (4914, 5143), 
                    (4916, 5141), (4920, 5137), (4921, 5136), (4928, 5129), (4929, 5128), (4932, 4999), 
                    (5113, 5124), (5114, 5123), (5115, 5122), (5116, 5121), (5939, 6181), (5948, 6172), 
                    (6010, 6107), (6021, 6097), (6116, 6160), (6362, 6369), (6422, 6458), (6423, 6457), 
                    (6476, 6537), (6479, 6494), (6482, 6491), (6634, 6839), (6641, 6647), (6664, 6688), 
                    (6706, 6721), (6871, 7327), (6933, 6946), (6934, 6945), (6952, 7056), (6953, 7055), 
                    (6958, 7048), (6960, 7046), (6961, 7045), (7068, 7085), (7070, 7083), (7072, 7081), 
                    (7073, 7080), (7074, 7079), (7090, 7132), (7106, 7113), (7137, 7154), (7168, 7181), 
                    (7192, 7240), (7389, 7981), (7563, 7809), (7576, 7793), (7577, 7792), (7578, 7791), 
                    (7579, 7790), (7580, 7789), (7581, 7788), (7584, 7786), (7585, 7785), (7587, 7783), 
                    (7590, 7780), (7591, 7779), (7597, 7778), (7598, 7777), (7599, 7776), (7602, 7775), 
                    (7603, 7774), (7604, 7773), (7605, 7772), (7606, 7771), (7607, 7770), (7608, 7769), 
                    (7609, 7768), (7611, 7767), (7612, 7766), (7614, 7765), (7615, 7764), (7616, 7763), 
                    (7617, 7762), (7618, 7761), (7619, 7760), (7620, 7759), (7621, 7758), (7629, 7750), 
                    (7630, 7749), (7632, 7747), (7633, 7746), (7641, 7740), (7644, 7737), (7645, 7736), 
                    (7646, 7735), (7648, 7733), (7649, 7732), (7650, 7731), (7651, 7730), (7653, 7729), 
                    (7654, 7728), (7655, 7727), (7656, 7726), (7657, 7725), (7658, 7724), (7665, 7720), 
                    (7666, 7719), (7667, 7718), (7668, 7717), (7672, 7713), (7673, 7712), (7677, 7707), 
                    (7679, 7705), (7681, 7701), (7682, 7700), (7917, 7929), (7918, 7928), (7935, 7950), 
                    (7936, 7949), (7937, 7948), (8019, 8025), (8147, 8580), (8161, 8188), (8162, 8187), 
                    (8206, 8218), (8208, 8215), (8226, 8513), (8232, 8508), (8241, 8499), (8259, 8478), 
                    (8281, 8454), (8284, 8451), (8305, 8412), (8320, 8400), (8763, 8998), (8765, 8996), 
                    (8766, 8995), (8767, 8994), (8772, 8985), (8774, 8983), (8775, 8982), (8778, 8979), 
                    (8779, 8978), (8787, 8974), (8788, 8973), (8789, 8972), (8791, 8970), (8795, 8966), 
                    (8797, 8961), (8799, 8959), (8800, 8958), (8802, 8956), (8806, 8951), (8808, 8949), 
                    (8809, 8948), (8810, 8947), (8811, 8946), (8816, 8941), (8819, 8935), (8823, 8932), 
                    (8824, 8931), (8828, 8927), (8829, 8926), (8830, 8925), (8831, 8924), (8832, 8923), 
                    (8834, 8919), (8835, 8918), (8837, 8916), (8853, 8901), (8854, 8900), (9487, 9548), 
                    (9576, 9699), (9579, 9696), (9581, 9694), (9587, 9602), (9607, 9689), (9873, 9999), 
                    (10204, 10313), (10206, 10216), (10223, 10276), (10240, 10258), (10381, 10481), 
                    (10385, 10477), (10395, 10467), (10396, 10466), (10405, 10421), (10406, 10420), 
                    (10440, 10451), (10441, 10450), (10442, 10449), (10485, 10503), (10486, 10502), 
                    (10487, 10501), (10488, 10500), (10489, 10499), (10515, 10680), (10516, 10679), 
                    (10604, 10673), (10605, 10672), (10606, 10671), (10607, 10670), (10608, 10669), 
                    (10612, 10644), (10616, 10640), (10617, 10639), (10618, 10638), (10619, 10637), 
                    (10620, 10636), (10621, 10635), (10622, 10634), (10699, 10710), (10700, 10709), 
                    (10712, 10795), (10713, 10794), (10714, 10793), (10717, 10790), (10718, 10789), 
                    (10719, 10788), (10720, 10787), (10721, 10786), (10723, 10784), (10724, 10783), 
                    (10725, 10782), (10728, 10779), (10730, 10778), (10731, 10777), (10732, 10776), 
                    (10733, 10775), (10734, 10774), (10735, 10773), (10736, 10772), (10737, 10771), 
                    (10738, 10770), (10739, 10769), (10740, 10768), (10741, 10767), (10742, 10766), 
                    (10745, 10757)]




