
"""

Predict structure of each domain and compare with PARIS duplex group

"""

import structure, icSHAPE, virus

def prepare_shape_file(shape_list, file_name):
    OUT = open(file_name, "w")
    for idx in range(len(shape_list)):
        if shape_list[idx] != "NULL":
            print >>OUT, "%d\t%s" % (idx+1, shape_list[idx])
    OUT.close()

def prepare_fasta_file(raw_seq, file_name):
    OUT = open(file_name, "w")
    print >>OUT, ">seq\n%s" % (raw_seq)
    OUT.close()

def dp_file_to_rainbow(input_dp_file, output_rb_file):
    IN = open(input_dp_file, "r")
    OUT = open(output_rb_file, "w")
    line = IN.readline()
    length = line.strip()
    print >>OUT, "#len:"+length
    line_idx = 0
    while line:
        if line_idx <= 1:
            pass
        else:
            left, right, score = line.strip().split()
            prob = 10**(-float(score))
            color = ""
            if prob >= 0.8:
                print >>OUT, "%s\t%s\t%s\t1\t1" % ( left, right, "#59CE83" )
            elif prob >= 0.3:
                print >>OUT, "%s\t%s\t%s\t1\t1" % ( left, right, "#67B2DD" )
            elif prob >= 0.1:
                print >>OUT, "%s\t%s\t%s\t1\t1" % ( left, right, "#FFDC6A" )
            elif prob >= 0.03:
                print >>OUT, "%s\t%s\t%s\t1\t1" % ( left, right, "#E0E0E0" )
            else:
                pass
        line_idx += 1
        line = IN.readline()
    IN.close()
    OUT.close()

def load_dp_file(input_dp_file):
    paring_list = []
    IN = open(input_dp_file)
    line = IN.readline()
    line_idx = 0
    while line:
        if line_idx <= 1:
            pass
        else:
            left, right, score = line.strip().split()
            prob = 10**(-float(score))
            if prob > 0.03:
                paring_list.append( (int(left), int(right), prob) )
        line_idx += 1
        line = IN.readline()
    IN.close()
    return paring_list


def seperate_fold_whole_virus(virus_domain, with_shape=True):
    seq_file = "/tmp/tmp_seq.fa"
    shape_file = "/tmp/tmp_shape.txt"
    ct_file = "/tmp/tmp_ct.ct"
    
    virus_len = virus_domain.length()
    paring_list = []
    for domain_idx in range(virus_domain.size()):
        print "==========> domain: %s <==========" % (domain_idx, )
        cur_seq = virus_domain.sequence(domain_idx)
        cur_shape = virus_domain.shape(domain_idx)
        cur_base_coor = virus_domain.at(domain_idx)[0] - 1
        
        if with_shape:
            dot_structure = structure.predictStructure(cur_seq, Shape=cur_shape, bp_constraint=[], mfe=True)
        else:
            dot_structure = structure.predictStructure(cur_seq, Shape=[], bp_constraint=[], mfe=True)
        bp_list = structure.dot2ct(dot_structure)
        for item in bp_list:
            paring_list.append( (item[0]+cur_base_coor, item[1]+cur_base_coor) )
    return paring_list

def rainbow_overlap_paris_stems(dg_list, stem_list, outFile_stem, outFile_paris, tolerant=20):
    import copy
    stem_list = copy.deepcopy(stem_list)
    PAIR = open(outFile_stem, 'w')
    PARIS = open(outFile_paris, 'w')
    print >>PAIR, "#len:10802"
    print >>PARIS, "#len:10802"
    delete_list = []
    for domain_idx in range(len(dg_list)):
        for dg in dg_list[domain_idx]:
            curColor='black'; curProb=0.0
            cur_delete_list = []
            for stem in stem_list:
                if dg[0]<stem[1]+tolerant and stem[0]<dg[1]+tolerant and dg[2]<stem[3]+tolerant and stem[2]<dg[3]+tolerant:
                    curColor = "red"
                    if stem in delete_list:
                        continue
                    delete_list.append( stem )
                    cur_delete_list.append( ( (stem[0]+stem[1])/2, (stem[2]+stem[3])/2, stem ) )
            print >>PARIS, "%s\t%s\t%s\t3\t1" % ((dg[0]+dg[1])/2, (dg[2]+dg[3])/2, curColor)
            for stem_pair in cur_delete_list:
                print >>PAIR, "%s\t%s\t%s\t3\t1" % (stem_pair[0], stem_pair[1], "red")
                #stem_list.remove(stem_pair[2])
    for stem in stem_list:
        if stem not in delete_list:
            print >>PAIR, "%s\t%s\t%s\t3\t1" % ((stem[0]+stem[1])/2, (stem[2]+stem[3])/2, "black")
    PAIR.close()
    PARIS.close()

def filter_stem_in_domain(domain_start, domain_end, in_stem_list):
    out_stem_list = []
    for stem in in_stem_list:
        assert( stem[0] < stem[1] < stem[2] < stem[3] )
        if domain_start <= stem[0] <= stem[3] <= domain_end:
            out_stem_list.append(stem)
    return out_stem_list

def shuffle_stem(domain_start, domain_end, stem_list):
    import random
    random_stem = []
    domain_wid = domain_end - domain_start
    for stem in stem_list:
        assert( domain_start <= stem[0] < stem[1] < stem[2] < stem[3] <= domain_end  )
        stem_len = stem[3] - stem[0] + 1
        max_dist = domain_wid - stem_len
        start_dist = stem[0] - domain_start
        rand_bias = random.randint(0,max_dist)
        random_stem.append( (stem[0]-start_dist+rand_bias, stem[1]-start_dist+rand_bias, stem[2]-start_dist+rand_bias, stem[3]-start_dist+rand_bias) )
        #print stem, domain_wid, stem_len, max_dist, rand_bias, random_stem[-1]
    return random_stem

def paris_stem_overlap(dg_list, stem_list, tolerant=0):
    common_dg_number = 0
    for dg in dg_list:
        for stem in stem_list:
            if dg[0]<stem[1]+tolerant and stem[0]<dg[1]+tolerant and dg[2]<stem[3]+tolerant and stem[2]<dg[3]+tolerant:
                common_dg_number += 1
                break
    return common_dg_number

def shuffle_domain_stems(in_stem_list, dg_list, domain_list, tolerant=0):
    random_shuffle_list = []
    index = 0
    for start, end in domain_list:
        print start, end
        random_shuffle_list.append([])
        stem_in_domain = filter_stem_in_domain(start, end, in_stem_list)
        for i in range(1000):
            random_stem = shuffle_stem(start, end, stem_in_domain)
            common_dg_nums = paris_stem_overlap(dg_list[index], random_stem, tolerant=tolerant)
            random_shuffle_list[-1].append(common_dg_nums)
        index += 1
    random_sum = [0]*1000
    for idx in range(1000):
        random_sum[idx] += sum( [ cur_list[idx] for cur_list in random_shuffle_list ] )
    return sorted(random_sum)

def get_common_dg_num(in_stem_list, dg_list, tolerant=0):
    num = 0
    for idx in range(len(dg_list)):
        num += paris_stem_overlap(dg_list[idx], in_stem_list, tolerant=tolerant)
    return num

def get_a_shuffle_stem_list(in_stem_list, domain_list):
    random_stem_list = []
    for start, end in domain_list:
        stem_in_domain = filter_stem_in_domain(start, end, in_stem_list)
        random_stem = shuffle_stem(start, end, stem_in_domain)
        random_stem_list += random_stem
    return random_stem_list

def read_paris_dg(in_tab_file, domain_list):
    dg_list = [ [] for i in range(len(domain_list)) ]
    IN = open(in_tab_file)
    line = IN.readline()
    while line:
        if line[0] == '>':
            data = line.strip().split()
            start_1, end_1, start_2, end_2 = int(data[4]), int(data[5]), int(data[8]), int(data[9])
            for idx in range(len(domain_list)):
                if domain_list[idx][0] <= start_1 < end_2 <= domain_list[idx][1]:
                    dg_list[idx].append((start_1, end_1, start_2, end_2))
        line = IN.readline()
    IN.close()
    return dg_list

def get_domain_list(domain_handle):
    return [ domain_handle.at(i) for i in range(domain_handle.size()) ]

def get_stem_list(stems_handle):
    return [ (stems_handle.at(idx).l_start, stems_handle.at(idx).l_end, stems_handle.at(idx).r_start, stems_handle.at(idx).r_end) for idx in range(stems_handle.size()) ]

def get_p_value(shuffle_value_list, true_value):
    shuffle_value_list.sort()
    tn = len(shuffle_value_list)
    for idx, cur_value in enumerate(shuffle_value_list):
        if cur_value >= true_value:
            return 1.0*(tn-idx-1)/tn
    return 0.0





##### 1. Load sequence and icSHAPE data

seq = icSHAPE.readSeq("/Users/lee/Desktop/Projects/Virus/Virus_Genome/final_virus_genome.fa")
shape = icSHAPE.loadicSHAPE("/Users/lee/Desktop/Projects/Virus/Virus_Genome/final_virus_icSHAPE.out")

##### 2. Read domain boundaries

domain_766 = virus.read_domain("/Users/lee/pCloud Drive/Projects/Virus/re-domain/766_domain.txt")
domain_766.load_seq(seq['AY632535'])
domain_766.load_shape(shape['AY632535'])

domain_59 = virus.read_domain("/Users/lee/pCloud Drive/Projects/Virus/re-domain/59_domain.txt")
domain_59.load_seq(seq['KU501215'])
domain_59.load_shape(shape['KU501215'])

##### 3. Predict seondary structure of each domain seperately

structure_59 = seperate_fold_whole_virus(domain_59, with_shape=True)
structure_766 = seperate_fold_whole_virus(domain_766, with_shape=True)

##### 4. base pairs clustering

structure_59.sort(key=lambda x: x[0])
stems_59 = virus.find_stem(structure_59, max_stem_gap=1, min_stem_len=5)

structure_766.sort(key=lambda x: x[0])
stems_766 = virus.find_stem(structure_766, max_stem_gap=1, min_stem_len=5)

##### 5. Write common or specific rainbow

# 59

domain_list = get_domain_list(domain_59)
dg_list = read_paris_dg("/tmp/59.tab", domain_list)
stem_list = get_stem_list(stems_59)

rainbow_overlap_paris_stems(dg_list, stem_list, "/tmp/virus_paris/59_ss.rainbow", "/tmp/virus_paris/59_paris.rainbow", tolerant=tolerant)
get_common_dg_num(stem_list, dg_list, tolerant=5)


# 766

domain_list = get_domain_list(domain_766)
dg_list = read_paris_dg("/tmp/virus_paris/766/766.tab", domain_list)
stem_list = get_stem_list(stems_766)

rainbow_overlap_paris_stems(dg_list, stem_list, "/tmp/virus_paris/766_ss.rainbow", "/tmp/virus_paris/766_paris.rainbow", tolerant=tolerant)
get_common_dg_num(stem_list, dg_list, tolerant=5)


##### 6. Visualization

"""

rainbow_plot 59_paris.rainbow 59_paris.pdf
rainbow_plot 59_ss.rainbow 59_ss.pdf 

rainbow_plot 766_paris.rainbow 766_paris.pdf
rainbow_plot 766_ss.rainbow 766_ss.pdf 

"""

##### 7. Shuffle and P-value

random_shuffle_list = shuffle_domain_stems(stem_list, dg_list, domain_list, tolerant=5)
print get_p_value(random_shuffle_list, get_common_dg_num(stem_list, dg_list, tolerant=5))



