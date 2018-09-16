
import tools, structure, visual


def readDG(inFile):
    dgList = []
    IN = open(inFile)
    line = IN.readline()
    while line:
        if line[0] == '>':
            data = line.strip().split()
            dgList.append((int(data[4]), int(data[5]), int(data[8]), int(data[9])))
        line = IN.readline()
    IN.close()
    return dgList

def uniq_list(inList):
    tmp_list = []
    for item in inList:
        if item not in tmp_list:
            tmp_list.append(item)
    return tmp_list

def short_long_dg(dgList):
    shortList = []
    longList = []
    for dg in dgList:
        if dg[2]-dg[1]>1000:
            longList.append(dg)
        else:
            shortList.append(dg)
    return shortList, longList

def classify_dg(dgList_59, dgList_766, tolerate=0):
    common_list = []
    for dg_59 in dgList_59:
        for dg_766 in dgList_766:
            s59_1, e59_1, s59_2, e59_2 = dg_59[:4]
            s766_1, e766_1, s766_2, e766_2 = dg_766[:4]
            if s766_1 > 1430:
                s766_1 += 12
            if e766_1 > 1430:
                e766_1 += 12
            if s766_2 > 1430:
                s766_2 += 12
            if e766_2 > 12:
                e766_2 += 12
            if s59_1<e766_1+tolerate and s766_1<e59_1+tolerate and s59_2<e766_2+tolerate and s766_2<e59_2+tolerate:
                common_list.append( (dg_59, dg_766) )
    uniq_59 = []
    uniq_766 = []
    common_59 = uniq_list([ item[0] for item in common_list ])
    common_766 = uniq_list([ item[1] for item in common_list ])
    for dg_59 in dgList_59:
        if dg_59 not in common_59:
            uniq_59.append(dg_59)
    for dg_766 in dgList_766:
        if dg_766 not in common_766:
            uniq_766.append(dg_766)
    return uniq_59, common_59, uniq_766, common_766


### Mac

DG_59 = readDG("/Users/lee/Desktop/Projects/Virus/Virus_Genome/PARIS/virus_paris/59/59.tab")
DG_766 = readDG("/Users/lee/Desktop/Projects/Virus/Virus_Genome/PARIS/virus_paris/766/766.tab")
Sequence = tools.readSeq("/Users/lee/Desktop/Projects/Virus/Virus_Genome/final_virus_genome.fa")


### Linux

DG_59 = readDG("/Share/home/zhangqf8/lipan/virus/PARIS/virus_paris/59/59.tab")
DG_766 = readDG("/Share/home/zhangqf8/lipan/virus/PARIS/virus_paris/766/766.tab")
Sequence = tools.readSeq("/Share/home/zhangqf8/lipan/virus/final_virus_genome.fa")

DG59_short, DG59_long = short_long_dg(DG_59)
DG766_short, DG766_long = short_long_dg(DG_766)


###########
# 统计59和766公共的和各自的dg数量
###########
key59, key766 = ['KU501215.1', 'AY632535.2']

uniq_59_short, common_59_short, uniq_766_short, common_766_short = classify_dg(DG59_short, DG766_short, tolerate=0)
print len(uniq_59_short), len(common_59_short), len(uniq_766_short), len(common_766_short)


uniq_59_long, common_59_long, uniq_766_long, common_766_long = classify_dg(DG59_long, DG766_long, tolerate=0)
print len(uniq_59_long), len(common_59_long), len(uniq_766_long), len(common_766_long)


###########
# 获得Lineage Specific Interaction

####################################################    Run alignment.py

# uniq_59_short, common_59_short, uniq_766_short, common_766_short = classify_dg(DG59_short, DG766_short, tolerate=0)
# uniq_59_long, common_59_long, uniq_766_long, common_766_long = classify_dg(DG59_long, DG766_long, tolerate=0)

###########


key59 = "KU501215.1"
key766 = "AY632535.2"
trans_list_59, muNum_59 = get_neighbor_strain(AlignCDS, key59, max_mutation=700)
trans_list_766, muNum_766 = get_neighbor_strain(AlignCDS, key766, max_mutation=700)


def count_ibp_num(ss, border):
    ctList = structure.dot2ct(ss)
    num = 0
    for bp in ctList:
        if bp[0] <= border < bp[1]:
            num += 1
    return num

def count_interaction_bp_num(seq_1, seq_2):
    seq_1 = seq_1.replace('-', 'I')
    seq_2 = seq_2.replace('-', 'I')
    ss = structure.bi_fold(seq_1, seq_2, local_pairing=True, mfe=True)
    left_border = len(seq_1)
    return count_ibp_num(ss, left_border)


def parse_dg(dg_list, Sequence, extend=20):
    import virus
    new_stem = []
    length = len(Sequence)
    for dg in dg_list:
        l_start, l_end, r_start, r_end = dg[0], dg[1], dg[2], dg[3]
        if not (l_start < l_end < r_start < r_end):
            print dg
            continue
        l_start = max(l_start-extend,1)
        r_end = min(r_end+extend,length)
        if r_start - l_end < extend*2:
            r_start = l_end = (l_end+r_start)/2
        else:
            l_end = min(l_end+extend,length)
            r_start = max(r_start-extend,1)
        
        seq_l = Sequence[l_start-1:l_end]
        seq_r = Sequence[r_start-1:r_end]
        full_seq = seq_l + "III" + seq_r
        
        try:
            ss = structure.bi_fold(seq_l, seq_r, local_pairing=True, mfe=True)
        except:
            continue
        
        stem = virus.find_stem(structure.dot2ct(ss), max_stem_gap=2, min_stem_len=8)
        break_point = len(seq_l) + 1
        dg_stem_list = []
        for idx in range(stem.size()):
            ls,le,rs,re = stem.at(idx).l_start, stem.at(idx).l_end, stem.at(idx).r_start, stem.at(idx).r_end
            if not ls<le<break_point<rs<re:
                continue
            abs_ls = ls + l_start - 1
            abs_le = le + l_start - 1
            abs_rs = rs - break_point - 3 + r_start + 1
            abs_re = re - break_point - 3 + r_start + 1
            
            sub_ss = ss[ls-1:le] + "..." + ss[rs-1:re]
            sub_seq = full_seq[ls-1:le] + "III" + full_seq[rs-1:re]
            border = le - ls + 1
            
            dg_stem_list.append((abs_ls, abs_le, abs_rs, abs_re, sub_seq, sub_ss, count_ibp_num(ss, border)))
        if dg_stem_list:
            new_stem.append( (dg_stem_list,dg) )
    return new_stem


stem_list_long_59 = parse_dg(uniq_59_long, Sequence[key59], extend=20)
stem_list_long_766 = parse_dg(uniq_766_long, Sequence[key766], extend=20)

stem_list_short_59 = parse_dg(uniq_59_short, Sequence[key59], extend=10)
stem_list_short_766 = parse_dg(uniq_766_short, Sequence[key766], extend=10)


def filter_strain_specific(dg_list_list, trans_list, key_name, AlignWhole, max_diff_num=10):
    def string_diff_num(str_1, str_2):
        assert(len(str_1) == len(str_2))
        num = 0
        for b1, b2 in zip(str_1, str_2):
            if b1 != b2 and '-' not in (b1, b2):
                num += 1
        return num
    newDG = []
    std_seq = AlignWhole[key_name]
    length = len(std_seq)
    
    for dg_list,raw_dg in dg_list_list:
        cur_new_dgList = []
        for dg in dg_list:
            diff_num = 0
            l_start, l_end, r_start, r_end = dg[0],dg[1],dg[2],dg[3]
            
            false_num = 0
            for trans_id in trans_list:
                str1 = AlignWhole[trans_id][l_start-1:l_end]
                str2 = std_seq[l_start-1:l_end]
                
                str3 = AlignWhole[trans_id][r_start-1:r_end]
                str4 = std_seq[r_start-1:r_end]
                
                diff_num += string_diff_num(str1, str2) + string_diff_num(str3, str4)
            
            if diff_num <= max_diff_num:
                cur_new_dgList.append( tuple(list(dg[:6])+[diff_num]) )
        if cur_new_dgList:
            #diff_num_list = [ it[-1] for it in cur_new_dgList ]
            newDG.append( [raw_dg, cur_new_dgList] )
    return newDG

key59 = "KU501215.1"
key766 = "AY632535.2"

DG_long_59 = filter_strain_specific(stem_list_long_59, trans_list_59, key59, AlignWhole, max_diff_num=20)
DG_long_766 = filter_strain_specific(stem_list_long_766, trans_list_766, key766, AlignWhole, max_diff_num=10)

DG_short_59 = filter_strain_specific(stem_list_short_59, trans_list_59, key59, AlignWhole, max_diff_num=20)
DG_short_766 = filter_strain_specific(stem_list_short_766, trans_list_766, key766, AlignWhole, max_diff_num=10)



def write_table(DG, outFile):
    OUT = open(outFile, 'w')
    print >>OUT, ""
    for raw_dg, stem_list in DG:
        stem_string = ""
        for stem in stem_list:
            stem_string += '%s-%s<==>%s-%s|%s|%s|%s\t' % (stem[0], stem[1], stem[2], stem[3], stem[4], stem[5], stem[6])
        stem_string = stem_string[:-1]
        print >>OUT, "%s\t%s\t%s\t%s\t%s" % ( raw_dg[0], raw_dg[1], raw_dg[2], raw_dg[3], stem_string )
    OUT.close()

write_table(DG_long_59, "/tmp/Zika/DG_long_59.txt")
write_table(DG_long_766, "/tmp/Zika/DG_long_766.txt")
write_table(DG_short_59, "/tmp/Zika/DG_short_59.txt")
write_table(DG_short_766, "/tmp/Zika/DG_short_766.txt")



