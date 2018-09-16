
import tools, structure, visual
reload(tools)


RED = "\033[31m"
GREEN = "\033[32m"
BLUE = "\033[34m"
DEF = "\033[39m"
GRAY = "\033[90m"
BLACK = "\033[30m"
BLOCK = "*"

dark_red = "\033[38;5;196m"
light_red = "\033[38;5;134m"
dark_blue = "\033[38;5;21m"
light_blue = "\033[38;5;33m"



raw_ids = """
KU740184|China|2016
KU820898|China|2016
KU501215.1|Puerto_Rico|2015
KU365778|Brazil|2015
KU312312|Suriname|2015
KU365777|Brazil|2015
KU365780|Brazil|2015
KU707826|Brazil|2015
KU365779|Brazil|2015
KU321639|Brazil|2015
KU509998|Haiti|2014
KU870645|USA|2016
KU501216|Guatemala|2015
KU501217|Guatemala|2015
KU497555|Brazil|2015
KU922960|Mexico|2016
KU922923|Mexico|2016
KU647676|Martinique|2015
KU820897|Colombia|2015
KU527068|Brazil|2015
KU866423|China|2016
KJ776791|French|2013
KX813683|Singapore|2016
KX827309|Singapore|2016
KU681081|Thailand|2014
KU955593|Cambodia|2010
KU681082|Philippines|2012
EU545988|Micronesia|2007
HQ234499|Malaysia|1966
KF383119|Senegal|2001
AY632535.2|Uganda|1947
KU955594|Uganda|1947
KU955591|Senegal|1984
KU955592|Senegal|1984
KU955595|Senegal|1984
"""
ids = raw_ids.strip().split('\n')
ids_african = ids[-6:]
ids_american = ids[:20]

"""
for trans_id in ids:
    if trans_id not in AlignCDS:
        print trans_id
KU527068|Brazil|2015
KJ776791|French|2013
HQ234499|Malaysia|1966
KF383119|Senegal|2001
"""

##############################
#   Some Important Functions
##############################


def clean_ID(inFastaList):
    cleanFastaList = {}
    for trans_id in inFastaList:
        cleanFastaList[trans_id.split('|')[0]] = inFastaList[trans_id]
    if len(inFastaList) != len(cleanFastaList):
        print "Warning: different length"
    return cleanFastaList

def get_sub_list(seq_dict, start, end):
    new_set = {}
    for trans_id in seq_dict:
        new_set[trans_id] = seq_dict[trans_id][start-1:end]
    return new_set

def filter_seq(seq_dict, std_strans_id="KU501215.1"):
    std_seq = seq_dict[std_strans_id]
    new_seq_dict = {std_strans_id: std_seq.replace('-', '')}
    for trans_id in seq_dict:
        if trans_id == std_strans_id:
            continue
        else:
            cur_seq = seq_dict[trans_id]
            assert(len(cur_seq) == len(std_seq))
            new_seq_dict[trans_id] = ""
            for idx,base in enumerate(list(std_seq)):
                if base != '-':
                    new_seq_dict[trans_id] += cur_seq[idx]
            assert( len(new_seq_dict[trans_id]) == len(new_seq_dict[std_strans_id]) )
    return new_seq_dict

def filter_shape(shape_dict, seq_dict, std_strans_id="KU501215.1"):
    assert(std_strans_id in shape_dict and std_strans_id in seq_dict)
    std_seq = seq_dict[std_strans_id]
    new_shape_dict = { std_strans_id: shape_dict[std_strans_id] }
    for trans_id in shape_dict:
        if trans_id == std_strans_id:
            continue
        else:
            cur_seq = seq_dict[trans_id]
            assert(len(cur_seq) == len(std_seq))
            new_shape_dict[trans_id] = []
            cur_index = 0
            for idx,base in enumerate(list(std_seq)):
                if base != '-':
                    if cur_seq[idx] == '-':
                        new_shape_dict[trans_id].append( "NULL" )
                    else:
                        new_shape_dict[trans_id].append( shape_dict[trans_id][cur_index] )
                        cur_index += 1
                else:
                    if cur_seq[idx] == '-':
                        pass
                    else:
                        cur_index += 1
    return new_shape_dict


def seq2codonList(raw_cds):
    assert(len(raw_cds) % 3 == 0)
    codon_list = []
    i = 0
    while i < len(raw_cds) - 1:
        codon_list.append(raw_cds[i:i+3])
        i += 3
    return codon_list

def generate_protein_mutation_profile(codonlist1, codonlist2, codonTable):
    assert(len(codonlist1) == len(codonlist2))
    profile = ""
    for c_1, c_2 in zip(codonlist1, codonlist2):
        if c_1 not in codonTable or c_2 not in codonTable:
            profile += "xxx"
        elif codonTable[c_1] == codonTable[c_2]:
            profile += "..."
        else:
            profile += "xxx"
    return profile

def generate_seq_mutation_profile(seq_1, seq_2):
    assert(len(seq_1) == len(seq_2))
    profile = ""
    for b_1, b_2 in zip(list(seq_1), list(seq_2)):
        if b_1 == b_2:
            profile += "."
        else:
            profile += "x"
    return profile

def generate_shape_profile(shape_list):
    shape_profile = ""
    for shape_value in shape_list:
        if shape_value == 'NULL':
            shape_profile += 'x'
        else:
            shape_value = float(shape_value)
            if shape_value < 0.3:
                shape_profile += '0'
            elif shape_value < 0.5:
                shape_profile += '1'
            elif shape_value < 0.7:
                shape_profile += '2'
            else:
                shape_profile += '3'
    return shape_profile

def generate_shape_profile_2(shape_list):
    shape_profile = ""
    for shape_value in shape_list:
        if shape_value == 'NULL':
            shape_profile += GRAY+BLOCK
        else:
            shape_value = float(shape_value)
            if shape_value < 0.3:
                shape_profile += BLACK+BLOCK
            elif shape_value < 0.5:
                shape_profile += BLUE+BLOCK
            elif shape_value < 0.7:
                shape_profile += GREEN+BLOCK
            else:
                shape_profile += RED+BLOCK
    shape_profile += DEF
    return shape_profile

def generate_shape_diff_profile(shape_list_1, shape_list_2):
    profile = ""
    for s_1, s_2 in zip(shape_list_1, shape_list_2):
        if 'NULL' not in (s_1, s_2):
            diff = float(s_1) - float(s_2)
            if diff < -0.7:
                profile += dark_blue+"*"
            elif diff < -0.3:
                profile += light_blue+"*"
            elif 0.3 < diff < 0.7:
                profile += light_red+"*"
            elif diff > 0.7:
                profile += dark_red+"*"
            else:
                profile += BLACK+"*"
    return profile


def generate_coordiation_profile(length, step=100, start=1):
    end = start + length
    profile = str(start)
    cur_coor = start
    blank_length = step-len(str(cur_coor))-1
    profile += " "*blank_length
    
    cur_coor=start+step-1
    while cur_coor < end+step:
        blank_length = step-len(str(cur_coor))
        profile += str(cur_coor) + " "*blank_length
        cur_coor += step
    return profile

def generate_mutliAlign_profile(Align, std_trans_id, start, end, transSet):
    profile = ""
    std_seq = Align[std_trans_id][start-1:end]
    max_len = len(std_trans_id)
    trans_ids = []
    for trans_id in transSet:
        if trans_id == std_trans_id:
            continue
        if len(trans_id) > max_len:
            max_len = len(trans_id)
        curSeq = Align[trans_id][start-1:end]
        assert( len(Align[trans_id]) == len(Align[std_trans_id]) )
        trans_ids.append(trans_id)
        profile += "%s"
        for idx,base in enumerate(list(curSeq)):
            if base == std_seq[idx]:
                profile += BLUE + base
            else:
                profile += RED + base
        profile += "\t"+trans_id+"\n"
    trans_ids.append(std_trans_id)
    profile += "%s"
    profile += GRAY + std_seq + "\t" + std_trans_id + "\n"
    profile = profile % tuple( item+" "*(max_len-len(item)+10) for item in trans_ids )
    return profile, max_len+10


def write_Align_Profile(Align, OUT, start, end, ids_american, ids_african):
    """
    Example: 
        write_Align_Profile(AlignCDS, outFile="/tmp/ss.txt", start=1050-107, end=1150-107, ids_american=ids_american, ids_african=ids_african)
    """
    ss_59 = structure.predictStructure(Align['KU501215.1'][start-1:end].replace('-','I'))
    ss_766 = structure.predictStructure(Align['AY632535.2'][start-1:end].replace('-','I'))
    ss_africa = structure.predictStructure(Align['KU955595'][start-1:end].replace('-','I'))
    
    #OUT = open(outFile, 'w')
    align_profile, max_len = generate_mutliAlign_profile(Align, 'KU501215.1', start, end, ids_american)
    print >>OUT, align_profile+" "*max_len+ss_59
    print >>OUT, " "*max_len+protein_mut_profile[start-1:end]
    print >>OUT, " "*max_len+generate_seq_mutation_profile(Align['KU501215.1'][start-1:end], Align['AY632535.2'][start-1:end])         #seq_mut_profile[start-1:end]
    print >>OUT, " "*max_len+generate_shape_diff_profile(Align['KU501215.1'][start-1:end], Align['AY632535.2'][start-1:end])
    print >>OUT, generate_mutliAlign_profile(Align, 'AY632535.2', start, end, ids_african)[0]+" "*max_len+ss_africa+"\n"+" "*max_len+ss_766
    print >>OUT, " "*max_len+generate_coordiation_profile(end-start+1, step=10)
    print >>OUT, " "*max_len+generate_coordiation_profile(end-start+1, step=10, start=start)
    OUT.close()


def get_nonsynonymouse_regions(protein_mut_profile, min_len=40, shrink=10):
    last_start = 0
    nonsynonymouse_regions = []
    for idx in range(len(protein_mut_profile)):
        if protein_mut_profile[idx] == 'x':
            if idx - last_start > min_len:
                if len(nonsynonymouse_regions) == 0:
                    if idx-shrink-last_start>min_len:
                        nonsynonymouse_regions.append((last_start+1, idx-shrink))
                else:
                    if (idx-shrink)-(last_start+shrink)>min_len:
                        nonsynonymouse_regions.append((last_start+shrink+1, idx-shrink))
            last_start = idx + 1
            continue
    if last_start != 0:
        if idx-(last_start+shrink)>min_len:
            nonsynonymouse_regions.append((last_start+shrink+1, idx))
    return nonsynonymouse_regions

def get_structure_change_dict(ShapeList_1, ShapeList_2, region_start, region_end, windowSize=5, step=2):
    shape_change_list = {}
    begin = region_start
    while begin < region_end:
        list_1 = ShapeList_1[begin-1:begin-1+windowSize]
        list_2 = ShapeList_2[begin-1:begin-1+windowSize]
        
        sum_diff = []
        for s_1, s_2 in zip(list_1, list_2):
            if s_1 != 'NULL' and s_2 != 'NULL':
                sum_diff.append( abs(float(s_1)-float(s_2)) )
        if len(sum_diff) > step/2:
            sum_diff = sum(sum_diff)/len(sum_diff)
            shape_change_list[begin] = round(sum_diff,3)
        
        begin += step
    return shape_change_list


def Structure_Change(nonsynonymouse_regions, AlignShapeCDS, windowSize=5, step=2):
    shape_change_dict = {}
    for region in nonsynonymouse_regions:
        myDict = get_structure_change_dict(AlignShapeCDS['KU501215.1'], AlignShapeCDS['AY632535.2'], region[0], region[1], windowSize=windowSize, step=step)
        for k in myDict:
            shape_change_dict[k] = myDict[k]
    shape_change_regions = sorted(shape_change_dict, key=lambda x: shape_change_dict[x], reverse=True)
    return shape_change_regions, shape_change_dict


def combine_sites(site_list, window_size=15):
    site_list = sorted(site_list)
    regions = []
    last_start = site_list[0]
    last_site = site_list[0]
    for idx in range(1, len(site_list)):
        if abs(site_list[idx]-last_site)<=window_size:
            last_site = site_list[idx]
        else:
            regions.append((last_start, last_site+15))
            last_start = last_site = site_list[idx]
    if last_site != site_list[-1]:
        regions.append((last_start, last_site+15))
    return regions

def write_mutation_matrix(Align, id_list, outFile):
    OUT = open(outFile, 'w')
    print >>OUT, " \t"+"\t".join(id_list)
    for id1 in id_list:
        print id1
        OUT.writelines(id1+"\t")
        for id2 in id_list:
            profile = generate_seq_mutation_profile(Align[id1], Align[id2])
            OUT.writelines("%s\t" % (profile.count('x'), ))
        OUT.writelines("\n")
    OUT.close()


def get_neighbor_strain(Align, std_trans_id, max_mutation):
    neighbor_strain_ids = []
    mutation_num_list = []
    for trans_id in Align:
        if trans_id == std_trans_id:
            continue
        assert( len(Align[std_trans_id]) == len(Align[trans_id]) )
        profile = generate_seq_mutation_profile(Align[std_trans_id], Align[trans_id])
        mutation_num = profile.count('x')
        mutation_num_list.append(mutation_num)
        if mutation_num <= max_mutation:
            neighbor_strain_ids.append(trans_id)
    return neighbor_strain_ids, sorted(mutation_num_list)



################################################################
#   尝试计数窗口中突变个数来判定变化区域
################################################################


######### MAC

key59, key766 = ['KU501215.1', 'AY632535.2']

Sequence = tools.readSeq("/Users/lee/Desktop/Projects/Virus/Virus_Genome/final_virus_genome.fa")
Shape = tools.loadicSHAPE("/Users/lee/Desktop/Projects/Virus/Virus_Genome/final_virus_icSHAPE.out")

del Shape[key59+'_vitro']; del Shape[key766+'_vitro'];
Align = tools.readSeq("/Users/lee/pCloud Drive/Virus_Structure/Covariation/Zika_358.afa")
#Align = tools.readSeq("/tmp/dengue.afa")
CodonTable = tools.loadCodon("/Users/lee/Desktop/Projects/Virus/Virus_Genome/Sequence/Codon_Table.txt")

######### Linux

key59, key766 = ['KU501215.1', 'AY632535.2']

Sequence = tools.readSeq("/Share/home/zhangqf8/lipan/virus/final_virus_genome.fa")
Shape = tools.loadicSHAPE("/Share/home/zhangqf8/lipan/virus/final_virus_icSHAPE.out")

del Shape[key59+'_vitro']; del Shape[key766+'_vitro'];
Align = tools.readSeq("/Share/home/zhangqf8/lipan/virus/select_virus_genome/Zika_358.afa")
CodonTable = tools.loadCodon("/Share/home/zhangqf8/lipan/virus/Codon_Table.txt")


########## left

align_genome_start = Align[key59].find( "GTTGTTGATCTGTGTG" ) + 1
align_genome_end = Align[key59].find("TGGGGAAATCCATGG" ) + 1 + len("TGGGGAAATCCATGG")

AlignPre = get_sub_list(Align, align_genome_start, align_genome_end)
AlignWhole = filter_seq(AlignPre, std_strans_id=key59)

seq_mut_profile = generate_seq_mutation_profile(AlignWhole[key59], AlignWhole[key766])

AlignShape = filter_shape(Shape, Align, std_strans_id=key59)


cds_start_59 = 107
cds_end_59 = 10375

align_cds_start = Align[key59].find( Sequence[key59][cds_start_59-1:cds_start_59+7] ) + 1
align_cds_end = Align[key59].find( Sequence[key59][cds_end_59-1:cds_end_59+15] ) + 1

AlignCDSPre = get_sub_list(Align, align_cds_start, align_cds_end)
AlignCDS = filter_seq(AlignCDSPre, std_strans_id=key59)

CodonList59 = seq2codonList(AlignCDS[key59])
CodonList766 = seq2codonList(AlignCDS[key59])

protein_mut_profile = generate_protein_mutation_profile(CodonList59, CodonList766, CodonTable)
full_protein_mut_profile = "."*(cds_start_59-1)+protein_mut_profile+"."*(10802-cds_end_59)


