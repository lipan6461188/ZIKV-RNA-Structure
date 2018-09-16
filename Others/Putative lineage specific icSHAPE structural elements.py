

####################################################    Run alignment.py


def diff_shape(shape_list_1, shape_list_2):
    assert(len(shape_list_1)==len(shape_list_2))
    diff_shape = []
    for shape_1, shape_2 in zip(shape_list_1, shape_list_2):
        if shape_1 != 'NULL' and shape_2 != 'NULL':
            diff_shape.append( round(abs(float(shape_1)-float(shape_2)),3) )
        else:
            diff_shape.append( 'NULL' )
    return diff_shape

def combine_regions(raw_regions):
    new_regions = [ [it[0], it[1]] for it in raw_regions ]
    #print new_regions
    new_regions.sort(key=lambda x: x[0])
    idx = 0
    while idx < len(new_regions)-1:
        if new_regions[idx+1][0] <= new_regions[idx][1]:
            new_regions[idx][1] = new_regions[idx+1][1]
            del new_regions[idx+1]
        else:
            idx += 1
    return new_regions

def scan_diff_region(shape_diff_list, windowSize=5, min_diff_points=3, cutoff=0.6):
    def shape_num(shape_list):
        return len([ it for it in shape_list if it!='NULL' and float(it)>=cutoff ])
    regions = []
    idx = 0
    while idx+windowSize < len(shape_diff_list):
        cur_values = shape_diff_list[idx:idx+windowSize]
        if shape_num(cur_values) >= min_diff_points:
            regions.append([idx+1,idx+windowSize])
        idx += 1
    return combine_regions(regions)


def scan_diff_region_ave(shape_diff_list, windowSize=5, cutoff=0.5):
    def mean_ave(shape_list):
        valid_shape = [abs(it) for it in shape_list if it!='NULL']
        if len(valid_shape) > 3:
            return sum(valid_shape)/len(valid_shape)
        return 0.0
    regions = []
    idx = 0
    while idx+windowSize < len(shape_diff_list):
        cur_values = shape_diff_list[idx:idx+windowSize]
        if mean_ave(cur_values) >= cutoff:
            regions.append([idx+1,idx+windowSize])
        idx += 1
    return combine_regions(regions)

def annotate_diff_regions(shape_diff_list, regions):
    def gds(start, end):
        nums = sorted([ float(it) for it in shape_diff_list[start-1:end] if it != 'NULL' ])[-3:]
        return round(sum(nums)/3,3)
    return [ (it[0], it[1], it[1]-it[0]+1, gds(it[0], it[1])) for it in  regions]

def filter_strain_specific(regions, trans_list_59, trans_list_766, AlignWhole, max_diff_num=10, extend=5):
    def string_diff_num(str_1, str_2):
        assert(len(str_1) == len(str_2))
        num = 0
        for b1, b2 in zip(str_1, str_2):
            if b1 != b2 and '-' not in (b1, b2):
                num += 1
        return num
    newRegions = []
    key59 = "KU501215.1"
    key766 = "AY632535.2"
    seq59 = AlignWhole[key59]
    seq766 = AlignWhole[key766]
    length = len(seq766)
    
    for region in regions:
        diff_num_59 = 0
        diff_num_766 = 0
        start = max(region[0]-extend,0)
        end = min(region[1]+extend,length)
        
        for trans_id in trans_list_766:
            str1 = AlignWhole[trans_id][start:end]
            str2 = seq766[start:end]
            #print str1, str2, string_diff_num(str1, str2)
            diff_num_766 += string_diff_num(str1, str2)
        for trans_id in trans_list_59:
            str1 = AlignWhole[trans_id][start:end]
            str2 = seq59[start:end]
            #print str1, str2, string_diff_num(str1, str2)
            diff_num_59 += string_diff_num(str1, str2)
        #print region, diff_num
        if diff_num_766+diff_num_59 <= max_diff_num:
            newRegions.append(list(region)+[diff_num_59, diff_num_766])
    return newRegions

def save_table(regions, outFile, annoTable, Align, AlignShape, full_protein_mut_profile):
    OUT = open(outFile, 'w')
    key59 = "KU501215.1"
    key766 = "AY632535.2"
    print >>OUT, "start\tend\tseq59\tseq766\tshape59\tshape766\tshape_diff_score\tlineage_mut_num_59\tlineage_mut_num_766\tseq_mut_profile\tprotein_mut_profile\tpriotein_name"
    for region in regions:
        start, end = region[0], region[1]
        region_name = ""
        for pName in annoTable:
            if annoTable[pName][0] < end and start < annoTable[pName][1]:
                region_name = pName
        seq59 = Align[key59][start-1:end]
        seq766 = Align[key766][start-1:end]
        print >>OUT, "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}".format(start, end,   
                                            seq59,
                                            seq766,
                                            ",".join([str(sh) for sh in AlignShape[key59][start-1:end]]),
                                            ",".join([str(sh) for sh in AlignShape[key766][start-1:end]]),
                                            region[3],
                                            region[4],
                                            region[5],
                                            generate_seq_mutation_profile(seq59, seq766),
                                            full_protein_mut_profile[start-1:end],
                                            region_name)
    OUT.close()

def plot_coor(region_clean, outFile):
    symbol = [0]*10802
    for region in sorted(region_clean,key=lambda x: x[3], reverse=True)[:10]+region_clean[-1:]:
        for idx in range(region[0], region[1]):
            symbol[idx] = 1
    symbol[0] = symbol[-1] = 1
    plt.bar(range(len(symbol)), symbol)
    plt.savefig(outFile)
    plt.close()


def read_annoTable(inFile):
    annoTable = {}
    for line in open(inFile).readlines():
        if line[0] != '#':
            data = line.strip().split()
            annoTable[data[3]] = (int(data[1]), int(data[2]))
    return annoTable



shape_diff = diff_shape(AlignShape['KU501215.1'], AlignShape['AY632535.2'])
regions = scan_diff_region(shape_diff, windowSize=5, min_diff_points=3, cutoff=0.6)

#regions = scan_diff_region_ave(shape_diff, windowSize=5, cutoff=0.5)

regions = annotate_diff_regions(shape_diff, regions)

trans_list_59, mutation_num_list = get_neighbor_strain(AlignCDS, 'KU501215.1', 700)
trans_list_766, mutation_num_list = get_neighbor_strain(AlignCDS, 'AY632535.2', 700)

region_clean = filter_strain_specific(regions, trans_list_59, trans_list_766, AlignWhole, max_diff_num=25, extend=5)
len(region_clean)

## Plot Top 10
xx = sorted(region_clean, key=lambda x: x[3], reverse=True)[:10]
plot_coor(xx, "/tmp/ss.pdf")


annoTable = read_annoTable("/Users/lee/Desktop/Projects/Virus/Virus_Genome/Figures/Figure_1/59_annotation.txt")
save_table(region_clean, "/tmp/Zika/structure_change.txt", annoTable, AlignWhole, AlignShape, full_protein_mut_profile)



