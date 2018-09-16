import os, sys, commands, re, numpy
import pandas as pd
#import matplotlib.pyplot as plt
#import seaborn as sns


################################
######## 加载icSHAPE和Sequence
################################

def readSeq(seq_name):
    # + 读取转录组序列
    transSeq = {}
    IN = open(seq_name)
    line = IN.readline()
    cur_trans = ''
    while line:
        if line[0] == '>':
            cur_trans = line[1:].split()[0]
            cur_trans = cur_trans.split(".")[0]
            transSeq[ cur_trans ] = ''
        else:
            transSeq[ cur_trans ] += line.strip()
        line = IN.readline()
    return transSeq

def loadicSHAPE(file_name):
    # + 读取icSHAPE数据
    # + 并用长度信息校正
    icSHAPE = {}
    IN = open(file_name)
    line = IN.readline()
    while line:
        arr = line.strip().split()
        trans_id = arr[0].split(".")[0]
        shape = arr[3:]
        for idx, shape_value in enumerate(shape):
            if shape_value != "NULL":
                shape[idx] = float(shape_value)
        icSHAPE[ trans_id ] = shape
        line = IN.readline()
    return icSHAPE

################################
######## 加载PARIS数据
################################

class ParisRead:
    def __init__(self, read_id, cigar, left_start, left_end, right_start, right_end):
        self.read_id = read_id
        self.cigar = cigar
        self.left_start = int(left_start)
        self.left_end = int(left_end)
        self.right_start = int(right_start)
        self.right_end = int(right_end)
        assert(self.left_start <= self.left_end)
        assert(self.right_start <= self.right_end)
        assert(self.left_end < self.right_start)
    def left_arm(self):
        return (self.left_start, self.left_end)
    def right_arm(self):
        return (self.right_start, self.right_end)
    def ovrehang(self):
        return self.right_start - self.left_end
    def seq(self, raw_seq):
        return (raw_seq[self.left_start-1:self.left_end], raw_seq[self.right_start-1:self.right_end])
    def left_len(self):
        return self.left_end - self.left_start + 1
    def right_len(self):
        return self.right_end - self.right_start + 1

class VirusParis:
    def __init__(self, virus_len):
        self.read_list = []
        self.virus_len = virus_len
        self.sorted = False
    def sort(self):
        self.read_list.sort(key=lambda x: [x.left_start, x.right_start])
        self.sorted = True
    def append(self, other):
        assert(isinstance(other, ParisRead))
        assert(other.right_end <= self.virus_len)
        self.read_list.append(other)
        self.sorted = False
    def build_interaction_block(self, start, end):
        assert(start>=1 and end <= self.virus_len)
        if not self.sorted:
            self.sort()
        block = numpy.array(0, dtype=float)
        block.resize((end-start+1, end-start+1))
        for read in self.read_list:
            if start <= read.left_end < read.right_start <= end:
                for idx in range(read.left_start, read.left_end+1):
                    if idx < start: continue
                    for idy in range(read.right_start, read.right_end+1):
                        if idy > end: continue
                        block[idx-start][idy-start] += 1
                        block[idy-start][idx-start] += 1
        return block
    def show(self, n=50):
        if not self.sorted:
            self.sort()
        for read_id in range(min(n, len(self.read_list))):
            read = self.read_list[read_id]
            print "%s\t%s\t%s\t%s\t%s\t%s" % ( read.read_id, read.cigar, read.left_start, read.left_end, read.right_start, read.right_end)
    def length(self):
        return self.virus_len
            
def read_PARIS(dg_file_name, virus_len):
    vp = VirusParis(virus_len)
    IN = open(dg_file_name)
    line = IN.readline()
    while line:
        data = line.strip().split()
        read = ParisRead(data[0], data[4], data[5], data[6], data[7], data[8])
        vp.append(read)
        line = IN.readline()
    return vp

# remove unpaired paris data from raw block
def verify_interaction_block_seq(block, sequence):
    assert(len(block) == len(sequence))
    for idx in range(len(block)):
        for idy in range(idx+1, len(block)):
            if sequence[idx]+sequence[idy] not in ('AT', 'TA', 'CG', 'GC', 'GT', 'TG', 'AU', 'UA', 'GU', 'UG'):
                block[idx][idy] = 0
                block[idy][idx] = 0

# remove hign icSHAPE score paris data from raw block
def verify_interaction_block_shape(block, shape, shape_cutoff=0.7):
    assert(len(block) == len(shape))
    for idx in range(len(block)):
        for idy in range(len(block)):
            if shape[idx] >= shape_cutoff or shape[idy] >= shape_cutoff:
                block[idx][idy] = block[idy][idx] = 0

def normalize_block(raw_block):
    print "Normalize Block Information..."
    positive_value = []
    for idx in range(len(raw_block)):
        for idy in range(idx+1, len(raw_block)):
            if raw_block[idx][idy] > 0:
                positive_value.append(raw_block[idx][idy])
    raw_positive_value = numpy.array(positive_value, dtype=float)
    raw_positive_value.sort()
    positive_value = raw_positive_value.copy()
    nums = len(positive_value)
    #print "Values: ", positive_value
    print "\tPositive Values Numbers: ", nums
    norm_start_pos = int(nums * 16.0 / 20)
    norm_end_pos = int(nums * 18.0 / 20)
    #print "start: %d; end: %d" % (norm_start_pos, norm_end_pos)
    normalized_factor = numpy.mean(positive_value[norm_start_pos:norm_end_pos])
    #print "Normalized Factor: %s" % (normalized_factor, )
    for idx in range(len(positive_value)):
        #print "1.0*%s/%s" % ( positive_value[idx], normalized_factor)
        positive_value[idx] = float(1.0*positive_value[idx] / normalized_factor)
    #print "positive_value After normalized: ", positive_value
    upper_bound = positive_value[ int( nums * 0.9 ) ]
    lower_bound = positive_value[ int( nums * 0.1 ) ]
    raw_upper_bound = raw_positive_value[ int( nums * 0.9 ) ]
    raw_lower_bound = raw_positive_value[ int( nums * 0.1 ) ]
    print "\tRaw --- Upper: %f; Lower: %f" % (raw_upper_bound, raw_lower_bound)
    print "\tNorm --- Upper: %f; Lower: %f" % (upper_bound, lower_bound)
    for idx in range(len(raw_block)):
        for idy in range(idx+1, len(raw_block)):
            if raw_block[idx][idy] > 0:
                raw_block[idx][idy] /= normalized_factor
                if raw_block[idx][idy] >= upper_bound:
                    raw_block[idx][idy] = 1
                elif raw_block[idx][idy] <= lower_bound:
                    raw_block[idx][idy] = 0
                else:
                    raw_block[idx][idy] = (raw_block[idx][idy]-lower_bound)/(upper_bound-lower_bound)
                raw_block[idx][idy] = raw_block[idx][idy]
                raw_block[idy][idx] = raw_block[idx][idy]

# save block data to file to readable by human
def save_block(block, sequence, file_name):
    assert(len(block) == len(sequence))
    OUT = open(file_name, "w")
    OUT.writelines(" \t")
    for idx in range(len(sequence)):
        OUT.writelines(sequence[idx]+"\t")
    OUT.writelines("\n")
    for idx in range(len(sequence)):
        OUT.writelines(sequence[idx]+"\t")
        for idy in range(len(sequence)):
            OUT.writelines(`round(block[idx][idy],3)`+"\t")
        OUT.writelines("\n")
    OUT.close()


################################
######## 加载Domain相关的函数和类
################################

class Domain:
    def __init__(self, start, end, is_linker=False):
        self.start = int(start)
        self.end = int(end)
        self.linker = is_linker
        assert(self.start <= self.end)
    def length(self):
        return self.end - self.start + 1
    def seq(self, raw_seq):
        return raw_seq[self.start-1:self.end]

class VirusDomain:
    def __init__(self):
        self.domain_list = []
        self.sorted = True
        self.seq = ""
        self.sh = []
        self.PARIS = []
        self.domain_without_linker = []
    def append(self, domain):
        assert(isinstance(domain, Domain))
        self.domain_list.append(domain)
        self.seq = ""
        self.sh = []
        self.PARIS = []
        self.sorted = False
    def valid(self):
        if not self.sorted:
            self.sort()
        for idx in range(len(self.domain_list)-1):
            if idx == 0:
                if self.domain_list[idx].start != 1:
                    return False
            if self.domain_list[idx].end != self.domain_list[idx+1].start - 1:
                return False
        return True
    def sort(self):
        self.domain_list.sort(key=lambda x: x.start)
        self.sorted = True
        self.domain_without_linker = [d for d in self.domain_list if d.linker==False]
    def size(self):
        if not self.sorted:
            self.sort()
        return len(self.domain_without_linker)
    def at(self, idx):
        assert(0 <= idx < self.size())
        return ( self.domain_without_linker[idx].start, self.domain_without_linker[idx].end )
    def show(self, OUT):
        for d in self.domain_list:
            if d.linker:
                print >>OUT,"%d => %d  linker" % (d.start, d.end)
            else:
                print >>OUT,"%d => %d  domain" % (d.start, d.end)
    def length(self):
        if not self.sorted:
            self.sort()
        return self.domain_list[-1].end
    def load_seq(self, raw_seq):
        assert(len(raw_seq) == self.length())
        self.seq = raw_seq
    def sequence(self, idx=None):
        assert(len(self.seq) == self.length())
        if idx is None:
            return self.seq
        else:
            assert(0 <= idx < self.size())
            d = self.domain_without_linker[idx]
            return self.seq[d.start-1:d.end]
    def load_shape(self, raw_shape):
        assert(len(raw_shape) == self.length())
        for item in raw_shape:
            if item != 'NULL':
                self.sh.append(float(item))
            else:
                self.sh.append('NULL')
    def shape(self, idx):
        assert(0 <= idx < self.size())
        assert(len(self.sh) == self.length())
        d = self.domain_without_linker[idx]
        return self.sh[d.start-1:d.end]
    def load_paris(self, raw_paris):
        assert( isinstance(raw_paris, VirusParis) )
        assert( raw_paris.length() == self.length() )
        self.PARIS = raw_paris
    def paris(self, idx):
        assert(0 <= idx < self.size())
        assert(len(self.sh) == self.length())
        start, end = self.at(idx)
        block = self.PARIS.build_interaction_block(start, end)
        verify_interaction_block_seq(block, self.sequence(idx))
        verify_interaction_block_shape(block, self.shape(idx), shape_cutoff=0.6)
        normalize_block(block)
        return block
    def global_coor_2_domain_coor(self, global_coor):
        assert( global_coor <= len(self.seq) )
        if not self.sorted:
            self.sort()
        for domain_idx, cur_domain in enumerate(self.domain_list):
            if cur_domain.start <= global_coor <=  cur_domain.end:
                return [domain_idx, global_coor-cur_domain.start+1]
        return None

def read_domain(file_name):
    vd = VirusDomain();
    IN = open(file_name)
    line = IN.readline()
    while line:
        start, end = line.strip().split()
        d = Domain(start, end)
        vd.append(d)
        line = IN.readline()
    if vd.valid():
        return vd
    else:
        vd.show(sys.stdout)
        raise Exception("Bad Domain Format")

################################
######## 寻找二级结构上的Stem
################################

def ct2dot(ct):
    assert(isinstance(ct, str))
    ctList = []
    Stack = []
    for idx, symbol in enumerate(list(ct)):
        if symbol in ('(', '[', '<', '{'):
            Stack.append(idx+1)
        elif symbol in (')', ']', '>', '}'):
            ctList.append( (Stack[-1], idx+1) )
            del Stack[-1]
        else:
            pass
    ctList.sort(key=lambda x:x[0])
    return ctList

class Stem:
    def __init__(self, l_start, l_end, r_start, r_end):
        self.l_start = int(l_start)
        self.l_end = int(l_end)
        self.r_start = int(r_start)
        self.r_end = int(r_end)

class VirusStem:
    def __init__(self):
        self.stem_list = []
        self.sorted = False
    def append(self, stem):
        assert( isinstance(stem, Stem) )
        self.stem_list.append(stem)
        self.sorted = False
    def sort(self):
        self.stem_list.sort(key=lambda x: [x.l_start, x.r_start])
        self.sorted = True
    def at(self, idx):
        assert( 0 <= idx < len(self.stem_list) )
        if not self.sorted:
            self.sort()
        return self.stem_list[idx]
    def valid(self):
        if not self.sorted:
            self.sort()
        return True
    def size(self):
        return len(self.stem_list)
    def length(self):
        if not self.sorted:
            self.sort()
        return self.stem_list[-1].r_end
    def show(self):
        if not self.sorted:
            self.sort()
        for stem in self.stem_list:
            print "%d-%d <==> %d-%d" % (stem.l_start, stem.l_end, stem.r_start, stem.r_end)
    def tag_read(self, read, base_coor=0):
        assert(isinstance(read, ParisRead))
        max_overlap = 0; max_overlap_id = 0
        read_left_start = read.left_start - base_coor; read_left_end = read.left_end - base_coor;
        read_right_start = read.right_start - base_coor; read_right_end = read.right_end - base_coor
        #print read_left_end, read_right_start, self.length
        if read_left_end <= 0 or read_right_start > self.length():
            return -1, 0
        for stem_id in range(len(self.stem_list)):
            cur_stem = self.stem_list[stem_id]
            if cur_stem.l_start <= read_left_end and read_left_start <= cur_stem.l_end and cur_stem.r_start <= read_right_end and read_right_start <= cur_stem.r_end:
                overlap = min(read_left_end, cur_stem.l_end) - max(read_left_start, cur_stem.l_start) + 1
                overlap += min(read_right_end, cur_stem.r_end) - max(read_right_start, cur_stem.r_start) + 1
                if overlap > max_overlap:
                    max_overlap = overlap
                    max_overlap_id = stem_id
        if max_overlap == 0:
            if read_right_start < self.length():
                return 0, 0
        return max_overlap_id + 1, max_overlap

def find_stem(ctList, max_stem_gap=3, min_stem_len=5):
    vs = VirusStem()
    l_start = 0; l_end = 0
    r_start = 0; r_end = 0
    for pair in ctList:
        if l_start == 0:
            l_start = l_end = pair[0]
            r_start = r_end = pair[1]
        elif abs(pair[0] - l_end) > max_stem_gap or abs(r_start - pair[1]) > max_stem_gap:
            #print l_start, l_end, r_start, r_end
            if l_end - l_start + 1 >= min_stem_len and r_end - r_start + 1 >= min_stem_len:
                stem = Stem(l_start, l_end, r_start, r_end)
                vs.append(stem)
            l_start = l_end = pair[0]
            r_start = r_end = pair[1]
        else:
            l_end = pair[0]
            r_start = pair[1]
    if l_end - l_start + 1 >= min_stem_len and r_end - r_start + 1 >= min_stem_len:
        stem = Stem(l_start, l_end, r_start, r_end)
        vs.append(stem)
    return vs


################################
######## 分类Sam中的Reads
################################


def get_match_region(globalBasePos, Cigar):
    """Parse Cigar Symbol
    Test:
        "3S3M2D3M2I3M3S":
        - - - 1 1 1 0 0 1 1 1 - - 1 1 1 - - - ref
        0 0 0 1 1 1 - - 1 1 1 0 0 1 1 1 0 0 0 que
        parseCigar(1, '3S3M2D3M2I3M3S') ([[1, 3], [6, 8], [9, 11]], [[4, 6], [7, 9], [12, 14]], 17)
        parseCigar(1, '3S3M2D3M2I3M3=3S')
        parseCigar(1, '3S3M2D3M1000N2I3M3=3S')
    """
    def remove(List, elem):
        return filter(lambda a: a != elem, List)
    numVec = [ int(a) for a in remove(re.split('[MIDNSHP=X]', Cigar), '') ]
    modVec = remove(re.split('\d', Cigar), '')
    End = globalBasePos - 1
    matchGlobalSeg = []
    matchLocalSeg = []
    SeqLen = 0
    for Len, Mode in zip(numVec, modVec):
        if Mode in ('M', '=', 'X'):
            matchGlobalSeg.append( [End+1, End + Len] )
            End += Len
            matchLocalSeg.append( [SeqLen+1, SeqLen+Len] )
        elif Mode in ('D', 'N'):
            End += Len
        if Mode in ('M', 'I', 'S', 'H', '=', 'X'):
            SeqLen += Len
    return matchGlobalSeg, matchLocalSeg, SeqLen

def tag_sam_reads(in_sam_file_name, virus_stem, out_sam_file_name, base_coor):
    #color_database = ["#4c72b0", "#55a868", "#c44e52", "#8172b2", "#ccb974", "#64b5cd", "#e294ea", "#11efaf" ,"#1368e5", "#eaa094", "#ea0c96"]
    # record reads numbers for each stem
    read_record = [0] * (virus_stem.size() + 1)
    IN = open(in_sam_file_name)
    OUT = open(out_sam_file_name, "w")
    line = IN.readline()
    while line:
        if line.startswith('@'):
            print >>OUT, line.strip()
            line = IN.readline()
            continue
        data = line.strip().split()
        matchGlobalSeg, matchLocalSeg, SeqLen = get_match_region( int(data[3]), data[5] )
        if len(matchGlobalSeg) == 2:
            #def __init__(self, read_id, cigar, left_start, left_end, right_start, right_end):
            paris_read = ParisRead(data[0], data[5], matchGlobalSeg[0][0], matchGlobalSeg[0][1], matchGlobalSeg[1][0], matchGlobalSeg[1][1])
            tag_id, max_overlap = virus_stem.tag_read(paris_read, base_coor)
            if tag_id != -1:
                read_record[tag_id] += 1
                print >>OUT, line.strip() + "\tDG:i:" + str(tag_id)
                #print "%s\t%s\t%s\t%s\t%s" % (data[0], data[3], data[5], max_overlap, tag_id)
        line = IN.readline()
    OUT.close()
    IN.close()
    # output information of read record
    #print "=================Reads Summary==================="
    print "Total Reads: %s" % (sum(read_record), )
    print "Unsuported Reads: %s\t%.2f%%" % (read_record[0], 100.0*read_record[0]/sum(read_record))
    for idx in range(1, virus_stem.size()):
        cur_stem = virus_stem.at(idx-1)
        print "Stem: %s (%s-%s <==> %s-%s) ---- %s\t%.2f%%" % (idx, cur_stem.l_start, cur_stem.l_end, cur_stem.r_start, cur_stem.r_end, read_record[idx], 100.0*read_record[idx]/sum(read_record))

def sam_unsupport_rate(in_sam_file_name, virus_stem, base_coor):
    read_record = [0] * (virus_stem.size() + 1)
    IN = open(in_sam_file_name)
    line = IN.readline()
    while line:
        if line.startswith('@'):
            line = IN.readline()
            continue
        data = line.strip().split()
        matchGlobalSeg, matchLocalSeg, SeqLen = get_match_region( int(data[3]), data[5] )
        if len(matchGlobalSeg) == 2:
            #def __init__(self, read_id, cigar, left_start, left_end, right_start, right_end):
            paris_read = ParisRead(data[0], data[5], matchGlobalSeg[0][0], matchGlobalSeg[0][1], matchGlobalSeg[1][0], matchGlobalSeg[1][1])
            tag_id, max_overlap = virus_stem.tag_read(paris_read, base_coor)
            if tag_id != -1:
                read_record[tag_id] += 1
        line = IN.readline()
    IN.close()
    # output information of read record
    return 1.0 * read_record[0]/sum(read_record)


################################
######## 预测二级结构的函数
################################

def prepare_VARNA_colormap_file(shape_list, file_name, base=0):
    OUT = open(file_name, "w")
    for init_idx in range(base):
        print >>OUT, "%d\t%s" % (init_idx+1, "0.0")
    for idx,shape in enumerate(shape_list):
        if shape != "NULL":
            print >>OUT, "%d\t%s" % (idx+1+base, shape)
        else:
            print >>OUT, "%d\t%s" % (idx+1+base, 0.0)
    OUT.close()

def build_RNAFold_PARIS_constraint(Paris_Block, file_name, slope=-1.0, intersect=0.0):
    assert( isinstance(Paris_Block, numpy.ndarray) )
    OUT = open(file_name, "w")
    for idx in range(len(Paris_Block)):
        for idy in range(idx+1, len(Paris_Block)):
            index = Paris_Block[idx][idy]
            if index > 0:
                #print index
                print >>OUT, "E %d %d 1 %f" % (idx+1, idy+1, intersect + slope * index)
    OUT.close()

def build_Fold_PARIS_constraint(Paris_Block, file_name, slope=-1.0, intersect=0.0):
    assert( isinstance(Paris_Block, numpy.ndarray) )
    OUT = open(file_name, "w")
    for idx in range(len(Paris_Block)):
        for idy in range(len(Paris_Block)):
            index = Paris_Block[idx][idy]
            if index > 0:
                OUT.writelines( "%.3f" % (intersect + slope * index, ) )
            else:
                OUT.writelines("0")
            if idy != len(Paris_Block) - 1:
                OUT.writelines("\t")
        OUT.writelines("\n")
    OUT.close
    
def build_SHAPE_constraint(shape_list, file_name):
    SHAPE = open(file_name, 'w')
    for idx in range(len(shape_list)):
        if shape_list[idx] != "NULL":
            print >>SHAPE, "%d\t%s" % (idx+1, shape_list[idx])
        else:
            pass
    SHAPE.close()

def predictStructure_RNAfold(Seq, Shape=None, Paris=None):
    """ Predict RNA Structure using RNAfold
    Seq: sequence
    Shape: shape
    The length of each seq must be equal to shape
    """
    randID = random.randint(10000,99999)
    tmp_fa_file = "/tmp/tmp_%s.fa" % (randID, )
    tmp_shape_file = "/tmp/tmp_%s.shape" % (randID, )
    tmp_constrain_file = "/tmp/tmp_%s.const" % (randID, )
    
    # prepare tmp.fa
    FA = open(tmp_fa_file, 'w')
    print >>FA, Seq
    FA.close()
    
    if Shape is not None:
        assert( len(Seq) == len(Shape) )
        build_SHAPE_constraint(Shape, tmp_shape_file)
    
    if Paris is not  None:
        assert( len(Seq) == len(Paris) )
        build_RNAFold_PARIS_constraint(Paris, tmp_constrain_file)
    
    Fold_CMD = "RNAfold --infile=%s --noPS"
    param_list = [ tmp_fa_file ]
    if Shape is not None:
        Fold_CMD += " --shapeMethod=\"Dm2.0b-0.6\" --shape=%s"
        param_list.append(tmp_shape_file)
    if Paris is not None:
        Fold_CMD += " --commands=%s"
        param_list.append(tmp_constrain_file)    
    
    CMD = Fold_CMD % tuple(param_list)
    print CMD
    return_code, predicted_structure = commands.getstatusoutput(CMD)
    predicted_structure = predicted_structure.strip().split()[1]
    
    if os.path.isfile(tmp_constrain_file):
        os.remove(tmp_constrain_file)
    if os.path.isfile(tmp_shape_file):
        os.remove(tmp_shape_file)
    os.remove(tmp_fa_file)
    
    return predicted_structure


def predictStructure_RNAstructure(Seq, Shape=None, Paris=None):
    """Predict RNA Structure using RNAfold
    Seq: sequence
    Shape: shape
    The length of each seq must be equal to shape
    """
    
    randID = random.randint(10000,99999)
    tmp_fa_file = "/tmp/tmp_%s.fa" % (randID, )
    tmp_shape_file = "/tmp/tmp_%s.shape" % (randID, )
    tmp_ct_file = "/tmp/tmp_%s.ct" % (randID, )
    tmp_constrain_file = "/tmp/tmp_%s.const" % (randID, )
    
    # prepare tmp.fa
    FA = open(tmp_fa_file, 'w')
    print >>FA, ">seq\n%s" % (Seq, )
    FA.close()
    
    if Shape is not None:
        assert( len(Seq) == len(Shape) )
        build_SHAPE_constraint(Shape, tmp_shape_file)
    
    if Paris is not  None:
        assert( len(Seq) == len(Paris) )
        build_Fold_PARIS_constraint(Paris, tmp_constrain_file)
    
    Fold_CMD = "/Users/lee/BItools/RNAStructure/exe/Fold-smp %s %s"
    param_list = [tmp_fa_file, tmp_ct_file]
    if Shape is not None:
        Fold_CMD += " --SHAPEintercept -0.6 --SHAPEslope 2.0  --SHAPE %s"
        param_list.append(tmp_shape_file)
    if Paris is not None:
        Fold_CMD += " --experimentalPairBonus %s -xo 0.0 -xs 1.0"
        param_list.append(tmp_constrain_file)    

    CMD = Fold_CMD % tuple(param_list)
    print CMD
    os.system(CMD)
    
    return_code, return_string = commands.getstatusoutput( "grep \"ENERGY\" %s | wc -l" % (tmp_ct_file, ) )
    structure_number = int( return_string.strip() )
    structure_list = []
    ct2dot_cmd = "ct2dot %s %d /dev/stdout"
    for idx in range(structure_number):
        return_code, return_string = commands.getstatusoutput( ct2dot_cmd % (tmp_ct_file, idx+1) )
        energy = float(return_string.split()[5])
        structure = return_string.split()[8]
        structure_list.append( (energy, structure) )
    os.remove(tmp_fa_file); os.remove(tmp_ct_file);
    if os.path.isfile(tmp_constrain_file):
        os.remove(tmp_constrain_file)
    if os.path.isfile(tmp_shape_file):
        os.remove(tmp_shape_file)
    return structure_list








