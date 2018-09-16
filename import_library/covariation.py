#-*- coding:utf-8 -*-

import os, sys, commands, random


"""
from covariation import *

######### 1. input files

seqdb_fa = "/Share/home/zhangqf8/lipan/virus/select_virus_genome/flavivirus-4256.fasta"
input_sequence = "GTTGTTGATCTGTGTGAATCAGACTGCGACAGTTCGAGTTTGAAGCGAAAGCTAGCAACAGTATCAACAGGTTTTATTTTGGATTTGGAAACGAGAGTTTCTGGTCATGAAAAACCCAAAAAAGAAATCCGGAGGATTCCGGATTGTCAATATGC"
input_ss = '...((((((((((((...(((((((.(((....))))))))))(((....))).)).)))).)))))).((((........(((.(((((((....))))))))))......)))).........((((((((....))))))))..........'

######### 2. do covariation in local host

# method is one of 'R-scape/R2R/both'
output = get_covariate_bp(input_sequence, input_ss, sequence_db_fasta=seqdb_fa, method='R-scape', clean=True, tmp_outdir="/tmp/", randID=0, OUTLOG=sys.stdout)

######### 3. do covariation in cluster

handle = submit_get_covariate_bp(input_sequence, input_ss, sequence_db_fasta=seqdb_fa, cluster='Z-BNODE', cpu=20, clean=True, tmp_outdir="/Share/home/zhangqf8/.covariation")

handle.show() # show information
handle.finish() # finish?
handle.wait() # Ctrl+C terminate has no effect
output = handle.get_output() # get covariation bps
handle.clean() # or del handle to delete tmp files

######### 4. visulization

# os.environ['VARNAProg'] = "/Users/lee/Documents/VARNAv3-93.jar"
RNAStructure_highlight(input_sequence, input_ss, hg_base_list=bp2blist(output['R2R_cov']), VARNAProg=os.environ['VARNAProg'])

"""



def make_single_strand_sto(seq, structure, outFile, title="testSeq"):
    assert(len(seq) == len(structure))
    OUT = open(outFile, 'w')
    print >>OUT, "# STOCKHOLM 1.0\n"
    seq_line_head = ("%s"+" "*10) % (title, )
    ss_line_head = "#=GC SS_cons"+" "*(len(seq_line_head)-12)
    print >>OUT, seq_line_head + seq
    print >>OUT, ss_line_head + structure
    print >>OUT, "\n//"
    OUT.close()

def sto2cm(sto_file, cm_file, OUTLOG=False):
    CMD = "cmbuild %s %s > /dev/null" % (cm_file, sto_file)
    if OUTLOG:
        assert( isinstance(OUTLOG, file) )
        print >>OUTLOG, CMD
    os.system(CMD)

def calibrate_cm(cm_file, cpu_num=0, OUTLOG=False):
    CMD = "cmcalibrate %s > /dev/null" % (cm_file, )
    if cpu_num > 0:
        CMD += "--cpu %s" % (cpu_num, )
    if OUTLOG:
        assert( isinstance(OUTLOG, file) )
        print >>OUTLOG, CMD
    os.system(CMD)

def search_cm(cm_file, sequence_db_fasta, output_sto, hmm=False, cm=True, OUTLOG=False, memory=128):
    if hmm and not cm:
        CMD = "cmsearch -A %s --toponly --notextw --hmmonly --smxsize %s %s %s > /dev/null" % (output_sto, memory, cm_file, sequence_db_fasta)
    elif cm and not hmm:
        CMD = "cmsearch -A %s --toponly --notextw --nohmm --smxsize %s %s %s > /dev/null" % (output_sto, memory, cm_file, sequence_db_fasta)
    elif hmm and cm:
        CMD = "cmsearch -A %s --toponly --notextw --smxsize %s %s %s > /dev/null" % (output_sto, memory, cm_file, sequence_db_fasta)
    else:
        raise Exception("Error: hmm and cm should have 1 True at least")
    if OUTLOG:
        isinstance(OUTLOG, file)
        print >>OUTLOG, CMD
    os.system(CMD)


def clean_sto(inFile, outFile, ref_chr_id, correct_id=False, OUTLOG=False):
    """
    inFile: input sto file
    outFile: output sto file
    correct_id: 
        if True -- using number to replace transID
        if False -- raw transID name
    """
    IN = open(inFile)
    OUT = open(outFile, "w")
    line = IN.readline()
    SS_cons = ""; RF = ""; Seq = {};
    SS_label = '#=GC SS_cons'; RF_label = '#=GC RF'
    while line:
        if line.startswith('//'):
            break
        elif line.strip() == "":
            OUT.writelines("\n")
            pass
        elif line.startswith(SS_label):
            data = line.strip().split()
            SS_cons += data[2]
        elif line.startswith(RF_label):
            data = line.strip().split()
            RF += data[2]
        elif line[0] == '#':
            if correct_id:
                if not line.startswith("#=GS"):
                    OUT.writelines(line)
            else:
                OUT.writelines(line)
        else:
            data = line.strip().split()
            Seq[data[0]] = Seq.get(data[0], "") + data[1]
        line = IN.readline()
    #OUT.writelines("\n")
    if ref_chr_id not in Seq:
        raise Exception(ref_chr_id+" not in sto file")
    ##### Find Max Length
    max_len = 0
    for chr_id in Seq: 
        if len(chr_id) > max_len:
            max_len = len(chr_id)
    if correct_id:
        max_len = len(ref_chr_id)
    ##### Filter
    ref_seq = Seq[ref_chr_id]
    index = 0
    for chr_id in Seq.keys():
        if chr_id == ref_chr_id:
            continue
        index += 1
        if correct_id:
            new_chr_id = str(index)
        else:
            new_chr_id = chr_id
        OUT.writelines(new_chr_id+" "*(10+max_len-len(new_chr_id)))
        for idx in range(len(ref_seq)):
            if ref_seq[idx] not in ('.', '-', '_'):
                OUT.writelines(Seq[chr_id][idx])
        OUT.writelines("\n")
    OUT.writelines(ref_chr_id+" "*(10+max_len-len(ref_chr_id)))
    for idx in range(len(ref_seq)):
        if ref_seq[idx] not in ('.', '-', '_'):
            OUT.writelines(ref_seq[idx])
    OUT.writelines("\n")
    ##### output SS_cons
    OUT.writelines(SS_label+ " "*(10+max_len-len(SS_label)))
    for idx in range(len(ref_seq)):
        if ref_seq[idx] not in ('.', '-', '_'):
            OUT.writelines(SS_cons[idx])
        else:
            if SS_cons[idx] in set("{}()<>"):
                if OUTLOG:
                    assert( isinstance(OUTLOG, file) )
                    print >>OUTLOG, "Warning: base pairing in refchr"
    OUT.writelines("\n")
    ##### output RF
    OUT.writelines(RF_label+ " "*(10+max_len-len(RF_label)))
    for idx in range(len(ref_seq)):
        if ref_seq[idx] not in ('.', '-', '_'):
            OUT.writelines(RF[idx])
    OUT.writelines("\n\n")
    OUT.writelines("//")
    OUT.close()

def add_single_seq_to_fa(in_seqdb_fa, input_sequence, out_seqdb_fa, title="my_seq", OUTLOG=False):
    import stat
    assert(isinstance(input_sequence, str))
    CMD_CP = "cp %s %s" % (in_seqdb_fa, out_seqdb_fa)
    if OUTLOG:
        assert( isinstance(OUTLOG, file) )
        print >>OUTLOG, CMD_CP
    os.system(CMD_CP)
    os.chmod(out_seqdb_fa, os.stat(out_seqdb_fa).st_mode | stat.S_IWRITE)
    OUT = open(out_seqdb_fa, "a")
    print >>OUT, "\n>%s\n%s\n" % (title, input_sequence)
    OUT.close()


def get_seq_id_from_sto(in_sto_file, rna_id_prefix, OUTLOG=False):
    import commands
    CMD = """grep "^%s" %s | cut -d " " -f 1""" % (rna_id_prefix, in_sto_file)
    if OUTLOG:
        assert( isinstance(OUTLOG, file) )
        print >>OUTLOG, CMD
    code, rna_id = commands.getstatusoutput(CMD)
    rna_id_list = list(set(rna_id.split()))
    if len(rna_id_list) >= 2:
        raise Exception("FATAL Error: Multiple rnaID in sto file")
    return rna_id_list[0]

def makedir(dir_name):
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)

def remove_postfix(in_file_name):
    arr = in_file_name.split('.')
    if len(arr) == 1:
        return in_file_name
    else:
        return ".".join(arr[:-1])

def pureFilename(full_file_name):
    return full_file_name.split('/')[-1]


def dot2bp(dot):
    assert(isinstance(dot, str))
    bp = []
    Stack = []
    for idx, symbol in enumerate(list(dot)):
        if symbol in ('(', '[', '<', '{'):
            Stack.append(idx+1)
        elif symbol in (')', ']', '>', '}'):
            bp.append( (Stack[-1], idx+1) )
            del Stack[-1]
        else:
            pass
    bp.sort(key=lambda x:x[0])
    return bp

def covariation_from_R2R_sto(inSto, OUTLOG=False):
    import commands
    covariat_bp_list = []; half_flips_list = []
    CMD_COV = """grep "#=GC cov_SS_cons" %s""" % (inSto, )
    CMD_SS = """grep "#=GC SS_cons" %s""" % (inSto, )
    if OUTLOG:
        assert( isinstance(OUTLOG, file) )
        print >>OUTLOG, CMD_COV
    code, COV_line = commands.getstatusoutput(CMD_COV)
    if OUTLOG:
        assert( isinstance(OUTLOG, file) )
        print >>OUTLOG, CMD_SS
    code, SS_line = commands.getstatusoutput(CMD_SS)
    if COV_line == "":
        raise Exception("Bad R2R sto file: no #=GC cov_SS_cons line")
    ss = SS_line.split()[2]
    bp_dict = dict(dot2bp(ss))
    items = COV_line.split()
    bps = 0
    for idx in range(len(items[2])):
        if items[2][idx] == '2' and idx+1 in bp_dict:
            covariat_bp_list.append( (idx+1, bp_dict[idx+1]) )
        elif items[2][idx] == '2' and idx+1 in bp_dict:
            half_flips_list.append( (idx+1, bp_dict[idx+1]) )
    return covariat_bp_list, half_flips_list


def read_covariation_sites(Rscape_out_file):
    covariat_bp = []
    IN = open(Rscape_out_file)
    line = IN.readline()
    while line:
        if line.startswith('*'):
            data = line.strip().split()
            left = int(data[1])
            right = int(data[2])
            covariat_bp.append((left, right))
        line = IN.readline()
    IN.close()
    covariat_bp.sort(key=lambda x:x[0])
    return covariat_bp



def Rscape_covariation(inSto, work_dir, OUTLOG=False):
    work_dir = work_dir.rstrip('/')
    CMD = "R-scape --outdir %s %s > /dev/null" % (work_dir, inSto)
    if OUTLOG:
        assert( isinstance(OUTLOG, file) )
        print >>OUTLOG, CMD
    os.system(CMD)
    Rscape_out = work_dir + "/" +remove_postfix(pureFilename(inSto)) + ".out"
    covartate_bp = read_covariation_sites(Rscape_out)
    return covartate_bp

# Rscape_covariation(inSto="/Share/home/zhangqf8/lipan/tmp/virus59/search_clean_sto_file/segment_0_search.clean.sto", work_dir="/tmp/R-scape")


def R2R_covariation(inSto, tmpSto, OUTLOG=False):
    CMD = "r2r --GSC-weighted-consensus %s %s 3 0.97 0.9 0.75 4 0.97 0.9 0.75 0.5 0.15" % (inSto, tmpSto)
    if OUTLOG:
        assert( isinstance(OUTLOG, file) )
        print >>OUTLOG, CMD
    os.system(CMD)
    return covariation_from_R2R_sto(tmpSto, OUTLOG=OUTLOG)

def remove_dir_file(dir_name, OUTLOG=False):
    if os.path.exists(dir_name):
        CMD = "rm -rf %s" % (dir_name, )
        if OUTLOG:
            assert( isinstance(OUTLOG, file) )
            print >>OUTLOG, CMD
        os.system(CMD)


def get_covariate_bp(input_sequence, input_ss, sequence_db_fasta, method='R-scape', clean=True, tmp_outdir="/tmp/", randID=0, memory=128, OUTLOG=sys.stdout):
    """
        Get covariate base pairs
    """
    assert(len(input_sequence) == len(input_ss))
    assert(method in ('R-scape', 'R2R', 'both'))
    assert(os.path.exists(sequence_db_fasta))
    if OUTLOG: assert(isinstance(OUTLOG, file))
    
    if randID == 0:
        randID = random.randint(10000,99999)
    ## Make some tmp files
    tmp_dir_name = tmp_outdir.rstrip('/') + "/tmp_dir_%s" % (randID, )
    
    tmp_init_sto_file = "%s/tmp.sto" % (tmp_dir_name, )
    tmp_cm_file = "%s/tmp.cm" % (tmp_dir_name, )
    tmp_search_sto_file = "%s/tmp_search.sto" % (tmp_dir_name, )
    tmp_seqdb_fa_file = "%s/tmp_seqdb.fasta" % (tmp_dir_name, )
    tmp_clean_sto_file = "%s/tmp_clean.sto" % (tmp_dir_name, )
    tmp_Rscape_dir = "%s/R-scape" % (tmp_dir_name, )
    tmp_R2R_sto_file = "%s/R2R.sto" % (tmp_dir_name, )
    
    my_seq_name = "my_seq_%s" % (randID, )
    
    output = {}
    
    ## 0. 
    makedir(tmp_dir_name)
    
    ## 1. make a cm file
    if OUTLOG:
        print >>OUTLOG, "Step 1 -- make a cm file"
    make_single_strand_sto(input_sequence, input_ss, tmp_init_sto_file, title=my_seq_name)
    sto2cm(tmp_init_sto_file, tmp_cm_file, OUTLOG=OUTLOG)
    calibrate_cm(tmp_cm_file, cpu_num=0, OUTLOG=OUTLOG)
    
    ## 2. search homologous sequence
    if OUTLOG:
        print >>OUTLOG, "Step 2 -- search homologous sequence"
    add_single_seq_to_fa(sequence_db_fasta, input_sequence, tmp_seqdb_fa_file, title=my_seq_name, OUTLOG=OUTLOG)
    search_cm(tmp_cm_file, tmp_seqdb_fa_file, tmp_search_sto_file, hmm=False, cm=True, memory=memory, OUTLOG=OUTLOG)
    
    ## 3. clean sto
    if OUTLOG:
        print >>OUTLOG, "Step 3 -- clean sto"
    search_chr_id = get_seq_id_from_sto(tmp_search_sto_file, my_seq_name, OUTLOG=OUTLOG)
    clean_sto(tmp_search_sto_file, tmp_clean_sto_file, search_chr_id, correct_id=True, OUTLOG=OUTLOG)
    
    ## 4. Covariation
    if OUTLOG:
        print >>OUTLOG, "Step 4 -- Covariation"
    if method == "R-scape":
        makedir(tmp_Rscape_dir)
        output['R-scape'] = Rscape_covariation(tmp_clean_sto_file, tmp_Rscape_dir, OUTLOG=OUTLOG)
    elif method == "R2R":
        output['R2R_cov'], output['R2R_flip'] = R2R_covariation(tmp_clean_sto_file, tmp_R2R_sto_file, OUTLOG=OUTLOG)
    else:
        makedir(tmp_Rscape_dir)
        output['R-scape'] = Rscape_covariation(tmp_clean_sto_file, tmp_Rscape_dir, OUTLOG=OUTLOG)
        output['R2R_cov'], output['R2R_flip'] = R2R_covariation(tmp_clean_sto_file, tmp_R2R_sto_file, OUTLOG=OUTLOG)
    
    ## 5. clean dir
    if clean:
        remove_dir_file(tmp_dir_name, OUTLOG=OUTLOG)
    
    return output




############################
###  bsub to clusters
############################

from cluster import *

class COV_SUB_HANDLE(JOB_HANDLE):
    def __init__(self, job_handle, script_file, output_file, log_file, error_file, get_function):
        self.job_handle = job_handle
        self.script = script_file
        self.output = output_file
        self.log_file = log_file
        self.error_file = error_file
        self.get_function = get_function
    def finish(self):
        return self.job_handle.has_finish()
    def wait(self):
        return self.job_handle.wait()
    def get_output(self):
        if not self.finish():
            print "Warning: not finished yet"
            return
        return self.get_function(self.output)
        """
        IN = open(self.output)
        co_bps = IN.readlines()
        IN.close()
        output = {}
        output['R-scape'], output['R2R_cov'], output['R2R_flip'] = eval(co_bps[0]), eval(co_bps[1]), eval(co_bps[2])
        return output
        """
    def __del__(self):
        import sys
        if not self.finish():
            self.job_handle.kill()
            print >>sys.stderr, "Warning: Job not finish"
        remove_dir_file(self.script)
        remove_dir_file(self.output)
        remove_dir_file(self.log_file)
        remove_dir_file(self.error_file)
    def clean(self):
        import sys
        if not self.finish():
            self.job_handle.kill()
            print >>sys.stderr, "Warning: Job not finish"
        remove_dir_file(self.script)
        remove_dir_file(self.output)
        remove_dir_file(self.log_file)
        remove_dir_file(self.error_file)
    def show(self):
        print "job_id: %s" % (self.job_handle.job_id, )
        print "script file: %s" % (self.script, )
        print "output file: %s" % (self.output, )
        print "log file: %s" % (self.log_file, )
        print "error file: %s" % (self.error_file, )

def submit_get_covariate_bp(input_sequence, input_ss, sequence_db_fasta, cluster='Z-ZQF', clean=True, tmp_outdir="/tmp/", randID=0, memory=128, cpu=20):
    def get_output(inFile):
        IN = open(inFile)
        co_bps = IN.readlines()
        IN.close()
        output = {}
        output['R-scape'], output['R2R_cov'], output['R2R_flip'] = eval(co_bps[0]), eval(co_bps[1]), eval(co_bps[2])
        return output
    import time
    assert(os.path.exists(sequence_db_fasta))
    
    covariation_home = os.environ['HOME']+'/.covariation'
    makedir(covariation_home)
    if randID == 0:
        randID = random.randint(10000,99999)
    
    tmp_script = covariation_home + "/%s.py" % (randID, )
    tmp_output = covariation_home + "/%s_output.py" % (randID, )
    tmp_log = covariation_home + "/%s.log" % (randID, )
    tmp_error = covariation_home + "/%s.error" % (randID, )
    
    Python_Script = """
#!/usr/bin/env python
from covariation import *
    
input_sequence = "%s"
input_ss = "%s"
sequence_db_fasta = "%s"
tmp_outdir = "%s"
clean = %s
randID = "%s"
OUTLOG = sys.stdout
    
output = get_covariate_bp(input_sequence, input_ss, sequence_db_fasta, method='both', tmp_outdir=tmp_outdir, clean=clean, randID=randID,  memory=%s, OUTLOG=sys.stdout)
    
OUT = open("%s", "w")
print >>OUT, output['R-scape']
print >>OUT, output['R2R_cov']
print >>OUT, output['R2R_flip']
OUT.close()
    
    """ % (input_sequence, input_ss, sequence_db_fasta, tmp_outdir, clean, randID, memory, tmp_output, get_output)
    
    print "=====>", randID, "<====="
    OUT = open(tmp_script, 'w')
    print >>OUT, Python_Script
    OUT.close()
    
    command = "python " + tmp_script
    my_job_handle = new_job(command, queue=cluster, cpu=cpu, job_name="cov_"+str(randID), log_file=tmp_log, error_file=tmp_error)
    my_job_handle.submit()
        
    return COV_SUB_HANDLE(my_job_handle, tmp_script, tmp_output, tmp_log, tmp_error)




def get_covariate_bp_ii(input_align_sequence_1, input_align_sequence_2, input_ss, sequence_db_fasta, method='R-scape', clean=True, tmp_outdir="/tmp/", randID=0, memory=128, OUTLOG=sys.stdout):
    """
        Get covariate base pairs
    """
    assert(len(input_align_sequence_1) == len(input_align_sequence_2) == len(input_ss))
    assert(method in ('R-scape', 'R2R', 'both'))
    assert(os.path.exists(sequence_db_fasta))
    if OUTLOG: assert(isinstance(OUTLOG, file))
    
    if randID == 0:
        randID = random.randint(10000,99999)
    ## Make some tmp files
    tmp_dir_name = tmp_outdir.rstrip('/') + "/tmp_dir_%s" % (randID, )
    
    tmp_init_sto_file = "%s/tmp.sto" % (tmp_dir_name, )
    tmp_cm_file = "%s/tmp.cm" % (tmp_dir_name, )
    tmp_search_sto_file = "%s/tmp_search.sto" % (tmp_dir_name, )
    tmp_seqdb_fa_file = "%s/tmp_seqdb.fasta" % (tmp_dir_name, )
    tmp_clean_sto_file_1 = "%s/tmp_clean_1.sto" % (tmp_dir_name, )
    tmp_clean_sto_file_2 = "%s/tmp_clean_2.sto" % (tmp_dir_name, )
    tmp_Rscape_dir_1 = "%s/R-scape-1" % (tmp_dir_name, )
    tmp_Rscape_dir_2 = "%s/R-scape-2" % (tmp_dir_name, )
    tmp_R2R_sto_file_1 = "%s/R2R_1.sto" % (tmp_dir_name, )
    tmp_R2R_sto_file_2 = "%s/R2R_2.sto" % (tmp_dir_name, )
    
    my_seq1_name = "my_seq_1_%s" % (randID, )
    my_seq2_name = "my_seq_2_%s" % (randID, )
    
    output = {1:{}, 2:{}}
    
    ## 0. 
    makedir(tmp_dir_name)
    
    ## 1. make a cm file
    if OUTLOG:
        print >>OUTLOG, "Step 1 -- make a cm file"
    build_sto_from_dyalign(input_align_sequence_1, input_align_sequence_2, input_ss, tmp_init_sto_file)
    #make_single_strand_sto(input_sequence, input_ss, tmp_init_sto_file, title=my_seq_name)
    sto2cm(tmp_init_sto_file, tmp_cm_file, OUTLOG=OUTLOG)
    calibrate_cm(tmp_cm_file, cpu_num=0, OUTLOG=OUTLOG)
    
    ## 2. search homologous sequence
    if OUTLOG:
        print >>OUTLOG, "Step 2 -- search homologous sequence"
    add_single_seq_to_fa(sequence_db_fasta, input_align_sequence_1.replace('-',''), tmp_seqdb_fa_file, title=my_seq1_name, OUTLOG=OUTLOG)
    add_single_seq_to_fa(tmp_seqdb_fa_file, input_align_sequence_2.replace('-',''), tmp_seqdb_fa_file, title=my_seq2_name, OUTLOG=OUTLOG)
    search_cm(tmp_cm_file, tmp_seqdb_fa_file, tmp_search_sto_file, hmm=False, cm=True, memory=memory, OUTLOG=OUTLOG)
    
    ## 3. clean sto
    if OUTLOG:
        print >>OUTLOG, "Step 3 -- clean sto"
    search_chr_id_1 = get_seq_id_from_sto(tmp_search_sto_file, my_seq1_name, OUTLOG=OUTLOG)
    search_chr_id_2 = get_seq_id_from_sto(tmp_search_sto_file, my_seq2_name, OUTLOG=OUTLOG)
    clean_sto(tmp_search_sto_file, tmp_clean_sto_file_1, search_chr_id_1, correct_id=True, OUTLOG=OUTLOG)
    clean_sto(tmp_search_sto_file, tmp_clean_sto_file_2, search_chr_id_2, correct_id=True, OUTLOG=OUTLOG)
    
    ## 4. Covariation
    if OUTLOG:
        print >>OUTLOG, "Step 4 -- Covariation"
    if method == "R-scape":
        makedir(tmp_Rscape_dir_1)
        makedir(tmp_Rscape_dir_2)
        output[1]['R-scape'] = Rscape_covariation(tmp_clean_sto_file_1, tmp_Rscape_dir_1, OUTLOG=OUTLOG)
        output[2]['R-scape'] = Rscape_covariation(tmp_clean_sto_file_2, tmp_Rscape_dir_2, OUTLOG=OUTLOG)
    elif method == "R2R":
        output[1]['R2R_cov'], output[1]['R2R_flip'] = R2R_covariation(tmp_clean_sto_file_1, tmp_R2R_sto_file_1, OUTLOG=OUTLOG)
        output[2]['R2R_cov'], output[2]['R2R_flip'] = R2R_covariation(tmp_clean_sto_file_2, tmp_R2R_sto_file_2, OUTLOG=OUTLOG)
    else:
        makedir(tmp_Rscape_dir_1)
        makedir(tmp_Rscape_dir_2)
        output[1]['R-scape'] = Rscape_covariation(tmp_clean_sto_file_1, tmp_Rscape_dir_1, OUTLOG=OUTLOG)
        output[2]['R-scape'] = Rscape_covariation(tmp_clean_sto_file_2, tmp_Rscape_dir_2, OUTLOG=OUTLOG)
        output[1]['R2R_cov'], output[1]['R2R_flip'] = R2R_covariation(tmp_clean_sto_file_1, tmp_R2R_sto_file_1, OUTLOG=OUTLOG)
        output[2]['R2R_cov'], output[2]['R2R_flip'] = R2R_covariation(tmp_clean_sto_file_2, tmp_R2R_sto_file_2, OUTLOG=OUTLOG)
    
    ## 5. clean dir
    if clean:
        remove_dir_file(tmp_dir_name, OUTLOG=OUTLOG)
    
    return output



def submit_get_covariate_bp_ii(input_align_sequence_1, input_align_sequence_2, input_ss, sequence_db_fasta, cluster='Z-ZQF', clean=True, tmp_outdir="/tmp/", randID=0, memory=128, cpu=20):
    def get_output(inFile):
        co_bps = open(inFile).readlines()
        output = {1:{}, 2:{}}
        output[1]['R-scape'] = eval(co_bps[0])
        output[2]['R-scape'] = eval(co_bps[1])
        output[1]['R2R_cov'] = eval(co_bps[2])
        output[1]['R2R_flip'] = eval(co_bps[3])
        output[2]['R2R_cov'] = eval(co_bps[4])
        output[2]['R2R_flip'] = eval(co_bps[5])
        return output
    import time
    assert(os.path.exists(sequence_db_fasta))
    
    covariation_home = os.environ['HOME']+'/.covariation'
    makedir(covariation_home)
    if randID == 0:
        randID = random.randint(10000,99999)
    
    tmp_script = covariation_home + "/%s.py" % (randID, )
    tmp_output = covariation_home + "/%s_output.py" % (randID, )
    tmp_log = covariation_home + "/%s.log" % (randID, )
    tmp_error = covariation_home + "/%s.error" % (randID, )
    
    Python_Script = """
#!/usr/bin/env python
from covariation import *

input_align_sequence_1 = "%s"
input_align_sequence_2 = "%s"
input_ss = "%s"
sequence_db_fasta = "%s"
tmp_outdir = "%s"
clean = %s
randID = "%s"
memory = %s
OUTLOG = sys.stdout

output = get_covariate_bp_ii(input_align_sequence_1, input_align_sequence_2, input_ss, sequence_db_fasta, method='both', clean=True, tmp_outdir=tmp_outdir, randID=randID, memory=memory, OUTLOG=sys.stdout)

OUT = open("%s", "w")
print >>OUT, output[1]['R-scape']
print >>OUT, output[2]['R-scape']
print >>OUT, output[1]['R2R_cov']
print >>OUT, output[1]['R2R_flip']
print >>OUT, output[2]['R2R_cov']
print >>OUT, output[2]['R2R_flip']
OUT.close()

    """ % (input_align_sequence_1, input_align_sequence_2, input_ss, sequence_db_fasta, tmp_outdir, clean, randID, memory, tmp_output)
    
    print "=====>", randID, "<====="
    OUT = open(tmp_script, 'w')
    print >>OUT, Python_Script
    OUT.close()
    
    command = "python " + tmp_script
    my_job_handle = new_job(command, queue=cluster, cpu=cpu, job_name="cov_"+str(randID), log_file=tmp_log, error_file=tmp_error)
    my_job_handle.submit()
    
    return COV_SUB_HANDLE(my_job_handle, tmp_script, tmp_output, tmp_log, tmp_error, get_output)

def build_sto_from_dyalign(align_seq_1, align_seq_2, align_ss, outFile):
    assert( len(align_seq_1) == len(align_seq_2) == len(align_ss) )
    OUT = open(outFile, 'w')
    print >>OUT, "# STOCKHOLM 1.0"
    print >>OUT, "test_1            %s" % (align_seq_1, )
    print >>OUT, "test_2            %s" % (align_seq_2, )
    print >>OUT, "#=GC SS_cons      %s" % (align_ss, )
    print >>OUT, "\n//"
    OUT.close()


############################
###  Visualization
############################

def RNAStructure_highlight(seq, ss, hg_base_list=[], VARNAProg="VARNAv3-93.jar"):
    """
    高亮其中的某些碱基
    Example: Plot_RNAStructure_highlight("AGCTGGGTTTCCCGATT", "....(((...)))....", hg_base_list=[1,2,3,4])
    """
    assert(len(seq) == len(ss))
    CMD = "java -cp "+VARNAProg+" fr.orsay.lri.varna.applications.VARNAcmd -sequenceDBN %s -structureDBN \"%s\" -backbone \"#FF0000\" " % (seq, ss)
    if hg_base_list:
        CMD += "-basesStyle1 \"label=#FF0000\" "
        hg_base_list = [`item` for item in hg_base_list]
        CMD += "-applyBasesStyle1on \"%s\" " % (",".join(hg_base_list), )
    print CMD

def bp2blist(bp):
    return [b[0] for b in bp] + [b[1] for b in bp]



