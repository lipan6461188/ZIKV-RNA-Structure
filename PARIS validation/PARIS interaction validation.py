
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
# count common and uniq dg number of 59 and 766
###########
key59, key766 = ['KU501215.1', 'AY632535.2']

uniq_59_short, common_59_short, uniq_766_short, common_766_short = classify_dg(DG59_short, DG766_short, tolerate=0)
print len(uniq_59_short), len(common_59_short), len(uniq_766_short), len(common_766_short)


uniq_59_long, common_59_long, uniq_766_long, common_766_long = classify_dg(DG59_long, DG766_long, tolerate=0)
print len(uniq_59_long), len(common_59_long), len(uniq_766_long), len(common_766_long)


###########
# Permutation long-range interaction
###########

def shuffle_single_interaction(genome_length, dg):
    import random
    dg_length = dg[3] - dg[0] + 1
    rand_start = random.randint(1, genome_length-dg_length)
    return (rand_start, dg[1]-dg[0]+rand_start, dg[2]-dg[0]+rand_start, dg[3]-dg[0]+rand_start)

def shuffle_interactions(genome_length, dgList):
    rand_dgList = []
    for dg in dgList:
        rand_dgList.append( shuffle_single_interaction(genome_length, dg) )
    return rand_dgList

def calc_fold_energy(dgList, sequence, expand):
    length = len(sequence)
    energy = []
    for dg in dgList:
        assert( dg[3] <= length )
        s_1, e_1, s_2, e_2 = max(1,dg[0]-expand), min(length, dg[1]+expand), max(1,dg[2]-expand), min(length, dg[3]+expand)
        seq_1 = sequence[s_1-1:e_1]
        seq_2 = sequence[s_2-1:e_2]
        try:
            mfe = structure.bi_fold(seq_1, seq_2, mfe=False)[1][0][0]
        except:
            continue
        energy.append(mfe)
    return energy


true_mfe_59 = calc_fold_energy(DG59_long, Sequence['KU501215.1'], expand=10)
true_mfe_766 = calc_fold_energy(DG766_long, Sequence['AY632535.2'], expand=10)


#### Strategy 1 -- Shuffle all interactions

rand_mean_list_59 = []
for i in range(100):
    print len(rand_mean_list_59)
    rand_DG59_long = shuffle_interactions(10802, DG59_long)
    rand_mfe = calc_fold_energy(rand_DG59_long, Sequence['KU501215.1'], expand=20)
    rand_mean_list_59.append( sum(rand_mfe)/len(rand_mfe) )

rand_mean_list_766 = []
for i in range(100):
    print len(rand_mean_list_766)
    rand_DG766_long = shuffle_interactions(10795, DG766_long)
    rand_mfe = calc_fold_energy(rand_DG766_long, Sequence['AY632535.2'], expand=20)
    rand_mean_list_766.append( sum(rand_mfe)/len(rand_mfe) )


#### Plot

rand_mean_list_766 = [-24.71573033707864, -24.661797752808983, -24.660674157303383, -24.63707865168541, -24.617977528089888, -24.60505617977528, -24.55505617977528, -24.541011235955054, -24.514044943820213, -24.291573033707877, -24.280337078651677, -24.253932584269663, -24.242696629213476, -24.240449438202248, -24.22977528089887, -24.215730337078643, -24.203370786516853, -24.198876404494385, -24.183707865168532, -24.17584269662923, -24.16516853932585, -24.091573033707856, -24.08988764044943, -24.08764044943819, -24.08370786516853, -24.068539325842707, -24.05842696629213, -24.052247191011226, -24.035955056179763, -24.023595505617973, -24.015730337078658, -24.000000000000014, -23.99213483146067, -23.98651685393257, -23.976404494382034, -23.974157303370788, -23.96573033707866, -23.9623595505618, -23.946629213483135, -23.924719101123603, -23.923595505617975, -23.91516853932583, -23.907865168539328, -23.90280898876404, -23.89101123595506, -23.890449438202257, -23.889887640449448, -23.887640449438198, -23.887078651685375, -23.869662921348322, -23.86348314606741, -23.859550561797754, -23.835955056179788, -23.814044943820228, -23.794943820224724, -23.783707865168523, -23.744943820224712, -23.744382022471925, -23.732584269662926, -23.7247191011236, -23.724719101123597, -23.708988764044953, -23.699999999999996, -23.69831460674157, -23.696629213483153, -23.692134831460674, -23.68764044943821, -23.68764044943821, -23.687640449438195, -23.661797752809, -23.661235955056178, -23.659550561797747, -23.64719101123596, -23.644382022471913, -23.635955056179775, -23.625842696629217, -23.617415730337083, -23.614606741573027, -23.598314606741564, -23.594943820224728, -23.582022471910104, -23.577528089887647, -23.57640449438203, -23.57134831460673, -23.5629213483146, -23.560674157303367, -23.55842696629214, -23.55842696629213, -23.557865168539323, -23.555617977528087, -23.55280898876405, -23.550561797752803, -23.546629213483143, -23.54606741573032, -23.541011235955054, -23.539325842696623, -23.538764044943818, -23.535955056179773, -23.53314606741573, -23.529213483146055, -23.527528089887632, -23.526966292134844, -23.525280898876396, -23.51797752808989, -23.51404494382022, -23.50112359550561, -23.49550561797754, -23.492696629213476, -23.49213483146066, -23.491011235955074, -23.49101123595505, -23.490449438202234, -23.473595505617972, -23.46966292134833, -23.464606741573032, -23.45786516853932, -23.452247191011242, -23.451123595505628, -23.41629213483147, -23.415730337078656, -23.411235955056195, -23.401685393258443, -23.3932584269663, -23.39101123595505, -23.387640449438198, -23.38314606741573, -23.372471910112367, -23.36460674157303, -23.359550561797747, -23.347752808988773, -23.34662921348314, -23.34269662921349, -23.338202247191024, -23.330898876404486, -23.33033707865168, -23.329213483146074, -23.318539325842686, -23.30112359550562, -23.299999999999994, -23.28651685393259, -23.27752808988763, -23.276966292134833, -23.276404494382017, -23.27415730337077, -23.273595505617973, -23.273033707865167, -23.267415730337074, -23.251123595505618, -23.251123595505597, -23.23033707865169, -23.19494382022472, -23.1808988764045, -23.15617977528089, -23.15112359550562, -23.150561797752804, -23.144382022471888, -23.132584269662917, -23.130337078651685, -23.11629213483147, -23.110674157303375, -23.11011235955055, -23.0870786516854, -23.0814606741573, -23.055617977528104, -23.05505617977527, -23.045505617977515, -23.03426966292136, -23.009550561797745, -23.003932584269663, -22.93539325842696, -22.90730337078651, -22.893820224719107, -22.868539325842683, -22.851685393258414, -22.762359550561808, -22.73651685393258, -22.67078651685393]
rand_mean_list_59 = [-24.407826086956515, -24.313043478260862, -24.120434782608704, -24.00565217391306, -23.973913043478277, -23.970869565217377, -23.952173913043477, -23.939999999999994, -23.9395652173913, -23.928695652173904, -23.774347826086967, -23.76739130434783, -23.75869565217391, -23.752173913043478, -23.729130434782594, -23.68739130434784, -23.665217391304353, -23.6613043478261, -23.63434782608695, -23.629130434782606, -23.57304347826087, -23.570869565217393, -23.55347826086957, -23.530434782608694, -23.52521739130436, -23.523913043478256, -23.519565217391303, -23.510869565217405, -23.50695652173914, -23.485217391304342, -23.481739130434775, -23.479565217391297, -23.47434782608697, -23.473913043478262, -23.47347826086956, -23.470434782608688, -23.460000000000008, -23.460000000000004, -23.458695652173922, -23.444782608695654, -23.421739130434798, -23.415652173913056, -23.414347826086956, -23.39999999999999, -23.385217391304366, -23.382608695652177, -23.375652173913043, -23.37130434782609, -23.362173913043485, -23.356086956521736, -23.353478260869565, -23.340000000000003, -23.326521739130428, -23.2913043478261, -23.288695652173914, -23.286956521739125, -23.283478260869565, -23.278260869565212, -23.273043478260885, -23.258260869565216, -23.253043478260864, -23.252608695652167, -23.247826086956522, -23.244782608695647, -23.242608695652155, -23.2304347826087, -23.226956521739123, -23.226086956521748, -23.215652173913043, -23.210434782608715, -23.20869565217393, -23.204347826086945, -23.187391304347816, -23.18521739130433, -23.181739130434774, -23.177826086956525, -23.17434782608695, -23.17043478260869, -23.165217391304356, -23.163478260869564, -23.160869565217393, -23.16000000000002, -23.159565217391307, -23.155217391304358, -23.14913043478261, -23.13999999999998, -23.136521739130462, -23.12260869565217, -23.111304347826085, -23.10565217391304, -23.10565217391304, -23.105652173913036, -23.102173913043465, -23.10173913043479, -23.09826086956522, -23.09434782608694, -23.09391304347826, -23.086086956521736, -23.084347826086944, -23.081739130434773, -23.078260869565216, -23.076086956521742, -23.05260869565217, -23.04173913043481, -23.034347826086975, -23.0304347826087, -23.02347826086957, -23.015652173913054, -23.013913043478254, -23.013478260869565, -23.010869565217376, -22.996956521739122, -22.995217391304344, -22.99173913043479, -22.9908695652174, -22.971739130434777, -22.96695652173913, -22.96652173913043, -22.95434782608695, -22.952608695652188, -22.95130434782609, -22.94913043478261, -22.943913043478272, -22.94391304347825, -22.916521739130435, -22.905217391304337, -22.884347826086948, -22.88260869565218, -22.864347826086963, -22.86347826086956, -22.850434782608694, -22.844782608695663, -22.83000000000001, -22.827391304347827, -22.822608695652175, -22.812608695652177, -22.811304347826088, -22.803478260869575, -22.793478260869577, -22.781304347826083, -22.78, -22.766521739130432, -22.762173913043476, -22.723478260869566, -22.687826086956512, -22.682173913043478, -22.67130434782608, -22.659565217391307, -22.64391304347827, -22.635652173913027, -22.596956521739127, -22.592608695652192, -22.562608695652184, -22.53347826086956, -22.509565217391312, -22.484347826086978, -22.446521739130436, -22.357826086956514, -22.34739130434784, -22.15, -22.11869565217391]
true_mfe_59 = -31.149
true_mfe_766 = -32.137


def shuffle_plot(input_data_list, true_score):
    table = pd.DataFrame( [(it, '766') for it in rand_mean_list_766], columns=['mfe', 'type_id'] )
    plt.figure(figsize=(3,10))
    #ax = sns.boxplot(x="type_id", y="mfe", data=table)
    ax = sns.violinplot(x="type_id", y="mfe", data=table)
    ax = sns.stripplot(x="type_id", y="mfe", data=table, jitter=True, color="#4C72B0", size=8)
    ax.plot(0, true_score, marker='o', markersize=3, color="red")
    plt.ylim(-35, -20)


shuffle_plot(rand_mean_list_766, true_mfe_766)
plt.savefig("/tmp/766.pdf")
plt.close()

shuffle_plot(rand_mean_list_59, true_mfe_59)
plt.savefig("/tmp/59.pdf")
plt.close()

####
#### Strategy 2 -- Shuffle each interaction each time
####

true_mean_list_59 = []
for dg in DG59_long:
    mfe_list = calc_fold_energy([dg], Sequence[key59], expand=10)
    true_mean_list_59.append( mfe_list[0] )

rand_mean_list_59 = []
for dg in DG59_long:
    rand_dg_list = [ shuffle_single_interaction(10802, dg) for idx in range(100) ]
    try:
        mfe_list = calc_fold_energy(rand_dg_list, Sequence[key59], expand=10)
    except:
        continue
    rand_mean_list_59.append(mfe_list)

OUT = open("/tmp/59.interaction.txt", 'w')
for idx in range(len(true_mfe_59)):
    OUT.writelines("%s\t%s\n" % (true_mfe_59[idx], sorted(rand_mean_list_59[idx], reverse=False)))

OUT.close()


p_values = []
for idx in range(len(true_mfe_59)):
    p_values.append( permutate_pValue(rand_mean_list_59[idx], true_mfe_59[idx], mode='low') )

p_values = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.06, 0.06, 0.06, 0.06, 0.06, 0.07, 0.07, 0.07, 0.07, 0.07, 0.07, 0.08, 0.08, 0.08, 0.08, 0.08, 0.09, 0.09, 0.09, 0.09, 0.09, 0.09, 0.09, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.11, 0.11, 0.11, 0.11, 0.11, 0.12, 0.12, 0.12, 0.12, 0.12, 0.13, 0.13, 0.13, 0.13, 0.14, 0.14, 0.14, 0.14, 0.15, 0.15, 0.15, 0.16, 0.16, 0.16, 0.16, 0.16, 0.16, 0.17, 0.17, 0.17, 0.18, 0.18, 0.18, 0.18, 0.18, 0.19, 0.19, 0.2, 0.2, 0.2, 0.21, 0.21, 0.21, 0.21, 0.21, 0.22, 0.23, 0.23, 0.24, 0.24, 0.24, 0.25, 0.25, 0.26, 0.26, 0.26, 0.28, 0.29, 0.29, 0.3, 0.3, 0.31, 0.31, 0.32, 0.32, 0.33, 0.33, 0.34, 0.34, 0.34, 0.36, 0.37, 0.38, 0.38, 0.4, 0.4, 0.4, 0.41, 0.41, 0.42, 0.42, 0.43, 0.44, 0.44, 0.44, 0.47, 0.48, 0.48, 0.64, 0.66, 0.68, 0.69, 0.72, 0.76, 0.81, 0.81, 0.84, 0.84, 0.85]
sns.boxplot(y=p_values)
plt.show()


rand_mean_list_766 = []
for dg in DG766_long:
    rand_dg_list = [ shuffle_single_interaction(10795, dg) for idx in range(100) ]
    try:
        mfe_list = calc_fold_energy(rand_dg_list, Sequence[key766], expand=10)
    except:
        continue
    rand_mean_list_766.append(mfe_list)


p_values = []
for idx in range(len(true_mfe_766)):
    p_values.append( permutate_pValue(rand_mean_list_766[idx], true_mfe_766[idx], mode='low') )


####
#### Strategy 3 -- Shuffle the interacted sequence
####

def permutate_dg(seq_1, seq_2, times=100):
    import random,structure
    true_mfe = structure.bi_fold(seq_1, seq_2, mfe=False)[1][0][0]
    random_energy_list = []
    while len(random_energy_list) < times:
        randSeq_1 = list(seq_1)
        random.shuffle(randSeq_1)
        randSeq_1 = "".join(randSeq_1)
        
        randSeq_2 = list(seq_2)
        random.shuffle(randSeq_2)
        randSeq_2 = "".join(randSeq_2)
        
        try:
            cur_mfe = structure.bi_fold(randSeq_1, randSeq_2, mfe=False)[1][0][0]
        except:
            print "Error Occurred"
            continue
        
        random_energy_list.append(cur_mfe)
    
    return true_mfe, random_energy_list

def permutate_single_dg(dg, sequence, expand=20):
    length = len(sequence)
    
    assert( dg[3] <= length )
    s_1, e_1, s_2, e_2 = max(1,dg[0]-expand), min(length, dg[1]+expand), max(1,dg[2]-expand), min(length, dg[3]+expand)
    
    seq_1 = sequence[s_1-1:e_1]
    seq_2 = sequence[s_2-1:e_2]
    
    true_mfe, random_energy_list = permutate_dg(seq_1, seq_2, times=100)
    
    return true_mfe, random_energy_list

def permutate_dgList(dg_list, sequence, expand=20):
    dg_mfe_list = []
    for dg in dg_list:
        dg_mfe_list.append( permutate_single_dg(dg, sequence, expand) )
    return dg_mfe_list

def save_permutate_mfe(permutate_list, outFile):
    OUT = open(outFile, 'w')
    for true_v,permutate_vList in permutate_list:
        print >>OUT, "%s\t%s\t%s" % (true_v, round(sum(permutate_vList)/len(permutate_vList),3), tools.permutate_pValue(permutate_vList,true_v))
    OUT.close()


permute_dg_59 = permutate_dgList(DG59_long, Sequence[key59], expand=20)
permute_dg_766 = permutate_dgList(DG766_long, Sequence[key766], expand=20)

save_permutate_mfe(permute_dg_59, "/tmp/59_dg.txt")
save_permutate_mfe(permute_dg_766, "/tmp/766_dg.txt")


"""
xx1=pd.read_csv("/tmp/59_dg.txt",sep="\t",header=None)
scipy.stats.ttest_rel(xx1.iloc[:,0], xx1.iloc[:,1])
"""

test_dg = permutate_single_dg(random.sample(DG59_long, 1)[0], Sequence[key59], expand=10)
test_dg = permutate_single_dg(random.sample(DG766_long,1)[0], Sequence[key766], expand=10)

print test_dg[0]
print sorted(test_dg[1])


###########
# Calculate free energy of common interactions between 59 and 766
###########

def convert766_to_59(dgList_766):
    cov_dgList_766 = []
    for dg in dgList_766:
        s_1, e_1, s_2, e_2 = dg
        s_1 = s_1+12 if s_1>1430 else s_1
        s_2 = s_2+12 if s_2>1430 else s_2
        e_1 = e_1+12 if e_1>1430 else e_1
        e_2 = e_2+12 if e_2>1430 else e_2
        cov_dgList_766.append((s_1,e_1,s_2,e_2))
    return cov_dgList_766

def convert59_to_766(dgList_59):
    cov_dgList_59 = []
    for dg in dgList_59:
        s_1, e_1, s_2, e_2 = dg
        s_1 = s_1-12 if s_1>1430 else s_1
        s_2 = s_2-12 if s_2>1430 else s_2
        e_1 = e_1-12 if e_1>1430 else e_1
        e_2 = e_2-12 if e_2>1430 else e_2
        assert( s_1 < e_1 < s_2 < e_2 )
        cov_dgList_59.append((s_1,e_1,s_2,e_2))
    return cov_dgList_59

def dg_overlap(dg_1, dg_2, tolerate=0):
    if dg_1[0]<dg_2[1]+tolerate and dg_2[0]<dg_1[1]+tolerate and dg_1[2]<dg_2[3]+tolerate and dg_2[2]<dg_1[3]+tolerate:
        return True
    else:
        return False

def combine_dg(dg_1, dg_2):
    return (min(dg_1[0], dg_2[0]), max(dg_1[1], dg_2[1]), min(dg_1[2], dg_2[2]), max(dg_1[3], dg_2[3]))

def combine_common_interaction(common_59_long, common_766_long):
    combine = common_59_long+convert766_to_59(common_766_long)
    combine.sort(key=lambda x: x[0])
    i = 0
    while i < len(combine)-1:
        if dg_overlap(combine[i], combine[i+1]):
            combine[i] = combine_dg(combine[i], combine[i+1])
            del combine[i+1]
        else:
            i += 1
    return combine


def write_mfe(mfe_1, mfe_2, outFile):
    OUT = open(outFile, 'w')
    for d_1, d_2 in zip(mfe_1, mfe_2):
        print >>OUT, "%s\t%s" % (d_1, d_2)
    OUT.close()


dg766_homo_59 = convert766_to_59(uniq_766_long)
dg59_homo_766 = convert59_to_766(uniq_59_long)


##### 766

dg766_mfe = calc_fold_energy(uniq_766_long, Sequence['AY632535.2'], expand=20)
dg766_homo59_mfe = calc_fold_energy(dg766_homo_59, Sequence['KU501215.1'], expand=20)
write_mfe(dg766_mfe, dg766_homo59_mfe, "/tmp/766.txt")


##### 59

dg59_mfe = calc_fold_energy(uniq_59_long, Sequence['KU501215.1'], expand=20)
dg59_homo766_mfe = calc_fold_energy(dg59_homo_766, Sequence['AY632535.2'], expand=20)
write_mfe(dg59_mfe, dg59_homo766_mfe, "/tmp/59.txt")


##### common

common_59Based = combine_common_interaction(common_59_long, common_766_long)
common_766Based = convert59_to_766(common_59Based)

mfe_1 = calc_fold_energy(common_59Based, Sequence['KU501215.1'], expand=20)
mfe_2 = calc_fold_energy(common_766Based, Sequence['AY632535.2'], expand=20)
write_mfe(mfe_1, mfe_2, "/tmp/common.txt")

print "P-Value: ", scipy.stats.ttest_rel(dg766_mfe, dg766_homo59_mfe)[1]
print "P-Value: ", scipy.stats.ttest_rel(dg59_mfe, dg59_homo766_mfe)[1]
print "P-Value: ", scipy.stats.ttest_rel(mfe_1, mfe_2)[1]


