
def readCluster(inFile):
    Cluster = []
    for line in open(inFile):
        data = line.strip().split()
        ls,le,rs,re = int(data[4]), int(data[5]), int(data[8]), int(data[9])
        support = int(data[10])
        Cluster.append( (ls, le, rs, re, support) )
    return Cluster

def readDG(inFile):
    DG = []
    for line in open(inFile):
        data = line.strip().split()
        ls,le,rs,re = int(data[5]), int(data[6]), int(data[12]), int(data[13])
        DG.append( (ls, le, rs, re) )
    return DG


def plot_read_span(reads):
    gap_len = [ it[2]-it[1] for it in reads ]
    gap_len.sort()
    gap_len = np.array(gap_len)
    gap_len = np.log10(gap_len)
    X = []
    Y = []
    for x in np.arange(0, 20, 0.01):
    #for x in range(0, 10900, 50):
        x = round(x, 2)
        y = len([ it for it in gap_len if it < x ])
        ratio = round(1.0*y/len(gap_len), 3)
        X.append(x)
        Y.append(ratio)
        if ratio == 1.0: break
    return X, Y

def classify_Cluster(cluster_vivo, cluster_vitro, min_diff=5, min_ratio=2, min_reads=2):
    higher_vivo = []
    higher_vitro = []
    common = []
    uniq_vivo = []
    uniq_vitro = []
    
    cluster_vivo = [ it for it in cluster_vivo if it[4]>=min_reads ]
    cluster_vitro = [ it for it in cluster_vitro if it[4]>=min_reads ]
    
    for dg_vivo in cluster_vivo:
        find = False
        for dg_vitro in cluster_vitro:
            if dg_vivo[0]<dg_vitro[1] and dg_vitro[0]<dg_vivo[1] and dg_vivo[2]<dg_vitro[3] and dg_vitro[2]<dg_vivo[3]:
                find = True
                if dg_vitro[4] - dg_vivo[4] >= min_diff and 1.0*dg_vitro[4]/dg_vivo[4] >= min_ratio:
                    higher_vitro.append( (dg_vivo, dg_vitro) )
                elif dg_vivo[4] - dg_vitro[4] >= min_diff and 1.0*dg_vivo[4]/dg_vitro[4] >= min_ratio:
                    higher_vivo.append( (dg_vivo, dg_vitro) )
                else:
                    common.append( (dg_vivo, dg_vitro) )
                break
        if not find:
            uniq_vivo.append(dg_vivo)
    
    for dg_vitro in cluster_vitro:
        find = False
        for dg_vivo in cluster_vivo:
            if dg_vivo[0]<dg_vitro[1] and dg_vitro[0]<dg_vivo[1] and dg_vivo[2]<dg_vitro[3] and dg_vitro[2]<dg_vivo[3]:
                find = True
        if not find:
            uniq_vitro.append(dg_vitro)
    
    return higher_vivo, higher_vitro, common, uniq_vivo, uniq_vitro

###### 1. Read reads

dg_59_vitro = readDG("/private/tmp/PARIS/PRVABC59_48_vitro.dg")
dg_59_vivo = readDG("/private/tmp/PARIS/PRVABC59_72_vivo.dg")

dg_766_vitro = readDG("/private/tmp/PARIS/MR766_72_vitro.dg")
dg_766_vivo = readDG("/private/tmp/PARIS/MR766_72_vivo.dg")

cluster_59_vitro = readCluster("/private/tmp/PARIS/PRVABC59_48_vitro.tab")
cluster_59_vivo = readCluster("/private/tmp/PARIS/PRVABC59_72_vivo.tab")

###### 2. Plot CDF

X_59_vitro, Y_59_vitro = plot_read_span(dg_59_vitro)
X_59_vivo, Y_59_vivo = plot_read_span(dg_59_vivo)
X_766_vitro, Y_766_vitro = plot_read_span(dg_766_vitro)
X_766_vivo, Y_766_vivo = plot_read_span(dg_766_vivo)

plt.figure(figsize=(10,5))
plt.subplot(1,2,1)
plt.plot(X_59_vitro, Y_59_vitro, '-', c='blue')
plt.plot(X_59_vivo, Y_59_vivo, '-', c='red')
plt.subplot(1,2,2)
plt.plot(X_766_vitro, Y_766_vitro, '-', c='blue')
plt.plot(X_766_vivo, Y_766_vivo, '-', c='red')
plt.savefig("/tmp/figure.pdf")
plt.show()

scipy.stats.ks_2samp( [ it[2]-it[1] for it in dg_59_vitro ], [ it[2]-it[1] for it in dg_59_vivo ] )
scipy.stats.ks_2samp( [ it[2]-it[1] for it in dg_766_vitro ], [ it[2]-it[1] for it in dg_766_vivo ] )

###### 3. Read matrix

def readMatrix(inFile):
    matrix = []
    for line in open(inFile):
        data = line.strip().split()
        data = [ int(it) for it in data ]
        matrix.append(data)
    return matrix

matrix_59_vivo = readMatrix("/private/tmp/PARIS/PRVABC59_72_vivo.matrix")
matrix_59_vitro = readMatrix("/private/tmp/PARIS/PRVABC59_48_vitro.matrix")

matrix_766_vivo = readMatrix("/private/tmp/PARIS/MR766_72_vivo.matrix")
matrix_766_vitro = readMatrix("/private/tmp/PARIS/MR766_72_vitro.matrix")


def shrink_matrix(raw_matrix, shrink_dim=400):
    s_matrix = []
    for i in range(shrink_dim):
        s_matrix.append( [0]*shrink_dim )
    
    raw_dim = len(raw_matrix)
    step = 1.0*raw_dim/shrink_dim
    
    for i in range(shrink_dim):
        print("index: "+str(i))
        row_start = int(step*i)
        row_end = min( int(step*(i+1)), raw_dim )
        for j in range(shrink_dim):
            col_start = int(step*j)
            col_end = min( int(step*(j+1)), raw_dim )
            count = 0
            num = 0
            #max_v = 0
            for ind_i in range(row_start, row_end):
                V_array = raw_matrix[ind_i]
                for ind_j in range(col_start, col_end):
                    cur_value = V_array[ind_j]
                    count += cur_value
                    num += 1
                    #if cur_value > max_v:
                    #    max_v = cur_value
            if num > 0:
                s_matrix[i][j] = 1.0*count/num
    
    return s_matrix

def diff_matrix(matrix_1, matrix_2):
    dim = len(matrix_1)
    d_matrix = []
    for i in range(dim):
        d_matrix.append( [0]*dim )
    
    for i in range(dim):
        for j in range(dim):
            d_matrix[i][j] = matrix_1[i][j] - matrix_2[i][j]
    
    return d_matrix

matrix_59_vivo_shrink = shrink_matrix(matrix_59_vivo, shrink_dim=300)
matrix_59_vitro_shrink = shrink_matrix(matrix_59_vitro, shrink_dim=300)
diff_59 = diff_matrix(matrix_59_vitro_shrink, matrix_59_vivo_shrink)

Red = [(0xFF/0xFF, 0x00/0xFF, 0x00/0xFF, it) for it in np.linspace(0,1,20)]
Red_cmap = sns.mpl.colors.ListedColormap(Red)

Blue = [(0x00/0xFF, 0x00/0xFF, 0xFF/0xFF, it) for it in np.linspace(0,1,20)]
Blue_cmap = sns.mpl.colors.ListedColormap(Blue)


plt.figure(figsize=(20,5))
plt.subplot(1,3,1)
sns.heatmap(matrix_59_vivo_shrink, vmin=0, vmax=10, cmap=Red_cmap) #sns.color_palette("Blues", 20))
plt.subplot(1,3,2)
sns.heatmap(matrix_59_vitro_shrink, vmin=0, vmax=10, cmap=Blue_cmap) #sns.color_palette("Blues", 20))
plt.subplot(1,3,3)
sns.heatmap(diff_59, vmin=-5, vmax=5, center=0, cmap=sns.color_palette("RdBu_r", 40))
plt.savefig("/tmp/test.pdf")
plt.show()


matrix_766_vivo_shrink = shrink_matrix(matrix_766_vivo, shrink_dim=300)
matrix_766_vitro_shrink = shrink_matrix(matrix_766_vitro, shrink_dim=300)
diff_766 = diff_matrix(matrix_766_vitro_shrink, matrix_766_vivo_shrink)

plt.figure(figsize=(20,5))
plt.subplot(1,3,1)
sns.heatmap(matrix_766_vivo_shrink, vmin=0, vmax=10, cmap=sns.color_palette("Blues", 20))
plt.subplot(1,3,2)
sns.heatmap(matrix_766_vitro_shrink, vmin=0, vmax=10, cmap=sns.color_palette("Blues", 20))
plt.subplot(1,3,3)
sns.heatmap(diff_766, vmin=-5, vmax=5, center=0, cmap=sns.color_palette("RdBu_r", 40))
plt.show()



def compare_heatmap(matrix_1, matrix_2, min_v=0, max_v=10):
    import matplotlib.pyplot as plt
    import matplotlib.patches as patches
    
    assert len(matrix_1) == len(matrix_2)
    dim = len(matrix_1)
    
    fig,ax = plt.subplots(1)
    for r in range(dim):
        print r
        Row1 = matrix_1[r]
        Row2 = matrix_2[r]
        for c in range(dim):
            v1 = min(Row1[c],max_v)
            v2 = min(Row2[c],max_v)
            r1 = 1.0*(v1-min_v)/(max_v-min_v)
            r2 = 1.0*(v2-min_v)/(max_v-min_v)
            color = (r1, r2, np.sqrt((1-r1)*(1-r2)))
            #print color
            if r1>0.1 or r2>0.1:
                rect = patches.Rectangle((c,dim-r),1,-1, linewidth=0, edgecolor=None, facecolor=color)
                ax.add_patch(rect)
        #if r == 50:
        #    break
    ax.set_xlim(0, dim)
    ax.set_ylim(0, dim)

def compare_heatmap_descrete(matrix_1, matrix_2, cutoff=1):
    import matplotlib.pyplot as plt
    import matplotlib.patches as patches
    
    assert len(matrix_1) == len(matrix_2)
    dim = len(matrix_1)
    
    fig,ax = plt.subplots(1)
    for r in range(dim):
        print r
        Row1 = matrix_1[r]
        Row2 = matrix_2[r]
        for c in range(dim):
            v1 = Row1[c]
            v2 = Row2[c]
            if v1>=cutoff and v2<cutoff:
                color = 'red'
            elif v1<cutoff and v2>=cutoff:
                color = 'green'
            elif v1>=cutoff and v2>=cutoff:
                color = 'blue'
            if v1>=cutoff or v2>=cutoff:
                rect = patches.Rectangle((c,dim-r),1,-1, linewidth=0, edgecolor=None, facecolor=color)
                ax.add_patch(rect)
        #if r == 50:
        #    break
    ax.set_xlim(0, dim)
    ax.set_ylim(0, dim)

compare_heatmap_descrete(matrix_59_vivo_shrink, matrix_59_vitro_shrink, cutoff=1.0)
plt.savefig("/tmp/test3.pdf")
plt.close()

















