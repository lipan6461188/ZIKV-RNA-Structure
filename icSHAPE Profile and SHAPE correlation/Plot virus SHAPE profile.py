import matplotlib.pyplot as plt
import numpy

def isfloat(value):
  try:
    float(value)
    return True
  except ValueError:
    return False

class SHAPE:
    def __init__(self, trans_id, rpkm, shape_array):
        self.trans_id = trans_id
        self.rpkm = float(rpkm)
        self.length = len(shape_array)
        self.shape_array = shape_array
    @staticmethod
    def replace(raw_shape_array, raw_item, new_item):
        replaced_shape_array = []
        for item in raw_shape_array:
            if item == raw_item:
                replaced_shape_array.append(new_item)
            else:
                replaced_shape_array.append(item)
        return replaced_shape_array
    @staticmethod
    def divide_median(raw_shape_array):
        import numpy
        normed_shape_array = []
        ## 1. find median
        median = numpy.median([ float(item) for item in raw_shape_array if isfloat(item) ])
        ## 2. each number divide this median
        for shape_score in raw_shape_array:
            if isfloat(shape_score):
                normed_shape_array.append( float(shape_score)/median )
            else:
                normed_shape_array.append( 'NULL' )
        return normed_shape_array
    @staticmethod
    def smooth(raw_shape_array, window_size=50, step=25):
        smoothed_shape_array = []
        start = 0
        while start + window_size < len(raw_shape_array):
            smoothed_reactivity = SHAPE.average_str_array(raw_shape_array[start:start+window_size], min_num=10)
            if smoothed_reactivity:
                smoothed_shape_array.append( smoothed_reactivity )
            else:
                smoothed_shape_array.append( 'NULL' )
            start += step
        # the end
        smoothed_reactivity = SHAPE.average_str_array(raw_shape_array[start:start+window_size], min_num=10)
        if smoothed_reactivity:
            smoothed_shape_array.append( smoothed_reactivity )
        else:
            smoothed_shape_array.append( 'NULL' )
        smoothed_shape_array += ["NULL"] * (len(raw_shape_array) - len(smoothed_shape_array))
        return smoothed_shape_array
    @staticmethod
    def average_str_array(str_array, min_num=10):
        Sum = 0
        Num = 0
        for item in str_array:
            if isfloat(item):
                Sum += float(item)
                Num += 1
        if Num >= min_num:
            return Sum / Num
        else:
            return None
    @staticmethod
    def diff_shape_list(shape_array_1, shape_array_2, ABS=False):
        assert( len(shape_array_1) == len(shape_array_2) )
        
        diff_array = []
        for i in range(len(shape_array_1)):
            if shape_array_1[i] != 'NULL' and shape_array_2[i] != 'NULL':
                if ABS:
                    diff_array.append( abs(float(shape_array_1[i]) - float(shape_array_2[i])) )
                else:
                    diff_array.append( float(shape_array_1[i]) - float(shape_array_2[i]) )
            else:
                diff_array.append( 'NULL' )
        
        return diff_array


class icSHAPE:
    def __init__(self, fileName):
        self.shape_list = {}
        IN = open(fileName)
        line = IN.readline()
        while line:
            if not line.startswith('#'):
                data = line.strip().split()
                self.shape_list[ data[0] ] = SHAPE(data[0], float(data[2]), data[3:])
            line = IN.readline()
        IN.close()
    def keys(self):
        return sorted(self.shape_list.keys())
    def value(self, trans_id):
        return self.shape_list[trans_id]

def barplot_shape_score(raw_shape_array):
    import numpy
    float_shape_array = []
    for item in raw_shape_array:
        if isfloat(item):
            float_shape_array.append(float(item) - 1)
        else:
            float_shape_array.append(0.0)
    float_shape_array = numpy.array(float_shape_array)
    #float_shape_array = float_shape_array - 1
    float_shape_array[0] = 1
    float_shape_array[-1] = 1
    plt.bar(range(len(float_shape_array)), float_shape_array, width=1, linewidth=0, color='black')

def lineplot_shape_score(raw_shape_array, start, end):
    import numpy
    float_shape_array = []
    for item in raw_shape_array:
        if isfloat(item):
            float_shape_array.append(item)
        else:
            float_shape_array.append(0.0)
    float_shape_array = numpy.array(float_shape_array)
    plt.plot( range(start, end), float_shape_array, '-' )


##############################
##### Read SHAPE
##############################

shape_set = icSHAPE("/Users/lee/pCloud Drive/Virus_Structure/icSHAPE/final_virus_icSHAPE.out")
key59, key766 = ['KU501215.1', 'AY632535.2']

virus_766 = shape_set.value(key766).shape_array
virus_59 = shape_set.value(key59).shape_array

virus_766_smoothed = SHAPE.smooth(virus_766, window_size=30, step=1)
virus_766_divide = SHAPE.divide_median(virus_766_smoothed)
virus_766_replace = SHAPE.replace(virus_766_divide, 'NULL', 0)

virus_59_smoothed = SHAPE.smooth(virus_59, window_size=30, step=1)
virus_59_divide = SHAPE.divide_median(virus_59_smoothed)
virus_59_replace = SHAPE.replace(virus_59_divide, 'NULL', 0)

##############################
##### Plot Profile
##############################

plt.figure(figsize=(15,6))
plt.subplot(2,1,1)
barplot_shape_score(virus_59_divide)
plt.ylim(-1,1)

plt.subplot(2,1,2)
barplot_shape_score(virus_766_divide)
plt.ylim(-1,1)

plt.show()


##############################
##### Scan structral conserved regions
##### AlignShape is from alignment.py
##############################

def scan_stable_region(shape_diff, window_size=30, max_diff_cutoff=0.1):
    def window_stable(shape_list):
        if shape_list.count('NULL') >= len(shape_list)/2:
            return False
        pure_list = [ float(it) for it in shape_list if it != 'NULL']
        if sum(pure_list)/len(pure_list) > max_diff_cutoff:
            return False
        return True
    my_list = [0]*len(shape_diff)
    start = 0
    while start+window_size < len(shape_diff):
        if window_stable(shape_diff[start:start+window_size]):
            for idx in range(start, start+window_size):
                my_list[idx] = 1
        start += 1
    return my_list

diff_shape = SHAPE.diff_shape_list(AlignShape[key766], AlignShape[key59], ABS=True)
my_list = scan_stable_region(diff_shape, window_size=30, max_diff_cutoff=0.15)
my_list[0] = my_list[-1] = 0.5
plt.bar(range(len(my_list)), my_list, width=1, linewidth=0, color='black')
plt.show()

