
#include <iostream>
#include <paris_sam.h>
#include <align.h>
using namespace std;

using uINT = unsigned;
using uLONG = unsigned long;
using uLONGLONG = unsigned long long;

using uLONG_matrix = vector< vector<uLONG> >;
using double_matrix = vector< vector<double> >;

using dh_array = vector<Duplex_hang>;
using uLONG_pair_pair = pair<pair<uLONG, uLONG>, pair<uLONG, uLONG>>;

const vector<string> color_base = {"#EFE892", "#94D6ED", "#9A94EA", "#48A8EF", "#11EFAF", "#CAE895", "#1368E5", "#F4E3AB" ,"#AB94EA", "#EAA094", "#4EF9F4", "#5ABA93","#C5CC2B","#11efaf"};


/***************************
 
 Code to find the interaction region from matrix data
 
 ***************************/

struct InterRegion {
    InterRegion(uLONG_pair_pair region, double max_cov): region(region),max_cov(max_cov){}
    
    uLONG_pair_pair region;
    uLONG max_cov = 0;
};

void scan_interaction(uLONG_matrix read_count_matrix,
                      vector<InterRegion> &interact_regions,
                      uINT min_dist, uINT min_window_size,
                      uINT max_window_size, uLONG percep_threshold,
                      uLONG extend_threshold);

void filter_interaction(uLONG_matrix read_count_matrix,
                        vector<InterRegion> &interact_regions,
                        const vector<InterRegion> &preserve_interact_regions,
                        const uLONG max_threshold);

template<typename T>
T matrix_block_max(const std::vector<std::vector<T>> &raw_matrix,
                   pair< typename std::vector<std::vector<T>>::size_type, typename std::vector<T>::size_type > &max_pos,
                   const paris_sam_type::uLONG x_lower,
                   const paris_sam_type::uLONG x_upper,
                   const paris_sam_type::uLONG y_lower,
                   const paris_sam_type::uLONG y_upper)
{
    using size_type = typename std::vector<std::vector<T>>::size_type;
    
    T block_max = 0.0;
    max_pos.first = max_pos.second = 0;
    if(x_lower > x_upper or y_lower > y_upper)
    {
        std::cerr << "x_lower: " << x_lower << " x_upper: " << x_upper << std::endl;
        std::cerr << "y_lower: " << y_lower << " y_upper: " << y_upper << std::endl;
        std::cerr << "Fatal Error: lower > upper " << std::endl;
        throw runtime_error("Bad Params");
    }else if(x_upper > raw_matrix.size() or y_upper > raw_matrix.size())
    {
        std::cerr << "Fatal Error: upper > size " << std::endl;
        throw runtime_error("Bad Params");
    }else{
        for(size_type idx = x_lower; idx < x_upper; idx++)
            for(size_type idy = y_lower; idy < y_upper; idy++)
                if(block_max < raw_matrix[idx][idy])
                {
                    max_pos.first = idx;
                    max_pos.second = idy;
                    block_max = raw_matrix[idx][idy];
                }
        
    }
    //cout << block_max << endl;
    return block_max;
}

template <typename T1, typename T2>
T1 min(T1 t1, T2 t2)
{
    return t1 > t2 ? t2 : t1;
}

/***************************
 
 Code of PARIS data analysis
 
 ***************************/

template <typename T>
void fill_matrix_with_aligned_dh(vector<vector<T>> &matrix,
                                 const vector<Duplex_hang> &dh_array,
                                 const Multi_align &ma,
                                 const string &chromosome_name)
{
    //clog << "Start to fill matrix..." << endl;
    using namespace std::placeholders;
    using size_type = typename vector<vector<T>>::size_type;
    
    if(not ma.has(chromosome_name))
    {
        cerr << "FATAL Error: " << chromosome_name << " not in multialign " << endl;
        throw runtime_error("Bad Params");
    }
    
    matrix.clear();
    init_matrix(matrix, ma.length());
    auto coor_covert = std::bind(&Multi_align::raw_coor_to_align_coor, _1, chromosome_name, _2);
    
    for(const Duplex_hang &dh: dh_array)
    {
        size_type start_1 = coor_covert(ma, dh.start_1-1);
        size_type end_1 = coor_covert(ma, dh.end_1-1);
        for(size_type x=start_1; x<end_1; x++)
        {
            size_type start_2 = coor_covert(ma, dh.start_2-1);
            size_type end_2 = coor_covert(ma, dh.end_2-1);
            for(size_type y=start_2; y<end_2; y++)
            {
                matrix[x][y]++;
                matrix[y][x]++;
            }
        }
    }
}

ostream & operator<<(ostream& OUT, const vector<InterRegion> inter_reg);

/***************************
 
 Code of Sequence data
 
 ***************************/

struct Codon_Table
{
    unordered_map<string, string> codon;
    
    void read_codon_table(const string &file_name)
    {
        ifstream IN(file_name, ifstream::in);
        if(not IN)
        {
            cerr << "FATAL Error: " << file_name << " can not readable" << endl;
            exit(-1);
        }
        string this_line;
        while(getline(IN, this_line))
        {
            vector<string> aa_codons = split(this_line, '=');
            if(aa_codons.size() != 2)
            {
                cerr << "A unexpected line: \n\t" << this_line << endl;
                continue;
            }
            vector<string> codons_of_this_aa;
            split(aa_codons.at(1), '|', codons_of_this_aa);
            for(string cur_codon: codons_of_this_aa)
            {
                std::replace (cur_codon.begin(), cur_codon.end(), 'U', 'T');
                codon[cur_codon] = aa_codons.at(0);
            }
        }
    }
    
    Codon_Table(const string &file_name){ read_codon_table(file_name); }
    
    void show()
    {
        for(auto iter=codon.cbegin(); iter!=codon.cend(); iter++)
            cout << iter->first << " ==> " << iter->second << endl;
    }
};

struct annotation
{
    string chromosome_name;
    vector< tuple<uLONG, uLONG, string> > exons;
    
    friend ostream & operator <<(ostream &, const annotation&);
};

struct Circot_Mutation
{
    double total_len;
    
    vector<double> mutations;
};

ostream & operator <<(ostream &OUT, const annotation& annot);
void read_annotation(annotation& annot, const string& anno_file_name);
void annotation_to_align_annotation(const Multi_align &ma,
                                    const annotation& raw_annot,
                                    annotation& align_annot);
Circot_Mutation get_aa_mutation(const Codon_Table &codon_table,
                                const string &seq_1,
                                const string   &seq_2,
                                long start,
                                long end);

/***************************
 
 Code of Domain
 
 ***************************/

struct domain
{
    uLONG start;
    uLONG end;
    
    double domain_index = -1;
    // double in_domain_reads;
    // double out_domain_reads;
    
    uLONG domain_size()const{ return end-start+1; }
    
    domain(uLONG start, uLONG end, double domain_index):start(start), end(end), domain_index(domain_index) {}
    domain()=default;
    
    //bool operator<(const domain &other){ return domain_index < other.domain_index; }
    //friend ostream& operator<<(ostream &, const domain &);
};

void read_domain(const string &file_name,
                 vector<domain> &domain_list,
                 uLONG multiplex_index=1);
ostream& operator<<(ostream &OUT, const vector<domain> &d);
void domain_coor_convert(vector<domain> &domain_list,
                         const Multi_align &ma,
                         const string &chromosome_name);











































