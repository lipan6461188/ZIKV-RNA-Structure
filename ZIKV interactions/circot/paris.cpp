
#include <iostream>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <string>
#include <cstring>


#include <ctime>
#include <random>
#include <tuple>
#include <cmath>
#include <memory>
#include <functional>
#include "paris.h"


using namespace std;


/***************************
  
 Code to find the interaction region from matrix data
 
 ***************************/

void scan_interaction(uLONG_matrix read_count_matrix,
                      vector<InterRegion> &interact_regions,
                      uINT min_dist, uINT min_window_size,
                      uINT max_window_size, uLONG percep_threshold,
                      uLONG extend_threshold)
{
    using size_type = uLONG_matrix::size_type;
    
    
    interact_regions.clear();
    uLONG chromosome_len = read_count_matrix.size();
    
    for(size_type x_idx=0; x_idx<chromosome_len; x_idx+=min_window_size)
    {
        //cout << "X Step " << x_idx << endl;
        for(size_type y_idx=x_idx+min_dist; y_idx<chromosome_len;y_idx+=min_window_size)
        {
            //cout << "Y Step " << y_idx << endl;
            
            size_type width = min_window_size;
            size_type height = min_window_size;
            size_type proper_x_upper = min(x_idx+width, chromosome_len);
            size_type proper_y_upper = min(y_idx+height, chromosome_len);
            pair<size_type, size_type> last_max_pos, max_pos;
            
            uLONG current_max_rc = matrix_block_max<uLONG>(read_count_matrix, last_max_pos, x_idx, proper_x_upper, y_idx, proper_y_upper);
            
            if( current_max_rc >= percep_threshold )
            {
                //cout << "Find One: " << x_idx << "\t" << y_idx << endl;
                // posit center
                x_idx = last_max_pos.first >= min_window_size / 2 ? last_max_pos.first - min_window_size / 2 : 0;
                y_idx = last_max_pos.second >= min_window_size / 2 ? last_max_pos.second - min_window_size / 2 : 0;
                proper_x_upper = min(x_idx+width, chromosome_len);
                proper_y_upper = min(y_idx+height, chromosome_len);
                uLONG max_rc_record = matrix_block_max<uLONG>(read_count_matrix, max_pos, x_idx, proper_x_upper, y_idx, proper_y_upper);
                
                while(max_pos != last_max_pos)
                {
                    // cout << x_idx << "\t" << y_idx << endl;
                    x_idx = last_max_pos.first >= min_window_size / 2 ? last_max_pos.first - min_window_size / 2 : 0;
                    y_idx = last_max_pos.second >= min_window_size / 2 ? last_max_pos.second - min_window_size / 2 : 0;
                    proper_x_upper = min(x_idx+width, chromosome_len);
                    proper_y_upper = min(y_idx+height, chromosome_len);
                    last_max_pos = max_pos;
                    max_rc_record = matrix_block_max<uLONG>(read_count_matrix, max_pos, x_idx, proper_x_upper, y_idx, proper_y_upper);
                }
                
                // cout << "Start to extend..." << endl;
                bool extend_up(true), extend_left(true), extend_right(true), extend_down(true);
                vector<uLONG> tmp_vec; uLONG max_rc;
                while( extend_up or extend_left or extend_right or extend_down )
                {
                    // extend left
                    //  if(extend_left)
                    //  {
                    tmp_vec.clear();
                    for(size_type idx=0;idx<min(height, chromosome_len-y_idx);idx++)
                        tmp_vec.push_back(read_count_matrix.at(x_idx).at(y_idx+idx));
                    max_rc = *max_element(tmp_vec.cbegin(), tmp_vec.cend());
                    if(max_rc>=extend_threshold and y_idx-x_idx>min_dist and x_idx>0)
                    {
                        x_idx--;
                        width++;
                    }else{
                        extend_left = false;
                    }
                    //  }
                    if(width >= max_window_size) break;
                    
                    // extend down
                    //  if(extend_down)
                    //  {
                    tmp_vec.clear();
                    for(size_type idx=0;idx<min(width, chromosome_len-x_idx);idx++)
                        tmp_vec.push_back(read_count_matrix.at(x_idx+idx).at(y_idx));
                    max_rc = *max_element(tmp_vec.cbegin(), tmp_vec.cend());
                    if(max_rc>=extend_threshold and y_idx-x_idx>min_dist and y_idx>0)
                    {
                        y_idx--;
                        height++;
                    }else{
                        extend_down = false;
                    }
                    //  }
                    if(height >= max_window_size) break;
                    
                    // extend right
                    //  if(extend_right)
                    //  {
                    tmp_vec.clear();
                    for(size_type idx=0;idx<min(height, chromosome_len-y_idx);idx++)
                        tmp_vec.push_back(read_count_matrix.at( min(x_idx+width-1,chromosome_len-1) ).at(y_idx+idx));
                    max_rc = *max_element(tmp_vec.cbegin(), tmp_vec.cend());
                    if(max_rc>=extend_threshold and x_idx+width<chromosome_len)
                    {
                        width++;
                    }else{
                        extend_right = false;
                    }
                    // }
                    if(width >= max_window_size) break;
                    
                    // extend upper
                    //  if(extend_up)
                    // {
                    tmp_vec.clear();
                    for(size_type idx=0;idx<min(width, chromosome_len-x_idx);idx++)
                        tmp_vec.push_back(read_count_matrix.at(x_idx+idx).at( min(y_idx+height-1,chromosome_len-1) ));
                    max_rc = *max_element(tmp_vec.cbegin(), tmp_vec.cend());
                    if(max_rc>=extend_threshold and y_idx+height<chromosome_len)
                    {
                        height++;
                    }else{
                        extend_up = false;
                    }
                    // }
                    if(height >= max_window_size) break;
                }
                //  cout << "extend finish..." << endl;
                pair<uLONG, uLONG> left(x_idx, min(x_idx+width, chromosome_len));
                pair<uLONG, uLONG> right(y_idx, min(y_idx+height, chromosome_len));
                interact_regions.push_back( InterRegion(make_pair(left, right), max_rc_record) );
                
                //  cout << x_idx << "-" << x_idx+width << "\t" << y_idx << "-" << y_idx+height << endl;
                
                //cout << "to zero..." << endl;
                for(auto x_iter=read_count_matrix.begin()+x_idx; x_iter<read_count_matrix.begin()+x_idx+min(width,chromosome_len-x_idx); x_iter++)
                    for(auto y_iter=x_iter->begin()+y_idx; y_iter<x_iter->begin()+y_idx+min(height, chromosome_len-y_idx); y_iter++)
                        *y_iter = 0;
                
                // again
                x_idx = 0;
                // cout << "end find..." << endl;
                break;
            }
        }
    }
}

void filter_interaction(uLONG_matrix read_count_matrix,
                        vector<InterRegion> &interact_regions,
                        const vector<InterRegion> &preserve_interact_regions,
                        const uLONG max_threshold)
{
    using size_type = uLONG_matrix::size_type;
    
    auto in_preserve_region = [&preserve_interact_regions](const InterRegion &inter_region)->bool
    {
        uLONG LS = inter_region.region.first.first;
        uLONG LE = inter_region.region.first.second;
        uLONG RS = inter_region.region.second.first;
        uLONG RE = inter_region.region.second.second;
        
        for(auto iter=preserve_interact_regions.cbegin(); iter!=preserve_interact_regions.cend(); iter++)
        {
            if( LS < iter->region.first.second and iter->region.first.first < LE and
               RS < iter->region.second.second and iter->region.second.first < RE )
                return true;
        }
        
        return false;
    };
    
    auto iter=interact_regions.begin();
    while(iter!=interact_regions.end())
    {
        if( in_preserve_region(*iter) )
        {
            iter++;
            continue;
        }
        
        uLONG LS = iter->region.first.first;
        uLONG LE = iter->region.first.second;
        uLONG RS = iter->region.second.first;
        uLONG RE = iter->region.second.second;
        

        pair<size_type, size_type> max_pos;
        uLONG max_rc = matrix_block_max<uLONG>(read_count_matrix, max_pos, LS, LE, RS, RE);
        
        if( max_rc >= max_threshold )
        {
            cout << "Remove(" << max_rc << ">=" << max_threshold <<
            "  " << LS << "-" << LE << " <==> " << RS << "-" << RE << endl;
            
            iter = interact_regions.erase(iter);
            continue;
        }
        
        iter++;
        
        /*
        
        bool remove_current = false;
        for(size_type idx=LS; idx<LE; idx++)
            for(size_type idy=RS; idy<RE; idy++)
                if( read_count_matrix.at(idx).at(idy) >= max_threshold)
                {
                    cout << "Remove(" << read_count_matrix.at(idx).at(idy) << ">=" << max_threshold <<
                    " in [" << idx << ", " << idy << "]) " << LS << "-" << LE << " <==> " << RS << "-" << RE << endl;
                    
                    iter = interact_regions.erase(iter);
                    idx = LE;
                    remove_current = true;
                    break;
                }
        
        if( not remove_current )
            iter++;
         */
    }
}

/***************************
 
 Code of PARIS data analysis
 
 ***************************/


ostream & operator<<(ostream& OUT, const vector<InterRegion> inter_reg)
{
    for(const InterRegion &reg: inter_reg)
    {
        OUT << reg.region.first.first << "-" << reg.region.first.second << " <===> " << reg.region.second.first << "-" << reg.region.second.second << "\t" << reg.max_cov << endl;
    }
    OUT << "\n";
    return  OUT;
}


/***************************
 
 Code of Sequence data
 
 ***************************/


ostream & operator <<(ostream &OUT, const annotation& annot)
{
    OUT << "#=========" << annot.chromosome_name << "=========\n";
    for(const auto &exon: annot.exons)
        OUT << get<0>(exon) << "-" << get<1>(exon) << "\t" << get<2>(exon) << endl;
    return OUT;
}

void read_annotation(annotation& annot, const string& anno_file_name)
{
    ifstream IN(anno_file_name, ifstream::in);
    if(not IN)
    {
        cerr << "FATAL Error: Bad file --- " << anno_file_name << " cannot be open" << endl;
        exit(EXIT_FAILURE);
    }
    annot.exons.clear();
    string this_line;
    while(getline(IN, this_line))
    {
        if(this_line.at(0) == '#')
            continue;
        istringstream ifStream(this_line);
        uLONG start, end;
        string chrom_name, func_anno;
        if(ifStream >> chrom_name >> start >> end >> func_anno)
        {
            if( (not annot.chromosome_name.empty()) and (annot.chromosome_name != chrom_name) )
            {
                cout << "FATAL Error: Different chromosome name: " << annot.chromosome_name << " and " << chrom_name << endl;
                exit(EXIT_FAILURE);
            }else
                annot.chromosome_name = chrom_name;
            annot.exons.push_back(make_tuple(start, end, func_anno));
        }
    }
    IN.close();
    sort(annot.exons.begin(), annot.exons.end(), [](const tuple<uLONG,uLONG,string> &exon_1, const tuple<uLONG,uLONG,string> &exon_2)->bool{ return get<0>(exon_1) < get<0>(exon_2)?true:false; });
}


void annotation_to_align_annotation(const Multi_align &ma, const annotation& raw_annot, annotation& align_annot)
{
    align_annot.exons.clear();
    align_annot.chromosome_name = raw_annot.chromosome_name;
    if(ma.has(raw_annot.chromosome_name))
    {
        string align_seq = ma.get_align(raw_annot.chromosome_name).align_seq;
        
        uLONG insert_count = 0;
        uLONG alpha_count = 0;
        uLONG next_start = get<0>(raw_annot.exons.at(0));
        uLONG next_end = get<1>(raw_annot.exons.at(0));
        uLONG align_start = 0;
        uLONG align_end = 0;
        
        uINT exon_idx = 0;
        for(size_t align_idx=0;align_idx<align_seq.size();align_idx++)
        {
            if(align_seq.at(align_idx) == '-')
            {
                insert_count++;
            }else{
                alpha_count++;
            }
            
            if( alpha_count == next_start )
            {
                align_start = alpha_count + insert_count;
            }
            if( alpha_count == next_end )
            {
                align_end = alpha_count + insert_count;
                if(exon_idx==0) align_start = 1;
                if(exon_idx==raw_annot.exons.size()-1) align_end = align_seq.size();
                align_annot.exons.push_back(make_tuple(align_start, align_end, get<2>(raw_annot.exons.at(exon_idx))));
                ++exon_idx;
                if(exon_idx==raw_annot.exons.size())
                    break;
                next_start = get<0>(raw_annot.exons.at(exon_idx));
                next_end = get<1>(raw_annot.exons.at(exon_idx));
            }
        }
        // get<0>(align_annot.exons.at(0)) = 1;
    }else{
        cerr << "FATAL ERROR: Annotation chromosome " << raw_annot.chromosome_name << " didn't appear in align file" << endl;
        exit(EXIT_FAILURE);
    }
}

Circot_Mutation get_aa_mutation(const Codon_Table &codon_table, const string &seq_1, const string   &seq_2, long start, long end)
{
    if( seq_1.size() != seq_2.size() )
    {
        cerr << "FATAL Error: Two Sequence have different length -- " << seq_1.size() << "\t" << seq_2.size() << endl;
        exit(-1);
    }
    
    Circot_Mutation cm;
    cm.total_len = seq_1.size();
    cm.mutations.clear();
    
    vector<string> protein_1, protein_2;
    
    string cur_codon_1, cur_codon_2;
    int gap_1=0, gap_2=0;
    for(string::size_type idx=start; idx<(unsigned)end; idx+=1)
    {
        if(seq_1.at(idx) != '-')
            cur_codon_1.push_back(seq_1.at(idx));
        else
            gap_1 ++;
        
        if(seq_2.at(idx) != '-')
            cur_codon_2.push_back(seq_2.at(idx));
        else
            gap_2 ++;
        
        if(gap_1 == 3)
        {
            protein_1.push_back( "NULL" );
            gap_1 = 0;
        }
        
        if(gap_2 == 3)
        {
            protein_2.push_back( "NULL" );
            gap_2 = 0;
        }
        
        if(cur_codon_1.size() == 3)
        {
            protein_1.push_back( codon_table.codon.at(cur_codon_1) );
            cur_codon_1.clear();
        }
        if(cur_codon_2.size() == 3)
        {
            protein_2.push_back( codon_table.codon.at(cur_codon_2) );
            cur_codon_2.clear();
        }
    }
    
    if(protein_1.size() != protein_2.size())
    {
        cerr << "FATAL Error: Different protein length..." << endl;
        exit(-1);
    }
    
    for(string::size_type idx=0; idx<protein_1.size(); idx++)
    {
        if(protein_1.at(idx) != protein_2.at(idx))
            for(string::size_type bias=0; bias<3; bias++)
                cm.mutations.push_back( start + 3*idx + bias );
    }
    
    return cm;
}

/***************************
 
 Code of Domain
 
 ***************************/


void read_domain(const string &file_name, vector<domain> &domain_list, uLONG multiplex_index)
{
    ifstream IN(file_name,ifstream::in);
    if(not IN)
    {
        cerr << "FATAL Error: " << file_name << " cannot be opened" << endl;
        exit(-1);
    }
    
    //vector<domain> domain_list;
    
    using size_type = vector<domain>::size_type;
    
    domain_list.clear();
    
    string this_line;
    while(getline(IN, this_line))
    {
        trim(this_line, ' ');
        trim(this_line, '\t');
        if(this_line.empty())
            continue;
        
        istringstream SS(this_line);
        domain d;
        SS >> d.start >> d.end;
        
        d.start *= multiplex_index;
        d.end *= multiplex_index;
        //d.end += multiplex_index;
        
        domain_list.push_back(d);
    }
    IN.close();
    
    domain_list.front().start = 1;
    for(size_type idx=0; idx<domain_list.size()-1; idx++)
        domain_list.at(idx).end += multiplex_index - 1;
    
    std::sort(domain_list.begin(), domain_list.end(), [](const domain &d_1, const domain &d_2)->bool{ return d_1.domain_index < d_2.domain_index ? true : false; });
}

ostream& operator<<(ostream &OUT, const vector<domain> &d)
{
    //vector<domain>::size_type i = 0;
    
    for(auto iter=d.cbegin(); iter!=d.cend(); iter++)
    {
        OUT << iter->start << "\t" << iter->end << "\n";
    }
    
    return OUT;
}

void domain_coor_convert(vector<domain> &domain_list, const Multi_align &ma, const string &chromosome_name)
{
    using namespace std::placeholders;
    using size_type = vector<domain>::size_type;
    
    if(not ma.has(chromosome_name))
    {
        cerr << "FATAL Error: " << chromosome_name << " not in multialign " << endl;
        exit(-1);
    }
    
    auto coor_covert = std::bind(&Multi_align::align_coor_to_raw_coor, _1, chromosome_name, _2);
    
    uLONG align_len = ma.length();
    uLONG chromosome_len = ma.get_align(chromosome_name).seq_length;
    uLONG bin_size = domain_list.back().end;
    
    double pre_step_size = 1.0 * align_len / bin_size;
    //clog << "Step Size: " << step_size << endl;
    
    for(size_type i=0; i<domain_list.size(); i++)
    {
        //uLONG start = (domain_list.at(i).start-1) * pre_step_size;
        uLONG end = (domain_list.at(i).end-1) * pre_step_size;
        
        //uLONG post_start = coor_covert(ma, start);
        uLONG post_end = coor_covert(ma, end);
        
        if(i == 0)
        {
            uLONG next_start = (domain_list.at(i+1).start-1) * pre_step_size;
            uLONG next_post_start = coor_covert(ma, next_start);
            
            domain_list.at(i).start = 1;
            domain_list.at(i).end = next_post_start - 1;
        }else if(i==domain_list.size()-1){
            domain_list.at(i).start = domain_list.at(i-1).end+1;
            domain_list.at(i).end = chromosome_len;
        }else{
            domain_list.at(i).start = domain_list.at(i-1).end+1;
            domain_list.at(i).end = post_end;
        }
    }
}







