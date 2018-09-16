
#include "circos.h"
#include "paris.h"
#include <paris_sam.h>
#include <align.h>
#include <cmath>
#include <functional>

#include <QApplication>
#include <QPainter>
#include <QPen>
#include <QPainterPath>
#include <QPdfWriter>
#include <QtMath>
#include <exception>
#include <QDebug>


void plot_circos(const vector<InterRegion> &inter_reg_59,
                 const vector<InterRegion> &inter_reg_766,
                 const Circot_Mutation &seq_mutations,
                 const Circot_Mutation &aa_mutations,
                 const annotation &protein_anno,
                 const vector<domain> &domain_list,
                 uLONG total_len);

void test_mutation_significant(const Circot_Mutation &mutation,
                               const vector<InterRegion> &inter_reg,
                               uLONG total_len,
                               uLONG times=1000);

void save_interaction_region_seq(const vector<InterRegion> &region,
                                 const string &ref_trans_id,
                                 const string &homo_trans_id,
                                 const Multi_align &ma,
                                 const string &file_name);

int main(int argc, char *argv[])
{
    QApplication app(argc, argv);
    
    
    /************
      PARIS data
     ************/
    
    Multi_align ma("/Users/lee/Desktop/Projects/Virus/Virus_Genome/final_align.stoch");
    
    string chromosome_name_59("KU501215.1");
    string chromosome_name_766("AY632535.2");
    string sam_file_name_59("/Users/lee/Desktop/Projects/Virus/Virus_Genome/PARIS/sam/59.sam");
    string sam_file_name_766("/Users/lee/Desktop/Projects/Virus/Virus_Genome/PARIS/sam/766.sam");

    uINT min_overhang = 5;
    vector<Duplex_hang> dh_array_59, dh_array_766;
    
    get_chromosome_hang(sam_file_name_59, dh_array_59, chromosome_name_59, min_overhang);
    get_chromosome_hang(sam_file_name_766, dh_array_766, chromosome_name_766, min_overhang);
    
    vector<vector<uLONG>> matrix_59, matrix_766;
    clog << "start to fill_matrix_with_aligned_dh..." << endl;
    fill_matrix_with_aligned_dh(matrix_59, dh_array_59, ma, chromosome_name_59);
    fill_matrix_with_aligned_dh(matrix_766, dh_array_766, ma, chromosome_name_766);
    
    cout << "59 matrix: " << matrix_59.size() << endl;
    cout << "766 matrix: " << matrix_766.size() << endl;
    
    vector<InterRegion> inter_reg_59, inter_reg_766;
    
    uINT min_dist = 1000;
    uINT min_window_size = 25;
    uINT max_window_size = 150;
    uLONG percep_threshold = 25;
    uLONG extend_threshold = 10;
    uLONG filter_threshold = 10;
    
    clog << "start to scan_interaction..." << endl;
    scan_interaction(matrix_59, inter_reg_59,
                     min_dist, min_window_size,
                     max_window_size, percep_threshold,
                     extend_threshold);
    scan_interaction(matrix_766, inter_reg_766,
                     min_dist, min_window_size,
                     max_window_size, percep_threshold,
                     extend_threshold);
    
    clog << "start to filter_interaction..." << endl;
    clog << "\t in Zika_766" << endl;
    filter_interaction(matrix_59,
                       inter_reg_766,
                       inter_reg_59,
                       filter_threshold);
    clog << "\t in Zika_59" << endl;
    filter_interaction(matrix_766,
                       inter_reg_59,
                       inter_reg_766,
                       filter_threshold);
    
    clog << "Interaction Numbers: \n\t59 -- " << inter_reg_59.size() << "\n\t766 -- " << inter_reg_766.size() << endl;
    
    cout << "Zika 59 long range interaction" << endl;
    cout << inter_reg_59 << endl;
    cout << "Zika 766 long range interaction" << endl;
    cout << inter_reg_766 << endl;
    
    /************
     Sequence Data
     ************/
    
    Codon_Table codon_table("/Users/lee/Desktop/Projects/Virus/Virus_Genome/Sequence/Codon_Table.txt");
    Circot_Mutation seq_mutation, aa_mutation;
    
    // 1. annotation
    annotation annot, align_annot;
    read_annotation(annot, "/Users/lee/Desktop/Projects/Virus/Virus_Genome/Figures/Figure_1/59_annotation.txt");
    annotation_to_align_annotation(ma, annot, align_annot);
    //cout << annot << "\n";
    //cout << align_annot << "\n";

    // 2. Statistic mutations( sequence and amino acid )
    string seq_59 = ma.get_align("KU501215.1").align_seq;
    string seq_766 = ma.get_align("AY632535.2").align_seq;
    seq_mutation.total_len = ma.length();
    for(int idx=0; idx<seq_mutation.total_len; idx++)
        if(seq_59[idx] != seq_766[idx])
            seq_mutation.mutations.push_back(idx);
    aa_mutation = get_aa_mutation(codon_table,
                                  ma.get_align("KU501215.1").align_seq,
                                  ma.get_align("AY632535.2").align_seq,
                                  get<0>(align_annot.exons.at(1))-1,
                                  get<1>(align_annot.exons.at(align_annot.exons.size()-2))
                                  );
    
    //3. domain
    const string domain_file("/Users/lee/pCloud Drive/Projects/Virus/re-domain/59_domain.txt");
    vector<domain> domain_list;
    read_domain(domain_file, domain_list, 1);
    //clog << domain_list << endl;
    domain_coor_convert(domain_list, ma, "KU501215.1");
    //clog << domain_list << endl;
    
    plot_circos(inter_reg_59, inter_reg_766, seq_mutation, aa_mutation, align_annot, domain_list, ma.length());

    /*
     // 相互作用区域的突变是否显著比较高：结论并不是
    test_mutation_significant(seq_mutation,
                              inter_reg_59,
                              ma.length(),
                              1000);
    */

    clog << "Zika 766 search interaction seq..." << endl;
    save_interaction_region_seq(inter_reg_766,
                                chromosome_name_766,
                                chromosome_name_59,
                                ma,
                                "766_region.txt");
    clog << "Zika 59 search interaction seq..." << endl;
    save_interaction_region_seq(inter_reg_59,
                                chromosome_name_59,
                                chromosome_name_766,
                                ma,
                                "59_region.txt");
    
     return 0;
}

void save_interaction_region_seq(const vector<InterRegion> &region,
                                 const string &ref_trans_id,
                                 const string &homo_trans_id,
                                 const Multi_align &ma,
                                 const string &file_name)
{
    using namespace std::placeholders;
    
    if(not ma.has(ref_trans_id))
    {
        cerr << "FATAL Error: " << ref_trans_id << " is not in the multi-alignments file" << endl;
        throw runtime_error("Bad Params");
    }
    
    ofstream OUT(file_name, ofstream::out);
    if(not OUT)
    {
        cerr << "FATAL Error: " << file_name << " cannot be opened" << endl;
        throw runtime_error("IO Error");
    }
    //const vector<string>& trans_list = ma.keys();
    const string &sequence = ma.get_align(ref_trans_id).seq;
    const string &homo_sequence = ma.get_align(homo_trans_id).seq;
    const long length = sequence.size();
    
    auto raw_to_align = std::bind(&Multi_align::raw_coor_to_align_coor, _1, ref_trans_id, _2);
    auto align_to_raw = std::bind(&Multi_align::align_coor_to_raw_coor, _1, homo_trans_id, _2);
    
    for(InterRegion inter_reg: region)
    {
        long LS = max(long(inter_reg.region.first.first)-5, 1L);
        long LE = min(long(inter_reg.region.first.second)+5, length);
        long RS = max(long(inter_reg.region.second.first)-5, 1L);
        long RE = min(long(inter_reg.region.second.second)+5, length);

        //cout << "First: " << LS << "-" << LE << "\t" << RS << "-" << RE << endl;
        
        long homo_LS = align_to_raw(ma, raw_to_align(ma, LS-1)+1)+1;
        long homo_LE = align_to_raw(ma, raw_to_align(ma, LE-1)+1);
        long homo_RS = align_to_raw(ma, raw_to_align(ma, RS-1)+1)+1;
        long homo_RE = align_to_raw(ma, raw_to_align(ma, RE-1)+1);

        //cout << "Second: " << homo_LS << "-" << homo_LE << "\t" << homo_RS << "-" << homo_RE << endl;
        
        OUT << ">" << LS << "-" << LE << "<==>" << RS << "-" << RE << "\t" << homo_LS << "-" << homo_LE << "<==>" << homo_RS << "-" << homo_RE << "\t" << inter_reg.max_cov << endl;
        OUT << sequence.substr(LS-1, LE-LS) << "III" << sequence.substr(RS-1, RE-RS) << "\n";
        OUT << homo_sequence.substr(homo_LS-1, homo_LE-homo_LS) << "III" << homo_sequence.substr(homo_RS-1, homo_RE-homo_RS) << "\n" << endl;
    }
    OUT.close();
}

// to see if interaction region has more sequence difference sites than other region
void test_mutation_significant(const Circot_Mutation &mutation,
                               const vector<InterRegion> &inter_reg,
                               uLONG total_len,
                               uLONG times)
{
    auto shuffle_region = [=](const vector<pair<uLONG, uLONG>> &regions)
    ->vector<pair<uLONG, uLONG>>
    {
        vector<pair<uLONG, uLONG>> shuffled_region;
        for(pair<uLONG, uLONG> reg: regions)
        {
            uLONG end = total_len - reg.second;
            uLONG len = reg.second - reg.first;
            
            uLONG rand_start = uLONG(rand() % end);
            shuffled_region.push_back( make_pair(rand_start, rand_start+len) );
            //cout << rand_start << "-" << rand_start+len << endl;
        }
        return shuffled_region;
    };
    
    auto count_mutate = [&mutation](const vector<pair<uLONG, uLONG>> &regions)->uLONG
    {
        uLONG mut_total = 0;
        for(pair<uLONG, uLONG> reg: regions)
        {
            uLONG start = reg.first;
            uLONG end = reg.second;
            
            for(uLONG mut: mutation.mutations)
                if( small_to_large(start, mut, end) )
                    mut_total++;
        }
        return mut_total;
    };
    
    vector<pair<uLONG, uLONG>> regions;
    for(auto ireg: inter_reg)
    {
        regions.push_back( make_pair(ireg.region.first.first, ireg.region.first.second) );
        regions.push_back( make_pair(ireg.region.second.first, ireg.region.second.second) );
    }
    
    vector<uLONG> mut_count;
    for(uLONG idx=0; idx<times; idx++)
    {
        vector<pair<uLONG, uLONG>> new_region = shuffle_region(regions);
        mut_count.push_back( count_mutate(new_region) );
        //return;
    }
    uLONG true_count = count_mutate(regions);
    sort(mut_count.begin(), mut_count.end());
    
    for(uLONG xx: mut_count)
        cout << xx << "\t";
    cout << "\n";
    
    cout << "true_count: -- " << true_count << endl;
}

vector<qreal> inter_norm_factor(const vector<InterRegion> &inter_reg, qreal min_cov=15, qreal max_cov=80)
{
    vector<qreal> cov_array;
    for(InterRegion ir: inter_reg)
        cov_array.push_back(ir.max_cov);
    //qreal max_cov = 0;
    //qreal min_cov = -1UL;
    for(auto iter=cov_array.begin(); iter!=cov_array.end(); iter++)
    {
        //cout << "score: " << *iter << " ==> " << log10( *iter ) << endl;
        *iter = log10( *iter );
        //max_cov = max( max_cov, *iter );
        //min_cov = min( min_cov, *iter );
    }
    min_cov = log10(min_cov);
    max_cov = log10(max_cov);
    for(auto iter=cov_array.begin(); iter!=cov_array.end(); iter++)
        *iter = min( (*iter - min_cov) / (max_cov - min_cov), 1.0);
    cout << "Coverage Range: " << min_cov << "\t" << max_cov << endl;
    return cov_array;
}

void plot_circos(const vector<InterRegion> &inter_reg_59,
                 const vector<InterRegion> &inter_reg_766,
                 const Circot_Mutation &seq_mutations,
                 const Circot_Mutation &aa_mutations,
                 const annotation &protein_anno,
                 const vector<domain> &domain_list,
                 uLONG total_len)
{
    QPdfWriter pdf_writer(QString::fromStdString("test.pdf"));
    pdf_writer.setPageSize( QPageSize( QSizeF(1000*1.5, 1000*1.5), QPageSize::Point ) );
    pdf_writer.setPageMargins(QMarginsF(0,0,0,0));
    
    QPainter painter( &pdf_writer );
    painter.setRenderHint( QPainter::Antialiasing );
    
    static QPen pen;
    static QColor color;
    static QBrush brush;
    
    double width = painter.device()->width();
    
    double elem_len = width/6.0;
    QPointF centre(elem_len*3, elem_len*3);
    qreal radius = elem_len*2;
    
    pen.setWidth(5);
    pen.setColor(Qt::black);
    painter.setPen(pen);
    
    color = QColor("#08F229");
    color.setAlpha(150);
    painter.setBrush(color);
    
    auto degree_convert = [](qreal raw_degree){ return 90 - raw_degree < 0 ? 360 + (90 - raw_degree) : 90 - raw_degree; };
    
    qreal min_cov=25;
    qreal max_cov=100;
    
    vector<qreal> norm_59 = inter_norm_factor(inter_reg_59, min_cov, max_cov);
    vector<qreal> norm_766 = inter_norm_factor(inter_reg_766, min_cov, max_cov);
    
    size_t idx = 0;
    for(InterRegion ir: inter_reg_59)
    {
        qreal S1 = degree_convert( 1.0 * ir.region.first.first / total_len * (360 - 30) );
        qreal S2 = degree_convert( 1.0 * ir.region.second.first / total_len * (360 - 30) );
        
        qreal E1 = degree_convert( 1.0 * ir.region.first.second / total_len * (360 - 30) );
        qreal E2 = degree_convert( 1.0 * ir.region.second.second / total_len * (360 - 30) );
        
        //qSqrt(ir.max_cov);
        //color.setAlpha( min( (1.0 * (ir.max_cov - 15) / 60) * 150 + 50, 255) );
        //painter.setBrush(color);
        
        qreal alpha = norm_59.at(idx) * 180 + 40;
        color.setAlpha( alpha );
        painter.setBrush(color);
        
        plot_circos_interaction_region(centre, radius,
                                       S1, E2,
                                       E1, S2,
                                       painter);
        idx++;
    }
    
    color = QColor("#FF06FF");
    color.setAlpha(150);
    painter.setBrush(color);
    idx = 0;
    for(InterRegion ir: inter_reg_766)
    {
        qreal S1 = degree_convert( 1.0 * ir.region.first.first / total_len * (360 - 30) );
        qreal S2 = degree_convert( 1.0 * ir.region.second.first / total_len * (360 - 30) );
        
        qreal E1 = degree_convert( 1.0 * ir.region.first.second / total_len * (360 - 30) );
        qreal E2 = degree_convert( 1.0 * ir.region.second.second / total_len * (360 - 30) );
        
        //color.setAlpha( min( (1.0 * (ir.max_cov - 15) / 60) * 150 + 50, 255) );
        //painter.setBrush(color);
        
        qreal alpha = norm_766.at(idx) * 180 + 30;
        color.setAlpha( alpha );
        painter.setBrush(color);
        
        plot_circos_interaction_region(centre, radius,
                                       S1, E2,
                                       E1, S2,
                                       painter);
        idx++;
    }
    
    color = QColor("black");
    color.setAlpha(255);
    pen.setColor(color);
    painter.setPen(pen);
    for(qreal mut: seq_mutations.mutations)
    {
        qreal angle = degree_convert( 1.0 * mut / seq_mutations.total_len * (360 - 30) );
        plot_circos_line(centre, radius, radius*1.05,
                          angle, painter);
    }
    
    color = QColor("#8172b2");
    color.setAlpha(255);
    pen.setColor(color);
    painter.setPen(pen);
    for(qreal mut: aa_mutations.mutations)
    {
        qreal angle = degree_convert( 1.0 * mut / aa_mutations.total_len * (360 - 30) );
        plot_circos_line(centre, radius*1.05, radius*1.10,
                         angle, painter);
    }


    // protein structure
    color = QColor("black");
    color.setAlpha(255);
    pen.setColor(color);
    //pen.setWidth(20);
    painter.setPen(pen);
    for(size_t idx=0; idx<protein_anno.exons.size(); idx++)
    {
        
        painter.setBrush( QColor(QString::fromStdString(color_base.at(idx%color_base.size()))) );
        qreal end = degree_convert( 1.0 * get<0>(protein_anno.exons.at(idx)) / total_len * (360 - 30) );
        qreal start = degree_convert( 1.0 * get<1>(protein_anno.exons.at(idx)) / total_len * (360 - 30) );
        
        // remember the order is always anti-closewise
        plot_circos_region(centre,
                           radius*1.10,
                           radius*1.15,
                           start, end,
                           painter);
    }
    
    // domain
    color = QColor("black");
    color.setAlpha(255);
    pen.setColor(color);
    pen.setWidth(40);
    painter.setPen(pen);
    painter.setBrush( Qt::lightGray );
    for(size_t idx=0; idx<domain_list.size(); idx++)
    {
        
        qreal end = degree_convert( 1.0 * domain_list.at(idx).start / total_len * (360 - 30) );
        qreal start = degree_convert( 1.0 * domain_list.at(idx).end / total_len * (360 - 30) );
        
        plot_circos_region(centre,
                           radius*1.15,
                           radius*1.20,
                           start, end,
                           painter);
    }
    
    //sticks
    color = QColor("black");
    color.setAlpha(255);
    pen.setColor(color);
    pen.setWidth(50);
    painter.setPen(pen);
    for(size_t idx=0; idx<total_len; idx+=100)
    {
        
        qreal raw_angle = 1.0 * idx / total_len * (360 - 30);
        qreal angle = degree_convert( raw_angle );
        
        if(idx % 1000 == 0)
        {
            painter.translate(centre);
            painter.rotate(raw_angle);
            painter.translate(-centre);
            
            painter.drawText(QRectF(centre.x()-3000, centre.y()-radius*2, 6000, radius*0.75), Qt::AlignHCenter | Qt::AlignBottom, QString::number(idx));
            
            painter.translate(centre);
            painter.rotate(-raw_angle);
            painter.translate(-centre);
            
            plot_circos_line(centre, radius*1.20, radius*1.25, angle, painter);
        }
        else
            plot_circos_line(centre, radius*1.20, radius*1.225, angle, painter);
        

    }
    
    pen.setColor("black");
    pen.setWidth(50);
    painter.setPen(pen);
    
    painter.drawArc( QRectF(centre.x()-radius, centre.y()-radius, 2*radius, 2*radius ), 120*16, 330*16 );
    painter.drawArc( QRectF(centre.x()-radius*1.05, centre.y()-radius*1.05, 2*radius*1.05, 2*radius*1.05 ), 120*16, 330*16 );
    painter.drawArc( QRectF(centre.x()-radius*1.10, centre.y()-radius*1.10, 2*radius*1.10, 2*radius*1.10 ), 120*16, 330*16 );
    painter.drawArc( QRectF(centre.x()-radius*1.15, centre.y()-radius*1.15, 2*radius*1.15, 2*radius*1.15 ), 120*16, 330*16 );
    
}
















