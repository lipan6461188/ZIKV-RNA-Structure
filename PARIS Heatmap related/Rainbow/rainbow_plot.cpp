#include <QApplication>
#include <QPainter>
#include <iostream>
#include <fstream>
#include <QPdfWriter>
#include <unordered_map>
#include <vector>
#include <string_split.h>
using namespace  std;

void print_usage(const string &argv_0)
{
    cerr << "Usage: " << argv_0 << " input_rainbow_file output_pdf_file\n"
         << "\t##what it is...\n"
         << "\tinput_rainbow_file: a file to label interacted bases:\n"
         << "\t\t#len:10675\n"
         << "\t\t11      69      blue    1       1\n"
         << "\t\t12      68      blue    1       1\n"
         << "\t\t[leftbase]      [rightbase]      [color]    [thickness 1-5]       [solid 1/ dashed 2]\n"
         << "\toutput_pdf_file: output pdf file\n\n"
         << "\tAuthor: Li Pan\n"
         << "\tEmail: hnsfyfyzlp@126.com\n"
         << "\tDate: 2017-8-6\n\n";
}

using uINT = unsigned int;
using uLONG = unsigned long;
using uLONGLONG = unsigned long long;
using small_code = unsigned char;

struct arcLine
{
    arcLine(uINT arcLeft, uINT arcRight, string color, small_code weight, small_code line_type):
        arcLeft(arcLeft), arcRight(arcRight), color(color), weight(weight),
        line_type(line_type){}
    uINT arcLeft;
    uINT arcRight;
    string color;
    small_code weight;
    small_code line_type;

    friend ostream&operator <<(ostream &, const arcLine &);
};

ostream&operator <<(ostream &OUT, const arcLine &arc)
{
    OUT << arc.arcLeft << '\t' << arc.arcRight << '\t' << arc.color << '\t' << static_cast<int>(arc.weight)
        << static_cast<int>(arc.line_type);
    return OUT;
}

void read_data(const string &data_file_name, vector<arcLine> &data, unordered_map<string, string> &head)
{
    data.clear();
    head.clear();

    ifstream IN(data_file_name, ifstream::in);
    if(not IN)
    {
        cerr << "FATAL Error: open file " << data_file_name << " Failed" << endl;
        exit(-1);
    }
    string this_line;
    uLONG len = 0;
    while(getline(IN, this_line))
    {
        if(this_line.empty())
            continue;
        if(this_line.at(0) == '#')
        {
            trim(this_line, '#');
            vector<string> headItems;
            split(this_line, ':', headItems);
            if(headItems.size() != 2)
            {
                cerr << "Warning: invalid head: " << this_line << " Skip it" << endl;
                continue;
            }
            trim(headItems.at(0), ' ');
            trim(headItems.at(1), ' ');
            head.insert(make_pair(headItems.at(0), headItems.at(1)));
        }else{
            if(len == 0)
            {
                auto lenInter = head.find("len");
                if( lenInter != head.end() )
                {
                    len = stoul(lenInter->second);
                }else{
                    cerr << "Fatal Error: #len head line not be found" << endl;
                    exit(-1);
                }
            }
            vector<string> arcItems;
            split(this_line, '\t', arcItems);
            if(arcItems.size() != 5)
            {
                cerr << "Warning: invalid arc line: " << this_line << " Skip it" << endl;
                continue;
            }
            uLONG arcLeft = stoul(arcItems.at(0));
            uLONG arcRight = stoul(arcItems.at(1));
            if(arcLeft >= arcRight or arcRight > len or arcLeft < 1)
            {
                cerr << "FATAL Error: invalid arc line: " << this_line << endl;
                exit(-1);
            }
            data.emplace_back( arcLeft, arcRight, arcItems.at(2), stoul(arcItems.at(3)), stoul(arcItems.at(4)) );
        }
    }
   // if(head.find("len"))
}

void draw_arc(const vector<arcLine> &data, uLONG length, const string &pdf_file_name)
{
    double suggest_width = max(1.25*length, 80.0);
    double suggest_height = max(0.78*length, 50.0);

    QPdfWriter pdf_writer(QString::fromStdString(pdf_file_name));
    pdf_writer.setPageSize( QPageSize( QSizeF(suggest_width, suggest_height), QPageSize::Point ) );
    pdf_writer.setPageMargins(QMarginsF(0,0,0,0));

    QPainter painter( &pdf_writer );
    painter.setRenderHint( QPainter::Antialiasing );
   // painter.setOpacity(0.5);

    double width = painter.device()->width();
    double height = painter.device()->height();

    cout << "total_width: " << width << " total_height:" << height << endl;

    double x_start = 100;
    double x_end = width - x_start;

    QLineF start_line(x_start, 0, x_start, height);
    QLineF end_line(x_end, 0, x_end, height);

    if(x_end - x_start < length*5)
    {
        cerr << "FATAL Error: Image is too small" << endl;
        exit(-1);
    }

    double base_width = double(x_end - x_start + 1) / length; //64

    QPen pen;

    QVector<qreal> dashes, solid;
    // dash line
    dashes << 1 << 4 << 1 << 4;
    // solid line
    solid << 1 << 1;
    for(const arcLine &arc: data)
    {
        double arc_width = (arc.arcRight - arc.arcLeft) * base_width;
        double x = x_start + (arc.arcLeft - 1) * base_width;
        double y = height - arc_width / 2 - 50;
       // cout << "x: " << x << " y: " << y << " width:" << arc_width << " right:" << x+arc_width << endl;

        QColor color(arc.color.c_str());
        //color.setAlphaF(0.5);
        if(arc.arcRight>10575)
            cout << "x: " << x << " y: " << y << " width:" << arc_width << " right:" << x+arc_width << endl;
        pen.setColor(color);
        pen.setWidthF(arc.weight);
        if(arc.line_type == 1)
            pen.setDashPattern(solid);
        else
            pen.setDashPattern(dashes);
        painter.setPen(pen);

        painter.drawArc( x,y,arc_width,arc_width,0*16,180*16 );
    }
    pen.setColor(Qt::black);
    pen.setDashPattern(solid);
    painter.setPen(pen);
    painter.drawLine(start_line);
    painter.drawLine(end_line);
}

int main( int argc, char *argv[] )
{

    if(argc != 3)
    {
        print_usage(argv[0]);
        exit(-1);
    }

    vector<arcLine> data;
    unordered_map<string, string> head;
    read_data(argv[1], data, head);
    draw_arc(data, stoul(head.at("len")), argv[2]);

    return 0;
}
