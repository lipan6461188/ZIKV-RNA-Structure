
#include <iostream>
#include <QApplication>

using namespace std;

template<typename T>
bool small_to_large(T t1, T t2, T t3)
{
    if( t1 <= t2 and t2 <= t3 )
        return true;
    else
        return false;
}

template<typename T>
bool small_to_large(T t1, T t2, T t3, T t4)
{
    if( t1 <= t2 and t2 <= t3 and t3 <= t4 )
        return true;
    else
        return false;
}

bool valid_degree(qreal raw_degree);
// if test_degree in the range of start_degree -> end_degree with close-wise direction
bool contains(qreal start_degree, qreal end_degree, qreal test_degree);
qreal point_dist(QPointF point_1, QPointF point_2);
// return degree span of start -> end with close-wise direction
qreal degree_span(qreal start_degree, qreal end_degree);
// Return degree of arccos with a ratio input
double arc_cos(double v);
qreal fix_precision(qreal raw_double, int precision);
bool valid_circos_band(qreal S1, qreal E1, qreal S2, qreal E2);
// return the cross point of vertical line from ref_point to line(point_1, point_2)
QPointF cross_point(QPointF point_1, QPointF point_2, QPointF ref_point);
// return degree of ref_point with reference of centre
qreal point_degree(QPointF ref_point, QPointF centre);


/************************

 Plot circos interaction
 
 ************************/
 
// start_1 -> start_2 -> end_2 -> end_1 -> start_1 -> close
void plot_circos_interaction_region(QPointF centre, qreal radius,
                                    qreal start_angle_1, qreal end_angle_1,
                                    qreal start_angle_2, qreal end_angle_2,
                                    QPainter &painter, bool debug=false);
void plot_circos_interaction(QPointF centre, qreal radius,
                             qreal start_angle, qreal end_angle,
                             QPainter &painter, bool debug=false);

/************************
 
 Plot circos surrounded lines
 
 ************************/

void plot_circos_line(QPointF centre, qreal inner, qreal outer,
                       qreal degree, QPainter &painter);

void plot_circos_lines(QPointF centre, qreal inner, qreal outer,
                      vector<qreal> degrees, QPainter &painter);

/************************
 
 Plot circos surrounded regions
 
 ************************/

void plot_circos_region(QPointF centre, qreal inner, qreal outer,
                        qreal start_degree, qreal end_degree,
                        QPainter &painter);






















