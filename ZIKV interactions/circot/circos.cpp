
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
#include <exception>

#include <QApplication>
#include <QPainter>
#include <QPen>
#include <QPainterPath>
#include <QPdfWriter>
#include <QtMath>
#include <QDebug>

#include "circos.h"

//using namespace std;
//using namespace  paris_sam_type;


bool valid_degree(qreal raw_degree)
{
    if( raw_degree >= 0 and raw_degree < 360 )
        return true;
    else
        return false;
}

// if test_degree in the range of start_degree -> end_degree with close-wise direction
bool contains(qreal start_degree, qreal end_degree, qreal test_degree)
{
    if( not valid_degree(start_degree) or not valid_degree(end_degree) or not valid_degree(test_degree) )
    {
        cerr << "Invalid Degree(from contains) -- " << start_degree << "\t" << end_degree << "\t" << test_degree << endl ;
        throw runtime_error("Bad Degree");
    }
    
    if (
        small_to_large( start_degree, test_degree, end_degree ) or
        small_to_large( test_degree, end_degree, start_degree ) or
        small_to_large( end_degree, start_degree, test_degree )
       )
        return true;
    else
        return false;
}

qreal point_dist(QPointF point_1, QPointF point_2)
{
    return qSqrt( qPow(point_1.x()-point_2.x(), 2) + qPow(point_1.y()-point_2.y(), 2) );
}

// return degree span of start -> end with close-wise direction
qreal degree_span(qreal start_degree, qreal end_degree)
{
    if( not valid_degree(start_degree) or not valid_degree(end_degree) )
    {
        cerr << "Invalid Degree(from degree_span) -- " << start_degree << "\t" << end_degree << endl ;
        throw runtime_error("Bad Degree");
    }
    
    qreal span = end_degree - start_degree;
    if(span < 0)
        span += 360;
    
    return span;
}

// Return degree with a radians input
double arc_cos(double v)
{
    return qAcos(v)/M_PI*180;
}

qreal fix_precision(qreal raw_double, int precision)
{
    qreal x = (long long)(raw_double * qPow(10, precision));
    return x / qPow(10, precision);
}


bool valid_circos_band(qreal S1, qreal E1, qreal S2, qreal E2)
{
    if( not valid_degree(S1) or not valid_degree(S2) or not valid_degree(E1) or not valid_degree(E2) )
    {
        cerr << "Invalid Degree(from valid_circos_band) -- " << S1 << "\t" << S2 << "\t" << E1 << "\t" << E2 << endl ;
        throw runtime_error("Bad Degree");
    }
    
    if  (
         small_to_large(S2, S1, E1, E2) or
         small_to_large(E2, S2, S1, E1) or
         small_to_large(E1, E2, S2, S1) or
         small_to_large(S1, E1, E2, S2) or
       
         small_to_large(E2, E1, S1, S2) or
         small_to_large(E1, S1, S2, E2) or
         small_to_large(S1, S2, E2, E1) or
         small_to_large(S2, E2, E1, S1)
        )
        return true;
    else
        return false;
}

// return the cross point of vertical line from ref_point to line(point_1, point_2)
QPointF cross_point(QPointF point_1, QPointF point_2, QPointF ref_point)
{
    
    qreal radius_1 = point_dist(point_1, ref_point);
    qreal radius_2 = point_dist(point_2, ref_point);

    if( qAbs(fix_precision(radius_1, 5)-fix_precision(radius_2, 5)) > 0.01 )
    {
        cerr << "FATAL Error: It seems that you provided bad centre or circle point: " <<radius_1 << "\t" << radius_2 << endl;
        
        throw runtime_error("Bad Cross Point");
    }
    
    qreal part_1 = ref_point.y() + ( point_1.x() - point_2.x() ) / ( point_1.y() - point_2.y() ) * ref_point.x();
    qreal part_2 = point_1.y() - ( point_1.y() - point_2.y() ) / ( point_1.x() - point_2.x() ) * point_1.x();
    qreal part_3 = ( point_1.y()-point_2.y() ) / ( point_1.x() - point_2.x() ) + ( point_1.x() - point_2.x() ) / ( point_1.y() - point_2.y() );
    
    qreal cross_x = (part_1 - part_2) / part_3;
    qreal cross_y = -(point_1.x()-point_2.x())/(point_1.y()-point_2.y())*cross_x + (ref_point.y()+(point_1.x()-point_2.x())/(point_1.y()-point_2.y())*ref_point.x());
    
    return QPointF(cross_x, cross_y );
}

qreal point_degree(QPointF ref_point, QPointF centre)
{
    qreal dist = point_dist(ref_point, centre);
    
    qreal ratio = fix_precision( (centre.y() - ref_point.y())/dist, 5 );
    double degree = arc_cos(ratio);
    //cout << "Ratio: " << ratio << endl;
    
    if(centre.x() < ref_point.x())
        degree = 360 - degree;
    
    degree = degree + 90;
    if(degree >= 360) degree -= 360;
    
    return degree;
}

// start_1 -> start_2 -> end_2 -> end_1 -> start_1 -> close
void plot_circos_interaction_region(QPointF centre, qreal radius,
                                    qreal start_angle_1, qreal end_angle_1,
                                    qreal start_angle_2, qreal end_angle_2,
                                    QPainter &painter, bool debug)
{
    // Check Valid
    if( not valid_circos_band(start_angle_1, end_angle_1, start_angle_2, end_angle_2) )
    {
        cerr << "FATAL Error: " << start_angle_1 << ", " << end_angle_1 << ", " << start_angle_2 << ", " << end_angle_2 << endl;
        throw runtime_error("Bad Circos Band");
    }
    
    QPointF start_point_1( centre.x()+radius*qCos(qDegreesToRadians(start_angle_1)),
                          centre.y()-radius*qSin(qDegreesToRadians(start_angle_1)) );
    QPointF start_point_2( centre.x()+radius*qCos(qDegreesToRadians(start_angle_2)),
                          centre.y()-radius*qSin(qDegreesToRadians(start_angle_2)) );
    
    QPointF end_point_1( centre.x()+radius*qCos(qDegreesToRadians(end_angle_1)),
                        centre.y()-radius*qSin(qDegreesToRadians(end_angle_1)) );
    QPointF end_point_2( centre.x()+radius*qCos(qDegreesToRadians(end_angle_2)),
                        centre.y()-radius*qSin(qDegreesToRadians(end_angle_2)) );
    
    if(debug)
    {
        painter.drawPoint(start_point_1);
        painter.drawPoint(start_point_2);
        painter.drawPoint(end_point_1);
        painter.drawPoint(end_point_2);
        //return;
    }
    
    QPointF cross_1 = cross_point(start_point_1, end_point_1, centre);
    QPointF cross_2 = cross_point(start_point_2, end_point_2, centre);
    
    if(debug)
    {
        painter.drawPoint(centre);
        painter.drawPoint(cross_1);
        painter.drawPoint(cross_2);
        //return;
    }
    
    QPointF ref_centre_1( centre.x() + 2.5*(cross_1.x() - centre.x()), centre.y() + 2.5*(cross_1.y() - centre.y()) );
    QPointF ref_centre_2( centre.x() + 2.5*(cross_2.x() - centre.x()), centre.y() + 2.5*(cross_2.y() - centre.y()) );
    
    if(debug)
    {
        painter.drawLine(ref_centre_2, end_point_2);
        painter.drawLine(ref_centre_2, start_point_2);

        QPen pen("black");
        pen.setWidth(10);
        painter.setPen(pen);
        
        painter.drawPoint(ref_centre_1);
        painter.drawPoint(ref_centre_2);
        //return;
    }

    qreal radius_1 = point_dist(ref_centre_1, start_point_1);
    qreal radius_2 = point_dist(ref_centre_2, start_point_2);
    
    qreal aux_start_angle_1 = point_degree(start_point_1, ref_centre_1);
    qreal aux_start_angle_2 = point_degree(start_point_2, ref_centre_2);
    qreal aux_end_angle_1 = point_degree(end_point_1, ref_centre_1);
    qreal aux_end_angle_2 = point_degree(end_point_2, ref_centre_2);
    
    qreal span;
    
    // end_1 -> start_1
    double begin_angle_1, angle_span_1;
    span = degree_span(aux_end_angle_1, aux_start_angle_1);
    if(span > 180)
    {
        begin_angle_1 = aux_end_angle_1;
        angle_span_1 = -(360 - span);
    }else{
        begin_angle_1 = aux_end_angle_1;
        angle_span_1 = span;
    }
    
    // start_2 -> end_2
    double begin_angle_2, angle_span_2;
    span = degree_span(aux_start_angle_2, aux_end_angle_2);
    if(span > 180)
    {
        begin_angle_2 = aux_start_angle_2;
        angle_span_2 = -(360 - span);
    }else{
        begin_angle_2 = aux_start_angle_2;
        angle_span_2 = span;
    }
    
    QRectF rect_1(ref_centre_1.x()-radius_1, ref_centre_1.y()-radius_1, 2*radius_1, 2*radius_1);
    QRectF rect_2(ref_centre_2.x()-radius_2, ref_centre_2.y()-radius_2, 2*radius_2, 2*radius_2);
    
    QRectF raw_rect( centre.x()-radius, centre.y()-radius, 2*radius, 2*radius );
    
    bool anticlose_wise = not contains(start_angle_1, start_angle_2, end_angle_1);
    span = degree_span(start_angle_1, start_angle_2);
    if(not anticlose_wise)
        span = -(360 - span);
    //cout << "start_angle_1 => start_angle_2: " << start_angle_1 << ", " << span << endl;
    
    // large span angle lead to line
    bool line_replace_arc = false;
    if(qAbs(angle_span_2) > 120 || qAbs(angle_span_1) > 120)
        line_replace_arc = true;
    
    QPainterPath path(start_point_1);
    path.arcTo(raw_rect, start_angle_1, span);
    if( line_replace_arc )
        path.lineTo(end_point_2);
    else
        path.arcTo(rect_2, begin_angle_2, angle_span_2);
    
    span = degree_span(end_angle_2, end_angle_1);
    if(not anticlose_wise)
        span = -(360 - span);
    path.arcTo(raw_rect, end_angle_2, span);
    
    if( line_replace_arc )
        path.lineTo(start_point_1);
    else
        path.arcTo(rect_1, begin_angle_1, angle_span_1);
    
    path.closeSubpath();
    
    painter.drawPath(path);
}


void plot_circos_interaction(QPointF centre, qreal radius,
                                    qreal start_angle, qreal end_angle,
                                    QPainter &painter, bool debug)
{
    QPointF start_point( centre.x()+radius*qCos(qDegreesToRadians(start_angle)),
                          centre.y()-radius*qSin(qDegreesToRadians(start_angle)) );
    QPointF end_point( centre.x()+radius*qCos(qDegreesToRadians(end_angle)),
                        centre.y()-radius*qSin(qDegreesToRadians(end_angle)) );
    
    if(debug)
    {
        painter.drawPoint(start_point);
        painter.drawPoint(end_point);
    }
    
    QPointF cross = cross_point(start_point, end_point, centre);
    
    if(debug)
    {
        painter.drawPoint(centre);
        painter.drawPoint(cross);
    }
    
    QPointF ref_centre( centre.x() + 2*(cross.x() - centre.x()), centre.y() + 2*(cross.y() - centre.y()) );
    
    if(debug)
    {
        painter.drawLine(ref_centre, end_point);

        QPen pen("black");
        pen.setWidth(10);
        painter.setPen(pen);
        
        painter.drawPoint(ref_centre);
    }
    
    qreal radius_1 = point_dist(ref_centre, start_point);
    
    qreal aux_start_angle = point_degree(start_point, ref_centre);
    qreal aux_end_angle = point_degree(end_point, ref_centre);
    
    qreal span;
    
    
    // start_2 -> end_2
    double begin_angle, angle_span;
    span = degree_span(aux_start_angle, aux_end_angle);
    if(span > 180)
    {
        begin_angle = aux_start_angle;
        angle_span = -(360 - span);
    }else{
        begin_angle = aux_start_angle;
        angle_span = span;
    }
    
    QRectF aux_rect(ref_centre.x()-radius_1, ref_centre.y()-radius_1, 2*radius_1, 2*radius_1);
    
    painter.drawArc( aux_rect, begin_angle*16, angle_span*16 );
}

void plot_circos_line(QPointF centre, qreal inner, qreal outer,
                       qreal degree, QPainter &painter)
{
    if(inner >= outer or inner <= 0)
    {
        cerr << "FATAL Error: invalid inner/outer parameters -- " << inner << "\t" << outer << endl;
        throw runtime_error("Bad Params");
    }
    if( not valid_degree(degree) )
    {
        cerr << "FATAL Error: invalid degree -- " << degree << endl;
        throw runtime_error("Bad Degree");
    }
    
    QPointF start_point( centre.x()+inner*qCos(qDegreesToRadians(degree)),
                        centre.y()-inner*qSin(qDegreesToRadians(degree)) );
    QPointF end_point( centre.x()+outer*qCos(qDegreesToRadians(degree)),
                      centre.y()-outer*qSin(qDegreesToRadians(degree)) );
    
    painter.drawLine(start_point, end_point);
}

void plot_circos_lines(QPointF centre, qreal inner, qreal outer,
                       vector<qreal> degrees, QPainter &painter)
{
    for(qreal degree: degrees)
        plot_circos_line(centre, inner, outer,
                         degree, painter);
}


void plot_circos_region(QPointF centre, qreal inner, qreal outer,
                        qreal start_degree, qreal end_degree,
                        QPainter &painter)
{
    if(inner >= outer or inner <= 0)
    {
        cerr << "FATAL Error: invalid inner/outer parameters -- " << inner << "\t" << outer << endl;
        throw runtime_error("Bad Params");
    }
    if( not valid_degree(start_degree) or not valid_degree(end_degree) )
    {
        cerr << "FATAL Error: invalid degree -- " << start_degree << "\t" << end_degree << endl;
        throw runtime_error("Bad Degree");
    }
    
    QRectF outRect( centre.x()-outer, centre.y()-outer, 2*outer, 2*outer );
    QRectF inRect( centre.x()-inner, centre.y()-inner, 2*inner, 2*inner );
    
    QPainterPath myPath;
    
    double start_point_x = centre.x() + inner*qSin( qDegreesToRadians(90-start_degree) );
    double start_point_y = centre.y() - inner*qCos( qDegreesToRadians(90-start_degree) );
    
    qreal span = end_degree - start_degree;
    if(span < 0)
        span += 360;
    
    myPath.moveTo(QPointF(start_point_x, start_point_y));
    myPath.arcTo(outRect, start_degree, span);
    
    myPath.arcTo(inRect, end_degree, -span);
    myPath.closeSubpath();
    painter.drawPath(myPath);
}

void test_plot_circos()
{
    QPdfWriter pdf_writer(QString::fromStdString("test.pdf"));
    pdf_writer.setPageSize( QPageSize( QSizeF(1000*1.5, 1000*1.5), QPageSize::Point ) );
    pdf_writer.setPageMargins(QMarginsF(0,0,0,0));
    
    QPainter painter( &pdf_writer );
    painter.setRenderHint( QPainter::Antialiasing );
    
    QPen pen;
    pen.setWidth(50);
    painter.setPen(pen);
    
    double width = painter.device()->width();
    
    double elem_len = width/6.0;
    QPointF centre(elem_len*3, elem_len*3);
    qreal radius = elem_len*2;
    
    painter.drawArc( QRectF(centre.x()-radius, centre.y()-radius, 2*radius, 2*radius ), 0, 360*16 );
    
    pen.setWidth(10);
    pen.setColor(Qt::red);
    painter.setPen(pen);
    QColor red("#c44e52");
    red.setAlpha(150);
    QBrush brush(red);
    painter.setBrush(brush);
    
    plot_circos_interaction_region(centre, radius,
                                   230, 340,
                                   231, 339,
                                   painter);
    plot_circos_interaction_region(centre, radius,
                                   10, 260,
                                   9, 261,
                                   painter);
    pen.setWidth(50);
    pen.setColor("#8172b2");
    painter.setPen(pen);
    plot_circos_interaction(centre, radius,
                            0, 270,
                            painter);
    plot_circos_interaction(centre, radius,
                            90, 210,
                            painter);
}



