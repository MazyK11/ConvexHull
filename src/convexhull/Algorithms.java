/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package convexhull;

import java.awt.geom.Path2D;
import java.awt.geom.Point2D;
import static java.lang.Math.abs;
import static java.lang.Math.sqrt;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Stack;
import java.util.TreeMap;
import java.util.concurrent.ThreadLocalRandom;

/**
 *
 * @author mazurd
 */
public class Algorithms {
    
    public static Path2D minimumAreaEnclosingBox(Point2D [] p){
        Path2D rectg = new Path2D.Double();
        Point2D minX = p[0];
        Point2D maxX = p[0];
        Point2D minY = p[0];
        Point2D maxY = p[0];
        // max a min
        for (Point2D pt : p) {
            if (pt.getX() < minX.getX()) {
                minX = pt;
            }
            if (pt.getX() > maxX.getX()) {
                maxX = pt;
            }
            if (pt.getY() < minY.getY()) {
                minY = pt;
            }
            if (pt.getY() > maxY.getY()) {
                maxY = pt;
            }
        }
        
        // délka mezi extrémy
        double d1 = distanceBetweenTwoPoints(minX,maxX);
        double d2 = distanceBetweenTwoPoints(minY,maxY);
        // fík a fíj       
        double fj = angle(minY.getX() - maxY.getX(),minY.getY() - maxY.getY(),
                0,minY.getY() - maxY.getY());       
        double fk = angle(minX.getX() - maxX.getX(),minX.getY() - maxX.getY(),
                minX.getX() - maxX.getX(),0);
        // plocha
        double l1 = abs(d1 * Math.cos((fj*180)/Math.PI));
        double l2 = abs(d2 * Math.cos((fk*180)/Math.PI));      
        double A =l1 * l2; 
        System.out.println(A);
        
        double Amin = A; double SumaT = 0;
        
        while(SumaT < Math.PI/4){
            double delty [] = deltaAngels(minX,maxX,minY,maxY);
            double delta = delty[0];
            for(int i = 0;i<4;i++){
                if (delta > delty[i]){
                    delta = delty[i];
                } 
            }
            
            
            
            double l1new = abs(d1 * Math.cos(((fj*180)/Math.PI) - delta));
            double l2new = abs(d1 * Math.cos(((fk*180)/Math.PI) - delta));
            double Anew = l1new * l2new; 
            
            SumaT = SumaT + delta;
        }
        
        
        
        
        
        return rectg;
    }
    
    private static double distanceBetweenTwoPoints(Point2D a, Point2D b){
        return sqrt((a.getX() - b.getX())*(a.getX() - b.getX()) + 
                (a.getY() - b.getY())*(a.getY() - b.getY()));
    }
    
    
    private static double [] deltaAngels(Point2D minX, Point2D maxX, Point2D minY,
            Point2D maxY){
        double [] d = new double[4];
        

        return d;
    }
    
    
    
    
    
    
    
    
    
    //method, which counts the dot product
    public static double dotProd(double ux,double uy,double vx,double vy){
        return ux*vx + uy*vy;
    }
    // method for counting length of vectors
    public static double len(double ux,double uy){
        return sqrt(ux*ux + uy*uy);
    }
    // method, which counts angle between two vectors 
    public static double angle(double ux, double uy, double vx, double vy){
        double skalarsoucin = dotProd(ux,uy,vx,vy);
        double ulen = len(ux,uy);
        double vlen = len(vx,vy);
        return Math.acos(skalarsoucin/(ulen*vlen));
    }
    
    
    public static Path2D jarvisScan(Point2D [] p, drawPanel d){
        Path2D poly = new Path2D.Double();
        LinkedList<Point2D> polyPoints =  new LinkedList<>();
        Point2D miny = p[0];
        // projdi pole proměnnou curPT což je Pont2D - v poli p 
        for (Point2D pt : p){
            if(pt.getY() < miny.getY()){
                miny = pt;
            }
        }
        // prvni bod do listu
        // 0 0 kvůli přímce(úhel vektory) následně vyhodím 
        // move to - kvůli path posuň se na ten první bod aby věděla odkud kreslit
        poly.moveTo(miny.getX(),miny.getY());
        
        polyPoints.add(new Point2D.Double(0,miny.getY()));
        polyPoints.add(miny);
        
        Point2D curPT = miny;
        Point2D prevpt = polyPoints.getFirst(); 
        while(true){
            double start = Double.MAX_VALUE;
            // mění se jen jednou vektor
            double ux = prevpt.getX() - curPT.getX();
            double uy = prevpt.getY() - curPT.getY();
            Point2D minPT = null;
            for(Point2D pp: p){
                if (pp == curPT || pp == prevpt){
                    continue;
                }
                // vektor body a current bod
                double vx = curPT.getX() - pp.getX();
                double vy = curPT.getY() - pp.getY();
                // výsledný úhel
                double min = angle(ux,uy,vx,vy);
                // pokud je výsledný úhel nižší něž dosavadný 
                if (start > min){
                    // start je min 
                    start = min;
                    // pamatuji minpt
                    minPT = pp;
                }
            }
            polyPoints.add(minPT);
            prevpt = curPT;
            curPT = minPT;
            if (polyPoints.getLast() == miny){
                break;
            }
        }
        d.recpt = new Point2D [polyPoints.size()-1];
        for (int i = 1;i < polyPoints.size();i++){
            poly.lineTo(polyPoints.get(i).getX(),polyPoints.get(i).getY());
            d.recpt[i-1] = polyPoints.get(i);
        }
        // line to si pamatuje co má spojit za body - přidáváme ještě jednou první bod
        poly.lineTo(polyPoints.get(1).getX(),polyPoints.get(1).getY());
        return poly;
    }
    public static Path2D quickHull(Point2D[] p) {
        // Path2D
        Path2D poly = new Path2D.Double();
        Point2D minX = p[0];
        Point2D maxX = p[0];
        // Max Min
        for (Point2D pt : p) {
            if (pt.getX() < minX.getX()) {
                minX = pt;
            }
            if (pt.getX() > maxX.getX()) {
                maxX = pt;
            }
        }

        poly.moveTo(maxX.getX(), maxX.getY());

        List<Point2D>[] s = determinant(maxX, minX, Arrays.asList(p));

        // volám obě půlky
        Point2D[] upperHull = qh(maxX, minX, s[0]);
        Point2D[] lowerHull = qh(minX, maxX, s[1]);

        for (Point2D pt : lowerHull) {
            poly.lineTo(pt.getX(), pt.getY());
        }
        poly.lineTo(minX.getX(), minX.getY());

        for (Point2D pt : upperHull) {
            poly.lineTo(pt.getX(), pt.getY());
        }

        poly.lineTo(maxX.getX(), maxX.getY());

        return poly;
    }

    private static Point2D[] qh(Point2D start, Point2D end, List<Point2D> points) {
        if (points.size() == 0) {
            return new Point2D[0];
        }
        double max = -1;
        Point2D maxPt = null;

        for (Point2D pt : points) {
            // vzdálenost bodů
            double dist = distanceFromLine(start, end, pt);
            if (dist > max) {
                max = dist;
                // pamatuji bod 
                maxPt = pt;
            }
        }
        
        // rozděl na půl
        List<Point2D>[] startPts = determinant(maxPt, start, points);
        List<Point2D>[] endPts = determinant(end, maxPt, points);
        // dokud to jde volej startHull 
        // přepisuji si body - ale o úroveň výše
        
        // pak se objeví po konci rekurzí
        Point2D[] startHull = qh(start, maxPt, startPts[1]);
        // pak volej endhull 
        Point2D[] endHull = qh(maxPt, end, endPts[1]);

        Point2D[] res = new Point2D[startHull.length + endHull.length + 1];
        int i = 0;
        for (Point2D pt : endHull) {
            res[i] = pt;
            i++;
        }
        // pokud se dostanu k tomu poslednímu a na žádné straně nic není, tak 
        //se uložím sem.
        res[i] = maxPt;
        i++;
        for (Point2D pt : startHull) {
            res[i] = pt;
            i++;
        }
        return res;
    }

    public static List<Point2D>[] determinant(Point2D start, Point2D end, List<Point2D> points) {
        // pole co uloží 2 listy
        List<Point2D>[] res = new List[2];
        res[0] = new LinkedList<>();
        res[1] = new LinkedList<>();
        // přidej do listu, podle strany na které leží od úsečky
        for (Point2D pt : points) {
            if (pt == start || pt == end) {
                continue;
            }
            double det = (end.getX() - start.getX()) * (pt.getY() - start.getY()) - (end.getY() - start.getY()) * (pt.getX() - start.getX());

            if (det > 0) {
                res[0].add(pt);
            }else if(det < 0){
                res[1].add(pt);
            }
        }
        return res;
    }

    public static double distanceFromLine(Point2D start, Point2D end, Point2D pt) {
        // bod od přímky vzdálenost
        double nom = abs((end.getY() - start.getY()) * pt.getX()
                - (end.getX() - start.getX()) * pt.getY()
                + (end.getX() * start.getY())
                - (end.getY() * start.getX()));

        double denom = sqrt((end.getY() - start.getY()) * (end.getY() - start.getY())
                + ((end.getX() - start.getX()) * (end.getX() - start.getX())));

        return nom / denom;

    }
    
    
        // Method, which does the Graham Scan algorithm 
    public static Path2D grahamScan(Point2D[] p) {
        // Path2D
        Path2D poly = new Path2D.Double();
        // minimum 
        Point2D miny = p[0];
        // find a point with minimum Y
        for (Point2D pt : p) {
            if (pt.getY() < miny.getY()) {
                miny = pt;
            }
        }
        // Point, which together with point miny representing X axis
        Point2D X;
        X = new Point2D.Double(0, miny.getY());
        // vector u
        double ux = X.getX() - miny.getX();
        double uy = X.getY() - miny.getY();
        // Creation of map, which automaticly sort values by its key
        Map<Double,Point2D> map;
        map = new TreeMap();
        // cycle through all points
        for (Point2D pp : p) {
            double a = 0;
            // if point pp is miny - set a to max value -> point will be sorted as last 
            if (pp == miny) {
                a = Double.MAX_VALUE;
            } else {
                // vector v - between point miny and point pp
                double vx = pp.getX() - miny.getX();
                double vy = pp.getY() - miny.getY();
                // calculates angle between vectors
                a = angle(ux, uy, vx, vy);
            }
            // add key - angle and value - point
            map.put(a, pp);
        }
        // stack
        Stack convexHull = new Stack();
        // Point miny and point with lowest angle  
        Point2D prevPt = miny;
        Point2D curPt = (Point2D) map.values().toArray()[0];
        // add to stack
        convexHull.push(prevPt);
        convexHull.push(curPt);
        // Cycle through all sorted points
        for (Point2D nextPt : map.values()) {
            if (nextPt == curPt) {
                continue;
            }
            // calculate determinant to find out, if nextPt is positioned left
            // or right against a line(prevPt, curPt)
            boolean d = deter(prevPt.getX(), prevPt.getY(), curPt.getX(), curPt.getY(),
                    nextPt.getX(), nextPt.getY());
            // if point is located left from line
            if (d == true) {
                // add the point to stack
                convexHull.push(nextPt);
                // change prevPt and curPt
                prevPt = curPt;
                curPt = nextPt;
            } else {
                // if not, remove current point from stack
                convexHull.pop();
                // cycle through stack to find out, if you need to remove more
                // points to make sure, taht convex hull is yet without "mistake"
                // if it is calculated for nextPt
                while (convexHull.size() > 1) {
                    // change curPt and prevPt
                    curPt = (Point2D) convexHull.pop();
                    prevPt = (Point2D) convexHull.peek();
                    // calculate determinant
                    d = deter(prevPt.getX(), prevPt.getY(), curPt.getX(), curPt.getY(),
                            nextPt.getX(), nextPt.getY());
                    // if point is located left from line
                    if (d == true) {
                        // add current point and nextPt
                        convexHull.push(curPt);
                        convexHull.push(nextPt);
                        // change prevPt and prevPt
                        prevPt = curPt;
                        curPt = nextPt;
                        // everything is as it should be, break through cycle
                        break;
                    }
                    
                }
            }
        }
        // move to minimum point
        poly.moveTo(miny.getX(), miny.getY());
        // Cycle through the stack and add points to Path2D
        while (convexHull.size() > 0) {
            Point2D w = (Point2D) convexHull.pop();
            poly.lineTo(w.getX(), w.getY());
        }
        // return Path2D
        return poly;
    }

    // method, which does determinant test
    public static boolean deter(double x1, double y1, double x2, double y2,
            double x, double y) {
        // calculate determinant
        double det = (x2 - x1) * (y - y1) - (y2 - y1) * (x - x1);
        // Point is located on one side of the line
        if (det > 0) {
            return false;
        } 
        // Point is located on the other side of the line
        return true;
    }
    //  Method, which generates set of points in defined shape
    public static Point2D [] generateSet(int npoints, double a, double b, String s){
        Point2D [] point =  new Point2D[npoints];
        Random rnd = new Random();
        double x;
        double y;
        // determine defined shape
        if (null != s)switch (s) {
            case "Rectangle":
                // generate points in a shape of rectangle
                for (int i = 0;i<npoints;i++){
                    x =  0.5 +ThreadLocalRandom.current().nextDouble(-1,1)*(a/2);
                    y =  0.5 +ThreadLocalRandom.current().nextDouble(-1,1)*(b/2);
                    point[i] = new Point2D.Double(x,y);
                }   break;
            case "Square":
                // generate points in a shape of square
                for (int i = 0;i<npoints;i++){
                    x = 0.5 +ThreadLocalRandom.current().nextDouble(-1,1)*(a/2);
                    y = 0.5 +ThreadLocalRandom.current().nextDouble(-1,1)*(a/2);
                    point[i] = new Point2D.Double(x,y);
                }   break;
            case "Circle":
                // generate points in a shape of circle
                for (int i = 0;i<npoints;i++){
                    double angle = rnd.nextDouble()*Math.PI*2;
                    x = 0.5 +(Math.cos(angle)*a*rnd.nextDouble());
                    y = 0.5 +(Math.sin(angle)*a*rnd.nextDouble());
                    point[i] = new Point2D.Double(x,y);
                }   break;
            case "Ellipse":
                // generate points in a shape of ellipse
                for (int i = 0;i<npoints;i++){
                    double angle = rnd.nextDouble()*Math.PI*2;
                    x = 0.5 +(Math.cos(angle)*a*rnd.nextDouble());
                    y = 0.5 +(Math.sin(angle)*(a/1.5)*rnd.nextDouble());
                    point[i] = new Point2D.Double(x,y);
                }   break;
            case "Cross":{
                // generate points in a shape of cross
                for (int i = 0;i<npoints;i++){
                    if(i*2 > npoints){
                        y = 0.5 +ThreadLocalRandom.current().nextDouble(-1,1)*(b);
                        x = 0.5 +(ThreadLocalRandom.current().nextDouble(-1,1)*(a))*(y-0.5);
                    }
                    else{
                        x = 0.5 +ThreadLocalRandom.current().nextDouble(-1,1)*(b);
                        y = 0.5 +(ThreadLocalRandom.current().nextDouble(-1,1)*(a))*(x-0.5);
                    }
                    point[i] = new Point2D.Double(x,y);
                }       break;
                }
            case "Triangle":{
                // generate points in a shape of triangle
                a = a/2;
                for (int i = 0;i<npoints;i++){
                    y = 0.5 +ThreadLocalRandom.current().nextDouble(-1,1)*(b/2);
                    x = 0.5 +(ThreadLocalRandom.current().nextDouble(-1,1)*(a/b))*(y-(b/2)-0.5);
                    point[i] = new Point2D.Double(x,y);
                }       break;
                }
            default:
                break;
        }
       return point; 
    }
}
