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
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;
import java.util.ListIterator;
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
    
    public static final double EPSILON = 0.0001;
    public enum OrientationEnum {CW,COLINEAR,CCW}
    
    public static Path2D minimumAreaEnclosingBox(Point2D [] p){
        Path2D rectg = new Path2D.Double();
        Point2D minX = p[0]; Point2D maxX = p[0];
        Point2D minY = p[0]; Point2D maxY = p[0];
        // pole pro body, které násedují po extrémních bodech - pro počítání
        // úhlů 
        Point2D [] nextPts = new Point2D [4];
        // max a min
        Point2D nextPt;
        int i = 1;
        for (Point2D pt : p) {
            
            if (i >= p.length){
                nextPt = p[0];
            }
            else {
                nextPt = p[i];
            }
            
            if (pt.getX() <= minX.getX()) {
                minX = pt;
                nextPts[0] = nextPt;
            }
            if (pt.getX() >= maxX.getX()) {
                maxX = pt;
                nextPts[1] = nextPt;
            }
            if (pt.getY() <= minY.getY()) {
                minY = pt;
                nextPts[2] = nextPt;
            }
            if (pt.getY() >= maxY.getY()) {
                maxY = pt;
                nextPts[3] = nextPt;
            }
            i++;
        }

        // fík a fíj       
        double fj = angle(minY.getX() - maxY.getX(),minY.getY() - maxY.getY(),
                0,minY.getY() - maxY.getY());       
        double fk = angle(minX.getX() - maxX.getX(),minX.getY() - maxX.getY(),
                minX.getX() - maxX.getX(),0);
 
//        fj = (fj*180)/Math.PI;
//        fk = (fk*180)/Math.PI; 

        // rohy obdélníku
        Point2D [] corners = new Point2D[4];
        corners[0] = new Point2D.Double(minX.getX(),minY.getY());
        corners[1] = new Point2D.Double(maxX.getX(),maxY.getY());
        corners[2] = new Point2D.Double(maxX.getX(),minY.getY());
        corners[3] = new Point2D.Double(minX.getX(),maxY.getY());
        // kopie rohů
        Point2D [] c = corners;
        
         // délka mezi extrémy  - strany obdélníku
        double l1 = maxX.getX() - minX.getX();
        double l2 = maxY.getY() - minY.getY();
        // plocha
        double A = l1 * l2;
//        System.out.println(A);
        
        double Amin = A; double SumaT = 0; double rot = 0;double Aprev = A;
        // algoritmus
        while(SumaT < Math.PI/4){
            // nalezení delta min
            double t = deltamin(minX,maxX,minY,maxY,nextPts,corners);
            double d = (t*180)/Math.PI;
            // počítat strany? jak zmenšit ten rotovaný obdélník?
            
            // otočení obdélníku
            corners = rotate(corners,t);
            
            System.out.println(t);
            // výpočet plochy? podle vzorce v prezentaci
            double sin = Math.sin(t);
            double cos = Math.cos(t);
            double Anew  = Aprev*(Math.pow(cos,2) + 
                    (Math.tan(fk) - Math.tan(fj))* Math.cos(t) *Math.sin(t) - 
                    Math.tan(fk)*Math.tan(fj)*(Math.pow(sin,2)));

            if (abs(Anew) < Amin){
                Amin = Anew;
                rot = rot + t;
            }       
            SumaT = SumaT + d;
            Aprev = Anew;
        }
        System.out.println(rot);
        // otoč rohy na výsledný obdélník
        c = rotate(c,rot);
   
        rectg.moveTo(c[0].getX(),c[0].getY());
        rectg.lineTo(c[2].getX(),c[2].getY());
        rectg.lineTo(c[1].getX(),c[1].getY());
        rectg.lineTo(c[3].getX(),c[3].getY());
        rectg.lineTo(c[0].getX(),c[0].getY());

        return rectg;
    }
    
    private static Point2D [] rotate(Point2D [] corners,double rotate){
        Point2D [] newPts = new Point2D[4];
        
        double centerX = (corners[0].getX() + corners[1].getX())/2;
        double centerY = (corners[0].getY() + corners[1].getY())/2;
        
        for(int i = 0;i<4;i++){
            double newX = centerX + (corners[i].getX()-centerX)*Math.cos(rotate) - 
                                    (corners[i].getY()-centerY)*Math.sin(rotate);
            double newY = centerY + (corners[i].getX()-centerX)*Math.sin(rotate) + 
                                    (corners[i].getY()-centerY)*Math.cos(rotate);
           newPts[i] = new Point2D.Double(newX,newY);
        }
        
        return newPts;
    }
    
//    private static double distanceBetweenTwoPoints(Point2D a, Point2D b){
//        return sqrt((a.getX() - b.getX())*(a.getX() - b.getX()) + 
//                (a.getY() - b.getY())*(a.getY() - b.getY()));
//    }
      
    private static double deltamin(Point2D minX, Point2D maxX, Point2D minY,
            Point2D maxY, Point2D [] nextPts, Point2D [] corners){
        double t = 0;
        Point2D [] maxMin = {minX,maxX,minY,maxY};
        
        double mint = Double.MAX_VALUE; 
        double ux; double uy; double vx; double vy;
//        System.out.println("uhly");
        for (int i = 0;i<nextPts.length;i++){
            
            ux = nextPts[i].getX()-maxMin[i].getX();
            uy = nextPts[i].getY()-maxMin[i].getY();
            vx = corners[i].getX()-maxMin[i].getX();
            vy = corners[i].getY()-maxMin[i].getY();

            t = angle(ux,uy,vx,vy);
//            System.out.println((t*180)/Math.PI);
            if (t < mint){
                mint = t;
            }
        }
//        System.out.println("konec");
        return mint;
    }
    
    public static OrientationEnum getOrientation(Point2D p1, Point2D p2, Point2D p3){
         double val = (p2.getY() - p1.getY()) * (p3.getX() - p2.getX()) - 
                  (p2.getX() - p1.getX()) * (p3.getY() - p2.getY()); 
         if (abs(val) < EPSILON){
             return OrientationEnum.COLINEAR;
         }else if (val > 0){
             return OrientationEnum.CW;
         }else {
             return OrientationEnum.CCW;
         }
    }
        
    public static Path2D sweepHull(Point2D [] points){
        Path2D hull;
        hull = new Path2D.Double();
        // parametry -> Double compare(souřadnice X)
        Arrays.sort(points, (Point2D p1, Point2D p2) -> Double.compare(p1.getX(), p2.getX()));
        
        List<Point2D> upperHull;
        List<Point2D> lowerHull;
        
        upperHull = new LinkedList<>();
        lowerHull = new LinkedList<>();
        // přidám prvky do konvexních obálek
        for (Point2D pt: points){
            upperHull.add(pt);
            lowerHull.add(pt);
        }
        // fixnu konvexitu - void jelikož používám iterator
        fixConvexity(upperHull, OrientationEnum.CW);
        fixConvexity(lowerHull, OrientationEnum.CCW);
        
        // přidám do path
        hull.moveTo(lowerHull.get(0).getX(), lowerHull.get(0).getY());
        
        for (Point2D pt : lowerHull){
            hull.lineTo(pt.getX(), pt.getY());
        }
        // otočím pole horní obálky, aby to bylo správně spojený
        Collections.reverse(upperHull);
        
        for (Point2D pt : upperHull){
            hull.lineTo(pt.getX(), pt.getY());
        }
        
        return hull;
    }
    
    public static void fixConvexity(List<Point2D> points,OrientationEnum orientation){
        if (points.size() < 3){
            // void jelikož nešahám do listu ale pomocí iterátoru ho měním
            return;
        }
        ListIterator<Point2D> iterator;
        // iterátor - pracuje s prvky (boxy v listu) - umí next
        iterator = points.listIterator();
        Point2D prev;
        Point2D cur;
        Point2D next;
        prev = iterator.next();
        cur = iterator.next();
        next = iterator.next();

        // jeď
        while (true){
            if (getOrientation(prev,cur,next) == orientation){
                prev = cur;
                cur = next;
                if (! iterator.hasNext()){
                    break;
                }
                next = iterator.next();
                continue;
            }
            // iterátor je v mezerách
            // next
            iterator.previous();
            // cur
            iterator.previous();
            // delete cur
            iterator.remove();
            
            // před prev
            iterator.previous();
            // ! kdyý už nemůže jít dál (prev)
            if (!iterator.hasPrevious()){
                // před prev
                iterator.next();
                // za cur 
                cur = iterator.next();
                if (! iterator.hasNext()){
                    break;
                }
                // za next
                next = iterator.next();   
            }
            else{
            // zpátky za prev 
            iterator.next();
            // nastavím cur a prev
            cur = iterator.previous();
            prev = iterator.previous();
            // vrátím se zpátky do nextu
            iterator.next();
            iterator.next();
            iterator.next();
            }
        }    
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
        
        polyPoints.add(new Point2D.Double(1,miny.getY()));
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
                if (start >= min){
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
        d.recpt = new Point2D [polyPoints.size()-2];
        for (int i = 2;i < polyPoints.size();i++){
            poly.lineTo(polyPoints.get(i).getX(),polyPoints.get(i).getY());
            d.recpt[i-2] = polyPoints.get(i);
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
