import java.lang.Object;
import org.apache.commons.math3.analysis.polynomials.PolynomialFunction;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.NonMonotonicSequenceException;
import org.apache.commons.math3.exception.NullArgumentException;
import org.apache.commons.math3.exception.NumberIsTooSmallException;
import org.apache.commons.math3.exception.util.LocalizedFormats;
import org.apache.commons.math3.util.FastMath;
import org.apache.commons.math3.util.MathArrays;
import org.apache.commons.math3.util.Precision;

int map = 1;

XML xml;

ArrayList<IrNode> coords = new ArrayList<IrNode>();
ArrayList<IrNode> inter1 = new ArrayList<IrNode>();
ArrayList<IrNode> inter2 = new ArrayList<IrNode>();
ArrayList<double[]> roadLens = new ArrayList<double []>();
CubicSpline cs_x;
CubicSpline cs_y;
double rLen;


void setup(){
  fullScreen();
  background(40);
  parse();
  interpolate();
  drawMap();
  fractalDim();
}

void parse(){
  print("Loading XML data...\n");
  xml = loadXML("map"+map+".osm");
  XML[] roads = xml.getChildren("way");
  XML[] nodes = xml.getChildren("node");
  print("Parsing and converting coordinates...\n");
  for(int i = 0;i<roads.length;i++){
    XML[] points = roads[i].getChildren("nd");
    for(int j = 0;j<points.length;j++)
    {
      String point = points[j].getString("ref");
      for(int k = 0;k<points.length;k++)
      {
        if(point.equals(nodes[k].getString("id"))){
          double [] irCo = convert(nodes[k].getDouble("lat"),nodes[k].getDouble("lon"));
          IrNode n = new IrNode(irCo[0],irCo[1]);
          coords.add(n);
        }
      }
    }
  }
  sortAr();
}

double[] convert(double lat, double lon){ 
  double deg2rad = Math.PI / 180;
  double phi = lat * deg2rad;      // convert latitude to radians
  double lam = lon * deg2rad;      // convert longitude to radians

  // Ireland
  double a = 6377340.189;      // OSI semi-major
  double b = 6356034.447;      // OSI semi-minor
  double e0 = 200000;          // OSI easting of false origin
  double n0 = 250000;          // OSI northing of false origin
  double f0 = 1.000035;        // OSI scale factor on central meridian
  double e2 = 0.00667054015;   // OSI eccentricity squared
  double lam0 = -0.13962634015954636615389526147909;   // OSI false east
  double phi0 = 0.93375114981696632365417456114141;    // OSI false north

  double af0 = a * f0;
  double bf0 = b * f0;
  // easting
  double slat2 = Math.sin(phi) * Math.sin(phi);
  double nu = af0 / (Math.sqrt(1 - (e2 * (slat2))));
  double rho = (nu * (1 - e2)) / (1 - (e2 * slat2));
  double eta2 = (nu / rho) - 1;
  double p = lam - lam0;
  double IV = nu * Math.cos(phi);
  double clat3 = Math.pow(Math.cos(phi), 3);
  double tlat2 = Math.tan(phi) * Math.tan(phi);
  double V = (nu / 6) * clat3 * ((nu / rho) - tlat2);
  double clat5 = Math.pow(Math.cos(phi), 5);
  double tlat4 = Math.pow(Math.tan(phi), 4);
  double VI = (nu / 120) * clat5 * ((5 - (18 * tlat2)) + tlat4 + (14 * eta2) - (58 * tlat2 * eta2));
  double east = e0 + (p * IV) + (Math.pow(p, 3) * V) + (Math.pow(p, 5) * VI);
  // northing
  double n = (af0 - bf0) / (af0 + bf0);
  double M = Marc(bf0, n, phi0, phi);
  double I = M + (n0);
  double II = (nu / 2) * Math.sin(phi) * Math.cos(phi);
  double III = ((nu / 24) * Math.sin(phi) * Math.pow(Math.cos(phi), 3)) * (5 - Math.pow(Math.tan(phi), 2) + (9 * eta2));
  double IIIA = ((nu / 720) * Math.sin(phi) * clat5) * (61 - (58 * tlat2) + tlat4);
  double north = I + ((p * p) * II) + (Math.pow(p, 4) * III) + (Math.pow(p, 6) * IIIA);

  //println("E = " + east + ", N = " + north); 
  double [] irCo = {east, north};
  return irCo;
}

private double Marc(double bf0, double n, double phi0, double phi)
{
  double Marc = bf0 * (((1 + n + ((5 / 4) * (n * n)) + ((5 / 4) * (n * n * n))) * (phi - phi0))
  - (((3 * n) + (3 * (n * n)) + ((21 / 8) * (n * n * n))) * (Math.sin(phi - phi0)) * (Math.cos(phi + phi0)))
  + ((((15 / 8) * (n * n)) + ((15 / 8) * (n * n * n))) * (Math.sin(2 * (phi - phi0))) * (Math.cos(2 * (phi + phi0))))
  - (((35 / 24) * (n * n * n)) * (Math.sin(3 * (phi - phi0))) * (Math.cos(3 * (phi + phi0)))));
  return (Marc);
}

void sortAr(){
  double [] easts = new double[coords.size()];
  
  for(int i = 0;i<coords.size();i++){
    easts[i] = coords.get(i).getEast();
  }
  for(int i = 0; i<coords.size();i++){
    for(int j = 0;j<coords.size();j++){
      if(easts[j]<easts[i]){
        double temp = easts[j];
        easts[j] = easts[i];
        easts[i] = temp;
        
        IrNode tempN = coords.get(j);
        coords.set(j,coords.get(i));
        coords.set(i,tempN);
      }
    }
  }
  
  
  for(int i = 0;i<coords.size();i++){
    for(int j = i;j<coords.size();j++){
      if((coords.get(i).getEast()==coords.get(j).getEast())&&(coords.get(i).getNor()==coords.get(j).getNor())&&(j!=i)){
        coords.remove(j);
      }
    }
  }
}

double dist(IrNode n1, IrNode n2){
  double dis;
  
  dis = sqrt(pow((float)(n2.getEast()-n1.getEast()),2.0)+pow((float)(n2.getNor()-n1.getNor()),2.0));
  
  return dis;
}

void interpolate(){
  print("Interpolating...\n");
  double [] x1 = new double[coords.size()];
  double [] y1 = new double[coords.size()];
    
  double [] x1t = new double[coords.size()];
  double [] y1t = new double[coords.size()];
  
  for(int i = 0;i<coords.size();i++){
    IrNode n = coords.get(i);
    x1[i] = n.getEast();
    y1[i] = n.getNor();
  }
  
  x1t[0] = 0.0;
  y1t[0] = 0.0;
  double dis = 0.0;
  for(int i = 0;i<coords.size()-1;i++){
    IrNode n1 = coords.get(i);
    IrNode n2 = coords.get(i+1);
    dis += dist(n1,n2);

    x1t[i+1] = dis;
    y1t[i+1] = dis;
  }
  rLen = dis;
  
  cs_x = new CubicSpline(x1t,x1);
  cs_y = new CubicSpline(y1t,y1);
}

void drawMap(){
  double eas;
  double nor;
  //print(coords.size()+"\n");
  double falsex = -coords.get(0).getEast();
  double falsey = coords.get(0).getNor();
  
  stroke(255,0,0);
  fill(255,0,0);
  for(int i = 0;i<coords.size();i++){
    eas = ((-coords.get(i).getEast() - falsex))+30;
    nor = (-(coords.get(i).getNor() - falsey))+600;
    
    ellipse((float)eas,(float)nor,20.0,20.0);
  }
  
  stroke(0,255,0);
  fill(0,0,255);
  
  for(double i = 0;i<pow(coords.size(),coords.size());i+=rLen/100.0){
    try{
      eas = ((-cs_x.interpolate(i) - falsex))+30;
      nor = (-(cs_y.interpolate(i) - falsey))+600;
    }
    catch(IllegalArgumentException e){
      break;
    }
    inter1.add(new IrNode(eas,nor));
    ellipse((float)eas,(float)nor,2.0,2.0);
  }
  
  stroke(0,0,255);
  fill(0,0,255);
  for(double i = 0;i<pow(coords.size(),coords.size());i+=rLen/10.0){
    try{
      eas = ((-cs_x.interpolate(i) - falsex))+30;
      nor = (-(cs_y.interpolate(i) - falsey))+600;
    }
    catch(IllegalArgumentException e){
      break;
    }
    inter2.add(new IrNode(eas,nor));
    ellipse((float)eas,(float)nor,10.0,10.0);
  }
}

void fractalDim(){
  print("Calculating fractal dimension...\n");
  double eas;
  double nor;
  
  double falsex = -coords.get(0).getEast();
  double falsey = coords.get(0).getNor();
  
  ArrayList<IrNode> r1 = new ArrayList<IrNode>();
  ArrayList<IrNode> r2 = new ArrayList<IrNode>();
  ArrayList<IrNode> r3 = new ArrayList<IrNode>();
  ArrayList<IrNode> r4 = new ArrayList<IrNode>();
  ArrayList<IrNode> r5 = new ArrayList<IrNode>();
  ArrayList<IrNode> r6 = new ArrayList<IrNode>();
  ArrayList<IrNode> r7 = new ArrayList<IrNode>();
  ArrayList<IrNode> r8 = new ArrayList<IrNode>();
  ArrayList<IrNode> r9 = new ArrayList<IrNode>();
  ArrayList<IrNode> r10 = new ArrayList<IrNode>();
  
  for(double i = 0;i<pow(coords.size(),coords.size());i+=rLen/1.0){
    try{
      eas = ((-cs_x.interpolate(i) - falsex))+1000;
      nor = (-(cs_y.interpolate(i) - falsey))+600;
      
      r1.add(new IrNode(eas,nor));
    }
    catch(IllegalArgumentException e){
      break;
    }
  }
  for(double i = 0;i<pow(coords.size(),coords.size());i+=rLen/2.0){
    try{
      eas = ((-cs_x.interpolate(i) - falsex))+1000;
      nor = (-(cs_y.interpolate(i) - falsey))+600;
      
      r2.add(new IrNode(eas,nor));
    }
    catch(IllegalArgumentException e){
      break;
    }
  }
  for(double i = 0;i<pow(coords.size(),coords.size());i+=rLen/4.0){
    try{
      eas = ((-cs_x.interpolate(i) - falsex))+1000;
      nor = (-(cs_y.interpolate(i) - falsey))+600;
      r3.add(new IrNode(eas,nor));
    }
    catch(IllegalArgumentException e){
      break;
    }
  }
  for(double i = 0;i<pow(coords.size(),coords.size());i+=rLen/8.0){
    try{
      eas = ((-cs_x.interpolate(i) - falsex))+1000;
      nor = (-(cs_y.interpolate(i) - falsey))+600;
      
      r4.add(new IrNode(eas,nor));
    }
    catch(IllegalArgumentException e){
      break;
    }
  }
  for(double i = 0;i<pow(coords.size(),coords.size());i+=rLen/16.0){
    try{
      eas = ((-cs_x.interpolate(i) - falsex))+1000;
      nor = (-(cs_y.interpolate(i) - falsey))+600;
      
      r5.add(new IrNode(eas,nor));
    }
    catch(IllegalArgumentException e){
      break;
    }
  }
  for(double i = 0;i<pow(coords.size(),coords.size());i+=rLen/32.0){
    try{
      eas = ((-cs_x.interpolate(i) - falsex))+1000;
      nor = (-(cs_y.interpolate(i) - falsey))+600;
      
      r6.add(new IrNode(eas,nor));
    }
    catch(IllegalArgumentException e){
      break;
    }
  }
  for(double i = 0;i<pow(coords.size(),coords.size());i+=rLen/64.0){
    try{
      eas = ((-cs_x.interpolate(i) - falsex))+1000;
      nor = (-(cs_y.interpolate(i) - falsey))+600;
      
      r7.add(new IrNode(eas,nor));
    }
    catch(IllegalArgumentException e){
      break;
    }
  }
  for(double i = 0;i<pow(coords.size(),coords.size());i+=rLen/128.0){
    try{
      eas = ((-cs_x.interpolate(i) - falsex))+1000;
      nor = (-(cs_y.interpolate(i) - falsey))+600;
      r8.add(new IrNode(eas,nor));
    }
    catch(IllegalArgumentException e){
      break;
    }
  }
  for(double i = 0;i<pow(coords.size(),coords.size());i+=rLen/256.0){
    try{
      eas = ((-cs_x.interpolate(i) - falsex))+1000;
      nor = (-(cs_y.interpolate(i) - falsey))+600;
      
      r9.add(new IrNode(eas,nor));
    }
    catch(IllegalArgumentException e){
      break;
    }
  }
  for(double i = 0;i<pow(coords.size(),coords.size());i+=rLen/512.0){
    try{
      eas = ((-cs_x.interpolate(i) - falsex))+1000;
      nor = (-(cs_y.interpolate(i) - falsey))+600;
      
      r10.add(new IrNode(eas,nor));
    }
    catch(IllegalArgumentException e){
      break;
    }
  }
  
  double dis1 = 0.0;
  double dis2 = 0.0;
  double dis3 = 0.0;
  double dis4 = 0.0;
  double dis5 = 0.0;
  double dis6 = 0.0;
  double dis7 = 0.0;
  double dis8 = 0.0;
  double dis9 = 0.0;
  double dis10 = 0.0;
  
  for(int i = 1; i<r1.size();i++){
    dis1 += dist(r1.get(i-1),r1.get(i));
  }
  
  for(int i = 1; i<r2.size();i++){
    dis2 += dist(r2.get(i-1),r2.get(i));
  }
  
  for(int i = 1; i<r3.size();i++){
    dis3 += dist(r3.get(i-1),r3.get(i));
  }
  
  for(int i = 1; i<r4.size();i++){
    dis4 += dist(r4.get(i-1),r4.get(i));
  }
  
  for(int i = 1; i<r5.size();i++){
    dis5 += dist(r5.get(i-1),r5.get(i));
  }
  
  for(int i = 1; i<r6.size();i++){
    dis6 += dist(r6.get(i-1),r6.get(i));
  }
  
  for(int i = 1; i<r7.size();i++){
    dis7 += dist(r7.get(i-1),r7.get(i));
  }
  
  for(int i = 1; i<r8.size();i++){
    dis8 += dist(r8.get(i-1),r8.get(i));
  }
  
  for(int i = 1; i<r9.size();i++){
    dis9 += dist(r9.get(i-1),r9.get(i));
  }
  
  for(int i = 1; i<r10.size();i++){
    dis10 += dist(r10.get(i-1),r10.get(i));
  }
  
  
  double xBar = (rLen/1.0+rLen/2.0+rLen/4.0+rLen/8.0+rLen/16.0+rLen/32.0+rLen/64.0+rLen/128.0+rLen/256.0+rLen/512.0)/10.0;
  double yBar = (dis1+dis2+dis3+dis4+dis5+dis6+dis7+dis8+dis9+dis10)/10.0;
  
  double slopeTop = (((rLen/1.0-xBar)*(dis1-yBar))+((rLen/2.0-xBar)*(dis2-yBar))+((rLen/4.0-xBar)*(dis3-yBar))+((rLen/8.0-xBar)*(dis4-yBar))+((rLen/16.0-xBar)*(dis5-yBar))+((rLen/32.0-xBar)*(dis6-yBar))+((rLen/64.0-xBar)*(dis7-yBar))+((rLen/128.0-xBar)*(dis8-yBar))+((rLen/256.0-xBar)*(dis9-yBar))+((rLen/512.0-xBar)*(dis10-yBar)));
  double slopeBottom = pow((float)(rLen/1.0-xBar),2.0)+pow((float)(rLen/2.0-xBar),2.0)+pow((float)(rLen/4.0-xBar),2.0)+pow((float)(rLen/8.0-xBar),2.0)+pow((float)(rLen/16.0-xBar),2.0)+pow((float)(rLen/32.0-xBar),2.0)+pow((float)(rLen/64.0-xBar),2.0)+pow((float)(rLen/128.0-xBar),2.0)+pow((float)(rLen/256.0-xBar),2.0)+pow((float)(rLen/512.0-xBar),2.0);
  double m = slopeTop/slopeBottom;
  print("Fractal dimension: "+(1-m) + "\n");
  fill(255,0,0);
  textSize(40);
  text("Fractal Dimension: "+(1-m),30,80);
  text("Legend: Origional map points",935,80);
  fill(0,255,0);
  text("        Interpolated at (origional length)/100",1000,120);
  fill(0,0,255);
  text("        Interpolated at (origional length)/10",1000,160);
  print("Complete");
}
