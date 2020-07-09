class IrNode{
  private double nor;
  private double east;
  IrNode(double east, double nor){
    this.east = east;
    this.nor = nor;
  }
  double getNor(){
    return nor;
  }
  double getEast(){
   return east; 
  }
  void setNor(double nor){
    this.nor = nor;
  }
  void setEast(double east){
    this.east = east;
  }
}
