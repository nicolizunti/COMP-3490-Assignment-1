class Triangle {
  Triangle(float[] V1, float[] V2, float[] V3) {  // does DEEP copy!!
    v1 = Arrays.copyOf(V1, V1.length); 
    v2 = Arrays.copyOf(V2, V2.length);
    v3 = Arrays.copyOf(V3, V3.length);
  }

  // position data. in 3D space
  float[] v1; // 3 triangle vertices
  float[] v2;
  float[] v3;

  // projected data. On the screen raster
  float[] pv1; // (p)rojected vertices
  float[] pv2;
  float[] pv3;

  // add other things as needed, like normals (face, vectors), edge vectors, colors, etc.  
}

Triangle[] sphereList;
Triangle[] rotatedList;

void setup() {
  ortho(-320, 320, 320, -320); // hard coded, 640x640 canvas, RHS
  resetMatrix();
  colorMode(RGB, 1.0f);

  sphereList = makeSphere(SPHERE_SIZE, 10);
  rotatedList = new Triangle[sphereList.length];
  announceSettings();
}

void settings() {
  size(640, 640, P3D); // hard coded 640x640 canvas
}

float theta = 0.0;
float delta = 0.01;
void draw() {
  clear();

  if (rotate)
  {
    theta += delta;
    while (theta > PI*2) theta -= PI*2;
  } 

  if (lineTest)
    lineTest();
  else
  {
    rotateSphere(sphereList, rotatedList, theta);
    drawSphere(rotatedList, lighting, shading);
  }
}

//////////////////////  MAIN PROGRAM
// creates a sphere made of triangles, centered on 0,0,0, with given radius
//
// also - 
// calculates the 3 edge vectors for each triangle
// calculates the face normal (unit length)
//
// HINT: first setup a loop to calculate all the points around the sphere,
//       store in an array then loop over those points and setup your triangles.
Triangle[] makeSphere(int radius, int divisions)
{
  //ALL variables
  int nPoints = divisions*divisions;
  int nTriang = (nPoints - divisions)*2;
  float[] xCoord = new float[nPoints];
  float[] yCoord = new float[nPoints];
  float[] zCoord = new float[nPoints];
  int coordPos = 0;
  Triangle[] returnTriang = new Triangle[nTriang];
  int triangPos = 0;
  float[] v1 = new float[3];
  float[] v2 = new float[3];
  float[] v3 = new float[3];
  float[] bottomPoint = new float[]{0, -1*radius, 0};
  float myPhi = 0, myTheta = 0;
  int phiCheck = 1;
  float phiIncrem = PI/divisions;
  float thetaIncrem = TWO_PI/divisions;
  
  //Calculation point coordinates
  for(int i = 0; i < divisions; i++){
    for(int j = 0; j < divisions; j++){
      myTheta = i*thetaIncrem;
      myPhi = j*phiIncrem;
      xCoord[coordPos] = radius*sin(myPhi)*sin(myTheta);
      yCoord[coordPos] = radius*cos(myPhi);
      zCoord[coordPos] = radius*sin(myPhi)*cos(myTheta);
      coordPos++;
    }
  }
  
  //Creating Triangles
  for(int k = 0; k < nPoints; k++){
    if(k == 10){
      v1[X] = xCoord[k%nPoints];
      v1[Y] = yCoord[k%nPoints];
      v1[Z] = zCoord[k%nPoints];
      v2 = Arrays.copyOf(bottomPoint, bottomPoint.length); 
      v3[X] = xCoord[(k+10)%nPoints];
      v3[Y] = yCoord[(k+10)%nPoints];
      v3[Z] = zCoord[(k+10)%nPoints];
      returnTriang[triangPos] = new Triangle(v1, v2, v3);
      triangPos++;
    }
    else{
      if(k != 1){ 
        v1[X] = xCoord[k%nPoints];
        v1[Y] = yCoord[k%nPoints];
        v1[Z] = zCoord[k%nPoints];
        v2[X] = xCoord[(k+1)%nPoints];
        v2[Y] = yCoord[(k+1)%nPoints];
        v2[Z] = zCoord[(k+1)%nPoints];
        v3[X] = xCoord[(k+10)%nPoints];
        v3[Y] = yCoord[(k+10)%nPoints];
        v3[Z] = zCoord[(k+10)%nPoints];
        returnTriang[triangPos] = new Triangle(v1, v2, v3);
        triangPos++;
      }
      v1[X] = xCoord[(k+11)%nPoints];
      v1[Y] = yCoord[(k+11)%nPoints];
      v1[Z] = zCoord[(k+11)%nPoints];
      v2[X] = xCoord[(k+10)%nPoints];
      v2[Y] = yCoord[(k+10)%nPoints];
      v2[Z] = zCoord[(k+10)%nPoints];
      v3[X] = xCoord[(k+1)%nPoints];
      v3[Y] = yCoord[(k+1)%nPoints];
      v3[Z] = zCoord[(k+1)%nPoints];
      returnTriang[triangPos] = new Triangle(v1, v2, v3);
      triangPos++;
    }
  }
  
  return returnTriang;
}

void tesselationPoints(int radius, int divisions){
  float phiIncrem = PI/divisions;
  float thetaIncrem = TWO_PI/divisions;
  float myX, myY;
  float myPhi = 0.0f, myTheta = 0.0f;
  
  for(int i = 0; i < 10; i++){//myTheta <= TWO_PI){
     for(int j = 0; j <= 10; j++){
       myPhi = j*phiIncrem;
       myTheta = i*thetaIncrem;
       myX = radius*sin(myPhi)*sin(myTheta);
       myY = radius*cos(myPhi);
       stroke(1,1,1);
       beginShape(POINTS);
         vertex(myX, myY);
       endShape();
     }
  }
  
  
}


// takes a new triangle, and calculates it's normals and edge vectors
Triangle setupTriangle(Triangle t)
{
  return t;
}

// This function draws the 2D, already projected triangle, on the raster
// - it culls degenerate or back-facing triangles
//
// - it calls fillTriangle to do the actual filling, and bresLine to
// make the triangle outline. 
//
// - implements the specified lighting model (using the global enum type)
// to calculate the vertex colors before calling fill triangle. Doesn't do shading
//
// - if needed, it draws the outline and normals (check global variables)
//
// HINT: write it first using the gl LINES/TRIANGLES calls, then replace
// those with your versions once it works.
void draw2DTriangle(Triangle t, Lighting lighting, Shading shading)
{
}

// uses a scanline algorithm to fill the 2D on-raster triangle
// - implements the specified shading algorithm to set color as specified
// in the global variable shading. Note that for NONE, this function simply
// returns without doing anything
// - uses POINTS to draw on the raster
void fillTriangle(Triangle t, Shading shading)
{
}

// given point p, normal n, eye location, light location, calculates phong
// - material represents ambient, diffuse, specular as defined at the top of the file
// - calculates the diffuse, specular, and multiplies it by the material and
// - fillcolor values
float[] phong(float[] p, float[] n, float[] eye, float[] light, 
  float[] material, float[] fillColor, float s)
{
  return new float[]{0, 0, 0};
}

// implements Bresenham's line algorithm
void bresLine(int fromX, int fromY, int toX, int toY) //<>//
{
  
  int myX = fromX, myY = fromY, sX ,sY;
  int dX = toX - fromX;
  int dY = toY - fromY;
  float slope, error = 0.5;
  
  if(dX == 0 && dY == 0){
    plotPoint(toX, toY);
  }
  else if(dY != 0 && abs(dX/dY) <= 1){ 
    
    if(dX != 0) sX = dX/abs(dX);
    else sX = 0;
    
    sY = dY/abs(dY);
    
    for(int i = 0; i < abs(dY); i++){
      slope = abs((float)dX/dY);
      myX = floor((fromX + i*slope*sX) + 0.5);
      myY = fromY + i*sY;
      plotPoint(myX,myY);
    }
  }
  else{ //(dX != 0 && abs(dY/dX) <= 1)
  
    if(dY != 0) sY = dY/abs(dY);
    else sY = 0;
    
    sX = dX/abs(dX);
    
    for(int i = 0; i < abs(dX); i++){
      slope = abs((float)dY/dX);
      myY = floor((fromY + i*slope*sY) + 0.5);
      myX = fromX + i*sX;
      plotPoint(myX,myY);
    }
  }
   //<>// //<>//
}

void plotPoint(int myX, int myY){
  beginShape(POINTS);
  vertex(myX, myY);
  endShape();
}
