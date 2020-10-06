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
  return new Triangle[0];
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
  
  if(dX == 0 && dY == 0){ //<>//
    plotPoint(toX, toY);
  }
  else if(dY != 0 && abs(dX/dY) <= 1){ 
    
    if(dX != 0) sX = dX/abs(dX);
    else sX = 0;
    
    sY = dY/abs(dY);
    /*slope = abs((float)dX/dY);
    
    while(myX < toX){
      plotPoint(myX,myY);
      myY += sY; 
      if(error >= 1){ 
        myX += sX;
        error = error -1;
      }*/
    
    for(int i = 0; i < abs(dY); i++){
      slope = abs((float)dX/dY);
      //myX = floor((fromX + i*slope*sX) + 0.5);
      error += slope;
      if(error >= 0.5){
        myX = myX + sX;
        error = error - 1;
      }
      myY = fromY + i*sY;
      plotPoint(myX,myY);
    }
  }
  else{ //(dX != 0 && abs(dY/dX) <= 1) //<>//
  
    if(dY != 0) sY = dY/abs(dY);
    else sY = 0;
    
    sX = dX/abs(dX);
    
    for(int i = 0; i < abs(dX); i++){
      slope = abs((float)dY/dX);
      //myY = floor((fromY + i*slope*sY) + 0.5);
      error += slope;
      if(error >= 0.5){
        myY = myY + sY;
        error = error - 1;
      }        
      myX = fromX + i*sX;
      plotPoint(myX,myY);
    }
  }
  
}

void plotPoint(int myX, int myY){
  beginShape(POINTS);
  vertex(myX, myY);
  endShape();
}
