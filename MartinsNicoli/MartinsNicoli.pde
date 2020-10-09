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
  float[] centroid;

  // projected data. On the screen raster
  float[] pv1; // (p)rojected vertices
  float[] pv2;
  float[] pv3;

  // add other things as needed, like normals (face, vectors), edge vectors, colors, etc.  
  float[] normal;
  float[] e1;
  float[] e2;
  float[] e3;
  
  float pnormal;
  float[] pe1;
  float[] pe2;
  float[] pe3;
  
  float[] normalV1, normalV2, normalV3;
  
  //Color value
  float[] colorV1, colorV2, colorV3;
}

Triangle[] sphereList;
Triangle[] rotatedList;

void setup() {
  ortho(-320, 320, 320, -320); // hard coded, 640x640 canvas, RHS
  resetMatrix();
  colorMode(RGB, 1.0f);

  sphereList = makeSphere(SPHERE_SIZE, 10);//10
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
  int nTriang = (divisions*(divisions-2))*2 + 2*divisions;
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
    if( (k-(divisions-1))%divisions == 0){
      v1[X] = xCoord[k%nPoints];
      v1[Y] = yCoord[k%nPoints];
      v1[Z] = zCoord[k%nPoints];
      v2 = Arrays.copyOf(bottomPoint, bottomPoint.length); 
      v3[X] = xCoord[(k+divisions)%nPoints];
      v3[Y] = yCoord[(k+divisions)%nPoints];
      v3[Z] = zCoord[(k+divisions)%nPoints];
      returnTriang[triangPos] = new Triangle(v1, v2, v3);
      setupTriangle(returnTriang[triangPos]);
      triangPos++;
    }
    else{
      if(k%divisions != 0){ 
        v1[X] = xCoord[k%nPoints];
        v1[Y] = yCoord[k%nPoints];
        v1[Z] = zCoord[k%nPoints];
        v2[X] = xCoord[(k+1)%nPoints];
        v2[Y] = yCoord[(k+1)%nPoints];
        v2[Z] = zCoord[(k+1)%nPoints];
        v3[X] = xCoord[(k+divisions+1)%nPoints];
        v3[Y] = yCoord[(k+divisions+1)%nPoints];
        v3[Z] = zCoord[(k+divisions+1)%nPoints];
        returnTriang[triangPos] = new Triangle(v1, v2, v3);
        setupTriangle(returnTriang[triangPos]);
        triangPos++;
      }
      v1[X] = xCoord[k%nPoints];
      v1[Y] = yCoord[k%nPoints];
      v1[Z] = zCoord[k%nPoints];
      v2[X] = xCoord[(k+divisions+1)%nPoints];
      v2[Y] = yCoord[(k+divisions+1)%nPoints];
      v2[Z] = zCoord[(k+divisions+1)%nPoints];
      v3[X] = xCoord[(k+divisions)%nPoints];
      v3[Y] = yCoord[(k+divisions)%nPoints];
      v3[Z] = zCoord[(k+divisions)%nPoints]; 
      returnTriang[triangPos] = new Triangle(v1, v2, v3);
      setupTriangle(returnTriang[triangPos]);
      triangPos++;
    }
  }
  return returnTriang;
}
  

// takes a new triangle, and calculates it's normals and edge vectors
Triangle setupTriangle(Triangle t)
{
  //edge vectors
  t.e1 = subtract(t.v2, t.v1); //v2 - v1
  t.e2 = subtract(t.v3, t.v2); //v3 - v2
  t.e3 = subtract(t.v1, t.v3); //v1 - v3
  
  t.normal = cross3(t.e1, t.e2);
  normalize(t.normal);
  
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
  
  //projection
  
  //"projected" edge vectors
  t.pe1 = new float[]{t.pv2[X] - t.pv1[X], t.pv2[Y] - t.pv1[Y]};//subtract(t.pv2, t.pv1); //v2 - v1
  t.pe2 = new float[]{t.pv3[X] - t.pv2[X], t.pv3[Y] - t.pv2[Y]}; //v3 - v2
  t.pe3 = new float[]{t.pv1[X] - t.pv3[X], t.pv1[Y] - t.pv3[Y]}; //v1 - v3
  
  t.pnormal = cross2(t.pe1, t.pe2);
  
  
  if(t.pnormal > 0){//dot(t.normal, new float[]{0,0,-1}) < 0){
    //println(t.pnormal);
    
    //define lighting color scheme
    
    lightColor(t, lighting);
    
    fillTriangle(t, shading);
    
    if(doOutline){
      stroke(OUTLINE_COLOR[R], OUTLINE_COLOR[G], OUTLINE_COLOR[B]);
      
      bresLine((int)t.pv1[X], (int)t.pv1[Y], (int)t.pv2[X], (int)t.pv2[Y]);
      bresLine((int)t.pv2[X], (int)t.pv2[Y], (int)t.pv3[X], (int)t.pv3[Y]);
      bresLine((int)t.pv3[X], (int)t.pv3[Y], (int)t.pv1[X], (int)t.pv1[Y]);
      
    }
  }
}

// uses a scanline algorithm to fill the 2D on-raster triangle
// - implements the specified shading algorithm to set color as specified
// in the global variable shading. Note that for NONE, this function simply
// returns without doing anything
// - uses POINTS to draw on the raster
void fillTriangle(Triangle t, Shading shading)
{
  if(shading != Shading.NONE){
    int xmin = (int)min(t.pv1[X], t.pv2[X], t.pv3[X]);
    int xmax = (int)max(t.pv1[X], t.pv2[X], t.pv3[X]);
    int ymin = (int)min(t.pv1[Y], t.pv2[Y], t.pv3[Y]);
    int ymax = (int)max(t.pv1[Y], t.pv2[Y], t.pv3[Y]);  
    float crossTri = abs(cross2(t.pe1,t.pe2)); 
    //area of the paralelograms from edge and bary vectors
    float a1, a2, a3;
    //Vectors of Barycentric Coordinates
    float[] baryE1, baryE2, baryE3;
    float[] avgColor;
        
    for(int i = ymin; i <= ymax; i++){
      for(int j = xmin; j <= xmax; j++){
        baryE1 = subtract(new float[]{j,i}, t.pv1);
        baryE2 = subtract(new float[]{j,i}, t.pv2);
        baryE3 = subtract(new float[]{j,i}, t.pv3);
        a1 = cross2(t.pe1, baryE1);
        a2 = cross2(t.pe2, baryE2);
        a3 = cross2(t.pe3, baryE3);
        if((a1 >= 0 && a2 >= 0 && a3 >= 0) || (a1 <= 0 && a2 <= 0 && a3 <= 0)){
          
          if(shading == Shading.FLAT){
              avgColor = new float[]{(t.colorV1[R] + t.colorV2[R] + t.colorV3[R])/3, 
                        (t.colorV1[R] + t.colorV2[R] + t.colorV3[R])/3,
                        (t.colorV1[R] + t.colorV2[R] + t.colorV3[R])/3};
              stroke(avgColor[R], avgColor[G], avgColor[B]);
          }
                              
            else if(shading == Shading.BARYCENTRIC)
              stroke((a1/crossTri), (a2/crossTri), (a3/crossTri));
              
            else if(shading == Shading.GOURAUD){
              stroke((t.colorV1[R]*a1+t.colorV2[R]*a2+t.colorV3[R]*a3)/crossTri, 
                        (t.colorV1[G]*a1+t.colorV2[G]*a2+t.colorV3[G]*a3)/crossTri, 
                        (t.colorV1[B]*a1+t.colorV2[B]*a2+t.colorV3[B]*a3)/crossTri); //c = v1*u+v2*v+v3*w
            }
            else{ //shading == Shading.PHONG
            stroke(0,0,0);
            }
              
              
            beginShape(POINTS);
            vertex(j, i);
            endShape();
            
        }
      }
    }
  }
  else stroke(0,0,0); //stateful
  
}

// given point p, normal n, eye location, light location, calculates phong
// - material represents ambient, diffuse, specular as defined at the top of the file
// - calculates the diffuse, specular, and multiplies it by the material and
// - fillcolor values
float[] phong(float[] p, float[] n, float[] eye, float[] light, 
  float[] material, float[] fillColor, float s)
{
  float[] viewV = subtract(eye, p);
  float[] lightV = subtract(light, p);
  float[] result = new float[3];
  normalize(viewV);
  normalize(lightV);
  float[] refV = subtract(lightV,multVS(2*dot(n, lightV),n));
  normalize(refV);
  float diffused = dot(lightV, n)*material[M_DIFFUSE];
  float specular = material[M_SPECULAR]*pow(dot(refV, viewV), s);
  float lightPhong = material[M_AMBIENT] + diffused + material[M_SPECULAR]*specular;
  result[R] = lightPhong*fillColor[R] ;
  result[G] = lightPhong*fillColor[G];
  result[B] = lightPhong*fillColor[B];
  
  return result;
}

// implements Bresenham's line algorithm
void bresLine(int fromX, int fromY, int toX, int toY)
{
  
  int myX = fromX, myY = fromY, sX ,sY;
  int dX = toX - fromX;
  int dY = toY - fromY;
  float slope;
  
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
  
}

void plotPoint(int myX, int myY){
  beginShape(POINTS);
  vertex(myX, myY);
  endShape();
}

void lightColor(Triangle t, Lighting lighting){
  if(lighting == Lighting.FLAT)
    t.colorV1 = t.colorV2 = t.colorV3 = FILL_COLOR;
  else if(lighting == Lighting.PHONG_FACE){
    float[] point = new float[]{(t.v1[X]+t.v2[X]+t.v3[X])/3,(t.v1[Y]+t.v2[Y]+t.v3[Y])/3,(t.v1[Z]+t.v2[Z]+t.v3[Z])/3};
    t.colorV1 = t.colorV2 = t.colorV3 = phong(point, t.normal, EYE, LIGHT, MATERIAL, FILL_COLOR, PHONG_SPECULAR);
  }
  else{ //if(lighting == Lighting.PHONG_VERTEX){
    //calculate normalV1,2,3
    t.normalV1 = subtract(t.v1, new float[]{0,0,0}); //<>//
    t.normalV2 = Arrays.copyOf(t.v2, t.v2.length);
    t.normalV3 = Arrays.copyOf(t.v3, t.v3.length);
    normalize(t.normalV1);
    normalize(t.normalV2);
    normalize(t.normalV3);
    
    //calculate colors
    t.colorV1 =  phong(t.v1, t.normalV1, EYE, LIGHT, MATERIAL, FILL_COLOR, PHONG_SPECULAR);
    t.colorV2 =  phong(t.v2, t.normalV2, EYE, LIGHT, MATERIAL, FILL_COLOR, PHONG_SPECULAR);
    t.colorV3 =  phong(t.v3, t.normalV3, EYE, LIGHT, MATERIAL, FILL_COLOR, PHONG_SPECULAR);
  }
  
}
