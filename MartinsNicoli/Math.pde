/// BASIC Math functions
// populate as needed, and add others you may need. I only needed these.
// HHIINNT:: use test cases and a test function to make sure you don't have a mistake!!!!!!
//   - I spent like 3 hours because I had a typo in my cross product :(

// the "2D cross product", as in class
float cross2(float[] e1, float[] e2)
{
  return (e1[X]*e2[Y] - e1[Y]*e2[X]);
}

float[] cross3(float[] a, float[] b)
{
  return new float[]{a[Y]*b[Z] - a[Z]*b[Y], a[Z]*b[X] - a[X]*b[Z], a[X]*b[Y] - a[Y]*b[X]};
}

// normalize v to length 1 in place
void normalize(float[] v)
{
  float lengthV = sqrt(v[X]*v[X] + v[Y]*v[Y] + v[Z]*v[Z]);
  for(int i = 0; i < v.length; i++){
    if(lengthV != 0.0)
      v[i] = v[i]/lengthV;
    else v[i] = 0.0;
  }
}

float dot(float[] v1, float[] v2)
{ 
  return (v1[X]*v2[X] + v1[Y]*v2[Y] + v1[Z]*v2[Z]);
}

float[] subtract(float[] v1, float v2[])
{  
  float[] result = new float[min(v1.length, v2.length)];
  for(int i = 0; i < result.length; i++){
    result[i] = v1[i] - v2[i];
  }
  
  return result;
}

float[] multVS(float scalar, float[] vector){
  return new float[]{scalar*vector[X], scalar*vector[Y], scalar*vector[Z]};
}
