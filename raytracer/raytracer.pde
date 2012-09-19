
// spheres information
private float[][] spherePositions = new float[][] {{-1,0,-11},{.3,-.5,-7},{.7,-.5,-4}};
private float[] sphereRadii = new float[] {.5,.5,.5};
private float[][] sphereReflectColors = new float[][] {{.8,.8,.1},{.1,.8,.1},{.1,.8,.8}};
private int sphereInd = 0;
private int numSpheres = 3;

// light information
private float[] lightPos = new float[] {9,9,9};
private float[] lightColor = new float[] {.8,.8,.8};
private float[] ambientLightColor = new float[] {0,0,0};

// camera inputs
private float[] eye = new float[3];
private float[] gaze;// = new float[3];
private float[] tilt;// = new float[3];
private float viewAngle = 0.0;
private float n, f;
private float nx, ny;

// things we have to calculate
private float[] u;// = new float[3];
private float[] v;// = new float[3];
private float[] w;// = new float[3];
private float l, r, b, t;
private float[] s;// = new float[3];
private float[] rayd;// = new float[3];
private float[] surfaceNormal;
private float intersectT;
private float tempIntersectT;
private float[] intersectionPoint;
private float[] colorAtPoint;
private float[] vectorToLight;


void setup() {
  size(800,600);
  colorMode(RGB,1.0);
  background(.5);
  //testMatrix();
  loadCameraAndScreen();
  
  calculateUVW();
  calculateTBRL();
  
  loadPixels();
  
  for (int i = 0; i < ny; i++)
  {
    for (int j = 0; j < nx; j++)
    {
      calculateS(i, j);
      calculateRayDirection();
      
      // find the right sphere
      boolean hitASphere = false;
      intersectT = 10000;
      for (int k = 0; k < numSpheres; k++)
      {
        if(rayIntersectSphere(k) && (tempIntersectT < intersectT))
        {
          hitASphere = true;
          sphereInd = k;
          intersectT = tempIntersectT;
        }
      }
      
      // found right sphere (if found at all), figure out the rest
      if (hitASphere)
      {
        calculateIntersectionPoint();
        calculateNormal();
        calculateToLight();
        calculateColorAtPoint();
        pixels[((int)ny - i - 1) * width + j] = color(colorAtPoint[0],colorAtPoint[1],colorAtPoint[2]);
      }
    }
  }
  updatePixels();
  // rejoice!!!
}


void printFloat4Vector(float[] a) {
  String out = "Vector: [";
  for (int i = 0; i < 4; i++) {
    out += "(" + a[i] + ")";
  }
  out += "]";
  
  //System.out.println(out);
}

void loadCameraAndScreen() {
  eye = new float[] {0,0,0};
  gaze = new float[] {0,0,-1};
  tilt = new float[] {0,1,0};
  viewAngle = 30;
  n = 2;
  f = 200;
  nx = 800;
  ny = 600;
}

void calculateUVW() {
  // w = -g / ||g||
  float lengthG = (float)Math.sqrt(Math.pow(gaze[0],2) + Math.pow(gaze[1],2) + Math.pow(gaze[2],2));
  w = new float[] {-gaze[0]/lengthG, -gaze[1]/lengthG, -gaze[2]/lengthG};
  
  //float[] txw = new float[] {tilt[2]*w[3] - tilt[3]*w[2], tilt[3]*w[1] - tilt[1]*w[3], tilt[1]*w[2] - tilt[2]*w[1]};
  float[] txw = new float[] {tilt[1]*w[2] - tilt[2]*w[1], tilt[2]*w[0] - tilt[0]*w[2], tilt[0]*w[1] - tilt[1]*w[0]};
  float lengthTXW = (float)Math.sqrt(Math.pow(txw[0],2) + Math.pow(txw[1],2) + Math.pow(txw[2],2));
  
  // u = t x w / ||t x w||
  u = new float[] {txw[0]/lengthTXW, txw[1]/lengthTXW, txw[2]/lengthTXW};
  
  // u = w x u
  //v = new float[] {w[2]*u[3] - w[3]*u[2], w[3]*u[1] - w[1]*u[3], w[1]*u[2] - w[2]*u[1]};
  v = new float[] {w[1]*u[2] - w[2]*u[1], w[2]*u[0] - w[0]*u[2], w[0]*u[1] - w[1]*u[0]};
  
  /*
  System.out.println("u: [" + u[0] + "," + u[1] + "," + u[2] + "]");
  System.out.println("v: [" + v[0] + "," + v[1] + "," + v[2] + "]");
  System.out.println("w: [" + w[0] + "," + w[1] + "," + w[2] + "]");
  */
}

void calculateTBRL() {
  t = (float)(Math.abs(n) * Math.tan( (viewAngle * (Math.PI / 180)) / 2));
  b = -t;
  r = t * nx / ny;
  l = -r;
  
  //System.out.println("t: " + t + "\nb: " + b + "\nr: " + r + "\nl: " + l);
}

void calculateS(int i, int j) {  // i is Y, j is X
  //System.out.println("i:" + i + "j:" + j);
  
  s = new float[3];
  s[0] = l + (r-l) * ((float)j + 0.5) / nx;
  s[1] = b + (t-b) * ((float)i + 0.5) / ny;
  s[2] = Math.abs(n);
  
  //System.out.println("S");
  //printFloat3(s);
}

void calculateRayDirection()
{
  rayd = new float[] {(s[0] * u[0]) + (s[1] * v[0]) - (s[2] * w[0]),
                      (s[0] * u[1]) + (s[1] * v[1]) - (s[2] * w[1]),
                      (s[0] * u[2]) + (s[1] * v[2]) - (s[2] * w[2])};
  //System.out.println("rayd");
  //printFloat3(rayd);
  
  // normalize rayd
  float lengthRayd = (float)Math.sqrt(rayd[0] * rayd[0] + rayd[1] * rayd[1] + rayd[2] * rayd[2]);
  //System.out.println("lengthRayd:" + lengthRayd);
  
  rayd[0] = rayd[0] / lengthRayd;
  rayd[1] = rayd[1] / lengthRayd;
  rayd[2] = rayd[2] / lengthRayd;
  
  //System.out.println("rayd");
  //printFloat3(rayd);
  
}

boolean rayIntersectSphere(int sInd)
{
  float[] OC = new float[] {spherePositions[sInd][0] - eye[0],
                            spherePositions[sInd][1] - eye[1],
                            spherePositions[sInd][2] - eye[2]};
  //System.out.println("OC");
  //printFloat3(OC);
                            
  float tca = OC[0] * rayd[0] + OC[1] * rayd[1] + OC[2] * rayd[2];
  //System.out.println("tca:" + tca);
  if (tca < 0)
    return false;
    
  float d2 = (float)Math.pow(Math.sqrt(OC[0] * OC[0] + OC[1] * OC[1] + OC[2] * OC[2]),2) - (tca * tca);
  float thc2 = (sphereRadii[sInd] * sphereRadii[sInd]) - d2;
  //System.out.println("thc2:" + thc2);
  if (thc2 < 0)
    return false;
  
  float thc = (float)Math.sqrt(thc2);
  float plusT = tca + thc;
  float minusT = tca - thc;  
  
  if (plusT < minusT)
    tempIntersectT = plusT;
  else
    tempIntersectT = minusT;  
  
  return true;
}

void calculateIntersectionPoint()
{
  intersectionPoint = new float[] {eye[0] + rayd[0] * intersectT,
                                   eye[1] + rayd[1] * intersectT,
                                   eye[2] + rayd[2] * intersectT};
}

void calculateNormal()
{
  surfaceNormal = new float[] {(intersectionPoint[0] - spherePositions[sphereInd][0]) / sphereRadii[sphereInd],
                               (intersectionPoint[1] - spherePositions[sphereInd][1]) / sphereRadii[sphereInd],
                               (intersectionPoint[2] - spherePositions[sphereInd][2]) / sphereRadii[sphereInd]};
}

void calculateToLight()
{
  vectorToLight = new float[] {lightPos[0] - intersectionPoint[0],
                               lightPos[1] - intersectionPoint[1],
                               lightPos[2] - intersectionPoint[2]};
  // normalize
  
  float lengthVector = (float)Math.sqrt(vectorToLight[0] * vectorToLight[0] +
                                      vectorToLight[1] * vectorToLight[1] +
                                      vectorToLight[2] * vectorToLight[2]);
                                      
  vectorToLight[0] = vectorToLight[0] / lengthVector;
  vectorToLight[1] = vectorToLight[1] / lengthVector;
  vectorToLight[2] = vectorToLight[2] / lengthVector;
  
}

void calculateColorAtPoint()
{
  float dotProduct = surfaceNormal[0] * vectorToLight[0] +
                     surfaceNormal[1] * vectorToLight[1] +
                     surfaceNormal[2] * vectorToLight[2];
  float maxDot = Math.max(0,dotProduct);
  colorAtPoint = new float[] {sphereReflectColors[sphereInd][0] * (ambientLightColor[0] + lightColor[0] * maxDot),
                              sphereReflectColors[sphereInd][1] * (ambientLightColor[1] + lightColor[1] * maxDot),
                              sphereReflectColors[sphereInd][2] * (ambientLightColor[2] + lightColor[2] * maxDot)};
}


// print a float[3]
void printFloat3(float[] inF)
{
  String out = "float[3] = ";
  for (int i = 0; i < 3; i++)
  {
    out += "(" + inF[i] + ")";
  }
  System.out.println(out);
}

// here begin methods associated with my Matrix class
void testMatrix() {
  float[][] a = {{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
  float[][] b = {{1,2,0,0},{0,1,2,0},{0,0,1,2},{2,0,0,1}};
  float[][] c = {{1},{1},{2},{1}};
  
  Matrix.print4x4(a);
  Matrix.print4x4(b);
  Matrix.print4x1(c);
  
  float[][] temp = Matrix.multiply4x4by4x4(a, b);
  Matrix.print4x4(temp);
  
  temp = Matrix.multiply4x4by4x1(a, c);
  Matrix.print4x1(temp);
  
  temp = Matrix.multiply4x4by4x1(b, c);
  Matrix.print4x1(temp);
  
}

private static class Matrix {
  public static float[][] multiply4x4by4x4(float[][] a, float[][] b) {
    // new temp
    float[][] temp = new float[4][4];
    // multiply the matrix entries together for each spot in temp
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
	temp[i][j] = a[i][0] * b[0][j] + a[i][1] * b[1][j] + a[i][2] * b[2][j] + a[i][3] * b[3][j];
      }
    }

    return temp; 
  }
  
  public static float[][] multiply4x4by4x1(float[][] a, float[][] b) {
    // new temp
    float[][] temp = new float[4][1];
    // multiply the matrix entries together for each spot in temp
    for (int i = 0; i < 4; i++) {
	temp[i][0] = a[i][0] * b[0][0] + a[i][1] * b[1][0] + a[i][2] * b[2][0] + a[i][3] * b[3][0];
    }

    return temp;
  }
  
  public static void print4x4(float[][] a) {
    String out = "4x4:";    
    for (int i = 0; i < 4; i++) {
      out += "\n[";
      for (int j = 0; j < 4; j++) {
	out += "(" + a[i][j] + ")";
      }
      out += "]";
    }
    out += "\n";
    
    System.out.println(out);
  }
  
  public static void print4x1(float[][] a) {
    String out = "4x1:";    
    for (int i = 0; i < 4; i++) {
      out += "\n[";
      out += "(" + a[i][0] + ")";
      out += "]";
    }
    out += "\n";
    
    System.out.println(out);    
  }
  
   
}
