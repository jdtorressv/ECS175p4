#include <GL/glut.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
#include <cmath>
using namespace std;

#define ERROR 1 
#define NORM 600.0 

/*****************************************************************************/
/* Classes                            		                             */
/*****************************************************************************/
 
class Vertex {
	public:
		double x;   
		double y;  
		double z; 
	        double normX;
		double normY; 
		double normZ;
		double red;
		double green;
		double blue; 
		Vertex(double x, double y, double z);
};

class Triangle {
	public: 
		int vertexA;
		int vertexB; 
		int vertexC;
		Triangle(int vertexA, int vertexB, int vertexC); 
};

class Normal {
        public:
                int VID;
                double xVal;
                double yVal;
                double zVal;
		Normal(int VID, double xVal, double yVal, double zVal); 
};

class Polyhedron {
	public: 
		vector<Vertex> vArr; 
		vector<Triangle> triArr; 
		vector<Normal> normArr; 
};	
class CoeffSet {
	public:
		float ka;
		float kd; 
		float ks;
	        CoeffSet(float ka, float kd, float ks); 	
};

/*****************************************************************************/
/* Globals                                                                   */
/*****************************************************************************/

int windowID, windowXY, windowXZ, windowYZ;
vector<Polyhedron> polyArr; 
float *PixelBufferXY, *PixelBufferXZ, *PixelBufferYZ;
float *PolygonBufferXY, *PolygonBufferXZ, *PolygonBufferYZ;

/*****************************************************************************/
/* Function Definitions                                                      */
/*****************************************************************************/

//Constructors
Vertex::Vertex(double x, double y, double z) : x(x), y(y), z(z) {}

Triangle::Triangle(int vertexA, int vertexB, int vertexC) : vertexA(vertexA), vertexB(vertexB), vertexC(vertexC) {}

Normal::Normal(int VID, double xVal, double yVal, double zVal) : VID(VID), xVal(xVal), yVal(yVal), zVal(zVal) {}

CoeffSet::CoeffSet(float ka, float kd, float ks) : ka(ka), kd(kd), ks(ks) {}

//Fill appropriate classes with input file specs
void populatePolyhedronInfo(vector<double> v)
{
        auto vpoint = v.begin();

        int polyTotal = (int)*(vpoint);

        for (int i = 0; i < polyTotal; i++) {
                Polyhedron poly;
                polyArr.push_back(poly);
                int vertices = (int)*(++vpoint);

                for (int j = 0; j < vertices; j++) {
                        double x = *(++vpoint)/NORM;
                        double y = *(++vpoint)/NORM;
                        double z = *(++vpoint)/NORM;
                        Vertex vertex(x, y, z);
                        polyArr.at(i).vArr.push_back(vertex);
                }

                int numTri = (int)*(++vpoint);

                for (int k = 0; k < numTri; k++) {
                        int first = (int)*(++vpoint);
                        int second = (int)*(++vpoint);
                        int third = (int)*(++vpoint);
                        Triangle triangle(first, second, third);
                        polyArr.at(i).triArr.push_back(triangle);
                }
        }
}
void swap(int &first, int &second, int &third) 
{
	int temp = first; 
	first = second;
        second = third; 
	third = temp; 	

}
void crossProduct(int index, int A, int B, int C)
{
			double Ux = polyArr.at(index).vArr.at(B).x - polyArr.at(index).vArr.at(A).x;
                        double Uy = polyArr.at(index).vArr.at(B).y - polyArr.at(index).vArr.at(A).y;
                        double Uz = polyArr.at(index).vArr.at(B).z - polyArr.at(index).vArr.at(A).z;

                        double Vx = polyArr.at(index).vArr.at(C).x - polyArr.at(index).vArr.at(B).x;
                        double Vy = polyArr.at(index).vArr.at(C).y - polyArr.at(index).vArr.at(B).y;
                        double Vz = polyArr.at(index).vArr.at(C).z - polyArr.at(index).vArr.at(B).z;

                        Vertex U(Ux, Uy, Uz);
                        Vertex V(Vx, Vy, Vz);
			
			double Nx = U.y*V.z - U.z*V.y;
			double Ny = U.z*V.x - U.x*V.z;
			double Nz = U.x*V.y - U.y*V.x;
			double mag = sqrt(pow(Nx, 2) + pow(Ny, 2) + pow(Nz, 2)); 

                        Normal surfNormal(B, Nx/mag, Ny/mag, Nz/mag);
                        polyArr.at(index).normArr.push_back(surfNormal);

}
void calculateNormals() 
{
	for (int i = 0; i < polyArr.size(); i++) {
		for (int j = 0; j < polyArr.at(i).triArr.size(); j++) {
			Triangle currTri = polyArr.at(i).triArr.at(j); 

			int A = currTri.vertexA;
			int B = currTri.vertexB; 
			int C = currTri.vertexC;
			
			crossProduct(i, A, B, C);  

			swap(A, B, C);

			crossProduct(i, A, B, C); 

			swap(A, B, C);
	
			crossProduct(i, A, B, C);
		}
	}	
}
void setNormals()
{
        for (int i = 0; i < polyArr.size(); i++) {
                for (int j = 0; j < polyArr.at(i).vArr.size(); j++) {
                        int count = 0;
                        for (int k = 0; k < polyArr.at(i).normArr.size(); k++) {
                                if (j == polyArr.at(i).normArr.at(k).VID) {
                                        polyArr.at(i).vArr.at(j).normX += polyArr.at(i).normArr.at(k).xVal;
                                        polyArr.at(i).vArr.at(j).normY += polyArr.at(i).normArr.at(k).yVal;
                                        polyArr.at(i).vArr.at(j).normY += polyArr.at(i).normArr.at(k).zVal;
                                        count++;
                                }
                        }
                        polyArr.at(i).vArr.at(j).normX /= (double)count; 
                        polyArr.at(i).vArr.at(j).normY /= (double)count;
                        polyArr.at(i).vArr.at(j).normZ /= (double)count; 
                }
        }
}
float ComputeIP(CoeffSet coeff, float Ia, float Il, float Vmag, float K, Vertex vecL, float Nx, float Ny, float Nz, Vertex vecR, Vertex vecV, int n)
{
	float inten = coeff.ka*Ia + (Il/(Vmag + K))*(coeff.kd*(vecL.x*Nx + vecL.y*Ny + vecL.z*Nz) + coeff.ks*pow((vecR.x*vecV.x + vecR.y*vecV.y + vecR.z*vecV.z), n));
	if (inten < 0)
		return 0;
	else
		return inten; 

}
void setIntensities(Vertex lightPt, Vertex viewPt, CoeffSet redC, CoeffSet greenC, CoeffSet blueC, float K, float Ia, float Il, int n)
{
	for (int i = 0; i < polyArr.size(); i++) {
		for (int j = 0; j < polyArr.at(i).vArr.size(); j++) {
			float x = polyArr.at(i).vArr.at(j).x;
			float y = polyArr.at(i).vArr.at(j).y;
			float z = polyArr.at(i).vArr.at(j).z; 
			
			//Calculate vectors: reflection r, light l, and viewing v 
			float Lx = x - lightPt.x; 
			float Ly = y - lightPt.y; 
			float Lz = z - lightPt.z; 
			float Lmag = sqrt(pow(Lx, 2) + pow(Ly, 2) + pow(Lz, 2)); 
			Vertex vecL(Lx/Lmag, Ly/Lmag, Lz/Lmag); 

			float Vx = viewPt.x - lightPt.x;
			float Vy = viewPt.y - lightPt.y;
			float Vz = viewPt.z - lightPt.z; 
			float Vmag = sqrt(pow(Vx, 2) + pow(Vy, 2) + pow(Vz, 2));
			Vertex vecV(Vx/Vmag, Vy/Vmag, Vz/Vmag); 

			float Nx = polyArr.at(i).vArr.at(j).normX; 
                        float Ny = polyArr.at(i).vArr.at(j).normY;
                        float Nz = polyArr.at(i).vArr.at(j).normZ;
			float NdotL = Nx*vecL.x + Ny*vecL.y + Nz*vecL.z; 
			Vertex vecR(2*NdotL*Nx - vecL.x, 2*NdotL*Ny - vecL.y, 2*NdotL*Nz - vecL.z);

			polyArr.at(i).vArr.at(j).red = ComputeIP(redC, Ia, Il, Vmag, K, vecL, Nx, Ny, Nz, vecR, vecV, n); 
		        polyArr.at(i).vArr.at(j).green = ComputeIP(greenC, Ia, Il, Vmag, K, vecL, Nx, Ny, Nz, vecR, vecV, n);
			polyArr.at(i).vArr.at(j).blue = ComputeIP(blueC, Ia, Il, Vmag, K, vecL, Nx, Ny, Nz, vecR, vecV, n);
			
			cout << "Polyhedron #" << i << ", vertex #" << j << ": (" << polyArr.at(i).vArr.at(j).red << ", " << polyArr.at(i).vArr.at(j).green << ", " << polyArr.at(i).vArr.at(j).blue << ")\n"; 	
		}
	}				
}
inline int roundOff(const double a) {return (int)(a+0.5);}
void makePix(int x, int y, int pid)
{

}
void copyBuffer(int pid)
{

}
void clearPixelBuffer()
{

}
void clearPolygonBuffer(int pid)
{

}
void lineDrawRaster()
{

}

//Inline function for mainMenu delegates all functionality to sub menus 
inline void mainMenu(int pid) {;}

//Initializes each of the subwindoes as well as the primary window 
void init()
{	
	if (glutGetWindow() == windowID) 
		glClearColor(1.0, 1.0, 1.0, 0.0); //Set color to white
	else 
        	glClearColor(0.0, 0.0, 0.0, 0.0); //Set color to black

        glMatrixMode(GL_PROJECTION);
}

//XY: All Z values are ignored 
void drawSceneXY()
{
	glClear(GL_COLOR_BUFFER_BIT); 
        glLoadIdentity();
	glBegin(GL_LINES);
		
		for (int i = 0; i < polyArr.size(); i++) { //For each polyhedron
			for (int j = 0; j < polyArr.at(i).triArr.size(); j++) { //For each each triangle in the polyhedron 
				int A = polyArr.at(i).triArr.at(j).vertexA;
				int B = polyArr.at(i).triArr.at(j).vertexB; 
				int C = polyArr.at(i).triArr.at(j).vertexC;
			
				
				float x1 = polyArr.at(i).vArr.at(A).x;
				float y1 = polyArr.at(i).vArr.at(A).y;

                                float x2 = polyArr.at(i).vArr.at(B).x;
                                float y2 = polyArr.at(i).vArr.at(B).y;

				float x3 = polyArr.at(i).vArr.at(C).x;
                                float y3 = polyArr.at(i).vArr.at(C).y;


				glColor3f(1.0, 1.0, 1.0); 
				glVertex2f(x1, y1);
				glVertex2f(x2, y2); 
                                glVertex2f(x2, y2);
                                glVertex2f(x3, y3);
                                glVertex2f(x3, y3);
                                glVertex2f(x1, y1);

			}
		}

	glEnd(); 
	glFlush(); 
}
//XZ: All Y values are ignored
void drawSceneXZ()
{
        glClear(GL_COLOR_BUFFER_BIT);
        glLoadIdentity();
	glBegin(GL_LINES);
		 
                for (int i = 0; i < polyArr.size(); i++) { //For each polyhedron
                        for (int j = 0; j < polyArr.at(i).triArr.size(); j++) { //For each each triangle in the polyhedron 
                                int A = polyArr.at(i).triArr.at(j).vertexA;
                                int B = polyArr.at(i).triArr.at(j).vertexB;
                                int C = polyArr.at(i).triArr.at(j).vertexC;


                                float x1 = polyArr.at(i).vArr.at(A).x;
                                float z1 = polyArr.at(i).vArr.at(A).z;

                                float x2 = polyArr.at(i).vArr.at(B).x;
                                float z2 = polyArr.at(i).vArr.at(B).z;

                                float x3 = polyArr.at(i).vArr.at(C).x;
                                float z3 = polyArr.at(i).vArr.at(C).z;


                                glColor3f(1.0, 1.0, 1.0);
                                glVertex2f(x1, z1);
                                glVertex2f(x2, z2);
                                glVertex2f(x2, z2);
                                glVertex2f(x3, z3);
                                glVertex2f(x3, z3);
                                glVertex2f(x1, z1);

                        }
                }


	glEnd(); 
	glFlush(); 
}
//XZ: All X values are ignored
void drawSceneYZ()
{
        glClear(GL_COLOR_BUFFER_BIT);
        glLoadIdentity();
	glBegin(GL_LINES);

                for (int i = 0; i < polyArr.size(); i++) { //For each polyhedron
                        for (int j = 0; j < polyArr.at(i).triArr.size(); j++) { //For each each triangle in the polyhedron 
                                int A = polyArr.at(i).triArr.at(j).vertexA;
                                int B = polyArr.at(i).triArr.at(j).vertexB;
                                int C = polyArr.at(i).triArr.at(j).vertexC;


                                float y1 = polyArr.at(i).vArr.at(A).y;
                                float z1 = polyArr.at(i).vArr.at(A).z;

                                float y2 = polyArr.at(i).vArr.at(B).y;
                                float z2 = polyArr.at(i).vArr.at(B).z;

                                float y3 = polyArr.at(i).vArr.at(C).y;
                                float z3 = polyArr.at(i).vArr.at(C).z;


                                glColor3f(1.0, 1.0, 1.0);
                                glVertex2f(y1, z1);
                                glVertex2f(y2, z2);
                                glVertex2f(y2, z2);
                                glVertex2f(y3, z3);
                                glVertex2f(y3, z3);
                                glVertex2f(y1, z1);

                        }
                }


	glEnd(); 
	glFlush(); 
}
//Set up main display 
void background()
{
        glClear(GL_COLOR_BUFFER_BIT);
        glLoadIdentity();
        glFlush();
}
//Master display function
void display()
{	
	glutSetWindow(windowID);
	background(); 
        glutSetWindow(windowXY);
        drawSceneXY();
        glutSetWindow(windowXZ);
        drawSceneXZ();
        glutSetWindow(windowYZ);
        drawSceneYZ();
}
//Provide current vertex information for specified polyhedron 
void vertexMenu(int pid)
{
	;
}
//Translation 
void translateMenu(int pid)
{
	glutPostRedisplay(); 
 	;	
}
//Scaling
void scaleMenu(int pid) 
{

        //glutPostRedisplay();
	;
}
void rotateMenu(int pid)
{
	
        //glutPostRedisplay();
	;

}
/*****************************************************************************/
/* Main                                                                      */
/*****************************************************************************/

int main(int argc, char** argv) 
{
	if (argc != 20) {
		cout << "Usage: p3 <Px Py Pz Fx Fy Fz kaR kdR ksR kaG kdG ksG kaB kdB ksB K Ia Il n> \n";
		exit(ERROR); 
	}	

	Vertex lightPt(atof(argv[1])/NORM, atof(argv[2])/NORM, atof(argv[3])/NORM);
	Vertex viewPt(atof(argv[4])/NORM, atof(argv[5])/NORM, atof(argv[6])/NORM);
	CoeffSet redC(atof(argv[7]), atof(argv[8]), atof(argv[9]));
        CoeffSet greenC(atof(argv[10]), atof(argv[11]), atof(argv[12]));
	CoeffSet blueC(atof(argv[13]), atof(argv[14]), atof(argv[15]));
	float K = atof(argv[16])/NORM; 
	float Ia = atof(argv[17]); 
	float Il = atof(argv[18]); 
	int n = atoi(argv[19]); 


	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
	glutInitWindowPosition(100, 100); 
	glutInitWindowSize(800, 800); 

	windowID = glutCreateWindow("Polyhedron Orthographic Projections: XY, XZ, YZ from left to right and top down");  
	init();

	int scale_menu, rotate_menu, translate_menu, vertex_menu; //For use in graphical menu

        vector<double> v;
        double num;
        fstream file;
        file.open("inputFile.txt");

        if (!file) {
                cerr << "Unable to open file\n";
                exit(ERROR);
        }
        while (file >> num)
                v.push_back(num); //Initial vector for all polygons
        file.close();
	
	populatePolyhedronInfo(v); 
	calculateNormals(); 
	setNormals(); 
	setIntensities(lightPt, viewPt, redC, greenC, blueC, K, Ia, Il, n); // Phong 

	
	//XY
	windowXY = glutCreateSubWindow(windowID, 25, 50, 320, 320);
	init();

	//XZ
	windowXZ = glutCreateSubWindow(windowID, 25, 450, 320, 320); 
	init(); 

	//YZ
	windowYZ = glutCreateSubWindow(windowID, 425, 450, 320, 320); 
	init();

	
	glutSetWindow(windowID); 

        // Offer the user opportunities to 3D transform! 
        translate_menu = glutCreateMenu(translateMenu);
                glutAddMenuEntry("Octahedron/Polyhedron 0", 0);
                glutAddMenuEntry("Tetrahedron/Polyhedron 1", 1);
                glutAddMenuEntry("Hexahedron/Polyhedron 2", 2);

        scale_menu = glutCreateMenu(scaleMenu);
                glutAddMenuEntry("Octahedron/Polyhedron 0", 0);
                glutAddMenuEntry("Tetrahedron/Polyhedron 1", 1);
                glutAddMenuEntry("Hexahedron/Polyhedron 2", 2);

        rotate_menu = glutCreateMenu(rotateMenu);
                glutAddMenuEntry("Octahedron/Polyhedron 0", 0);      
      		glutAddMenuEntry("Tetrahedron/Polyhedron 1", 1);
                glutAddMenuEntry("Hexahedron/Polyhedron 2", 2);
	
	vertex_menu = glutCreateMenu(vertexMenu);
                glutAddMenuEntry("Octahedron/Polyhedron 0", 0);
                glutAddMenuEntry("Tetrahedron/Polyhedron 1", 1);
                glutAddMenuEntry("Hexahedron/Polyhedron 2", 2);

        glutCreateMenu(mainMenu);
                glutAddSubMenu("Translate", translate_menu);
                glutAddSubMenu("Scale", scale_menu);
                glutAddSubMenu("Rotate", rotate_menu);
		glutAddSubMenu("Vertex Dump", vertex_menu); 
        glutAttachMenu(GLUT_RIGHT_BUTTON);


        glutDisplayFunc(display);


	glutMainLoop(); 

	return 0; 
}
