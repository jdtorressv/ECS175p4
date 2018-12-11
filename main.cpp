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
#define NORM 500.0

/*****************************************************************************/
/* Classes                            		                                   */
/*****************************************************************************/
class Point {
	public:
		float x;
		float y;
	  Point(float x, float y);
};
class BezSpline {
	public:
		 vector<Point> pointArr;
};
class BSpline {
	public:
		vector<Point> pointArr;
		int k;
};

/*****************************************************************************/
/* Globals                                                                   */
/*****************************************************************************/

//global variables
BezSpline bez1, bez2;
BSpline b1, b2;
int windowID, windowBez1, windowBez2, windowB1, windowB2;

/*****************************************************************************/
/* Function Definitions                                                      */
/*****************************************************************************/

//Constructors
Point::Point(float x, float y) : x(x), y(y) {}

void populateCurves(vector<float> v)
{
		auto vpoint = v.begin();
		int bezTotal = (int)*vpoint;
		if (bezTotal != 2) {
				cout << "This program only supports two Bezier curve instances\n";
				exit(ERROR);
		}
		int cntrlTotal = (int)*(++vpoint);
		for (int i = 0; i < cntrlTotal; i++) {
				Point pt(*(++vpoint) / NORM, *(++vpoint) / NORM);
				bez1.pointArr.push_back(pt);
		}
		cntrlTotal = (int)*(++vpoint);
		for (int i = 0; i < cntrlTotal; i++) {
				Point pt(*(++vpoint)/ NORM, *(++vpoint) / NORM);
				bez2.pointArr.push_back(pt);
		}
		int bTotal = (int)*(++vpoint);
		if (bezTotal != 2) {
				cout << "This program only supports two B-Spline curve instances\n";
				exit(ERROR);
		}
		b1.k = (int)*(++vpoint);
		cntrlTotal = (int)*(++vpoint);
		for (int i = 0; i < cntrlTotal; i++) {
				Point pt(*(++vpoint)/ NORM, *(++vpoint) / NORM);
				b1.pointArr.push_back(pt);
		}
		b2.k = (int)*(++vpoint);
		cntrlTotal = (int)*(++vpoint);
		for (int i = 0; i < cntrlTotal; i++) {
				Point pt(*(++vpoint)/ NORM, *(++vpoint) / NORM);
				b2.pointArr.push_back(pt);
		}

		cout << "Bez1 has points:\n";
		for (int i = 0; i < bez1.pointArr.size(); i++)
				cout << "(" << bez1.pointArr.at(i).x << ", " << bez1.pointArr.at(i).y << ")\n";
		cout << "Bez2 has points:\n";
		for (int i = 0; i < bez2.pointArr.size(); i++)
				cout << "(" << bez2.pointArr.at(i).x << ", " << bez2.pointArr.at(i).y << ")\n";
		cout << "B1 has k: " << b1.k << " and points:\n";
		for (int i = 0; i < b1.pointArr.size(); i++)
				cout << "(" << b1.pointArr.at(i).x << ", " << b1.pointArr.at(i).y << ")\n";
		cout << "B2 has k: " << b2.k << " and points:\n";
		for (int i = 0; i < b2.pointArr.size(); i++)
				cout << "(" << b2.pointArr.at(i).x << ", " << b2.pointArr.at(i).y << ")\n";
}

//Initializes each of the subwindows as well as the primary window
void init()
{
	if (glutGetWindow() == windowID)
		glClearColor(1.0, 1.0, 1.0, 0.0); //Set color to white
	else
        	glClearColor(0.0, 0.0, 0.0, 0.0); //Set color to black

        glMatrixMode(GL_PROJECTION);
}

//Draws the lines as specified by the lines via vertex pairs in the input file; XY: All Z values are ignored
void drawSceneBez1()
{
		glClear(GL_COLOR_BUFFER_BIT);
  	glLoadIdentity();
		glBegin(GL_LINES);

			glColor3f(.04, .15, 1);
			glVertex2f(-1,0);
			glVertex2f(1,0);
			glVertex2f(0,-1);
			glVertex2f(0,1);

		/*	for (int i = 0; i < lArr.size(); i++) { //For each polyhedron
					for (int j = 1; j < lArr.at(i).size(); j = j+2) { //For each each line in the polyhedron
					//First vArr entry appears as [6,0,100,200,200,100,200,0,300,200,200,300,200,100,200,400,100,200,0]
					//First lArr entry appears as [12,1,2,1,3,2,4,3,4,1,5,2,5,3,5,4,5,1,6,2,6,3,6,4,6]
							float x1 = vArr.at(i).at((lArr.at(i).at(j) - 1)*3+1);    //X of first point
            	float y1 = vArr.at(i).at((lArr.at(i).at(j) - 1)*3+2);    //Y of first point
            	float x2 = vArr.at(i).at((lArr.at(i).at(j+1) - 1)*3+1);   //X of second point
            	float y2 = vArr.at(i).at((lArr.at(i).at(j+1) - 1)*3+2);   //Y of second point

							glColor3f(1.0, 1.0, 1.0);
							glVertex2f(x1, y1);
							glVertex2f(x2, y2);
			}
		}
*/

	glEnd();
	glFlush();
}
//Draws the lines as specified by the lines via vertex pairs in the input file; XZ: All Y values are ignored
void drawSceneBez2()
{
        glClear(GL_COLOR_BUFFER_BIT);
        glLoadIdentity();
	glBegin(GL_LINES);

					glColor3f(.04, .15, 1);
					glVertex2f(-1,0);
					glVertex2f(1,0);
					glVertex2f(0,-1);
					glVertex2f(0,1);
					/*
        	for (int i = 0; i < lArr.size(); i++) { //For each polyhedron
                	for (int j = 1; j < lArr.at(i).size(); j = j+2) { //For each each line in the polyhedron
                        	//First vArr entry appears as [6,0,100,200,200,100,200,0,300,200,200,300,200,100,200,400,100,200,0]
                        	//First lArr entry appears as [12,1,2,1,3,2,4,3,4,1,5,2,5,3,5,4,5,1,6,2,6,3,6,4,6]
                        	float x1 = vArr.at(i).at((lArr.at(i).at(j) - 1)*3+1);    //X of first point
                        	float z1 = vArr.at(i).at((lArr.at(i).at(j) - 1)*3+3);    //Z of first point

                        	float x2 = vArr.at(i).at((lArr.at(i).at(j+1) - 1)*3+1);   //X of second point
                        	float z2 = vArr.at(i).at((lArr.at(i).at(j+1) - 1)*3+3);   //Z of second point

                                glColor3f(1.0, 1.0, 1.0);
                                glVertex2f(x1, z1);
                                glVertex2f(x2, z2);
                	}
        	}
*/
	glEnd();
	glFlush();
}
//Draws the lines as specified by the lines via vertex pairs in the input file; XZ: All Y values are ignored
void drawSceneB1()
{
        glClear(GL_COLOR_BUFFER_BIT);
        glLoadIdentity();
				glBegin(GL_LINES);

					glColor3f(.04, .15, 1);
					glVertex2f(-1,0);
					glVertex2f(1,0);
					glVertex2f(0,-1);
					glVertex2f(0,1);
					/*
        	for (int i = 0; i < lArr.size(); i++) { //For each polyhedron
                	for (int j = 1; j < lArr.at(i).size(); j = j+2) { //For each each line in the polyhedron
				//First vArr entry appears as [6,0,100,200,200,100,200,0,300,200,200,300,200,100,200,400,100,200,0]
                        	//First lArr entry appears as [12,1,2,1,3,2,4,3,4,1,5,2,5,3,5,4,5,1,6,2,6,3,6,4,6]
                        	float y1 = vArr.at(i).at((lArr.at(i).at(j) - 1)*3+2);    //Y of first point
                        	float z1 = vArr.at(i).at((lArr.at(i).at(j) - 1)*3+3);    //Z of first point

                        	float y2 = vArr.at(i).at((lArr.at(i).at(j+1) - 1)*3+2);   //Y of second point
                        	float z2 = vArr.at(i).at((lArr.at(i).at(j+1) - 1)*3+3);   //Z of second point

				glColor3f(1.0, 1.0, 1.0);
                                glVertex2f(y1, z1);
                                glVertex2f(y2, z2);
                	}
        	}

  */

	glEnd();
	glFlush();
}

void drawSceneB2()
{
        glClear(GL_COLOR_BUFFER_BIT);
        glLoadIdentity();
				glBegin(GL_LINES);

					glColor3f(.04, .15, 1);
					glVertex2f(-1,0);
					glVertex2f(1,0);
					glVertex2f(0,-1);
					glVertex2f(0,1);

					/*
        	for (int i = 0; i < lArr.size(); i++) { //For each polyhedron
                	for (int j = 1; j < lArr.at(i).size(); j = j+2) { //For each each line in the polyhedron
				//First vArr entry appears as [6,0,100,200,200,100,200,0,300,200,200,300,200,100,200,400,100,200,0]
                        	//First lArr entry appears as [12,1,2,1,3,2,4,3,4,1,5,2,5,3,5,4,5,1,6,2,6,3,6,4,6]
                        	float y1 = vArr.at(i).at((lArr.at(i).at(j) - 1)*3+2);    //Y of first point
                        	float z1 = vArr.at(i).at((lArr.at(i).at(j) - 1)*3+3);    //Z of first point

                        	float y2 = vArr.at(i).at((lArr.at(i).at(j+1) - 1)*3+2);   //Y of second point
                        	float z2 = vArr.at(i).at((lArr.at(i).at(j+1) - 1)*3+3);   //Z of second point

				glColor3f(1.0, 1.0, 1.0);
                                glVertex2f(y1, z1);
                                glVertex2f(y2, z2);
                	}
        	}
*/
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
        glutSetWindow(windowBez1);
        drawSceneBez1();
        glutSetWindow(windowBez2);
        drawSceneBez2();
        glutSetWindow(windowB1);
        drawSceneB1();
				glutSetWindow(windowB2);
				drawSceneB2();
}
/*
//Provide current vertex information for specified polyhedron
void vertexMenu(int pid)
{
        switch (pid)
        {
                case 0: cout << "Present vertices of the octahedron or polhedron 0:\n";
                        break;
                case 1: cout << "Present vertices of the tetrahedron or polyhedron 1:\n";
                        break;
                case 2: cout << "Present vertices of the hexahedron or polyhedron 2:\n";
                        break;
        }
        for (int i = 1; i < vArr.at(pid).size(); i+=3)
                cout << "(" << vArr.at(pid).at(i)*NORM << ", " << vArr.at(pid).at(i+1)*NORM << ", " << vArr.at(pid).at(i+2)*NORM << ")\n";

}
//Translation
void translateMenu(int pid)
{
        float x, y, z;
        int vertices = vArr.at(pid).at(0);
	cout << "Number of vertices is " << vertices << endl;
        cout << "Please enter the x, y, and z translation values:\n";
        cin >> x >> y >> z;
	//Norm is used to normalize world coordinates
	cout << "After normalization, you entered " << x/NORM << ", " << y/NORM << ", " << z/NORM << endl;
        for (int i = 0; i < vertices; i++) {
                vArr.at(pid).at(1+3*i) += x/NORM;
                vArr.at(pid).at(2+3*i) += y/NORM;
		vArr.at(pid).at(3+3*i) += z/NORM;

        }

	glutPostRedisplay();

	//Write changes back to file
	ofstream file;
        file.open(fileName, std::ofstream::out | std::ofstream::trunc);
        if (!file) {
                cerr << "Unable to open file\n";
                exit(ERROR);   // call system to stop
        }

        file << vArr.size() << '\n';
        for (int i = 0; i < vArr.size(); i++) {
		file << '\n' << vArr.at(i).at(0) << '\n';
                for (int j = 1; j < vArr.at(i).size(); j+=3)
                        file << vArr.at(i).at(j)*NORM << " " << vArr.at(i).at(j+1)*NORM << " " << vArr.at(i).at(j+2)*NORM << '\n';
		file << lArr.at(i).at(0) << '\n';
		for (int k = 1; k < lArr.at(i).size(); k+=2)
			file << lArr.at(i).at(k) << " " << lArr.at(i).at(k+1) << '\n';
        }
}
//Scaling
void scaleMenu(int pid)
{
        double scale;
        int vertices = vArr.at(pid).at(0);
        double xSum = 0;
        double ySum = 0;
	double zSum = 0;
        double centX, centY, centZ;

        cout << "Please enter the magnitude you'd like to scale by:\n";
        cin >> scale;

	// Calculate centroids
        for (int i = 0; i < vertices; i++)
                xSum += vArr.at(pid).at(1+i*3);

        for (int i = 0; i < vertices; i++)
                ySum += vArr.at(pid).at(2+i*3);

	for (int i = 0; i < vertices; i++)
		zSum += vArr.at(pid).at(3+i*3);

	centX = xSum / (double)vertices;
	centY = ySum / (double)vertices;
	centZ = zSum / (double)vertices;

        for (int i = 0; i < vertices; i++) {
                vArr.at(pid).at(1+i*3) = scale*(vArr.at(pid).at(1+i*3) - centX) + centX;
                vArr.at(pid).at(2+i*3) = scale*(vArr.at(pid).at(2+i*3) - centY) + centY;
		vArr.at(pid).at(3+i*3) = scale*(vArr.at(pid).at(3+i*3) - centZ) + centZ;
        }

        glutPostRedisplay();

	//Write changes back to file
        ofstream file;
        file.open(fileName, std::ofstream::out | std::ofstream::trunc);
        if (!file) {
                cerr << "Unable to open file\n";int scale_menu, rotate_menu, translate_menu, vertex_menu; //For use in graphical menu
                exit(ERROR);   // call system to stop
        }

        file << vArr.size() << '\n';
        for (int i = 0; i < vArr.size(); i++) {
                file << '\n' << vArr.at(i).at(0) << '\n';
                for (int j = 1; j < vArr.at(i).size(); j+=3)
                        file << vArr.at(i).at(j)*NORM << " " << vArr.at(i).at(j+1)*NORM << " " << vArr.at(i).at(j+2)*NORM << '\n';
                file << lArr.at(i).at(0) << '\n';
                for (int k = 1; k < lArr.at(i).size(); k+=2)
                        file << lArr.at(i).at(k) << " " << lArr.at(i).at(k+1) << '\n';
        }
}
void rotateMenu(int pid)
{
        float alpha;
	rotate = true; //Set rotate flag to true
        int vertices = vArr.at(pid).at(0);
	cout << "Please enter the x1, y1, z1, x2, y2, and z2 values to define an axis of rotation, followed by the angle of rotation in radians\n";
	cin >> rx1 >> ry1 >> rz1 >> rx2 >> ry2 >> rz2 >> alpha;
	rx1 /= NORM;
	ry1 /= NORM;
	rz1 /= NORM;
	rx2 /= NORM;
	ry2 /= NORM;
        rz2 /= NORM;
	float mag = sqrt(pow((rx2 - rx1), 2) + pow((ry2 -ry1), 2) + pow((rz2 - rz1), 2));
	float a = (rx2 - rx1) / mag;
	float b = (ry2 -ry1) / mag;
	float c = (rz2 - rz1) / mag;
	float d = sqrt(pow(b, 2) + pow(c, 2));

	// Phases are as outlined in the text, i.e. R(theta) = T^-1 * Rx^-1(alpha) * Ry^-1(beta) * Rz(theta) * Ry(beta) * Rx(alpha) * T
	for (int i = 0; i < vertices; i++) {
		float oldX, oldY, oldZ, tempX, tempY, tempZ;
		//Phase 1/7
		oldX = vArr.at(pid).at(1+i*3) - rx1;
		oldY = vArr.at(pid).at(2+i*3) - ry1;
		oldZ = vArr.at(pid).at(3+i*3) - rz1;

		//Phase 2/7
		tempX = oldX;
		tempY = (c/d)*oldY + (-1*b/d)*oldZ;
		tempZ = (b/d)*oldY + (c/d)*oldZ;

		//Phase 3/7
		oldX = d*tempX - a*tempZ;
		oldY = tempY;
		oldZ = a*tempX + d*tempZ;

		//Phase 4/7
		tempX = oldX*cos(alpha) - oldY*sin(alpha);
		tempY = oldX*sin(alpha) + oldY*cos(alpha);
		tempZ = oldZ;

		//Phase 5/7
		oldX = d*tempX + a*tempZ;
		oldY = tempY;
		oldZ = d*tempZ - a*tempX;

		//Phase 6/7
		tempX = oldX;
		tempY = (c/d)*oldY + (b/d)*oldZ;
		tempZ = (c/d)*oldZ - (b/d)*oldY;

		//Phase 7/7
                vArr.at(pid).at(1+i*3) = tempX + rx1;
                vArr.at(pid).at(2+i*3) = tempY + ry1;
                vArr.at(pid).at(3+i*3) = tempZ + rz1;
	}

        glutPostRedisplay();

	//Write changes back to file
        ofstream file;
        file.open(fileName, std::ofstream::out | std::ofstream::trunc);
        if (!file) {
                cerr << "Unable to open file\n";
                exit(ERROR);   // call system to stop
        }

        file << vArr.size() << '\n';
        for (int i = 0; i < vArr.size(); i++) {
                file << '\n' << vArr.at(i).at(0) << '\n';
                for (int j = 1; j < vArr.at(i).size(); j+=3)
                        file << vArr.at(i).at(j)*NORM << " " << vArr.at(i).at(j+1)*NORM << " " << vArr.at(i).at(j+2)*NORM << '\n';
                file << lArr.at(i).at(0) << '\n';
                for (int k = 1; k < lArr.at(i).size(); k+=2)
                        file << lArr.at(i).at(k) << " " << lArr.at(i).at(k+1) << '\n';
        }
}*/
void mainMenu(int pid)
{




}

int main(int argc, char** argv)
{
	if (argc != 1) {
		cout << "Usage: p4\n";
		exit(ERROR);
	}

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
	glutInitWindowPosition(200, 100);
	glutInitWindowSize(900, 900);

	windowID = glutCreateWindow("Bezier (Left) and B-Spine (Right) Curves");
	init();

        vector<float> v;
        float num;
				fstream file;
        file.open("inputFile.txt");
        if (!file) {
                cerr << "Unable to open file\n";
                exit(ERROR);
        }
        while (file >> num)
                v.push_back(num); //Initial vector for all polygons
        file.close();

				populateCurves(v);

	//Bez1
	windowBez1 = glutCreateSubWindow(windowID, 50, 75, 350, 350);
	init();

	//Bez2
	windowBez2 = glutCreateSubWindow(windowID, 50, 475, 350, 350);
	init();

	//B1
	windowB1 = glutCreateSubWindow(windowID, 450, 475, 350, 350);
	init();

	//B2
	windowB2 = glutCreateSubWindow(windowID, 450, 75, 350, 350);
	init();

				glutSetWindow(windowBez1);
        glutCreateMenu(mainMenu);
								glutAddMenuEntry("Insert Control Point", 0);
								glutAddMenuEntry("Delete Control Point", 1);
								glutAddMenuEntry("Modify Control Point", 2);
        glutAttachMenu(GLUT_RIGHT_BUTTON);

			  glutSetWindow(windowBez2);
				glutCreateMenu(mainMenu);
								glutAddMenuEntry("Insert Control Point", 3);
								glutAddMenuEntry("Delete Control Point", 4);
								glutAddMenuEntry("Modify Control Point", 5);
				glutAttachMenu(GLUT_RIGHT_BUTTON);

				glutSetWindow(windowB1);
				glutCreateMenu(mainMenu);
								glutAddMenuEntry("Insert Control Point", 6);
								glutAddMenuEntry("Delete Control Point", 7);
								glutAddMenuEntry("Modify Control Point", 8);
				glutAttachMenu(GLUT_RIGHT_BUTTON);

				glutSetWindow(windowB2);
				glutCreateMenu(mainMenu);
								glutAddMenuEntry("Insert Control Point", 9);
								glutAddMenuEntry("Delete Control Point", 10);
								glutAddMenuEntry("Modify Control Point", 11);
				glutAttachMenu(GLUT_RIGHT_BUTTON);


        glutDisplayFunc(display);


	glutMainLoop();

	return 0;
}
