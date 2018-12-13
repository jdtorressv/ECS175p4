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
	  Point(float xVal, float yVal);
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
float bez1Res, bez2Res, b1Res, b2Res;
int windowID, windowBez1, windowBez2, windowB1, windowB2;

/*****************************************************************************/
/* Function Definitions                                                      */
/*****************************************************************************/

//Constructors
//Point::Point(float x, float y) : x(x), y(y) {}
Point::Point(float xVal, float yVal) : x(xVal), y(yVal) {}

void populateCurves()
{
		vector<float> v;
		float num;
		fstream file;
		file.open("inputFile.txt");
		if (!file) {
						cerr << "Unable to open file!\n";
						exit(ERROR);
		}
		while (file >> num)
						v.push_back(num);
		file.close();

		auto vpoint = v.begin();
		int bezTotal = (int)*vpoint;
		if (bezTotal != 2) {
				cout << "This program only supports two Bezier curve instances\n";
				exit(ERROR);
		}
		int cntrlTotal = (int)*(++vpoint);
		for (int i = 0; i < cntrlTotal; i++) {
				float xVal = *(++vpoint) / NORM;
				float yVal = *(++vpoint) / NORM;
				Point pt(xVal, yVal);
				bez1.pointArr.push_back(pt);
		}
		cntrlTotal = (int)*(++vpoint);
		for (int i = 0; i < cntrlTotal; i++) {
				float xVal = *(++vpoint) / NORM;
				float yVal = *(++vpoint) / NORM;
				Point pt(xVal, yVal);
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
				float xVal = *(++vpoint) / NORM;
				float yVal = *(++vpoint) / NORM;
				Point pt(xVal, yVal);
				b1.pointArr.push_back(pt);
		}
		b2.k = (int)*(++vpoint);
		cntrlTotal = (int)*(++vpoint);
		for (int i = 0; i < cntrlTotal; i++) {
				float xVal = *(++vpoint) / NORM;
				float yVal = *(++vpoint) / NORM;
				Point pt(xVal, yVal);
				b2.pointArr.push_back(pt);
		}
}
Point getNextCastelPt(int s, int j, float t, BezSpline spline)
{
		if (s == 0)
				return spline.pointArr.at(j);
		Point p1 = getNextCastelPt(s - 1, j, t, spline);
		Point p2 = getNextCastelPt(s - 1, j + 1, t, spline);
		Point nextPt(((1-t) * p1.x + t * p2.x), ((1-t) * p1.y + t * p2.y));
		return nextPt;
}
void casteljau(BezSpline spline, float resolution)
{
		glBegin(GL_POINTS);

				Point currPt(0, 0);
				for (float t = 0; t <= 1.0; t += resolution)
				{
					currPt = getNextCastelPt(spline.pointArr.size() - 1, 0, t, spline);
					glVertex2f(currPt.x, currPt.y);
				}

		glEnd();
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

				glColor3f(.96, .38, .53);
				for (int i = 0; i < bez1.pointArr.size()-1; i++) {
						glVertex2f(bez1.pointArr.at(i).x, bez1.pointArr.at(i).y);
						glVertex2f(bez1.pointArr.at(i+1).x, bez1.pointArr.at(i+1).y);
				}

		glEnd();
		for (int i = 0; i < bez1.pointArr.size(); i++) {
				cout << bez1.pointArr.at(i).x << " " << bez1.pointArr.at(i).y << endl;
		}
		casteljau(bez1, bez1Res);
		glFlush();
}
//
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

		glColor3f(.96, .38, .53);
		for (int i = 0; i < bez2.pointArr.size()-1; i++) {
				glVertex2f(bez2.pointArr.at(i).x, bez2.pointArr.at(i).y);
				glVertex2f(bez2.pointArr.at(i+1).x, bez2.pointArr.at(i+1).y);
		}

		glEnd();
		casteljau(bez2, bez2Res);
		glFlush();
}
//
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

		glColor3f(.96, .38, .53);
		for (int i = 0; i < b1.pointArr.size()-1; i++) {
				glVertex2f(b1.pointArr.at(i).x, b1.pointArr.at(i).y);
				glVertex2f(b1.pointArr.at(i+1).x, b1.pointArr.at(i+1).y);
		}
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

		glColor3f(.96, .38, .53);
		for (int i = 0; i < b2.pointArr.size()-1; i++) {
				glVertex2f(b2.pointArr.at(i).x, b2.pointArr.at(i).y);
				glVertex2f(b2.pointArr.at(i+1).x, b2.pointArr.at(i+1).y);
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
		oldZ = d*tempZ - a*tempX;	for (int i = 0; i < v.size(); i++)
						cout << v.at(i) << endl;

		//Phase 6/7
		tempX = oldX;
		tempY = (c/d)*oldY + (b/d)*oldZ;
		tempZ = (c/d)*oldZ - (b/d)*oldY;

		//Phase 7/7
                vArr.at(pid).at(1+i*3) = tempX + rx1;
                vArr.at(pid).at(2+i*3) = tempY + ry1;
                vArr.at(pid).at(3+i*3) = tempZ + rz1;
	}
	for (int i = 0; i < v.size(); i++)
						cout << v.at(i) << endl;
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
void insertBezPt(BezSpline& spline)
{
			float xVal;
			float yVal;
			cout << "Please enter the x and y coordinates for the control point you would like to add\n";
			cin >> xVal >> yVal;
			if (fabs(xVal) > NORM || fabs(yVal) > NORM) {
					cout << "Values entered must be in the range [-500.0, 500.0]\n";
					exit(ERROR);
			}
			Point pt(xVal / NORM, yVal / NORM);
			spline.pointArr.push_back(pt);
			glutSetWindow(windowID);
			glutPostRedisplay();
}
void deleteBezPt(BezSpline& spline)
{
			int index;
			cout << "This curve is governed by the following control points:\n\n";
			for (int i = 0; i < spline.pointArr.size(); i++) {
					cout << i+1 <<"	(" << spline.pointArr.at(i).x * NORM << ", " << spline.pointArr.at(i).y * NORM << ")\n";
			}
			cout << "\nPlease enter the number of the point you'd like removed:\n";
			cin >> index;
			index--;
			vector<Point> temp;
			for (int i = 0; i < spline.pointArr.size(); i++) {
					if (i != index)
					temp.push_back(spline.pointArr.at(i));
			}
			spline.pointArr.clear();
			for (int i = 0; i < temp.size(); i++) {
					spline.pointArr.push_back(temp.at(i));
			}
			glutSetWindow(windowID);
			glutPostRedisplay();

}
void modifyBezPt(BezSpline& spline)
{
			int index;
			float xVal, yVal;
			cout << "This curve is governed by the following control points:\n\n";
			for (int i = 0; i < spline.pointArr.size(); i++) {
					cout << i+1 <<"	(" << spline.pointArr.at(i).x * NORM << ", " << spline.pointArr.at(i).y * NORM << ")\n";
			}
			cout << "\nPlease enter the number of the point you'd like to modify, followed by the new x and y coordinates:\n";
			cin >> index >> xVal >> yVal;
			index--;
			spline.pointArr.at(index).x = xVal / NORM;
			spline.pointArr.at(index).y = yVal / NORM;
			glutSetWindow(windowID);
			glutPostRedisplay();
}

void mainMenu(int pid)
{
		switch(pid)
		{
					case 0: // Insert Bez1
							insertBezPt(bez1);
							break;
					case 1: // Delete Bez1
							deleteBezPt(bez1);
							break;
					case 2: // Modify Bez1
							modifyBezPt(bez1);
							break;
					case 3: // Insert Bez2
							insertBezPt(bez2);
							break;
					case 4: // Delete Bez2
							deleteBezPt(bez2);
							break;
					case 5: // Modify Bez2
							modifyBezPt(bez2);
							break;
					case 6: // Insert B1
							//
							break;
					case 7: // Delete B1
							//
							break;
					case 8: // Modify B1
							//
							break;
					case 9: // Insert B2
							//
							break;
					case 10: // Delete B2
							//
							break;
					case 11: // Modify B2
							//
							break;
		}
}
/*****************************************************************************/
/* Main Function                                                             */
/*****************************************************************************/
int main(int argc, char** argv)
{
	if (argc != 5) {
		cout << "Usage: p4 <bez1Res, bez2Rez, b1Res, b2Res>\n";
		exit(ERROR);
	}

	bez1Res = atof(argv[1]);
	bez2Res = atof(argv[2]);
	b1Res = atof(argv[3]);
	b2Res = atof(argv[4]);

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
	glutInitWindowPosition(200, 100);
	glutInitWindowSize(900, 900);

	windowID = glutCreateWindow("Bezier (Left) and B-Spine (Right) Curves");
	init();

	populateCurves();

	//Bez1
	windowBez1 = glutCreateSubWindow(windowID, 50, 75, 350, 350);
	init();

	//Bez2
	windowBez2 = glutCreateSubWindow(windowID, 50, 475, 350, 350);
	init();

	//B1
	windowB1 = glutCreateSubWindow(windowID, 450, 75, 350, 350);
	init();

	//B2
	windowB2 = glutCreateSubWindow(windowID, 450, 475, 350, 350);
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

	glutSetWindow(windowID);
  glutDisplayFunc(display);

	glutMainLoop();

	return 0;
}
