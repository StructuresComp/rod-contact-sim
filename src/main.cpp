/**
 * simDER
 * simDER stands for "[sim]plified [D]iscrete [E]lastic [R]ods"
 * Dec 2017
 * This code is based on previous iterations.
 * */

//This line is for mac
//#include <GLUT/glut.h>

//This is for linux
#include <GL/glut.h>

#include <iostream>
#include <fstream>
#include <ctime>
#include "eigenIncludes.h"

// Rod and stepper are included in the world
#include "world.h"
#include "setInput.h"


world myWorld;
int NPTS;
ofstream simfile;
ofstream configfile;
int rods;

clock_t start;
clock_t finish;
double time_taken;

static void Key(unsigned char key, int x, int y) {
    switch (key) // ESCAPE to quit
    {
        case 27:
            exit(0);
    }
}

/* Initialize OpenGL Graphics */
void initGL() {
    glClearColor(0.7f, 0.7f, 0.7f, 0.0f); // Set background color to black and opaque
    glClearDepth(10.0f);                   // Set background depth to farthest
    //glEnable(GL_DEPTH_TEST);   // Enable depth testing for z-culling
    //glDepthFunc(GL_LEQUAL);    // Set the type of depth-test
    glShadeModel(GL_SMOOTH);   // Enable smooth shading
    //glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);  // Nice perspective corrections

    glLoadIdentity();
//	gluLookAt(0.05, 0.05, 0.1, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
//	gluLookAt(0.00, -0.3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0);

// gluLookAt(0.05, 0.05, 0.1, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);

    gluLookAt(0.00, 0.00, 0.3, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0);

    glPushMatrix();

    //glMatrixMode(GL_MODELVIEW);
}

void display(void) {
    while (myWorld.simulationRunning() > 0) {
        //  Clear screen and Z-buffer
        glClear(GL_COLOR_BUFFER_BIT);

        // draw axis
        double axisLen = 1;
        glLineWidth(0.5);

        glBegin(GL_LINES);
        glColor3f(1.0, 0.0, 0.0);
        glVertex3f(0.0, 0.0, 0.0);
        glVertex3f(axisLen, 0.0, 0.0);

        glColor3f(0.0, 1.0, 0.0);
        glVertex3f(0.0, 0.0, 0.0);
        glVertex3f(0.0, axisLen, 0.0);

        glColor3f(0.0, 0.0, 1.0);
        glVertex3f(0.0, 0.0, 0.0);
        glVertex3f(0.0, 0.0, axisLen);
        glEnd();

        //draw a line
        glColor3f(0.1, 0.1, 0.1);
        glLineWidth(3.0);

        for (int i = 0; i < rods; i++) {
            glBegin(GL_LINES);
            for (int j = 0; j < NPTS - 1; j++) {
                glVertex3f(myWorld.getScaledCoordinate(i, 4 * j), myWorld.getScaledCoordinate(i, 4 * j + 1),
                           myWorld.getScaledCoordinate(i, 4 * j + 2));
                glVertex3f(myWorld.getScaledCoordinate(i, 4 * (j + 1)), myWorld.getScaledCoordinate(i, 4 * (j + 1) + 1),
                           myWorld.getScaledCoordinate(i, 4 * (j + 1) + 2));
            }
            glEnd();
        }

        glFlush();

        // Update step
        start = clock();
        myWorld.updateTimeStep();
        finish = clock();
        time_taken = double(finish - start) / double(CLOCKS_PER_SEC);
        myWorld.CoutData(simfile, configfile, time_taken);
    }
    exit(1);
}

int main(int argc, char *argv[]) {
    setInput inputData;
    inputData = setInput();
    inputData.LoadOptions(argv[1]);
    inputData.LoadOptions(argc, argv);

    myWorld = world(inputData);
    myWorld.setRodStepper();
    rods = inputData.GetIntOpt("numFlagella");


    myWorld.OpenFile(simfile, configfile);

    bool render = myWorld.isRender();
    if (render) // if OpenGL visualization is on
    {
        NPTS = myWorld.numPoints();

        glutInit(&argc, argv);
        glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
        glutInitWindowSize(1000, 1000);
        glutInitWindowPosition(100, 100);
        glutCreateWindow("simDER");
        initGL();
        glutKeyboardFunc(Key);
        glutDisplayFunc(display);
        glutMainLoop();
    } else {
        while (myWorld.simulationRunning() > 0) {
            start = clock();
            myWorld.updateTimeStep(); // update time step
            finish = clock();
            time_taken = double(finish - start) / double(CLOCKS_PER_SEC);
            myWorld.CoutData(simfile, configfile, time_taken); // write data to file
        }
    }

    // Close (if necessary) the data file
    myWorld.CloseFile(simfile, configfile);

    return 0;
}
