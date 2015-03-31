// 
// Torbert, 3.12.2009
// 
// OpenGL Demo, 3-D Example
// 

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <GL/glut.h>
#include <time.h>
#include <ctype.h>
#define N 3000//N is the max amount of bodies
#define G .005
///////////////////////////////
double pi=3.1415926;

int    w=640,h=480;
int currentindex=0;
double rho,phi,theta,up=1.0;
// the center of the screen
double xc,yc,zc;
double xe,ye,ze;
int collision=0;
int centered=0;
///////////////////////////////
typedef struct
{
    double mass;
    double radius;
    double x;
    double y;
    double z;
    double xv;
    double yv;
    double zv;
    double r;
    double g;
    double b;
    int flag;

} Body;

Body bodies[N];
double t=0.0;
int previousx;
int previousy;
double dt=.01;
double dist(double x1, double y1, double z1, double x2, double y2, double z2)
{
    double xpart=x1-x2;
    double ypart=y1-y2;
    double zpart=z1-z2;
    return sqrt( xpart*xpart + ypart*ypart + zpart*zpart);
}
void display(void)
{
    //    double t;

    // clear the screen

    glClear(GL_COLOR_BUFFER_BIT);    
    glColor3f(0.0,0.0,0.0); 
    //draw lines for stuff
    glBegin(GL_LINES);
    glVertex3f(0.0f+xc, 0.0f+yc, 0.0f+zc);
    glVertex3f(50.0f+xc, 0.0f+yc, 0.0f+zc);

    glVertex3f(0.0f+xc, 0.0f+yc, 0.0f+zc);
    glVertex3f(0.0f+xc, 10.0f+yc, 0.0f+zc);

    glVertex3f(0.0f+xc, 0.0f+yc, 0.0f+zc);
    glVertex3f(0.0f+xc, 0.0f+yc, 10.0f+zc);
    glEnd();
    glRasterPos3f(50.0f+xc, 0.0f+yc, 0.0f+zc);
    glutBitmapString(GLUT_BITMAP_HELVETICA_18,"x");
    glRasterPos3f(40.0f+xc, -2.0f+yc, 0.0f+zc);
    glutBitmapString(GLUT_BITMAP_TIMES_ROMAN_24,"w");
    glRasterPos3f(10.0f+xc, -2.0f+yc, 0.0f+zc);
    glutBitmapString(GLUT_BITMAP_TIMES_ROMAN_24,"s");

    glRasterPos3f(0.0f+xc, 50.0f+yc, 0.0f+zc);
    glutBitmapString(GLUT_BITMAP_HELVETICA_18,"y"); 
    glRasterPos3f(-2.0f+xc, 40.0f+yc, 0.0f+zc);
    glutBitmapString(GLUT_BITMAP_TIMES_ROMAN_24,"e");
    glRasterPos3f(-2.0f+xc, 10.0f+yc, 0.0f+zc);
    glutBitmapString(GLUT_BITMAP_TIMES_ROMAN_24,"d");

    glRasterPos3f(0.0f+xc, 0.0f+yc, 50.0f+zc);
    glutBitmapString(GLUT_BITMAP_HELVETICA_18,"z"); 
    glRasterPos3f(0.0f+xc, -2.0f+yc, 40.0f+zc);
    glutBitmapString(GLUT_BITMAP_TIMES_ROMAN_24,"r");
    glRasterPos3f(0.0f+xc, -2.0f+yc, 10.0f+zc);
    glutBitmapString(GLUT_BITMAP_TIMES_ROMAN_24,"f");

    int i;
    for(i=0;i<currentindex;i++)
    {
        int q; 
        double fx=0;
        double fy=0;
        double fz=0;
        for(q=0; q<currentindex;q++)
        {

            if(q!=i)
            {
                double d=dist(bodies[i].x,bodies[i].y, bodies[i].z,bodies[q].x, bodies[q].y,bodies[q].z);
                double top=(G*bodies[q].mass);
                double magnitude=top/(d);
                fx-=magnitude*((bodies[i].x-bodies[q].x)/d);
                fy-=magnitude*((bodies[i].y-bodies[q].y)/d);
                fz-=magnitude*((bodies[i].z-bodies[q].z)/d);
                if(d< bodies[i].radius*7/8+bodies[q].radius*7/8 && bodies[i].flag==0 && collision==0)
                {
                    if(bodies[q].mass<bodies[i].mass)
                        bodies[q].flag=1;
                    else
                        bodies[i].flag=1;
                }
            }

        }
        //get correct direction

        bodies[i].xv+=fx*dt;
        bodies[i].yv+=fy*dt;
        bodies[i].zv+=fz*dt; 
        bodies[i].x+=bodies[i].xv*dt;
        bodies[i].y+=bodies[i].yv*dt;
        bodies[i].z+=bodies[i].zv*dt;
        int m=(int)bodies[i].mass;
        //   glColor3f((double)(m%10000)/10000.0,(double)(m%1000)/1000.0,(double)(m%100)/100.0);
        glColor3f(bodies[i].r,bodies[i].g,bodies[i].b);
        glMatrixMode(GL_MODELVIEW);
        glPushMatrix();
        glTranslatef(bodies[i].x,bodies[i].y,bodies[i].z);
        glutWireSphere(bodies[i].radius, 10, 10);
        glPopMatrix();
    }
    int flagged=0;
    if( collision==0)
#pragma omp for private(i)
        for(i=0;i<currentindex;i++)
        {
            if(bodies[i].flag==1)
            {
                //     printf("%d\n", bodies[i].flag);
                flagged++;
                int q;
                for(q=0; q<currentindex;q++)
                {
                    if(q!=i)
                    {
                        double d=dist(bodies[i].x,bodies[i].y, bodies[i].z,bodies[q].x, bodies[q].y,bodies[q].z);
                        if(d< bodies[i].radius*7/8+bodies[q].radius*7/8 && bodies[q].flag==0)
                        {
                            bodies[q].mass+=bodies[i].mass;
                            bodies[q].radius+= pow((4/3*pi*bodies[i].radius*bodies[i].radius+4/3*pi*bodies[q].radius*bodies[q].radius)/(4/3),1/10);
                            bodies[q].xv+=bodies[i].xv*bodies[i].mass/bodies[q].mass;
                            bodies[q].yv+=bodies[i].yv*bodies[i].mass/bodies[q].mass;
                            bodies[q].zv+=bodies[i].zv*bodies[i].mass/bodies[q].mass;
                            bodies[q].r=(bodies[q].r+bodies[i].r)/2;
                            bodies[q].g=(bodies[q].g+bodies[i].g)/2;
                            bodies[q].b=(bodies[q].b+bodies[i].b)/2;
                            bodies[i].r=(bodies[q].r+bodies[i].r)/2;
                            bodies[i].g=(bodies[q].g+bodies[i].g)/2;
                            bodies[i].b=(bodies[q].b+bodies[i].b)/2;
                        }


                    }
                }
            }
        }
    int count=0;
    if( collision==0)
        for(i=0; i<currentindex; i++)
        {
            if(bodies[i].flag==0){
                bodies[count].mass	=bodies[i].mass;
                bodies[count].radius	=bodies[i].radius;
                bodies[count].x	=bodies[i].x;
                bodies[count].y	=bodies[i].y;
                bodies[count].z	=bodies[i].z;
                bodies[count].xv	=bodies[i].xv;
                bodies[count].yv	=bodies[i].yv;
                bodies[count].zv	=bodies[i].zv;
                bodies[count].r	=bodies[i].r;
                bodies[count].g	=bodies[i].g;
                bodies[count].b	=bodies[i].b;
                bodies[count].flag=0;
                count++;
            }
        }
    currentindex-=flagged;




    //     glutWireSphere(.5, 20, 16);
    t+=0.001;
    // printf("%f\n",t);
    glutSwapBuffers();
}
void look()
{
    xe=xc+rho*sin(theta)*sin(phi); // y
    ye=yc+rho*cos(theta);          // z
    ze=zc+rho*sin(theta)*cos(phi); // x
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(xe,ye,ze, xc,yc,zc, 0.0,up,0.0);
}
void idle(void)
{

    if(up>0.0)
    {
        if(theta+0.001>pi)
        {
            up*=-1.0;
            theta=-pi;
        }
        if(theta-0.001<0)
        {
            up*=-1.0;
        }
    }
    else{
        if(theta+0.001>0.0)
            up*=-1.0;
        if(theta-0.001<-pi)
        {
            theta=pi;
        }
    }
    if(centered==1)
    {
        double topx=0;
        double topy=0;
        double topz=0;
        double total=0;
        int in;
        for(in=0;in<currentindex;in++)
        {
            double x=bodies[in].x; 
            double y=bodies[in].y; 
            double z=bodies[in].z;
            double m=bodies[in].mass;
            topx+=x*m;
            topy+=y*m;
            topz+=z*m;
            total+=m;
        }
        xc=topx/total;
        yc=topy/total;
        zc=topz/total;
    }

    look();
    glutPostRedisplay();
}
void mouse(int button,int state,int xscr,int yscr)
{
    if(button==0 && state==0)
    {
        previousx=xscr;
        previousy=yscr;
    }
}
void motion(int xscr,int yscr)
{
    int deltax=xscr-previousx;
    int deltay=yscr-previousy;
    previousx=xscr;
    previousy=yscr;
    theta += (double)deltay/100.0;
    phi   += (double)deltax/100.0;
}
void mouse_wheel(int wheel,int direction,int xscr,int yscr)
{
    //   printf("%d\n", direction);
    if(direction==-1)
        rho+=2;
    else
        rho-=2;

    if(rho<0.1)
        rho=0.1;
}

int move_val = 10;
void keyfunc(unsigned char key,int xscr,int yscr)
{
    if(key=='q')
    {
        exit(0);
    }
    if(key=='w')
    {
        xc+=move_val;
    }
    if(key=='s')
    {
        xc-=move_val;
    }
    if(key=='e')
    {
       yc+=move_val; 
    }
    if(key=='d')
    {
        yc-=move_val; 
    }
    if(key=='r')
    {
        zc+=move_val; 
    }
    if(key=='f')
    {
        zc-=move_val; 
    }
    if(key==' ')
    {
        xc-=xe/5;
        yc-=ye/5;
        zc-=ze/5;
    }
    if(key=='c')
    {
        centered=(!centered); 
    }
}
void reshape(int wscr,int hscr)
{
    GLfloat aspect_ratio;

    w=wscr;
    h=hscr;
    aspect_ratio=(GLfloat)w/(GLfloat)h;
    glViewport(0,0,(GLsizei)w,(GLsizei)h);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(60.0,aspect_ratio,0.1,100000.0);

    look();
}

void makeBody(double m, double r, double xpos, double ypos, double zpos,double i, double j, double k,double rc, double gc, double bc)
{
    if(currentindex==N)
    {
        return;
    }
    int index=currentindex;
    bodies[index].mass=m;
    bodies[index].radius=r;
    bodies[index].x=xpos;
    bodies[index].y=ypos;
    bodies[index].z=zpos;
    bodies[index].xv=i;
    bodies[index].yv=j;
    bodies[index].zv=k;
    bodies[index].r=rc;
    bodies[index].g=gc;
    bodies[index].b=bc;

    bodies[index].flag=0;
    currentindex++;
}
void randombody()
{
    int index;
    double mx=0;
    double my=0;
    double mz=0;
    for(index=0;index<N;index++)
    {
        if(index<N-1)
        {
            double m=1;

            m=rand()%5000;
            double r=m/1000;
            double xpos=(rand()%100)-50;
            //   double xpos=0;
            double ypos=(rand()%100)-50;
            double zpos=(rand()%100)-50;
            //   double zpos=0;
            double i=(rand()%20)-10;
            double j=(rand()%20)-10;
            double k=(rand()%20)-10;
            //   double i=0;
            //   double j=0;
            //   double k=0;
            mx+=m*i;
            my+=m*j;
            mz+=m*k;
            double rc=(double)(rand()%256)/256.0;
            double g=(double)(rand()%256)/256.0;
            double b=(double)(rand()%256)/256.0;

            makeBody(m,r,xpos,ypos,zpos,i,j,k,rc,g,b); 
        }
        else
        {
            double m=1;

            m=rand()%500+500;
            double r=m/1000;
            double xpos=(rand()%50)-25;
            //   double xpos=0;
            double ypos=(rand()%50)-25;
            double zpos=(rand()%50)-25;
            //   double zpos=0;
            double i=mx/m;
            double j=my/m;
            double k=mz/m;

            double rc=(double)(rand()%256)/256.0;
            double g=(double)(rand()%256)/256.0;
            double b=(double)(rand()%256)/256.0; 
            makeBody(m,r,xpos,ypos,zpos,i,j,k,rc,g,b); 
        }
    }
}
void singlebody()
{
    int factor=2000;
    double m=rand()%factor;
    double r=m/factor;
    double xpos=(rand()%500)-250;
    //   double xpos=0;
    double ypos=(rand()%500)-250;
    //   double zpos=(rand()%100)-50;
    double zpos=0;
    double i=(rand()%40)-20;
    double j=(rand()%40)-20;
    double k=0;

    //   double k=(rand()%20)-10;

    double rc=(double)(rand()%256)/256.0;
    double g=(double)(rand()%256)/256.0;
    double b=(double)(rand()%256)/256.0;
    makeBody(m,r,xpos,ypos,zpos,i,j,k,rc,g,b); 

}
void asteroid()
{
    double m=rand()%50;
    double r=m/10;
    double xpos=(double)(rand()%400)-200;
    //   double xpos=0;
    double ypos=(rand()%400)-200;
    //   double zpos=(rand()%100)-50;
    double zpos=0;
    double i=(rand()%50)-25;
    double j=(rand()%50)-25;
    double k=0;

    //   double k=(rand()%20)-10;

    double rc=(double)(rand()%256)/256.0;
    double g=(double)(rand()%256)/256.0;
    double b=(double)(rand()%256)/256.0;
    makeBody(m,r,xpos,ypos,zpos,i,j,k,rc,g,b); 
}
void ringsystem(double x, double y,double z, double xv, double yv, double zv, double r, int planets)
{
    double sunmass=100000;
    double planetmass=100;
    makeBody(sunmass,40, x, y, z, xv, yv, zv,1,222.0/255.0,0);
    int e;
    double ov=10*sqrt((sunmass*sunmass*G)/((planetmass+sunmass)*r));
    printf("%f %f\n",ov, sunmass*G);
    double angle=2*pi/planets;
    for(e=1;e<=planets;e++)
    {
        makeBody(planetmass   ,5 ,x+cos(angle*e)*r, y+sin(angle*e)*r, z, xv+cos(angle*e - pi/2)*ov, yv+sin(angle*e - pi/2)*ov, zv,0,0,1);

    }
}




int main(int argc,char* argv[])
{  
    int x;
    srand(time(NULL));
    ringsystem(0,0,  0,0,0,0,100,1);
    ringsystem(0,300,  0,0,0,0,100,1);


    rho=1000.1;
    phi=0.0;
    theta=pi/2.0;
    xc=yc=zc=0.0;

    glutInit(&argc,argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
    glutInitWindowSize(w,h);
    glutInitWindowPosition(100,50);
    glutCreateWindow("OpenGL Demo");

    glClearColor(1.0,1.0,1.0,0.0);
    glShadeModel(GL_SMOOTH);

    glutDisplayFunc(display);
    glutIdleFunc(idle);
    glutMouseFunc(mouse);
    glutMotionFunc(motion);
    glutMouseWheelFunc(mouse_wheel);
    glutKeyboardFunc(keyfunc);
    glutReshapeFunc(reshape);

    glutMainLoop();

    return 0;
}
