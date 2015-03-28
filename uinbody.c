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
#define N 500 //N is the max amount of bodies
#define G .005
///////////////////////////////
double pi=3.1415926;

int    w=640,h=480;
int currentindex=0;
double rho,phi,theta,up=1.0;
double xc,yc,zc,xe,ye,ze;
int num=N;
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
// void setplanet(int index, double m, double r, double xpos, double ypos, double zpos,double i, double j, double k)
// {
//   bodies[index].mass=m;
//   bodies[index].radius=r;
//   bodies[index].x=xpos;
//   bodies[index].y=ypos;
//   bodies[index].z=zpos;
//   bodies[index].xv=i;
//   bodies[index].yv=j;
//   bodies[index].zv=k;
//   bodies[index].flag=0;
// }
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
  glBegin(GL_LINES);
  glVertex3f(0.0f+xc, 0.0f+yc, 0.0f+zc);
  glVertex3f(10.0f+xc, 0.0f+yc, 0.0f+zc);

  glVertex3f(0.0f+xc, 0.0f+yc, 0.0f+zc);
  glVertex3f(0.0f+xc, 10.0f+yc, 0.0f+zc);
  
  glVertex3f(0.0f+xc, 0.0f+yc, 0.0f+zc);
  glVertex3f(0.0f+xc, 0.0f+yc, 10.0f+zc);
  glEnd();
  glRasterPos3f(10.0f+xc, 0.0f+yc, 0.0f+zc);
  glutBitmapString(GLUT_BITMAP_HELVETICA_18,"x");
  glRasterPos3f(7.0f+xc, -2.0f+yc, 0.0f+zc);
  glutBitmapString(GLUT_BITMAP_TIMES_ROMAN_10,"w");
  glRasterPos3f(2.0f+xc, -2.0f+yc, 0.0f+zc);
  glutBitmapString(GLUT_BITMAP_TIMES_ROMAN_10,"s");
  
  glRasterPos3f(0.0f+xc, 10.0f+yc, 0.0f+zc);
  glutBitmapString(GLUT_BITMAP_HELVETICA_18,"y"); 
  glRasterPos3f(-2.0f+xc, 7.0f+yc, 0.0f+zc);
  glutBitmapString(GLUT_BITMAP_TIMES_ROMAN_10,"e");
  glRasterPos3f(-2.0f+xc, 2.0f+yc, 0.0f+zc);
  glutBitmapString(GLUT_BITMAP_TIMES_ROMAN_10,"d");
  
  glRasterPos3f(0.0f+xc, 0.0f+yc, 10.0f+zc);
  glutBitmapString(GLUT_BITMAP_HELVETICA_18,"z"); 
  glRasterPos3f(0.0f+xc, -2.0f+yc, 7.0f+zc);
  glutBitmapString(GLUT_BITMAP_TIMES_ROMAN_10,"r");
  glRasterPos3f(0.0f+xc, -2.0f+yc, 2.0f+zc);
  glutBitmapString(GLUT_BITMAP_TIMES_ROMAN_10,"f");
  
int i;
for(i=0;i<num;i++)
{
//   Body p=/*bodies[i]*/;

//   printf("%f, %f, %f\n",bodies[i].x,bodies  for(i=0;i<num;i++)

  //force(bodies[i],i);
    int q; 
  double fx=0;
  double fy=0;
  double fz=0;
  for(q=0; q<num;q++)
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
  for(i=0;i<num;i++)
{
  if(bodies[i].flag==1)
  {
//     printf("%d\n", bodies[i].flag);
    flagged++;
    int q;
      for(q=0; q<num;q++)
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
// printf("flagged: %d\n",flagged);
// Body temp[num-flagged];
int count=0;
if( collision==0)
  for(i=0;i<num;i++)
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
num-=flagged;


    

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
if(centered==0)
{
  double topx=0;
  double topy=0;
  double topz=0;
  double total=0;
//   double cx=0;
//   double cy=0;
//   double cz=0;
  int in;
  for(in=0;in<num;in++)
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
//   if(deltay<0)
//     up=1
// //   int current=(int)(theta/pi);
// //   if(theta>1)
// //   phi-=(current*phi);
// //   else
// //     phi+=(current*phi);
// //   if(theta+(double)deltay/100.0<current && theta>current)
// //     deltax*=-1;
// //   else if(theta+(double)deltay/100.0>current && theta<current)
// //     deltax*=-1;
    theta += (double)deltay/100.0;
      phi   += (double)deltax/100.0;
//   theta=abs(theta);
//      printf("%f, %f\n", phi, theta);
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
void keyfunc(unsigned char key,int xscr,int yscr)
{
	if(key=='q')
	{
		exit(0);
	}
	if(key=='w')
	{
	  xc+=3;
	}
	if(key=='s')
	{
	  xc-=3;
	}
	if(key=='e')
	{
	 yc+=3; 
	}
	if(key=='d')
	{
	 yc-=3; 
	}
	if(key=='r')
	{
	 zc+=3; 
	}
	if(key=='f')
	{
	 zc-=3; 
	}
	if(key==' ')
	{
	 xc-=xe/5;
	 yc-=ye/5;
	 zc-=ze/5;
	}
	if(key='c')
	{
	 centered=(centered+1)%2; 
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
  
  double m=rand()%250;
  double r=m/20;
  double xpos=(rand()%200)-100;
//   double xpos=0;
  double ypos=(rand()%200)-100;
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
     double ov=20*sqrt((sunmass*sunmass*G)/((planetmass+sunmass)*r));
     printf("%f %f\n",ov, sunmass*G);
  double angle=2*pi/planets;
  for(e=1;e<=planets;e++)
  {
/*  printf("%f, %f %f\n",angle ,x+cos(angle)*r,y+sin(angle)*r);*/
// printf("%f, %f %f\n",angle ,xv+cos(angle)*ov,yv+sin(angle)*ov);
  makeBody(planetmass   ,5 ,x+cos(angle*e)*r, y+sin(angle*e)*r, z, xv+cos(angle*e - pi/2)*ov, yv+sin(angle*e - pi/2)*ov, zv,0,0,1);

  }
}


int mygeti(int *result)
{
char c, buff [ 13 ]; /* signed 32-bit value, extra room for '\n' and '\0' */
return fgets(buff, sizeof buff, stdin) && !isspace(*buff) &&
sscanf(buff, "%d%c", result, &c) == 2 && (c == '\n' || c == '\0');
}



int read()
{
  int x;
      do {
fputs("Enter value: ", stdout);
fflush(stdout);

} while ( !mygeti(&x) );
return x;
}

int main(int argc,char* argv[])
{  
  int totalbodies=0;
  int x;
srand(time(NULL));
int numringsystems;

printf("Input ring systems? (y/n)\n");
  char c=getchar();
  if(c=='y')
  {
printf("Number of ring systems: \n");
numringsystems=read();
totalbodies+=numringsystems;


int ringx[numringsystems];
int ringy[numringsystems];
int ringradius[numringsystems];
int ringxv[numringsystems];
int ringyv[numringsystems];
int ringzv[numringsystems];
int planetsinrings[numringsystems];

for(x=0;x<numringsystems;x++)
{
  
printf("Enter x position of system %d: \n",x+1);
double xpos=(double)read();

printf("Enter y position of system %d: \n",x+1);
double ypos=(double)read();

printf("Enter z position of system %d: \n",x+1);
double zpos=(double)read();

printf("Enter x velocity of system %d: \n",x+1);
double xv=(double)read();

printf("Enter y velocity of system %d: \n",x+1);
double yv=(double)read();

printf("Enter z velocity of system %d: \n",x+1);
double zv=(double)read();

printf("Enter planet distance from the sun of system %d: \n",x+1);
double r=(double)read();

printf("Enter number of planets in system %d: \n",x+1);
int numplanets=read();

ringsystem(xpos, ypos, zpos,  xv,  yv,  zv,  r, numplanets);

//  printf("value = %d\n", numringsystems);
totalbodies+=numplanets;
}
}
//      ringsystem(-200,-200,0,10,0,0,200,8);
//      ringsystem(200,200,  0,0,-10,0,200,8);
//      ringsystem(0,0,  0,0,0,0,100,10);
//   makeBody(10000,50, 0, 0, 0, 0, 0, 0,1,222.0/255.0,0);
//   makeBody(100   ,1 ,100, 0, 0, 0, 7, 0,0,0,1);
//   makeBody(1    ,.5 ,102, 0, 0, 0, 6, 0,50.0/255.0,50.0/255.0,50.0/255.0);
  printf("Input single planets? (y/n):\n");
  c=getchar();
  if(c=='y')
  {
  printf("Enter number of single planets: \n");
  int planets=read();
  
  for(x=0;x<planets;x++)
  {
    singlebody();
   totalbodies++; 
  }
  }
  printf("Input asteroids? (y/n):\n");
  c=getchar();
  if(c=='y')
  {
     printf("Enter number of asteroids: \n");
     int asteroids=read();
  for(x=0;x<asteroids;x++)
  {
    
    asteroid();
   totalbodies++; 
  }

  }
  num=totalbodies;

//   makeBody(3,1    ,.5,61, 0, 0, 0, 4, 0);
//   makeBody(4,1000, 8, 100, 0, 0, 0, 2, 0);

// randombody();
  
  
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
