#ifndef IGLOO_BLUEPRINT_H
#define IGLOO_BLUEPRINT_H
#include <stdbool.h>
#include <stdlib.h>
enum { IGLOO_SIDE=11, IGLOO_HEIGHT=8, IGLOO_MAX_CELLS=512 };
#define IGLOO_CELL_SIZE .5f
/* y=0 is a flush foundation; y=1 starts above the support plane. */
static bool igloo_cell(int x,int y,int z,bool cutaway){
 if(x<0||x>=IGLOO_SIDE||z<0||z>=IGLOO_SIDE||y<0||y>=IGLOO_HEIGHT)return false;
 x-=5;z-=5;if(cutaway&&z<0)return false;
 if(y==0)return x*x+z*z<=25;
 if(y>=6)return abs(x)<=1&&abs(z)<=1&&(x||z);
 if(!x&&!z)return false;
 int radius[]={5,5,5,4,3,2};int r=radius[y];
 bool shell=x*x+z*z<=r*r && x*x+z*z>(r-2)*(r-2);
 if(y==5)shell=abs(x)<=2&&abs(z)<=2&&(x||z);
 if(z<0&&abs(x)<=1&&y<=3)return false;
 return shell;
}
#endif
