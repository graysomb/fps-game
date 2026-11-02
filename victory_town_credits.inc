#ifndef VICTORY_TOWN_CREDITS_INCLUDED
#define VICTORY_TOWN_CREDITS_INCLUDED
/* Presentation-only credits, drawn after scene composition and bloom.
 * The authored timeline also drives music playback, pause and seeking. */
typedef struct { const char *heading,*lines[3]; } TownCredit;
static const TownCredit townCredits[]={
 {"Game by",{"Michael",NULL,NULL}},
 {"Music by",{"Michael",NULL,NULL}},
 {"Slop by",{"chat",NULL,NULL}},
 {"Physics engine adapted from",{"Tim McGraw","Gram-Schmidt voxel constraints for","real-time destructible soft bodies."}},
 {"Play testing by",{"Puhstinkums, Aman, Guillaume,","Shane, Galen, Beck.",NULL}}
};
static void town_credits_draw(float authoredTime,int width,int height){
 const int count=(int)(sizeof(townCredits)/sizeof(townCredits[0]));
 float slot=town_duration()/count;
 float t=fmodf(fmaxf(authoredTime,0),town_duration());
 int index=(int)(t/slot);if(index>=count)index=count-1;
 float local=t-index*slot;
 /* Fade over about 1.1 real seconds at the soundtrack playback rate. */
 float alpha=Clamp(fminf(local,slot-local)/.4f,0,1);
 if(alpha<=0)return;
 float scale=fminf(width/1280.f,height/720.f);
 int x=(int)(48*scale),bottom=height-(int)(42*scale);
 const TownCredit *credit=&townCredits[index];
 int lines=0;while(lines<3&&credit->lines[lines])lines++;
 int lineHeight=(int)(35*scale),titleSize=(int)(20*scale),textSize=(int)(28*scale);
 int top=bottom-lines*lineHeight-(int)(38*scale);
 DrawRectangleGradientV(0,top-(int)(35*scale),width,height-top+(int)(35*scale),Fade(BLACK,0),Fade((Color){7,14,20,255},.87f*alpha));
 DrawRectangle(x,top,(int)(38*scale),(int)(3*scale),Fade((Color){227,188,116,255},alpha));
 DrawText(credit->heading,x+(int)(51*scale),top-(int)(8*scale),titleSize,Fade((Color){227,200,151,255},alpha));
 for(int i=0;i<lines;i++){
  int y=top+(int)(29*scale)+i*lineHeight;
  DrawText(credit->lines[i],x+1,y+2,textSize,Fade(BLACK,.75f*alpha));
  DrawText(credit->lines[i],x,y,textSize,Fade((Color){247,242,230,255},alpha));
 }
}

#endif
