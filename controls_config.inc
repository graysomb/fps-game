#include <stddef.h>
typedef struct {const char *name,*label;int context,key[2],pad,mouse;} CtDefinition;
#define CT_DEF(id,label,ctx,k1,k2,p,m) {#id,label,ctx,{k1,k2},p,m},
static const CtDefinition ctDefinitions[]={CT_ACTIONS(CT_DEF)};
#undef CT_DEF
typedef struct {const char *name,*label;size_t offset;float min,max,step;} CtField;
#define CTF(n,l,lo,hi,step) {#n,l,offsetof(CtTuning,n),lo,hi,step}
static const CtField ctFields[]={CTF(x,"Horizontal sensitivity",.1f,5,.1f),CTF(y,"Vertical sensitivity",.1f,5,.1f),CTF(zoom,"Zoom sensitivity",.1f,3,.1f),CTF(accel,"Acceleration",0,5,.1f),CTF(brake,"Stopping deceleration",.1f,5,.1f),CTF(inner,"Inner deadzone",0,.5f,.01f),CTF(outer,"Outer deadzone",0,.3f,.01f),CTF(curve,"Response curve (0 Classic, 1 Linear, 2 Expo)",0,2,1),CTF(trigger,"Trigger threshold",.05f,.95f,.05f),CTF(assist,"Aim assist strength",0,1,.1f),CTF(smooth,"Mouse smoothing (seconds)",0,.2f,.01f),CTF(moveInner,"Move stick inner deadzone",0,.5f,.01f),CTF(moveOuter,"Move stick outer deadzone",0,.3f,.01f)};
#undef CTF
static char ctError[192];
static float *ct_value(CtTuning *t,int field){return (float*)((unsigned char*)t+ctFields[field].offset);}
static CtProfile *ct_profile(int seat){if(seat<0||seat>=4)seat=0;return &ctLibrary.profiles[ctLibrary.assigned[seat]];}
static void ct_default_profile(CtProfile *p,int layout,const char *name){
    memset(p,0,sizeof(*p));p->layout=layout;snprintf(p->name,sizeof(p->name),"%s",name);
    for(int d=0;d<3;d++){p->tuning[d]=(CtTuning){.x=1,.y=1,.zoom=1,.accel=d==1?0:1,.brake=1,.inner=.15f,.moveInner=.15f,.trigger=.5f,.assist=1};p->aimToggle[d]=1;p->combined[d]=d==2;}
    for(int a=0;a<CT_COUNT;a++){
        int key=ctDefinitions[a].key[layout];if(key)p->bindings[0][a][0]=(CtBinding){1,key,0};
        if(ctDefinitions[a].mouse>=0)p->bindings[1][a][0]=(CtBinding){2,ctDefinitions[a].mouse,0};
        if(ctDefinitions[a].pad)p->bindings[2][a][0]=(CtBinding){4,ctDefinitions[a].pad,0};
    }
    p->bindings[0][CT_UNDO][0].mods=layout?0:1;p->bindings[0][CT_REDO][0].mods=1;
    p->bindings[1][CT_SWITCH][0]=(CtBinding){3,1,0};p->bindings[1][CT_SWITCH][1]=(CtBinding){3,-1,0};
    if(!layout){p->bindings[1][CT_DIST_UP][0]=(CtBinding){3,1,0};p->bindings[1][CT_DIST_DOWN][0]=(CtBinding){3,-1,0};}
}
static void ct_defaults(CtLibrary *l){memset(l,0,sizeof(*l));l->version=1;l->count=4;for(int i=0;i<4;i++){char name[48];snprintf(name,sizeof(name),"Player %d",i+1);ct_default_profile(&l->profiles[i],i==1,name);l->assigned[i]=i;}}
static bool ct_binding_valid(CtBinding b,int device){
    if(b.mods<0||b.mods>7)return false;if(!b.kind)return b.code==0;
    if(device==0)return b.kind==1&&b.code>0&&b.code<512;
    if(device==1)return (b.kind==2&&b.code>=0&&b.code<7)||(b.kind==3&&(b.code==1||b.code==-1));
    return b.kind==4&&b.code>0&&b.code<=GAMEPAD_BUTTON_RIGHT_THUMB;
}
static bool ct_discrete(int a){return a!=CT_FIRE&&a!=CT_TETHER&&a!=CT_BUILD&&a!=CT_REVIVE&&a!=CT_INTERACT&&a!=CT_PLACE&&a!=CT_ERASE&&a!=CT_FAST&&a!=CT_UP&&a!=CT_DOWN&&a>CT_LOOK_DOWN;}
static bool ct_validate(const CtLibrary *l){
    if(l->count<1||l->count>CT_PROFILES)return false;
    for(int i=0;i<4;i++)if(l->assigned[i]<0||l->assigned[i]>=l->count)return false;
    for(int i=0;i<l->count;i++){
        const CtProfile *p=&l->profiles[i];if(!p->name[0]||!memchr(p->name,0,sizeof(p->name))||strchr(p->name,'\n')||strchr(p->name,'\r')||p->layout<0||p->layout>1)return false;
        for(int d=0;d<3;d++){
            CtTuning t=p->tuning[d];for(int f=0;f<13;f++){float value=*ct_value(&t,f);if(!isfinite(value)||value<ctFields[f].min||value>ctFields[f].max)return false;}
            if(t.inner+t.outer>=.95f||t.moveInner+t.moveOuter>=.95f||t.invertX<0||t.invertX>1||t.invertY<0||t.invertY>1||t.direct<0||t.direct>1||t.swap<0||t.swap>1)return false;
            if(p->aimToggle[d]<0||p->aimToggle[d]>1||p->tetherToggle[d]<0||p->tetherToggle[d]>1||p->combined[d]<0||p->combined[d]>1)return false;
            for(int a=0;a<CT_COUNT;a++)for(int b=0;b<2;b++)if(!ct_binding_valid(p->bindings[d][a][b],d)||(p->bindings[d][a][b].kind==3&&!ct_discrete(a)))return false;
        }
    }return true;
}
static bool ct_save_file(const char *path,const CtLibrary *l){
    if(!ct_validate(l)){snprintf(ctError,sizeof(ctError),"Invalid control settings");return false;}
    char tmp[600];if(snprintf(tmp,sizeof(tmp),"%s.tmp",path)>=(int)sizeof(tmp))return false;FILE *f=fopen(tmp,"w");if(!f){snprintf(ctError,sizeof(ctError),"Cannot write %s",path);return false;}
    fprintf(f,"CONTROLS 1 %d\nSEATS %d %d %d %d\n",l->count,l->assigned[0],l->assigned[1],l->assigned[2],l->assigned[3]);
    for(int i=0;i<l->count;i++){const CtProfile *p=&l->profiles[i];fprintf(f,"PROFILE %d %d\nNAME %s\n",i,p->layout,p->name);
        for(int d=0;d<3;d++){CtTuning t=p->tuning[d];fprintf(f,"OPTIONS %d %d %d %d %d %d %d %d\n",d,p->aimToggle[d],p->tetherToggle[d],p->combined[d],t.invertX,t.invertY,t.direct,t.swap);
            for(int n=0;n<13;n++)fprintf(f,"TUNE %d %s %.9g\n",d,ctFields[n].name,*ct_value(&t,n));
            for(int a=0;a<CT_COUNT;a++)for(int b=0;b<2;b++){CtBinding v=p->bindings[d][a][b];fprintf(f,"BIND %d %s %d %d %d %d\n",d,ctDefinitions[a].name,b,v.kind,v.code,v.mods);}}
        fputs("ENDPROFILE\n",f);
    }fputs("END\n",f);bool ok=!ferror(f);if(fclose(f))ok=false;
#ifdef _WIN32
    if(ok)ok=MoveFileExA(tmp,path,MOVEFILE_REPLACE_EXISTING|MOVEFILE_WRITE_THROUGH)!=0;
#else
    if(ok)ok=rename(tmp,path)==0;
#endif
    if(!ok){remove(tmp);snprintf(ctError,sizeof(ctError),"Save failed; previous controls kept");}else ctError[0]=0;return ok;
}
static bool ct_read_file(const char *path,CtLibrary *out){
    FILE *f=fopen(path,"r");if(!f)return false;CtLibrary *l=calloc(1,sizeof(*l));if(!l){fclose(f);return false;}ct_defaults(l);
    char line[512],name[64];int version,count;bool ok=fgets(line,sizeof(line),f)&&sscanf(line,"CONTROLS %d %d",&version,&count)==2&&version==1&&count>0&&count<=CT_PROFILES;
    if(ok){l->count=count;for(int i=4;i<count;i++)ct_default_profile(&l->profiles[i],0,"Imported");}
    int profile=-1,seen=0;bool end=false;
    while(ok&&fgets(line,sizeof(line),f)){
        if(!strchr(line,'\n')&&!feof(f)){ok=false;break;}line[strcspn(line,"\r\n")]=0;
        if(!strcmp(line,"END")){end=true;break;}
        int d,a,b,c,e,g,h;float value;
        if(!strncmp(line,"SEATS ",6)){ok=sscanf(line+6,"%d %d %d %d",&l->assigned[0],&l->assigned[1],&l->assigned[2],&l->assigned[3])==4;continue;}
        if(!strncmp(line,"PROFILE ",8)){ok=sscanf(line+8,"%d %d",&a,&b)==2&&a==seen&&a<count&&b>=0&&b<2;if(ok){profile=a;ct_default_profile(&l->profiles[a],b,"Imported");seen++;}continue;}
        if(!strcmp(line,"ENDPROFILE")){profile=-1;continue;}if(profile<0){ok=false;break;}CtProfile *p=&l->profiles[profile];
        if(!strncmp(line,"NAME ",5)){ok=strlen(line+5)>0&&strlen(line+5)<CT_NAME;if(ok)strcpy(p->name,line+5);}
        else if(!strncmp(line,"OPTIONS ",8)){ok=sscanf(line+8,"%d %d %d %d %d %d %d %d",&d,&a,&b,&c,&e,&g,&h,&version)==8&&d>=0&&d<3;if(ok){p->aimToggle[d]=a;p->tetherToggle[d]=b;p->combined[d]=c;p->tuning[d].invertX=e;p->tuning[d].invertY=g;p->tuning[d].direct=h;p->tuning[d].swap=version;}}
        else if(!strncmp(line,"TUNE ",5)){ok=sscanf(line+5,"%d %63s %f",&d,name,&value)==3&&d>=0&&d<3;if(ok)for(int n=0;n<13;n++)if(!strcmp(name,ctFields[n].name))*ct_value(&p->tuning[d],n)=value;}
        else if(!strncmp(line,"BIND ",5)){ok=sscanf(line+5,"%d %63s %d %d %d %d",&d,name,&a,&b,&c,&e)==6&&d>=0&&d<3&&a>=0&&a<2;if(ok)for(int n=0;n<CT_COUNT;n++)if(!strcmp(name,ctDefinitions[n].name))p->bindings[d][n][a]=(CtBinding){b,c,e};}
        else if(line[0]&&line[0]!='#')ok=false;
    }
    ok=ok&&end&&seen==count&&ct_validate(l);fclose(f);if(ok)*out=*l;free(l);return ok;
}
static void ct_load(void){
    ct_defaults(&ctLibrary);ctReady=true;
    if(FileExists("controls.cfg")){if(!ct_read_file("controls.cfg",&ctLibrary))snprintf(ctError,sizeof(ctError),"controls.cfg could not be loaded; using safe defaults");return;}
    FILE *f=fopen("couch_settings.cfg","r");if(!f)return;
    for(int i=0;i<4;i++){float sensitivity,deadzone;int invert;if(fscanf(f,"%f %f %d",&sensitivity,&deadzone,&invert)!=3)break;
        if(isfinite(sensitivity)&&sensitivity>=.5f&&sensitivity<=2&&isfinite(deadzone)&&deadzone>=.05f&&deadzone<=.3f&&(invert==0||invert==1)){
            CtProfile *p=&ctLibrary.profiles[i];p->tuning[1].x=p->tuning[1].y=p->tuning[2].x=p->tuning[2].y=sensitivity;p->tuning[2].inner=p->tuning[2].moveInner=deadzone;p->tuning[1].invertY=p->tuning[2].invertY=invert;}}
    fclose(f);ct_save_file("controls.cfg",&ctLibrary);
}
