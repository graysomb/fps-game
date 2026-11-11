/* Menu edits a draft; only Apply changes active controls and the saved file. */
static CtLibrary ctDraft;
static int ctTab,ctRow,ctScroll,ctEditSeat,ctDevice,ctBindingSlot,ctCapture=-1,ctConflict=-1,ctConflictChoice;
static CtBinding ctCandidate;
static bool ctCaptureArmed,ctTyping,ctDeleteConfirm;
static int ctNameAction,ctCharacter;
static char ctNameBuffer[CT_NAME];
static Vector2 ctPreviewAim,ctPreviewVelocity,ctPreviewSmooth;
static float ctMenuRepeat;
static const char *ctTabs[]={"Profiles","Bindings","Mouse","Controller","Keyboard Look"};
static const char *ctDevices[]={"Keyboard","Mouse","Controller"};
static CtProfile *ct_edit_profile(void){return &ctDraft.profiles[ctDraft.assigned[ctEditSeat]];}
static void ct_open(void){if(ctMenuOpen)return;ct_cancel_all();debugConsoleOpen=false;ctDraft=ctLibrary;ctMenuOpen=true;ctCapture=-1;ctConflict=-1;ctTyping=false;ctDeleteConfirm=false;ctRow=ctScroll=0;ctPreviewAim=ctPreviewVelocity=ctPreviewSmooth=(Vector2){0};EnableCursor();mouseCaptured=false;mouseButtonsArmed=false;}
static void ct_close(bool apply){if(apply){if(!ct_save_file("controls.cfg",&ctDraft))return;ctLibrary=ctDraft;}ctMenuOpen=false;couchSettings=false;ct_cancel_all();}
static bool ct_equal_binding(CtBinding a,CtBinding b){return a.kind&&a.kind==b.kind&&a.code==b.code&&a.mods==b.mods;}
static bool ct_overlap(int a,int b){int ca=ctDefinitions[a].context,cb=ctDefinitions[b].context;return (ca&cb)||((ca==CT_GAME&&cb==CT_LAB)||(ca==CT_LAB&&cb==CT_GAME))||(ca==CT_COMMON&&cb!=CT_PREP)||(cb==CT_COMMON&&ca!=CT_PREP);}
static int ct_conflict(CtProfile *p,int action,CtBinding binding){for(int a=0;a<CT_COUNT;a++)if(a!=action&&ct_overlap(a,action))for(int b=0;b<2;b++)if(ct_equal_binding(binding,p->bindings[ctDevice][a][b]))return a;return -1;}
static void ct_assign_binding(int mode){CtProfile *p=ct_edit_profile();CtBinding old=p->bindings[ctDevice][ctCapture][ctBindingSlot];
    int matches=0;for(int a=0;a<CT_COUNT;a++)if(a!=ctCapture&&ct_overlap(a,ctCapture))for(int b=0;b<2;b++)matches+=ct_equal_binding(ctCandidate,p->bindings[ctDevice][a][b]);
    if(mode==0&&matches>1){snprintf(ctError,sizeof(ctError),"Multiple conflicts: choose Replace or Cancel");return;}
    for(int a=0;a<CT_COUNT;a++)if(a!=ctCapture&&ct_overlap(a,ctCapture))for(int b=0;b<2;b++)if(ct_equal_binding(ctCandidate,p->bindings[ctDevice][a][b]))p->bindings[ctDevice][a][b]=mode==0?old:(CtBinding){0};
    p->bindings[ctDevice][ctCapture][ctBindingSlot]=ctCandidate;ctConflict=ctCapture=-1;
}
static void ct_capture_candidate(CtBinding b){
    if((b.kind==1&&(b.code==KEY_ESCAPE||b.code==KEY_PAUSE))||(b.kind==4&&b.code==GAMEPAD_BUTTON_MIDDLE_RIGHT)){snprintf(ctError,sizeof(ctError),"Escape, Pause and Start remain recovery controls");return;}
    if(b.kind==3&&!ct_discrete(ctCapture)){snprintf(ctError,sizeof(ctError),"This action needs a held button, not the wheel");return;}
    ctCandidate=b;ctConflict=ct_conflict(ct_edit_profile(),ctCapture,b);ctConflictChoice=0;if(ctConflict<0)ct_assign_binding(1);
}
static void ct_name_begin(int action){ctTyping=true;ctNameAction=action;ctCharacter=0;snprintf(ctNameBuffer,sizeof(ctNameBuffer),"%s",action==2?ct_edit_profile()->name:"");}
static void ct_name_finish(void){
    if(!ctNameBuffer[0])return;for(int i=0;i<ctDraft.count;i++)if((ctNameAction!=2||i!=ctDraft.assigned[ctEditSeat])&&!strcmp(ctDraft.profiles[i].name,ctNameBuffer)){snprintf(ctError,sizeof(ctError),"That profile name already exists");return;}
    if(ctNameAction==2)snprintf(ct_edit_profile()->name,CT_NAME,"%s",ctNameBuffer);
    else if(ctDraft.count<CT_PROFILES){int n=ctDraft.count++;if(ctNameAction==1)ctDraft.profiles[n]=*ct_edit_profile();else ct_default_profile(&ctDraft.profiles[n],0,ctNameBuffer);strcpy(ctDraft.profiles[n].name,ctNameBuffer);ctDraft.assigned[ctEditSeat]=n;}
    else {snprintf(ctError,sizeof(ctError),"Maximum 16 profiles; delete an unused profile first");return;}ctTyping=false;
}
static int ct_tuning_rows(int device,int *fields){
    static const int mouse[]={0,1,2,3,10,100,101,104,105};
    static const int pad[]={0,1,2,3,4,5,6,11,12,7,8,9,102,103,100,101,104,105,106};
    static const int key[]={0,1,2,3,4,102,100,101,104,105,106};
    const int *list=device==1?mouse:device==2?pad:key;int count=device==1?9:device==2?19:11;memcpy(fields,list,count*sizeof(int));fields[count++]=200;return count;
}
static int *ct_option(CtProfile *p,int d,int field){switch(field){case 100:return &p->tuning[d].invertX;case 101:return &p->tuning[d].invertY;case 102:return &p->tuning[d].direct;case 103:return &p->tuning[d].swap;case 104:return &p->aimToggle[d];case 105:return &p->tetherToggle[d];default:return &p->combined[d];}}
static const char *ct_option_name(int field){switch(field){case 100:return "Invert horizontal";case 101:return "Invert vertical";case 102:return "Immediate turning (no ramp)";case 103:return "Swap movement / look sticks";case 104:return "Aim: toggle (off = hold)";case 105:return "Tether: toggle (off = hold)";default:return "Combined tap switch / hold interact";}}
static bool ct_shared_keys(void){for(int s=0;s<4;s++)if(s!=ctEditSeat){CtProfile *p=&ctDraft.profiles[ctDraft.assigned[s]],*q=ct_edit_profile();for(int a=0;a<CT_COUNT;a++)for(int b=0;b<CT_COUNT;b++)if(ct_overlap(a,b)&&ct_equal_binding(p->bindings[0][a][0],q->bindings[0][b][0]))return true;}return false;}
static void ct_preview(int x,int y,int w,int h,float dt){
    CtProfile *p=ct_edit_profile();int d=ctTab==2?1:ctTab==3?2:0;CtTuning *t=&p->tuning[d];Vector2 raw={0},shaped={0};int device=couch_menu_device();
    if(device>=0){raw=(Vector2){ctRaw.axes[device][GAMEPAD_AXIS_RIGHT_X],ctRaw.axes[device][GAMEPAD_AXIS_RIGHT_Y]};shaped=(Vector2){ct_axis(raw.x,t),ct_axis(raw.y,t)};}
    DrawRectangle(x,y,w,h,(Color){10,16,23,255});DrawRectangleLines(x,y,w,h,(Color){65,99,116,255});DrawText("AIM PREVIEW (no gameplay)",x+12,y+10,16,SKYBLUE);
    Vector2 mid={x+w*.5f,y+h*.52f};float radius=fminf(w*.3f,h*.28f);DrawCircleLines(mid.x,mid.y,radius,GRAY);DrawCircleLines(mid.x,mid.y,radius*t->inner,DARKGRAY);
    DrawCircleV(Vector2Add(mid,Vector2Scale(raw,radius)),5,GRAY);DrawCircleV(Vector2Add(mid,Vector2Scale(shaped,radius)),5,SKYBLUE);
    bool preview=CheckCollisionPointRec(GetMousePosition(),(Rectangle){x,y,w,h})&&IsMouseButtonDown(MOUSE_BUTTON_LEFT);
    if(d==1&&preview){float gain=1+t->accel*fminf(Vector2Length(ctRaw.delta)/fmaxf(dt,.001f)/1000,2);Vector2 delta=ctRaw.delta;if(t->smooth>0){ctPreviewSmooth=Vector2Lerp(ctPreviewSmooth,Vector2Scale(delta,1/fmaxf(dt,.001f)),1-expf(-dt/t->smooth));delta=Vector2Scale(ctPreviewSmooth,dt);}ctPreviewAim.x+=delta.x*.15f*t->x*gain*(t->invertX?-1:1);ctPreviewAim.y+=delta.y*.15f*t->y*gain*(t->invertY?-1:1);}
    if(d==0){int actions[]={CT_LOOK_RIGHT,CT_LOOK_LEFT,CT_LOOK_DOWN,CT_LOOK_UP};bool down[4]={0};for(int a=0;a<4;a++)for(int b=0;b<2;b++)down[a]|=ct_binding_down(p->bindings[0][actions[a]][b],ctEditSeat,-1);shaped=(Vector2){down[0]-down[1],down[2]-down[3]};}
    if(d!=1){ctPreviewAim.x+=ct_turn(shaped.x*(t->invertX?-1:1),&ctPreviewVelocity.x,90*t->x,400*t->x*t->accel,400*t->brake,dt,t->direct||t->accel<=0);ctPreviewAim.y+=ct_turn(shaped.y*(t->invertY?-1:1),&ctPreviewVelocity.y,90*t->y,400*t->y*t->accel,400*t->brake,dt,t->direct||t->accel<=0);}
    ctPreviewAim.x=Clamp(ctPreviewAim.x,-w*.4f,w*.4f);ctPreviewAim.y=Clamp(ctPreviewAim.y,-h*.3f,h*.3f);
    Vector2 at=Vector2Add(mid,ctPreviewAim);DrawLine(at.x-8,at.y,at.x+8,at.y,GOLD);DrawLine(at.x,at.y-8,at.x,at.y+8,GOLD);
    DrawText(d==1?"Drag here to test mouse":d==0?"Use your keyboard look bindings":"Right stick: raw gray / filtered blue",x+12,y+h-26,14,LIGHTGRAY);
}
static void ct_menu_text(const char *text,int x,int y,int size,int width,Color color){char shown[512];snprintf(shown,sizeof(shown),"%s",text);int n=(int)strlen(shown);bool truncated=false;while(n>3&&MeasureText(shown,size)>width){shown[--n]=0;truncated=true;}if(truncated&&n>3){shown[n-3]='.';shown[n-2]='.';shown[n-1]='.';}DrawText(shown,x,y,size,color);}
static bool ct_menu_frame(void){
    if(!ctMenuOpen)return false;float dt=GetFrameTime();ctMenuRepeat-=dt;int device=couch_menu_device();
    int vertical=(IsKeyPressed(KEY_DOWN)?1:0)-(IsKeyPressed(KEY_UP)?1:0),horizontal=(IsKeyPressed(KEY_RIGHT)?1:0)-(IsKeyPressed(KEY_LEFT)?1:0);
    if(device>=0){vertical+=(IsGamepadButtonPressed(device,GAMEPAD_BUTTON_LEFT_FACE_DOWN)?1:0)-(IsGamepadButtonPressed(device,GAMEPAD_BUTTON_LEFT_FACE_UP)?1:0);horizontal+=(IsGamepadButtonPressed(device,GAMEPAD_BUTTON_LEFT_FACE_RIGHT)?1:0)-(IsGamepadButtonPressed(device,GAMEPAD_BUTTON_LEFT_FACE_LEFT)?1:0);float axis=GetGamepadAxisMovement(device,GAMEPAD_AXIS_LEFT_Y);if(fabsf(axis)>.6f&&ctMenuRepeat<=0){vertical+=axis>0?1:-1;ctMenuRepeat=.18f;}}
    bool confirm=IsKeyPressed(KEY_ENTER)||(device>=0&&IsGamepadButtonPressed(device,GAMEPAD_BUTTON_RIGHT_FACE_DOWN));bool back=IsKeyPressed(KEY_ESCAPE)||(device>=0&&IsGamepadButtonPressed(device,GAMEPAD_BUTTON_RIGHT_FACE_RIGHT));
    bool click=IsMouseButtonPressed(MOUSE_BUTTON_LEFT);Vector2 mouse=GetMousePosition();int width=GetScreenWidth(),height=GetScreenHeight();bool compact=width<900;int listWidth=compact?width-60:(int)(width*.62f);int left=30,top=compact?150:170,rowH=compact?26:32,visible=(height-top-115)/rowH;if(visible<3)visible=3;
    if(ctTyping){
        int c;while((c=GetCharPressed())>0)if(c>=32&&c<127&&strlen(ctNameBuffer)<CT_NAME-1){size_t n=strlen(ctNameBuffer);ctNameBuffer[n]=(char)c;ctNameBuffer[n+1]=0;}
        if(IsKeyPressed(KEY_BACKSPACE)||(device>=0&&IsGamepadButtonPressed(device,GAMEPAD_BUTTON_RIGHT_FACE_LEFT))){size_t n=strlen(ctNameBuffer);if(n)ctNameBuffer[n-1]=0;}
        const char *alphabet="ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789 _-";int len=(int)strlen(alphabet);ctCharacter=(ctCharacter+horizontal+vertical*8+len*2)%len;
        if(device>=0&&IsGamepadButtonPressed(device,GAMEPAD_BUTTON_RIGHT_FACE_DOWN)&&strlen(ctNameBuffer)<CT_NAME-1){size_t n=strlen(ctNameBuffer);ctNameBuffer[n]=alphabet[ctCharacter];ctNameBuffer[n+1]=0;}
        if(IsKeyPressed(KEY_ENTER)||(device>=0&&IsGamepadButtonPressed(device,GAMEPAD_BUTTON_MIDDLE_RIGHT)))ct_name_finish();if(back)ctTyping=false;
        BeginDrawing();ClearBackground((Color){16,21,30,255});DrawText("PROFILE NAME",40,40,30,RAYWHITE);DrawText(ctNameBuffer,40,110,28,GOLD);DrawText(TextFormat("Controller character: [ %c ]",alphabet[ctCharacter]),40,180,24,SKYBLUE);DrawText("Type a name | D-pad chooses character | A adds | X deletes | Start saves",40,240,18,LIGHTGRAY);DrawText(ctError,40,290,18,ORANGE);EndDrawing();return true;
    }
    if(ctCapture>=0){
        if(IsKeyPressed(KEY_ESCAPE)||(back&&(ctDevice!=2||ctConflict>=0))){ctCapture=ctConflict=-1;}
        else if(ctConflict>=0){ctConflictChoice=(ctConflictChoice+horizontal+3)%3;if(confirm){if(ctConflictChoice==2)ctCapture=ctConflict=-1;else ct_assign_binding(ctConflictChoice);}}
        else {
            bool any=false;if(ctDevice==0){for(int k=1;k<512;k++)any|=ctRaw.keys[k];}else if(ctDevice==1){for(int k=0;k<7;k++)any|=ctRaw.mouse[k];}else if(device>=0){for(int k=1;k<32;k++)any|=ctRaw.pad[device][k];}
            if(!any)ctCaptureArmed=true;
            if(ctCaptureArmed){int mods=(ctRaw.keys[KEY_LEFT_CONTROL]||ctRaw.keys[KEY_RIGHT_CONTROL]?1:0)|(ctRaw.keys[KEY_LEFT_SHIFT]||ctRaw.keys[KEY_RIGHT_SHIFT]?2:0)|(ctRaw.keys[KEY_LEFT_ALT]||ctRaw.keys[KEY_RIGHT_ALT]?4:0);
                if(ctDevice==0){for(int k=1;k<512;k++)if(IsKeyPressed(k)&&k!=KEY_ESCAPE&&k!=KEY_LEFT_CONTROL&&k!=KEY_RIGHT_CONTROL&&k!=KEY_LEFT_SHIFT&&k!=KEY_RIGHT_SHIFT&&k!=KEY_LEFT_ALT&&k!=KEY_RIGHT_ALT){ct_capture_candidate((CtBinding){1,k,mods});break;}if(ctCapture>=0&&ctConflict<0)for(int k=KEY_LEFT_SHIFT;k<=KEY_RIGHT_SUPER;k++)if(IsKeyReleased(k)&&(k==KEY_LEFT_SHIFT||k==KEY_RIGHT_SHIFT||k==KEY_LEFT_CONTROL||k==KEY_RIGHT_CONTROL||k==KEY_LEFT_ALT||k==KEY_RIGHT_ALT)){ct_capture_candidate((CtBinding){1,k,0});break;}}
                if(ctDevice==1){for(int k=0;k<7;k++)if(IsMouseButtonPressed(k)){ct_capture_candidate((CtBinding){2,k,0});break;}if(ctCapture>=0&&ctRaw.wheel)ct_capture_candidate((CtBinding){3,ctRaw.wheel>0?1:-1,0});}
                if(ctDevice==2&&device>=0)for(int k=1;k<=GAMEPAD_BUTTON_RIGHT_THUMB;k++)if(IsGamepadButtonPressed(device,k)){ct_capture_candidate((CtBinding){4,k,0});break;}
            }
        }
        BeginDrawing();ClearBackground((Color){16,21,30,255});DrawText(ctConflict>=0?"BINDING CONFLICT":"PRESS A CONTROL",40,40,30,RAYWHITE);
        if(ctCapture>=0)DrawText(ctDefinitions[ctCapture].label,40,110,24,GOLD);
        if(ctConflict>=0){DrawText(TextFormat("%s is used by %s",ct_binding_name(ctCandidate),ctDefinitions[ctConflict].label),40,155,20,LIGHTGRAY);DrawText(TextFormat("%s Swap     %s Replace     %s Cancel",ctConflictChoice==0?">":" ",ctConflictChoice==1?">":" ",ctConflictChoice==2?">":" "),40,210,24,SKYBLUE);}else DrawText(TextFormat("Device: %s | Escape cancels",ctDevices[ctDevice]),40,170,20,LIGHTGRAY);
        DrawText(ctError,40,270,18,ORANGE);EndDrawing();return true;
    }
    if(ctDeleteConfirm){if(confirm&&ctDraft.count>1){int n=ctDraft.assigned[ctEditSeat];memmove(&ctDraft.profiles[n],&ctDraft.profiles[n+1],(ctDraft.count-n-1)*sizeof(CtProfile));ctDraft.count--;for(int s=0;s<4;s++){if(ctDraft.assigned[s]==n)ctDraft.assigned[s]=0;else if(ctDraft.assigned[s]>n)ctDraft.assigned[s]--;}ctDeleteConfirm=false;}if(back)ctDeleteConfirm=false;BeginDrawing();ClearBackground((Color){16,21,30,255});DrawText("Delete this profile? Assigned seats will use the first profile.",30,100,22,RAYWHITE);DrawText("Enter / A: delete    Escape / B: cancel",30,160,22,GOLD);EndDrawing();return true;}
    int tabStep=(IsKeyPressed(KEY_TAB)?1:0);if(device>=0)tabStep+=IsGamepadButtonPressed(device,GAMEPAD_BUTTON_RIGHT_TRIGGER_1)-IsGamepadButtonPressed(device,GAMEPAD_BUTTON_LEFT_TRIGGER_1);
    if(tabStep){ctTab=(ctTab+tabStep+5)%5;ctRow=ctScroll=0;}
    for(int t=0;t<5;t++)if(click&&CheckCollisionPointRec(mouse,(Rectangle){left+t*(width-60)/5.f,82,(width-60)/5.f,38})){ctTab=t;ctRow=ctScroll=0;click=false;}
    if(back){ct_close(false);return true;}
    CtProfile *p=ct_edit_profile();int fields[32],count=ctTab==0?9:ctTab==1?CT_COUNT+3:ct_tuning_rows(ctTab==2?1:ctTab==3?2:0,fields);
    ctRow=(ctRow+vertical+count)%count;if(ctRaw.wheel&&mouse.x<width*.65f)ctRow=Clamp(ctRow-(int)ctRaw.wheel,0,count-1);
    if(ctRow<ctScroll)ctScroll=ctRow;if(ctRow>=ctScroll+visible)ctScroll=ctRow-visible+1;
    if(click&&mouse.x>=left&&mouse.x<left+listWidth&&mouse.y>=top&&mouse.y<top+visible*rowH){int row=ctScroll+(int)((mouse.y-top)/rowH);if(row<count){ctRow=row;confirm=true;}}
    if(ctTab==0){if(ctRow==0&&horizontal)ctEditSeat=(ctEditSeat+horizontal+4)%4;if(ctRow==1&&(horizontal||confirm))ctDraft.assigned[ctEditSeat]=(ctDraft.assigned[ctEditSeat]+(horizontal?horizontal:1)+ctDraft.count)%ctDraft.count;
        if(confirm){if(ctRow==2)ct_name_begin(0);if(ctRow==3)ct_name_begin(1);if(ctRow==4)ct_name_begin(2);if(ctRow==5&&ctDraft.count>1)ctDeleteConfirm=true;if(ctRow==6||ctRow==7){char name[CT_NAME];strcpy(name,p->name);ct_default_profile(p,ctRow==7,name);}if(ctRow==8){CtLibrary loaded;if(ct_read_file("controls.cfg",&loaded)){ctDraft=loaded;ctError[0]=0;}else snprintf(ctError,sizeof(ctError),"Could not load controls.cfg; current draft kept");}}
    }else if(ctTab==1){if(ctRow==0&&(horizontal||confirm))ctDevice=(ctDevice+(horizontal?horizontal:1)+3)%3;else if(ctRow==1&&(horizontal||confirm))ctBindingSlot=1-ctBindingSlot;else if(ctRow==2&&confirm){CtProfile defaults;ct_default_profile(&defaults,p->layout,p->name);memcpy(p->bindings[ctDevice],defaults.bindings[ctDevice],sizeof(p->bindings[ctDevice]));}else if(ctRow>=3){int action=ctRow-3;if(IsKeyPressed(KEY_DELETE)||(device>=0&&IsGamepadButtonPressed(device,GAMEPAD_BUTTON_RIGHT_FACE_LEFT)))p->bindings[ctDevice][action][ctBindingSlot]=(CtBinding){0};if(confirm){ctCapture=action;ctCaptureArmed=false;ctError[0]=0;}}}
    else {int d=ctTab==2?1:ctTab==3?2:0,field=fields[ctRow];if(horizontal||confirm){int direction=horizontal?horizontal:1;if(field==200){CtProfile defaults;ct_default_profile(&defaults,p->layout,p->name);p->tuning[d]=defaults.tuning[d];p->aimToggle[d]=defaults.aimToggle[d];p->tetherToggle[d]=defaults.tetherToggle[d];p->combined[d]=defaults.combined[d];}else if(field<100){float *value=ct_value(&p->tuning[d],field);*value=Clamp(*value+direction*ctFields[field].step,ctFields[field].min,ctFields[field].max);}else{int *value=ct_option(p,d,field);*value=!*value;}}}
    bool apply=IsKeyPressed(KEY_F5)||(device>=0&&IsGamepadButtonPressed(device,GAMEPAD_BUTTON_MIDDLE_RIGHT));
    Rectangle applyBox={30,height-64,140,36},cancelBox={180,height-64,140,36},resetBox={330,height-64,170,36};
    if(click&&CheckCollisionPointRec(mouse,applyBox))apply=true;if(click&&CheckCollisionPointRec(mouse,cancelBox)){ct_close(false);return true;}
    if((click&&CheckCollisionPointRec(mouse,resetBox))||IsKeyPressed(KEY_F6)){CtProfile defaults;ct_default_profile(&defaults,p->layout,p->name);if(ctTab>=2){int d=ctTab==2?1:ctTab==3?2:0;p->tuning[d]=defaults.tuning[d];p->aimToggle[d]=defaults.aimToggle[d];p->tetherToggle[d]=defaults.tetherToggle[d];p->combined[d]=defaults.combined[d];}else if(ctTab==1)memcpy(p->bindings[ctDevice],defaults.bindings[ctDevice],sizeof(p->bindings[ctDevice]));}
    if(apply){ct_close(true);if(!ctMenuOpen)return true;}
    BeginDrawing();ClearBackground((Color){16,21,30,255});DrawText("CONTROLS & SETTINGS",30,24,30,RAYWHITE);DrawText(TextFormat("Editing P%d: %s | Mouse belongs to %s",ctEditSeat+1,ct_edit_profile()->name,mouseSeat>=0?TextFormat("P%d",mouseSeat+1):"first keyboard seat"),30,59,17,SKYBLUE);
    for(int t=0;t<5;t++){int x=left+t*(width-60)/5;DrawRectangle(x,82,(width-60)/5-4,38,t==ctTab?(Color){40,77,94,255}:(Color){25,34,43,255});ct_menu_text(compact&&t==4?"Keyboard":ctTabs[t],x+8,92,compact?13:17,(width-60)/5-18,t==ctTab?GOLD:LIGHTGRAY);}
    DrawText(ctTab==1?"Enter/A: bind | Delete/X: clear | Left/right: device or slot":"Arrows / D-pad: select and adjust | Tab / LB / RB: tabs",30,134,compact?12:16,LIGHTGRAY);
    p=ct_edit_profile();for(int n=ctScroll;n<count&&n<ctScroll+visible;n++){const char *label;
        if(ctTab==0){switch(n){case 0:label=TextFormat("Local player: %d",ctEditSeat+1);break;case 1:label=TextFormat("Assigned profile: %s",p->name);break;case 2:label="Create profile";break;case 3:label="Duplicate profile";break;case 4:label="Rename profile";break;case 5:label="Delete profile";break;case 6:label="Reset profile: first keyboard defaults";break;case 7:label="Reset profile: second keyboard defaults";break;default:label="Reload saved controls";}}
        else if(ctTab==1){if(n==0)label=TextFormat("Device: %s",ctDevices[ctDevice]);else if(n==1)label=TextFormat("Binding slot: %s",ctBindingSlot?"Secondary":"Primary");else if(n==2)label="Reset this device's bindings";else{int a=n-3;const char *group=ctDefinitions[a].context==CT_EDIT?"Editor":ctDefinitions[a].context==CT_LAB?"Lab":ctDefinitions[a].context==CT_PREP?"Prep":"Play";label=TextFormat("[%s] %s: %s",group,ctDefinitions[a].label,ct_binding_name(p->bindings[ctDevice][a][ctBindingSlot]));}}
        else {int d=ctTab==2?1:ctTab==3?2:0,f=fields[n];label=f==200?"Reset this tab":f<100?TextFormat("%s: %.2f",ctFields[f].label,*ct_value(&p->tuning[d],f)):TextFormat("%s: %s",ct_option_name(f),*ct_option(p,d,f)?"On":"Off");}
        int y=top+(n-ctScroll)*rowH;if(n==ctRow)DrawRectangle(left-6,y-3,listWidth,rowH-2,(Color){40,77,94,255});ct_menu_text(label,left,y+3,compact?14:17,listWidth-12,n==ctRow?GOLD:RAYWHITE);
    }
    if(!compact&&ctTab>=2)ct_preview(width*.67f,170,width*.30f,height-300,dt);
    else if(!compact) {DrawText("Profiles are local to this machine.",width*.67f,180,15,LIGHTGRAY);DrawText("Apply saves to controls.cfg.",width*.67f,207,15,LIGHTGRAY);DrawText("Escape / Start recovery stays available.",width*.67f,234,14,LIGHTGRAY);if(ctTab==1&&ctDevice==0&&ct_shared_keys())DrawText("Shared keyboard overlap:\ncheck the other player's bindings.",width*.67f,275,15,ORANGE);}
    if(count>visible)DrawText(TextFormat("%d-%d of %d | scroll / D-pad",ctScroll+1,(int)fminf(count,ctScroll+visible),count),30,height-113,compact?12:14,GRAY);
    DrawText(ctError,30,height-91,16,ORANGE);DrawRectangleRec(applyBox,(Color){40,77,94,255});DrawText("Apply (F5)",40,height-56,18,WHITE);DrawRectangleRec(cancelBox,(Color){40,50,64,255});DrawText("Cancel / Back",190,height-56,18,WHITE);DrawRectangleRec(resetBox,(Color){40,50,64,255});DrawText("Reset tab (F6)",340,height-56,18,WHITE);EndDrawing();return true;
}
