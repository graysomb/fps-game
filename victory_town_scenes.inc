/* Shared montage metadata: renderer, controls and capture manifest use this table. */
/* Authored animation seconds remain unchanged; playback stretches them to the song.
 * MP3 duration measured from its decoded sample count (44,100 Hz). */
#define TOWN_MUSIC_PATH "assets/audio/pretty-good.mp3"
#define TOWN_MUSIC_SECONDS 176.5355f
static float town_playback_scale(void);
enum TownSceneId { SC_POWER,SC_HOMES,SC_MARKET,SC_GARDEN,SC_PLAYGROUND,SC_SUPPER,SC_OBSERVATORY,SC_WIDE,SC_COUNT };
typedef struct {const char *id,*title;float duration;Camera3D camera;bool oneWay;} TownSceneDef;
static const TownSceneDef townScenes[SC_COUNT]={
 {"power","Restoring the town's power",8,{{-16,3.1f,-33},{-21,1.3f,-28},{0,1,0},57,0},true},
 {"homes","Homes in the Foundry",8,{{23.5f,3.5f,10.8f},{21.2f,1.25f,14.8f},{0,1,0},63,0},false},
 {"market","Market beneath the blast shelter",8,{{-16.5f,2.6f,-4},{-21.5f,1.1f,0},{0,1,0},56,0},true},
 {"garden-workshop","Garden and workshop",8,{{21.5f,5.4f,13.3f},{15.7f,.8f,9.5f},{0,1,0},60,0},false},
 {"playground","A war machine becomes a playground",8,{{29,3.1f,-27},{32,1.1f,-31},{0,1,0},56,0},false},
 {"supper","Supper on the platform",8,{{4.1f,3.6f,-4.8f},{0,2.65f,0},{0,1,0},48,0},false},
 {"observatory","The rooftop observatory",8,{{27,8,16},{21,7.1f,22},{0,1,0},58,0},false},
 {"town-under-stars","A place to call home",8,{{39,17,-38},{0,6,0},{0,1,0},57,0},false}
};
static float town_duration(void){float t=0;for(int i=0;i<SC_COUNT;i++)t+=townScenes[i].duration;return t;}
static float town_playback_scale(void){return TOWN_MUSIC_SECONDS/town_duration();}
static int town_timeline(float time,float *local){float t=fmodf(fmaxf(time,0),town_duration());for(int i=0;i<SC_COUNT;i++){if(t<townScenes[i].duration){*local=t;return i;}t-=townScenes[i].duration;}*local=0;return 0;}
enum StoryRole { R_WHEEL,R_CABLE,R_READER,R_VENDOR,R_BUYER,R_SWING,R_CLIMB,R_PEEK,R_ASTRONOMER,R_CHART };
typedef struct {int role;CharacterPose pose;Vector3 floor,hands[2],shoulders[2];float scale;} StoryActor;
typedef struct {StoryActor actors[4];int count;float power,wheel,swingAngle,telescopeAngle;Vector3 swingSeat,bread,scrap;int breadOwner,scrapOwner;Vector4 meteor;float meteorGlow;} TownStoryFrame;
