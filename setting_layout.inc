/* Decorative map data: coordinates survive voxel-array compaction. */
#define SETTING_DETAILS_MAX 256
typedef struct {
    Vector3 center; float width,height,yaw; int style,label;
    int cells[6][3]; Color colors[6]; bool visible;
} SettingDetail;
static SettingDetail settingDetails[SETTING_DETAILS_MAX];
static int settingDetailCount,settingMap,settingPendingMap;
static bool settingLayoutOnly;
static uint64_t settingCheckedGeneration,settingLayoutVersion=1;
static bool settingVisualsEnabled=true;
static void setting_begin(int map) {
    settingMap=map;settingDetailCount=0;settingCheckedGeneration=0;++settingLayoutVersion;
}
static void setting_register_box(int x0,int x1,int y0,int y1,int z0,int z1,Color expected) {
    // Keep load-bearing orange supports immediately recognizable.
    if(expected.r>200 && expected.g>100 && expected.b<80)return;
    float sx=(x1-x0+1)*VOXEL_SIZE,sy=(y1-y0+1)*VOXEL_SIZE,sz=(z1-z0+1)*VOXEL_SIZE;
    if(sy<1 || fmaxf(sx,sz)<1.4f)return;
    for(int side=0;side<2 && settingDetailCount<SETTING_DETAILS_MAX;side++) {
        SettingDetail *d=&settingDetails[settingDetailCount];*d=(SettingDetail){0};
        bool front=sx>=sz;int face=front?(side?z1:z0):(side?x1:x0);
        d->center=(Vector3){(x0+x1+1)*VOXEL_SIZE*.5f,(y0+y1+1)*VOXEL_SIZE*.5f,(z0+z1+1)*VOXEL_SIZE*.5f};
        d->width=fminf((front?sx:sz)-.15f,5.5f);d->height=fminf(sy-.12f,2.4f);
        if(front){d->center.z=(face+(side?1:0))*VOXEL_SIZE+(side?.018f:-.018f);d->yaw=side?0:180;}
        else {d->center.x=(face+(side?1:0))*VOXEL_SIZE+(side?.018f:-.018f);d->yaw=side?90:270;}
        d->style=settingDetailCount%4;d->label=(settingDetailCount/2)%8;
        for(int n=0;n<6;n++) {
            float u=(n%3-1)*d->width*.499f,v=(n/3?1:-1)*d->height*.499f;
            Vector3 at=Vector3RotateByQuaternion((Vector3){u,v,0},QuaternionFromAxisAngle((Vector3){0,1,0},d->yaw*DEG2RAD));at=v_add(at,d->center);
            int x=(int)floorf(at.x/VOXEL_SIZE),y=(int)floorf(at.y/VOXEL_SIZE),z=(int)floorf(at.z/VOXEL_SIZE);
            if(front)z=face;else x=face;
            d->cells[n][0]=x;d->cells[n][1]=y;d->cells[n][2]=z;
            d->colors[n]=expected;
        }
        settingDetailCount++;
    }
}
