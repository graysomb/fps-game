#import <Foundation/Foundation.h>
#import <Metal/Metal.h>
#include "adaptive_physics.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef struct {float pos[4],prev[4],pred[4],vel[4];} Particle;
typedef struct {uint32_t corner[8];float pos[4],vel[4];int32_t flags[4];float bmin[4],bmax[4];uint32_t life[4];} Voxel;
static void moments(const Particle *p,uint32_t n,double out[7]){
    memset(out,0,7*sizeof(*out));for(uint32_t i=0;i<n;i++)if(p[i].prev[3]>0){double m=1.0/p[i].prev[3];out[0]+=m;for(int d=0;d<3;d++)out[1+d]+=m*p[i].vel[d];
        out[4]+=m*(p[i].pos[1]*p[i].vel[2]-p[i].pos[2]*p[i].vel[1]);out[5]+=m*(p[i].pos[2]*p[i].vel[0]-p[i].pos[0]*p[i].vel[2]);out[6]+=m*(p[i].pos[0]*p[i].vel[1]-p[i].pos[1]*p[i].vel[0]);}
}
static void same_moments(const double *a,const double *b){for(int i=0;i<7;i++){double err=fabs(a[i]-b[i])/fmax(1,fabs(a[i]));if(err>(i?1e-4:1e-5)){fprintf(stderr,"moment %d error %.9g before %.9g after %.9g\n",i,err,a[i],b[i]);abort();}}}
static void finish(id<MTLCommandBuffer> command,id<MTLComputeCommandEncoder> encoder){[encoder endEncoding];[command commit];[command waitUntilCompleted];if(command.status!=MTLCommandBufferStatusCompleted){fprintf(stderr,"%s\n",command.error.localizedDescription.UTF8String);abort();}}
static void fixture(id<MTLDevice> device,id<MTLCommandQueue> queue,bool duplicates,bool two_arenas) {
    enum{SIDE=8,N=SIDE*SIDE*SIDE,GRID_P=(SIDE+1)*(SIDE+1)*(SIDE+1)};
    const int P=duplicates?8*N:GRID_P;
    AdaptiveInput input[N];AdaptiveParticleInput info[P];memset(info,0,sizeof(info));
    id<MTLBuffer> particles=[device newBufferWithLength:P*sizeof(Particle) options:MTLResourceStorageModeShared];
    id<MTLBuffer> voxels=[device newBufferWithLength:N*sizeof(Voxel) options:MTLResourceStorageModeShared];
    id<MTLBuffer> sims=[device newBufferWithLength:P*4 options:MTLResourceStorageModeShared];
    id<MTLBuffer> control=[device newBufferWithLength:16 options:MTLResourceStorageModeShared];
    Particle *p=particles.contents;Voxel *v=voxels.contents;uint32_t *ids=sims.contents,*pc=control.contents;
    pc[0]=pc[1]=pc[2]=P;pc[3]=0;
    for(int i=0;i<P;i++)ids[i]=i;
    for(int z=0;z<=SIDE;z++)for(int y=0;y<=SIDE;y++)for(int x=0;x<=SIDE;x++){
        int i=x+(SIDE+1)*(y+(SIDE+1)*z);ids[i]=i;
        p[i]=(Particle){{x,y,z,.5},{x,y,z,1},{x,y,z,1},{1-.2f*y,2+.2f*x,.3f,0}};
        info[i].xyz[0]=x;info[i].xyz[1]=y;info[i].xyz[2]=z;info[i].inverse_mass=info[i].base_inverse_mass=1;info[i].simulated=1;
    }
    for(int z=0;z<SIDE;z++)for(int y=0;y<SIDE;y++)for(int x=0;x<SIDE;x++){
        int i=x+SIDE*(y+SIDE*z);int ox=x+(two_arenas?(x<4?1:5):0);uint32_t arena=two_arenas?(x<4?101:205):0;
        input[i]=(AdaptiveInput){{ox,y,z},arena,0,0};v[i]=(Voxel){0};v[i].pos[3]=v[i].vel[3]=1;v[i].flags[0]=v[i].flags[3]=1;v[i].life[2]=63;
        for(int k=0;k<8;k++){
            int px=ox+(k&1),py=y+((k>>1)&1),pz=z+((k>>2)&1);
            int id=duplicates?8*i+k:px+(SIDE+1)*(py+(SIDE+1)*pz);
            if(duplicates){p[id]=(Particle){{px,py,pz,.5},{px,py,pz,1},{px,py,pz,1},{1-.2f*py,2+.2f*px,.3f,0}};
                info[id]=(AdaptiveParticleInput){0};info[id].xyz[0]=px;info[id].xyz[1]=py;info[id].xyz[2]=pz;
                info[id].inverse_mass=info[id].base_inverse_mass=1;info[id].simulated=1;info[id].arena=arena;}
            v[i].corner[k]=id;info[id].incident[info[id].count++]=8*i+k;
        }
    }
    char error[4096]={0};AdaptivePhysics *s=adaptive_physics_create((void *)device,input,N,info,P,error,sizeof(error));if(!s){fprintf(stderr,"%s\n",error);abort();}
    FpsGpuUniforms u={0};u.voxel_count=N;u.particle_count=u.sim_count=P;u.sor=1;u.vgs_epsilon=1e-6;u.vgs_alpha=.1;u.voxel_size=1;u.dt=1.0/120;
    double before[7],after[7];moments(p,P,before);
    id<MTLCommandBuffer> cb=[queue commandBuffer];id<MTLComputeCommandEncoder> enc=[cb computeCommandEncoder];
    assert(adaptive_physics_build(s,(void *)enc,(void *)particles,(void *)voxels,(void *)sims,(void *)control,&u));
    assert(adaptive_physics_stage(s,(void *)enc,(void *)particles,(void *)voxels,(void *)sims,(void *)control,&u,ADAPTIVE_STAGE_REFRESH));finish(cb,enc);
    uint32_t leaves,generation,err;adaptive_physics_diagnostics(s,&leaves,&generation,&err);assert(!err&&leaves<N&&pc[1]<(uint32_t)P);
    moments(p,P,after);same_moments(before,after);assert(adaptive_physics_mass_error(s,(void *)particles)<1e-5);
    uint32_t last_leaves=leaves,last_generation=generation;
    Particle *snapshot=malloc((size_t)P*sizeof(*snapshot));assert(snapshot);
    for(int repeat=0;repeat<3;repeat++){
        memcpy(snapshot,p,(size_t)P*sizeof(*snapshot));
        int cell=two_arenas?2+SIDE*(3+repeat+SIDE*3):(2+repeat)+SIDE*(3+SIDE*3);v[cell].flags[0]=0;v[cell].life[2]=0;
        moments(p,P,before);cb=[queue commandBuffer];enc=[cb computeCommandEncoder];
        assert(adaptive_physics_stage(s,(void *)enc,(void *)particles,(void *)voxels,(void *)sims,(void *)control,&u,ADAPTIVE_STAGE_REFRESH));
        assert(adaptive_physics_stage(s,(void *)enc,(void *)particles,(void *)voxels,(void *)sims,(void *)control,&u,ADAPTIVE_STAGE_FRACTURE));
        finish(cb,enc);moments(p,P,after);same_moments(before,after);
        if(two_arenas)for(int i=0;i<P;i++)if(info[i].arena==205)assert(memcmp(snapshot+i,p+i,sizeof(*p))==0);
        adaptive_physics_diagnostics(s,&leaves,&generation,&err);assert(!err&&leaves>=last_leaves&&generation>=last_generation);last_leaves=leaves;last_generation=generation;
        assert(adaptive_physics_mass_error(s,(void *)particles)<1e-5);
    }
    /* At-rest shape/interface passes must not introduce a seam or NaNs. */
    cb=[queue commandBuffer];enc=[cb computeCommandEncoder];
    for(int it=0;it<4;it++){
        assert(adaptive_physics_stage(s,(void *)enc,(void *)particles,(void *)voxels,(void *)sims,(void *)control,&u,ADAPTIVE_STAGE_SHAPE));
        for(int j=0;j<2;j++)assert(adaptive_physics_stage(s,(void *)enc,(void *)particles,(void *)voxels,(void *)sims,(void *)control,&u,ADAPTIVE_STAGE_ATTACH));
    }
    finish(cb,enc);for(int i=0;i<P;i++)for(int d=0;d<3;d++)assert(isfinite(p[i].pred[d])&&fabs(p[i].pred[d]-info[i].xyz[d])<.05);
    free(snapshot);
    adaptive_physics_destroy(s);[particles release];[voxels release];[sims release];[control release];
}
int main(void){@autoreleasepool{
    id<MTLDevice> device=MTLCreateSystemDefaultDevice();if(!device)return 77;
    id<MTLCommandQueue> queue=[device newCommandQueue];
    fixture(device,queue,false,false);fixture(device,queue,true,false);fixture(device,queue,true,true);
    [queue release];[device release];
    puts("adaptive physics: shared and coincident IDs, GPU refinement, mass, momentum and rest-state constraints passed");return 0;
}}
