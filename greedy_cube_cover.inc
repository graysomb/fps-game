/* Deterministic greedy cover of occupied integer cells by disjoint cubes.
 * This file deliberately has no dependency on the voxel solver so the cover
 * can be exercised by the headless invariant test. */
typedef struct { int x,y,z,type,source; } GreedyCell;
typedef struct { int x,y,z,size,type,first_cell,cell_count; } GreedyCube;
typedef struct {
    GreedyCube *cubes; int cube_count,cube_capacity;
    int *cell_cube; int cell_count;
    double milliseconds;
} GreedyCubeCover;
typedef struct { int x,y,z,type; signed char dx,dy,dz; } GreedyCandidate;

static int greedy_cell_compare(const void *aa,const void *bb) {
    const GreedyCell *a=aa,*b=bb;
    if(a->type!=b->type)return a->type<b->type?-1:1;
    if(a->x!=b->x)return a->x<b->x?-1:1;
    if(a->y!=b->y)return a->y<b->y?-1:1;
    if(a->z!=b->z)return a->z<b->z?-1:1;
    return a->source<b->source?-1:a->source>b->source;
}
static int greedy_candidate_compare(const void *aa,const void *bb) {
    const GreedyCandidate *a=aa,*b=bb;
    if(a->x!=b->x)return a->x<b->x?-1:1;if(a->y!=b->y)return a->y<b->y?-1:1;
    if(a->z!=b->z)return a->z<b->z?-1:1;if(a->type!=b->type)return a->type<b->type?-1:1;if(a->dx!=b->dx)return a->dx<b->dx?-1:1;
    if(a->dy!=b->dy)return a->dy<b->dy?-1:1;return a->dz<b->dz?-1:a->dz>b->dz;
}
static int greedy_find_cell(const GreedyCell *cells,int count,int type,int x,int y,int z) {
    GreedyCell key={x,y,z,type,0};int lo=0,hi=count;
    while(lo<hi){int m=lo+(hi-lo)/2,c=greedy_cell_compare(&cells[m],&key);
        if(c<0)lo=m+1;else hi=m;}
    return lo<count&&cells[lo].type==type&&cells[lo].x==x&&cells[lo].y==y&&cells[lo].z==z?lo:-1;
}
static bool greedy_present(const GreedyCell *cells,const unsigned char *remaining,int count,
                           int type,int x,int y,int z) {
    int i=greedy_find_cell(cells,count,type,x,y,z);return i>=0&&remaining[i];
}
static bool greedy_push_candidate(GreedyCandidate **items,int *count,int *capacity,GreedyCandidate c){
    if(*count>=*capacity){int next=*capacity?*capacity*2:256;void *p=realloc(*items,(size_t)next*sizeof(**items));
        if(!p)return false;*items=p;*capacity=next;}(*items)[(*count)++]=c;return true;
}
static bool greedy_add_vertex_candidates(GreedyCandidate **items,int *count,int *capacity,
        const GreedyCell *cells,const unsigned char *remaining,int n,int type,int vx,int vy,int vz){
    bool occupied=false,empty=false;
    for(int dz=-1;dz<=0;dz++)for(int dy=-1;dy<=0;dy++)for(int dx=-1;dx<=0;dx++)
        if(greedy_present(cells,remaining,n,type,vx+dx,vy+dy,vz+dz))occupied=true;else empty=true;
    if(!occupied||!empty)return true;
    for(int dz=-1;dz<=1;dz+=2)for(int dy=-1;dy<=1;dy+=2)for(int dx=-1;dx<=1;dx+=2){
        int x=vx+(dx<0?-1:0),y=vy+(dy<0?-1:0),z=vz+(dz<0?-1:0);
        if(greedy_present(cells,remaining,n,type,x,y,z)&&
           !greedy_push_candidate(items,count,capacity,(GreedyCandidate){vx,vy,vz,type,dx,dy,dz}))return false;
    }return true;
}
static int greedy_largest_size(GreedyCandidate c,int type,const GreedyCell *cells,
                               const unsigned char *remaining,int count){
    int size=0;
    for(;;){int ns=size+1;bool complete=true;
        /* Three disjoint pieces of the newly added shell. */
        for(int b=0;b<ns&&complete;b++)for(int a=0;a<ns&&complete;a++){
            int u[3]={ns-1,a,b};int xyz[3];
            xyz[0]=c.x+(c.dx>0?u[0]:-u[0]-1);xyz[1]=c.y+(c.dy>0?u[1]:-u[1]-1);xyz[2]=c.z+(c.dz>0?u[2]:-u[2]-1);
            if(!greedy_present(cells,remaining,count,type,xyz[0],xyz[1],xyz[2]))complete=false;
        }
        for(int b=0;b<ns&&complete;b++)for(int a=0;a<size&&complete;a++){
            int u[3]={a,ns-1,b};int xyz[3];
            xyz[0]=c.x+(c.dx>0?u[0]:-u[0]-1);xyz[1]=c.y+(c.dy>0?u[1]:-u[1]-1);xyz[2]=c.z+(c.dz>0?u[2]:-u[2]-1);
            if(!greedy_present(cells,remaining,count,type,xyz[0],xyz[1],xyz[2]))complete=false;
        }
        for(int b=0;b<size&&complete;b++)for(int a=0;a<size&&complete;a++){
            int u[3]={a,b,ns-1};int xyz[3];
            xyz[0]=c.x+(c.dx>0?u[0]:-u[0]-1);xyz[1]=c.y+(c.dy>0?u[1]:-u[1]-1);xyz[2]=c.z+(c.dz>0?u[2]:-u[2]-1);
            if(!greedy_present(cells,remaining,count,type,xyz[0],xyz[1],xyz[2]))complete=false;
        }
        if(!complete)return size;size=ns;
    }
}
static void greedy_candidate_min(GreedyCandidate c,int size,int *x,int *y,int *z){
    *x=c.x+(c.dx<0?-size:0);*y=c.y+(c.dy<0?-size:0);*z=c.z+(c.dz<0?-size:0);
}
static bool greedy_cube_better(int size,int x,int y,int z,GreedyCandidate c,
                               int best,int bx,int by,int bz,GreedyCandidate bc){
    if(size!=best)return size>best;if(x!=bx)return x<bx;if(y!=by)return y<by;if(z!=bz)return z<bz;
    return greedy_candidate_compare(&c,&bc)<0;
}
static void greedy_cube_cover_free(GreedyCubeCover *cover){
    if(!cover)return;free(cover->cubes);free(cover->cell_cube);memset(cover,0,sizeof(*cover));
}
static bool greedy_cube_cover_build(const GreedyCell *input,int input_count,GreedyCubeCover *out){
    if(!out||input_count<0)return false;memset(out,0,sizeof(*out));if(input_count==0)return true;
    GreedyCell *cells=malloc((size_t)input_count*sizeof(*cells));unsigned char *remaining=malloc((size_t)input_count);
    GreedyCandidate *frontier=NULL;int frontier_count=0,frontier_capacity=0;
    if(!cells||!remaining)goto fail;memcpy(cells,input,(size_t)input_count*sizeof(*cells));qsort(cells,input_count,sizeof(*cells),greedy_cell_compare);
    for(int i=1;i<input_count;i++)if(cells[i].type==cells[i-1].type&&cells[i].x==cells[i-1].x&&cells[i].y==cells[i-1].y&&cells[i].z==cells[i-1].z)goto fail;
    memset(remaining,1,(size_t)input_count);out->cell_cube=malloc((size_t)input_count*sizeof(*out->cell_cube));if(!out->cell_cube)goto fail;
    for(int i=0;i<input_count;i++)out->cell_cube[i]=-1;out->cell_count=input_count;
    for(int i=0;i<input_count;i++)for(int dz=0;dz<=1;dz++)for(int dy=0;dy<=1;dy++)for(int dx=0;dx<=1;dx++)
        if(!greedy_add_vertex_candidates(&frontier,&frontier_count,&frontier_capacity,cells,remaining,input_count,cells[i].type,cells[i].x+dx,cells[i].y+dy,cells[i].z+dz))goto fail;
    int left=input_count;
    while(left){
        qsort(frontier,frontier_count,sizeof(*frontier),greedy_candidate_compare);int w=0;
        for(int i=0;i<frontier_count;i++)if(i==0||greedy_candidate_compare(&frontier[i],&frontier[i-1]))frontier[w++]=frontier[i];
        frontier_count=w;int best=0,bx=0,by=0,bz=0,btype=0;GreedyCandidate bc={0};
        for(int i=0;i<frontier_count;i++){
            GreedyCandidate c=frontier[i];int fx=c.x+(c.dx<0?-1:0),fy=c.y+(c.dy<0?-1:0),fz=c.z+(c.dz<0?-1:0);
            int ci=greedy_find_cell(cells,input_count,c.type,fx,fy,fz);
            if(ci<0)continue;int size=greedy_largest_size(c,cells[ci].type,cells,remaining,input_count),x,y,z;greedy_candidate_min(c,size,&x,&y,&z);
            if(size>0&&greedy_cube_better(size,x,y,z,c,best,bx,by,bz,bc)){best=size;bx=x;by=y;bz=z;btype=cells[ci].type;bc=c;}
        }
        if(best==0){int i=0;while(i<input_count&&!remaining[i])i++;if(i==input_count)break;best=1;bx=cells[i].x;by=cells[i].y;bz=cells[i].z;btype=cells[i].type;}
        if(out->cube_count>=out->cube_capacity){int cap=out->cube_capacity?out->cube_capacity*2:32;void *p=realloc(out->cubes,(size_t)cap*sizeof(*out->cubes));if(!p)goto fail;out->cubes=p;out->cube_capacity=cap;}
        int cube_index=out->cube_count++;GreedyCube *cube=&out->cubes[cube_index];*cube=(GreedyCube){bx,by,bz,best,btype,-1,0};
        for(int z=bz;z<bz+best;z++)for(int y=by;y<by+best;y++)for(int x=bx;x<bx+best;x++){
            int i=greedy_find_cell(cells,input_count,btype,x,y,z);if(i<0||!remaining[i])goto fail;remaining[i]=0;left--;out->cell_cube[cells[i].source]=cube_index;if(cube->first_cell<0)cube->first_cell=cells[i].source;cube->cell_count++;
        }
        for(int z=bz;z<=bz+best;z++)for(int y=by;y<=by+best;y++)for(int x=bx;x<=bx+best;x++)
            if(x==bx||x==bx+best||y==by||y==by+best||z==bz||z==bz+best)
                if(!greedy_add_vertex_candidates(&frontier,&frontier_count,&frontier_capacity,cells,remaining,input_count,btype,x,y,z))goto fail;
    }
    free(cells);free(remaining);free(frontier);return true;
fail: free(cells);free(remaining);free(frontier);greedy_cube_cover_free(out);return false;
}
