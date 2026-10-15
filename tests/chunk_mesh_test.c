#define ACTIVATION_PROFILE
#define main fps_application_main
#include "../fps_ray.c"
#undef main
#include <assert.h>
/* Compare the two real shaders at identical cube positions, before simulation
 * can move or deform the activated representation. */
static void check_voxel_appearance(void) {
    Color color={150,175,200,255};
    const Vector3 offsets[]={{6,3,6},{-6,3,-6},{6,0,0},{-6,0,0},{0,6,0},{0,-6,0},{0,0,6},{0,0,-6}};
    RenderTexture2D target=LoadRenderTexture(640,360);
    for(int view=0;view<8;view++) {
        clear_all_chunks();clear_world_voxels();
        for(int x=-3;x<3;x++)for(int y=1;y<4;y++)for(int z=-1;z<2;z++)
            add_static_voxel_at_grid(x,y,z,color,0);
        Camera3D cam={v_add((Vector3){0,1.25f,.25f},offsets[view]),{0,1.25f,.25f},
            view==4 || view==5?(Vector3){0,0,1}:(Vector3){0,1,0},4,CAMERA_ORTHOGRAPHIC};
        Image images[2];
        for(int dynamic=0;dynamic<2;dynamic++) {
            if(dynamic) {
                // Tether activation is budgeted; activate the remaining groups
                // without advancing physics until the whole fixture is dynamic.
                for(int batch=0;batch<54;batch++) {
                    int next=-1;for(int i=0;i<voxel_count;i++)if(!voxels[i].simulate){next=i;break;}
                    if(next<0)break;
                    assert(activate_static_voxel_for_tether(next,0,ACTIVATION_TETHER_BELIEF));
                }
                for(int i=0;i<voxel_count;i++)assert(voxels[i].simulate);
            }
            voxel_render_frame++;BeginTextureMode(target);ClearBackground((Color){30,35,40,255});BeginMode3D(cam);
            rlDisableBackfaceCulling();
            DrawVoxels(cam);
            rlEnableBackfaceCulling();EndMode3D();EndTextureMode();
            images[dynamic]=LoadImageFromTexture(target.texture);ImageFlipVertical(&images[dynamic]);
            if(view==0)ExportImage(images[dynamic],dynamic?"artifacts/voxel-edges-dynamic.png":"artifacts/voxel-edges-static.png");
        }
        Color *a=LoadImageColors(images[0]),*b=LoadImageColors(images[1]);int different=0;
        for(int i=0;i<640*360;i++)if(abs(a[i].r-b[i].r)>3 || abs(a[i].g-b[i].g)>3 || abs(a[i].b-b[i].b)>3)different++;
        fprintf(stderr,"Voxel appearance view %d: %.4f%% differing pixels\n",view,100.0*different/(640*360));
        // Triangle rasterization can disagree at a handful of silhouette pixels.
        assert(different<640*360/100);
        UnloadImageColors(a);UnloadImageColors(b);UnloadImage(images[0]);UnloadImage(images[1]);
    }
    UnloadRenderTexture(target);
}
static void check_mesh_faces(void) {
    static unsigned char seen[MAX_VOXELS];memset(seen,0,sizeof(seen));
    for(int c=0;c<voxel_chunk_count;c++) {
        VoxelChunk *chunk=&voxel_chunks[c];rebuild_chunk(chunk);Mesh *m=&chunk->mesh;
        for(int v=0;v<m->vertexCount;v+=4) {
            Vector3 p[4];for(int k=0;k<4;k++)p[k]=(Vector3){m->vertices[(v+k)*3]/VOXEL_SIZE,m->vertices[(v+k)*3+1]/VOXEL_SIZE,m->vertices[(v+k)*3+2]/VOXEL_SIZE};
            Vector3 n={m->normals[v*3],m->normals[v*3+1],m->normals[v*3+2]};
            int axis=n.x?0:n.y?1:2;bool positive=n.x+n.y+n.z>0;
            int face=axis*2+!positive;
            Vector3 du=v_sub(p[1],p[0]),dv=v_sub(p[3],p[0]);
            int width=(int)lroundf(v_length(du)),height=(int)lroundf(v_length(dv));assert(width>0 && height>0);
            // Every merged face repeats the edge shader once per voxel.
            // This also covers negative coordinates and restored chunk meshes.
            for(int k=0;k<4;k++) {
                assert(m->texcoords[(v+k)*2]==((k==1 || k==2)?width:0));
                assert(m->texcoords[(v+k)*2+1]==(k>=2?height:0));
            }
            du=v_mul(du,1.0f/width);dv=v_mul(dv,1.0f/height);
            for(int j=0;j<height;j++)for(int i=0;i<width;i++) {
                Vector3 cell=v_add(p[0],v_add(v_mul(du,i+.5f),v_mul(dv,j+.5f)));cell=v_sub(cell,v_mul(n,.5f));
                int x=(int)floorf(cell.x),y=(int)floorf(cell.y),z=(int)floorf(cell.z);
                int idx=table_get_static_only(x,y,z);assert(idx>=0);
                assert(table_get_static_only(x+(int)n.x,y+(int)n.y,z+(int)n.z)<0);
                assert(!(seen[idx]&(1<<face)));seen[idx]|=1<<face;
            }
        }
    }
    for(int i=0;i<voxel_count;i++)if(!voxels[i].simulate)for(int face=0;face<6;face++) {
        bool exposed=table_get_static_only(voxels[i].gx+face_offsets[face][0],voxels[i].gy+face_offsets[face][1],voxels[i].gz+face_offsets[face][2])<0;
        if(((seen[i]>>face)&1)!=exposed)fprintf(stderr,"Face mismatch idx=%d grid=%d,%d,%d face=%d seen=%u expected=%d chunks=%d\n",i,voxels[i].gx,voxels[i].gy,voxels[i].gz,face,seen[i],exposed,voxel_chunk_count);
        assert(((seen[i]>>face)&1)==exposed);
    }
}
int main(void) {
    SetLoggingEnabled(false);SetTraceLogLevel(LOG_WARNING);SetConfigFlags(FLAG_WINDOW_HIDDEN);InitWindow(640,360,"Chunk mesh regression");
    assert(grid_to_chunk_coord(-1)==-1 && grid_to_chunk_coord(-31)==-1 && grid_to_chunk_coord(-32)==-1 && grid_to_chunk_coord(-33)==-2);
    clear_world_voxels();clear_all_chunks();
    /* Full-width mask, a cavity, negative coordinates, and all three chunk seams. */
    for(int z=-1;z<=1;z++)for(int y=30;y<=33;y++)for(int x=-33;x<=33;x++)add_static_voxel_at_grid(x,y,z,GRAY,0);
    int hole=table_get_static_only(0,31,0);remove_voxel_index(hole);check_mesh_faces();
    /* A distant mesh must keep its GPU allocation through local restore and removal. */
    add_static_voxel_at_grid(200,10,200,BLUE,0);check_mesh_faces();
    VoxelChunk *distant=get_or_create_chunk(grid_to_chunk_coord(200),grid_to_chunk_coord(10),grid_to_chunk_coord(200));unsigned int vao=distant->mesh.vaoId;
    int idx=table_get_static_only(31,31,0);Voxel snapshot=voxels[idx];assert(recycle_queue_push(&snapshot));remove_voxel_index(idx);
    assert(!distant->dirty);check_mesh_faces();assert(distant->mesh.vaoId==vao);
    recycleFrameCounter+=RECYCLE_STATIC_RESTORE_DELAY;meshDirty=false;recycle_dead_voxels(1.0f/60);
    assert(regeneration.active && table_get_static_only(31,31,0)<0);
    for(int f=0;f<REGEN_PREVIEW_FRAMES;f++)recycle_dead_voxels(1.0f/60);
    assert(table_get_static_only(31,31,0)>=0 && !meshDirty && !distant->dirty);check_mesh_faces();assert(distant->mesh.vaoId==vao);
    idx=table_get_static_only(31,31,0);voxels[idx].activationCooldownFrames=0;
    assert(rip_single_static_voxel(idx,0)>=0 && !meshDirty && !distant->dirty);
    check_mesh_faces();assert(distant->mesh.vaoId==vao);
    /* Sparse deterministic geometry catches holes and incorrect greedy merges. */
    clear_world_voxels();clear_all_chunks();unsigned int seed=91;
    for(int k=0;k<300;k++){seed=seed*1664525u+1013904223u;int x=(seed>>8)%70-35;seed=seed*1664525u+1013904223u;int y=(seed>>8)%36;seed=seed*1664525u+1013904223u;int z=(seed>>8)%70-35;add_static_voxel_at_grid(x,y,z,GRAY,0);}
    check_mesh_faces();
    /* Tether activation crosses a negative chunk seam without invalidating the
     * unrelated chunk. Compare every remaining static face after the removal. */
    clear_all_chunks();clear_world_voxels();activePlayers=0;
    for(int x=-34;x<=-30;x++)for(int y=10;y<=12;y++)add_static_voxel_at_grid(x,y,0,GRAY,0);
    add_static_voxel_at_grid(200,10,200,BLUE,0);check_mesh_faces();
    distant=get_or_create_chunk(grid_to_chunk_coord(200),grid_to_chunk_coord(10),grid_to_chunk_coord(200));vao=distant->mesh.vaoId;
    meshDirty=false;memset(&activationProfile,0,sizeof(activationProfile));
    assert(activate_static_voxel_for_tether(table_get_static_only(-32,11,0),0,ACTIVATION_TETHER_BELIEF));
    assert(!meshDirty && !distant->dirty && activationProfile.columnScans==0);
    assert(staticBeliefBatchDepth==0 && recycleQueueCount>0);
    check_mesh_faces();assert(distant->mesh.vaoId==vao);
    check_voxel_appearance();
    clear_all_chunks();clear_world_voxels();CloseWindow();puts("Chunk mesh regression passed");return 0;
}
