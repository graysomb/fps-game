#define main fps_application_main
#include "../fps_ray.c"
#undef main
#include <assert.h>
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
    recycleFrameCounter+=RECYCLE_STATIC_RESTORE_DELAY;meshDirty=false;recycle_dead_voxels();
    assert(table_get_static_only(31,31,0)>=0 && !meshDirty && !distant->dirty);check_mesh_faces();assert(distant->mesh.vaoId==vao);
    idx=table_get_static_only(31,31,0);voxels[idx].activationCooldownFrames=0;
    assert(rip_single_static_voxel(idx,0)>=0 && !meshDirty && !distant->dirty);
    check_mesh_faces();assert(distant->mesh.vaoId==vao);
    /* Sparse deterministic geometry catches holes and incorrect greedy merges. */
    clear_world_voxels();clear_all_chunks();unsigned int seed=91;
    for(int k=0;k<300;k++){seed=seed*1664525u+1013904223u;int x=(seed>>8)%70-35;seed=seed*1664525u+1013904223u;int y=(seed>>8)%36;seed=seed*1664525u+1013904223u;int z=(seed>>8)%70-35;add_static_voxel_at_grid(x,y,z,GRAY,0);}
    check_mesh_faces();clear_all_chunks();clear_world_voxels();CloseWindow();puts("Chunk mesh regression passed");return 0;
}
