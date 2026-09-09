#include <metal_stdlib>
using namespace metal;
#include "adaptive_types.inc"
#define AD_FN inline
#define AD_MAX(a,b) max(float(a),float(b))
#define AD_CELL(i) cells[i]
#define AD_RANGE(i) ranges[i]
#define AD_LINK(i) links[i]
#define AD_ID(i) ids[i]
#define AD_CRANGE(i) ranges[i]
#define AD_CONTACT(i) contacts[i]
#define AD_CID(i) ids[i]
#define AD_LINK_PARAMS device const AdCell *cells,device const AdLink *links,device const AdRange *ranges,device const int *ids
#define AD_CONTACT_PARAMS device const AdCell *cells,device const AdContact *contacts,device const AdRange *ranges,device const int *ids
#include "adaptive_kernel.inc"
struct AdUniforms { int count,mode;float unit,padding; };
kernel void adaptive_pipeline(device const AdCell *cells [[buffer(0)]],
    device AdCell *out [[buffer(1)]],device const AdLink *links [[buffer(2)]],
    device const AdRange *ranges [[buffer(3)]],device const int *ids [[buffer(4)]],
    device const AdContact *contacts [[buffer(5)]],device const AdRange *cr [[buffer(6)]],
    device const int *ci [[buffer(7)]],constant AdUniforms &u [[buffer(8)]],uint index [[thread_position_in_grid]]) {
    if(index>=uint(u.count))return;
    if(u.mode==0)out[index]=ad_shape_cell(cells[index],u.unit);
    else if(u.mode==1)out[index]=ad_solve_links(int(index),cells,links,ranges,ids);
    else out[index]=ad_solve_contacts(int(index),cells,contacts,cr,ci);
}
