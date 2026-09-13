struct AdVec { float x,y,z,w; };
struct AdCell { int grid[4]; uint meta[4]; AdVec x[8],v[8],old[8]; };
struct AdLink { int a,b,axis,side; float relaxation,padding[3]; };
struct AdRange { int begin,count; };
struct AdContact { int a,b,padding[2]; AdVec normal; float wa[8],wb[8]; AdVec pa,pb; };
