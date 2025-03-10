
/* --------------------------------------------------------------------------
 * Project: SpinWalk
 * File: phantom_ply.cpp
 *
 * Author   : Ali Aghaeifar <ali.aghaeifar@tuebingen.mpg.de>
 * Date     : 20.01.2022
 * Descrip  : 
 * -------------------------------------------------------------------------- */
#include <filesystem>
#include <random>
#include <execution>
#include <limits>
#include <ranges>
// boost includes
#include <boost/log/trivial.hpp> 

#include <highfive/highfive.hpp>
#include "barkeep.h"
#include "happly.h"
#include "phantom_ply.h"

// -------------------------------------------------------------------------- //
namespace phantom
{
// Helper functions (since C++20 doesn’t natively provide these)
Vec3 cross(const Vec3& a, const Vec3& b) {
    return {a.y * b.z - a.z * b.y,
            a.z * b.x - a.x * b.z,
            a.x * b.y - a.y * b.x};
}

float dot(const Vec3& a, const Vec3& b) { return a.x * b.x + a.y * b.y + a.z * b.z;}

struct BVHNode {
    AABB bounds;
    std::vector<Triangle> triangles;
    std::unique_ptr<BVHNode> left, right;
    bool isLeaf() const { return !left && !right; }
};

AABB computeAABB(std::vector<Triangle>& tris, size_t start, size_t end){
    AABB aabb{{std::numeric_limits<double>::max()   , std::numeric_limits<double>::max()   , std::numeric_limits<double>::max()   },
              {std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest()}};
    for (size_t i = start; i < end; ++i) {
        aabb.expand(tris[i].v0);
        aabb.expand(tris[i].v1);
        aabb.expand(tris[i].v2);
    }
    return aabb;
}

int partitionTriangles(std::vector<Triangle>& tris, size_t start, size_t end, const AABB& bounds) {
    // which axis to split along?
    Vec3 extent = bounds.max - bounds.min;
    int axis = 0;
    if (extent.y > extent.x) axis = 1;
    if (extent.z > std::max(extent.x, extent.y)) axis = 2;
    // sort triangles based on centroid along axis
    auto compareCentroid = [axis](const Triangle& a, const Triangle& b) {
        auto getCoord = [](const Vec3& v, int ax) { return (ax == 0) ? v.x : (ax == 1) ? v.y : v.z; };
        double ca = (getCoord(a.v0, axis) + getCoord(a.v1, axis) + getCoord(a.v2, axis)) / 3.0f;
        double cb = (getCoord(b.v0, axis) + getCoord(b.v1, axis) + getCoord(b.v2, axis)) / 3.0f;
        return ca < cb;
    };

    std::sort(tris.begin() + start, tris.begin() + end, compareCentroid);
    return static_cast<int>(start + (end - start) / 2);
}

std::unique_ptr<BVHNode> buildBVH(std::vector<Triangle>& tris, int start, int end) {
    auto node = std::make_unique<BVHNode>();
    node->bounds = computeAABB(tris, start, end);
    if (end - start <= 4) {
        node->triangles.assign(tris.begin() + start, tris.begin() + end);
        return node;
    }
    // sorts triangles based on centroid along the longest axis and returns the median index
    int mid = partitionTriangles(tris, start, end, node->bounds); // Median split
    node->left  = buildBVH(tris, start, mid);
    node->right = buildBVH(tris, mid, end);
    return node;
}


// Simple ray-triangle intersection (Möller-Trumbore)
bool rayIntersectsTriangle(const Vec3& origin, const Vec3& dir, const Triangle& tri) {
    constexpr float EPSILON = 1e-6f;
    Vec3 e1 = tri.v1 - tri.v0;
    Vec3 e2 = tri.v2 - tri.v0;
    Vec3 h = cross(dir, e2);
    float a = dot(e1, h);

    if (std::abs(a) < EPSILON) return false; // Parallel to triangle

    float f = 1.0f / a;
    Vec3 s = origin - tri.v0;
    float u = f * dot(s, h);
    if (u < 0.0f || u > 1.0f) return false;

    Vec3 q = cross(s, e1);
    float v = f * dot(dir, q);
    if (v < 0.0f || u + v > 1.0f) return false;

    float t = f * dot(e2, q);
    return t >= 0.0f; // Hit if t is positive
}

// Recursive ray traversal
void intersectRayNode(const BVHNode *node, const Vec3 &origin, const Vec3 &dir, std::vector<const Triangle *> &hits)
{
    if (!node) return;
    if (!node->bounds.intersectsRay(origin))
        return;

    if (node->isLeaf()){
        for (const auto &tri : node->triangles)
            if (rayIntersectsTriangle(origin, dir, tri))
                hits.push_back(&tri);
    }else{
        intersectRayNode(node->left.get(), origin, dir, hits);
        intersectRayNode(node->right.get(), origin, dir, hits);
    }
}


bool is_inside_particle(const std::array<double, 3> &ray_origin, BVHNode *root) 
{
    Vec3 origin{ray_origin[0], ray_origin[1], ray_origin[2]};
    Vec3 dir{1., 0., 0.}; // do not change direction! 
    std::vector<const Triangle*> hits;
    intersectRayNode(root, origin, dir, hits);
    
    return hits.size() % 2;
}

ply::ply()
{

}

ply::ply(float fov_um, size_t resolution, float dChi, float Y, std::string filename_in, std::string filename_out)
: phantom_base(fov_um, resolution, dChi, Y, 0, 0, filename_out)
{
    filename_ply = filename_in;
}

ply::~ply()
{
}

bool ply::generate_mask_fieldmap()
{
    BOOST_LOG_TRIVIAL(info) << "Reading mesh..." << std::endl;
    auto start = std::chrono::high_resolution_clock::now();

    happly::PLYData plyIn(filename_ply);
    
    // Get mesh-style data from the object
    std::vector<std::array<double, 3>> vPos = plyIn.getVertexPositions();
    std::vector<std::vector<size_t>> fInd   = plyIn.getFaceIndices<size_t>();
    for(const auto &v : fInd)
        if(v.size() != 3) {
            BOOST_LOG_TRIVIAL(error) << "Only triangular mesh is supported!";
            return false;
        }

    // copy vertices to own data structure + convert to um 
    std::vector<Triangle> triangles;
    for(const auto &f : fInd){
        Triangle tri;
        tri.v0.x=vPos[f[0]][0] * 1e3; tri.v0.y=vPos[f[0]][1] * 1e3; tri.v0.z=vPos[f[0]][2] * 1e3;
        tri.v1.x=vPos[f[1]][0] * 1e3; tri.v1.y=vPos[f[1]][1] * 1e3; tri.v1.z=vPos[f[1]][2] * 1e3;
        tri.v2.x=vPos[f[2]][0] * 1e3; tri.v2.y=vPos[f[2]][1] * 1e3; tri.v2.z=vPos[f[2]][2] * 1e3;
        triangles.push_back(tri);
    }
    // center the mesh
    boundry = computeAABB(triangles, 0, triangles.size());
    BOOST_LOG_TRIVIAL(info) << "mesh FoV min (um): " << boundry.min.x << " " << boundry.min.y << " " << boundry.min.z << std::endl;
    BOOST_LOG_TRIVIAL(info) << "mesh FoV max (um): " << boundry.max.x << " " << boundry.max.y << " " << boundry.max.z << std::endl;
    Vec3 shift = (boundry.max + boundry.min)/2.;
    for(auto &tri : triangles){
        tri.v0 = tri.v0 + m_fov/2. - shift;
        tri.v1 = tri.v1 + m_fov/2. - shift;
        tri.v2 = tri.v2 + m_fov/2. - shift;
    }
    // update boundry
    boundry = computeAABB(triangles, 0, triangles.size());
    // Bounding Volume Hierarchy (BVH)
    BOOST_LOG_TRIVIAL(info) << "Building Bounding Volume Hierarchy..." << std::endl;
    std::unique_ptr<BVHNode> root = buildBVH(triangles, 0, triangles.size());
    BOOST_LOG_TRIVIAL(info) << "Shifted mesh FoV min (um): " << root->bounds.min.x << " " << root->bounds.min.y << " " << root->bounds.min.z << std::endl;
    BOOST_LOG_TRIVIAL(info) << "Shifted mesh FoV max (um): " << root->bounds.max.x << " " << root->bounds.max.y << " " << root->bounds.max.z << std::endl;

    size_t res1 = m_resolution;
    size_t res2 = res1 * res1;
    size_t res3 = res1 * res2;
    m_calc_fieldmap = false;
    m_fieldmap.resize(res3, 0);
    m_mask.resize(res3, 0);
    m_volume_fraction = 0;
    m_Y = -1;

    size_t c = 0;
    auto bar = barkeep::ProgressBar(&c, {.total = res1, .message = "Progress", .style = barkeep::ProgressBarStyle::Rich,});
    
    std::vector<uint32_t> v(res1*res1);
    std::generate(std::execution::seq, v.begin(), v.end(), [n = 0] () mutable { return n++; }); 

    for(int32_t pz=0; pz<res1; pz++, c++){
    std::for_each(std::execution::par_unseq, v.begin(), v.end(), [&](int idx) {
        int32_t py = idx / res1;
        int32_t px = idx % res1;
        std::array<double, 3> point;
        size_t p = px*res2 + py*res1 + pz;  
        // convert to mm      
        point[0] = m_grid[3*p+0]; 
        point[1] = m_grid[3*p+1]; 
        point[2] = m_grid[3*p+2];
        // point lies outside volume?
        if(point[0] < root->bounds.min.x || point[0] > root->bounds.max.x || point[1] < root->bounds.min.y || point[1] > root->bounds.max.y || point[2] < root->bounds.min.z || point[2] > root->bounds.max.z)
            return;

        m_mask[p] = is_inside_particle(point, root.get());
    });
    }
    bar->done(); 
    auto end = std::chrono::high_resolution_clock::now();
    BOOST_LOG_TRIVIAL(info) << "Mesh converted successfully! " << "Elapsed Time: " << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() << " s\n";

    return true;
}


std::ostream& operator<<(std::ostream& os, const ply& obj)
{
    os << static_cast<const phantom_base&>(obj);
    return os;
}

// -------------------------------------------------------------------------- //

bool ply::run(bool write_to_disk)
{
    BOOST_LOG_TRIVIAL(info) << *this << std::endl;
    if(create_grid() == false)
        return false;
    if(generate_mask_fieldmap() == false)
        return false;
    if(write_to_disk)
        if(save() == false)
            return false;
    return true;
}

}