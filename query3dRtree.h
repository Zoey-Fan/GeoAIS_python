/*
 * date: 2024-04-26
 * author: Mengyu Ma@National University of Defense Technology
 *         Yifan Zhang@National University of Defense Technology
 * e-mail: mamengyu10@nudt.edu.cn
 *         zhangyifan.000@nudt.edu.cn
 * description:  3d index query function.
 */
#ifndef QUERY3DRTREE_H
#define QUERY3DRTREE_H

#include <functional>
#include "RTree.h"
#include <functional>

#define RTREE_TEMPLATE template <class DATATYPE, class ELEMTYPE, int NUMDIMS, class ELEMTYPEREAL, int TMAXNODES, int TMINNODES>
#define RTREE_QUAL RTree<DATATYPE, ELEMTYPE, NUMDIMS, ELEMTYPEREAL, TMAXNODES, TMINNODES>
#define MAXH 1000000000
#define MINH -1000000000

// This function checks for intersection between a 3D segment and elements within a R-tree.
// It takes the start and end points of the segment as arrays of three elements (ELEMTYPE),
// and a callback function that is called when an intersection is detected.
RTREE_TEMPLATE
bool RTREE_QUAL::Intersect3d(const ELEMTYPE s_start[3], const ELEMTYPE s_end[3],
                             std::function<bool(const DATATYPE &, const ELEMTYPE[3], const ELEMTYPE[3], ELEMTYPE &)> callback,ELEMTYPE& outDist) const
{
    Seg3d seg;
    for (int axis = 0; axis < NUMDIMS; ++axis)
    {
        seg.s_start[axis] = s_start[axis];
        seg.s_end[axis] = s_end[axis];
    }
    bool IntersectFlag = false; // Initialize a flag to track whether an intersection has been found.

    // Recursively check for intersections starting from the root of the 3D R-tree.
    Intersect3d(m_root, &seg, IntersectFlag, callback,outDist);
    //std::cout<<outDist<<std::endl;

    return IntersectFlag;
}

// This function retrieves the height at a specific 3D location using a 3D R-tree structure.
// It takes the x and y coordinates of the location and a callback function that processes the data.
// The function returns the height value at the given location.
RTREE_TEMPLATE
ELEMTYPE RTREE_QUAL::Getheight3d(const ELEMTYPE x, const ELEMTYPE y,
                                 std::function<ELEMTYPE(const DATATYPE &, const ELEMTYPE, const ELEMTYPE)> callback) const
{
    double Height = MINH - 1;
    // Recursively search for the height at the given x, y coordinates starting from the root of the 3D R-tree.
    Getheight3d(m_root, x, y, Height, callback);
    return Height;
}

// This function checks for overlap between two 2D line segments.
// It takes the coordinates of the endpoints of the first segment (x1, y1, x2, y2)
// and the coordinates of the endpoints of the second segment (tx1, ty1, tx2, ty2).
// The function returns true if the segments overlap, otherwise false.
RTREE_TEMPLATE
inline bool RTREE_QUAL::OverlapSeg2d(ELEMTYPE x1, ELEMTYPE y1, ELEMTYPE x2, ELEMTYPE y2,
                                     ELEMTYPE tx1, ELEMTYPE ty1, ELEMTYPE tx2, ELEMTYPE ty2) const
{
    if ((x1 > x2 ? x1 : x2) < (tx1 < tx2 ? tx1 : tx2) ||
        (y1 > y2 ? y1 : y2) < (ty1 < ty2 ? ty1 : ty2) ||
        (tx1 > tx2 ? tx1 : tx2) < (x1 < x2 ? x1 : x2) ||
        (ty1 > ty2 ? ty1 : ty2) < (y1 < y2 ? y1 : y2))
        return false;
    // Check for overlap using the cross product method.
    if ((((x1 - tx1) * (ty2 - ty1) - (y1 - ty1) * (tx2 - tx1)) *
         ((x2 - tx1) * (ty2 - ty1) - (y2 - ty1) * (tx2 - tx1))) > 0 ||
        (((tx1 - x1) * (y2 - y1) - (ty1 - y1) * (x2 - x1)) *
         ((tx2 - x1) * (y2 - y1) - (ty2 - y1) * (x2 - x1))) > 0)
        return false;

    // If none of the above conditions are met, the segments overlap.
    return true;
}

// This function checks for an overlap between a 2D line segment and a rectangle.
// It takes the coordinates of the endpoints of the line segment (x1, y1, x2, y2)
// and the minimum and maximum x and y coordinates of the rectangle (minx, miny, maxx, maxy).
// The function returns true if there is an overlap, otherwise false.

RTREE_TEMPLATE
inline bool RTREE_QUAL::OverlapSegRect2d(ELEMTYPE x1, ELEMTYPE y1, ELEMTYPE x2, ELEMTYPE y2,
                                         ELEMTYPE minx, ELEMTYPE miny, ELEMTYPE maxx, ELEMTYPE maxy) const
{
    // Check if any endpoint of the segment is inside the rectangle.
    if ((x1 >= minx && x1 <= maxx && y1 >= miny && y1 <= maxy) ||
        (x2 >= minx && x2 <= maxx && y2 >= miny && y2 <= maxy))
        return true;

    // Check for overlap between the segment and the rectangle using the OverlapSeg2d function.
    // This checks for proper overlap, not just touching at the corners.
    return OverlapSeg2d(x1, y1, x2, y2, minx, miny, maxx, maxy) ||
           OverlapSeg2d(x1, y1, x2, y2, minx, maxy, maxx, miny);
}

// This function checks for overlap between a 3D segment and a 3D rectangle.
// It takes pointers to a 3D segment (a_segA) and a 3D rectangle (a_rectB).
// The function returns true if there is an overlap, otherwise false.
RTREE_TEMPLATE
bool RTREE_QUAL::Overlap3d(Seg3d *a_segA, Rect *a_rectB) const
{
    ASSERT(a_segA && a_rectB); // Ensure that the pointers to the segment and rectangle are valid.

    // Check for overlap between the segment's XY plane projection and the rectangle's XY plane projection.
    if (!(OverlapSegRect2d(a_segA->s_start[0], a_segA->s_start[1],
                           a_segA->s_end[0], a_segA->s_end[1],
                           a_rectB->m_min[0], a_rectB->m_min[1],
                           a_rectB->m_max[0], a_rectB->m_max[1])))
        return false;
    // Check for overlap between the segment's XZ plane projection and the rectangle's XZ plane projection.
    if (!(OverlapSegRect2d(a_segA->s_start[0], a_segA->s_start[2],
                           a_segA->s_end[0], a_segA->s_end[2],
                           a_rectB->m_min[0], a_rectB->m_min[2],
                           a_rectB->m_max[0], a_rectB->m_max[2])))
        return false;
    // Check for overlap between the segment's YZ plane projection and the rectangle's YZ plane projection.
    if (!(OverlapSegRect2d(a_segA->s_start[1], a_segA->s_start[2],
                           a_segA->s_end[1], a_segA->s_end[2],
                           a_rectB->m_min[1], a_rectB->m_min[2],
                           a_rectB->m_max[1], a_rectB->m_max[2])))
        return false;
    // If all projections overlap, then the 3D segment and rectangle overlap.
    return true;
}

// This function checks for intersection between a 3D segment and the nodes of an R-tree.
RTREE_TEMPLATE
bool RTREE_QUAL::Intersect3d(Node *a_node, Seg3d *a_seg, bool &IntersectFlag,
                             std::function<bool(const DATATYPE &, const ELEMTYPE[NUMDIMS], const ELEMTYPE[NUMDIMS], ELEMTYPE &)> callback,ELEMTYPE& outDist) const
{
    
    // Ensure that the node and segment pointers are valid.
    ASSERT(a_node);
    ASSERT(a_node->m_level >= 0);
    ASSERT(a_seg);
    // If the current node is an internal node (not a leaf), recursively check its children.
    if (a_node->IsInternalNode())
    {

        // Check for overlap between the segment and the bounding rectangle of the current branch.
        for (int index = 0; index < a_node->m_count; ++index)
        {
            if (Overlap3d(a_seg, &a_node->m_branch[index].m_rect))
            {
                if (!Intersect3d(a_node->m_branch[index].m_child, a_seg, IntersectFlag, callback,outDist))
                {
                    return false; // If an intersection is found, stop the search.
                }
            }
        }
    }
    else
    {
        // If the current node is a leaf node, check for intersection with the data it contains.
        for (int index = 0; index < a_node->m_count; ++index)
        {
            // Check for overlap between the segment and the bounding rectangle of the current branch.
            if (Overlap3d(a_seg, &a_node->m_branch[index].m_rect))
            {
                // Retrieve the data associated with the current branch.
                DATATYPE &id = a_node->m_branch[index].m_data;
                // Variable to store the distance from the start point to the intersection point
                //double dise =0.0;
                // Call the callback function with the data, start, and end points of the segment.
                if (callback && !callback(id, a_seg->s_start, a_seg->s_end, outDist))
                {
                    //std::cout << "交点距离: " << outDist << std::endl;
                    IntersectFlag = true; // Set the intersection flag if the callback returns false (indicating an intersection).
                    return false;         // Stop the search.
                }
            }
        }
    }
    // If no intersections were found, continue the search.
    return true;
}

// This function retrieves the maximum height at a specific 2D location (x, y) within an 3D R-tree.
RTREE_TEMPLATE
bool RTREE_QUAL::Getheight3d(Node *a_node, ELEMTYPE x, ELEMTYPE y, ELEMTYPE &Height,
                             std::function<ELEMTYPE(const DATATYPE &, const ELEMTYPE, const ELEMTYPE)> callback) const
{
    ASSERT(a_node);
    ASSERT(a_node->m_level >= 0);
    bool found = false;

    if (a_node->IsInternalNode())
    {

        // If the current node is an internal node (not a leaf), recursively check its children.
        for (int index = 0; index < a_node->m_count; ++index)
        {
            // Check if the location (x, y) is within the bounding rectangle of the current branch.
            if (!(x < (&a_node->m_branch[index].m_rect)->m_min[0]) &&
                !(y < (&a_node->m_branch[index].m_rect)->m_min[1]) &&
                !(x > (&a_node->m_branch[index].m_rect)->m_max[0]) &&
                !(y > (&a_node->m_branch[index].m_rect)->m_max[1]))
            {
                if (!Getheight3d(a_node->m_branch[index].m_child, x, y, Height, callback))
                {
                    found = true;
                }
            }
        }
    }
    else
    {
        // std::cout << "333" << " ";
        //  If the current node is a leaf node, check for the maximum height at the location (x, y).
        for (int index = 0; index < a_node->m_count; ++index)
        {
            if (!(x < (&a_node->m_branch[index].m_rect)->m_min[0]) &&
                !(y < (&a_node->m_branch[index].m_rect)->m_min[1]) &&
                !(x > (&a_node->m_branch[index].m_rect)->m_max[0]) &&
                !(y > (&a_node->m_branch[index].m_rect)->m_max[1]))
            {
                // Retrieve the data associated with the current branch.
                DATATYPE &id = a_node->m_branch[index].m_data;
                {

                    ELEMTYPE h = callback(id, x, y); // h is the actural height of xy
                    // std::cout << "h= " << h << "  ";
                    //  If the height returned by the callback is greater than the current maximum, update the maximum height.
                    if (h > Height)
                    {
                        Height = h;
                        found = true; // Stop the search
                    }
                }
            }
        }
    }
    return found;
}

using namespace std;

using DATATYPE = std::tuple<double, double, double, double, double, double, double, double, double>;

struct Rect3d
{
    Rect3d() {}
    Rect3d(double minX, double minY, double minZ, double maxX, double maxY, double maxZ)
    {
        min[0] = minX;
        min[1] = minY;
        min[2] = minZ;
        max[0] = maxX;
        max[1] = maxY;
        max[2] = maxZ;
    }
    double min[3];
    double max[3];
};
// various operations
class Vec3d
{
public:
    double x, y, z;
    // Calculate the dot product of this vector with another vector b.
    double dot(const Vec3d &b)
    {
        return Vec3d::x * b.x + Vec3d::y * b.y + Vec3d::z * b.z;
    }
    // Calculate the cross product of this vector with another vector b.
    Vec3d cross(const Vec3d &b)
    {
        return Vec3d(
            Vec3d::y * b.z - Vec3d::z * b.y,
            Vec3d::z * b.x - Vec3d::x * b.z,
            Vec3d::x * b.y - Vec3d::y * b.x);
    }


    // 计算向量的欧氏长度（模长）
    double length() const
    {
        return sqrt(x * x + y * y + z * z);
    }

    // Normalize this vector to make its length equal to 1.
    Vec3d normalize()
    {
        const double s = 1.0f / sqrtf(Vec3d::x * Vec3d::x + Vec3d::y * Vec3d::y + Vec3d::z * Vec3d::z);
        return Vec3d(Vec3d::x * s, Vec3d::y * s, Vec3d::z * s);
    }
    Vec3d operator+(const Vec3d &b)
    {
        return Vec3d(
            Vec3d::x + b.x,
            Vec3d::y + b.y,
            Vec3d::z + b.z);
    }
    Vec3d operator+=(const Vec3d &b)
    {
        *this = Vec3d::operator+(b);
        return *this;
    }
    Vec3d operator-(const Vec3d &b)
    {
        return Vec3d(
            Vec3d::x - b.x,
            Vec3d::y - b.y,
            Vec3d::z - b.z);
    }
    Vec3d operator-=(const Vec3d &b)
    {
        *this = Vec3d::operator-(b);
        return *this;
    }
    Vec3d operator*(const Vec3d &b)
    {
        return Vec3d(
            Vec3d::x * b.x,
            Vec3d::y * b.y,
            Vec3d::z * b.z);
    }
    Vec3d operator*=(const Vec3d &b)
    {
        *this = Vec3d::operator*(b);
        return *this;
    }
    Vec3d operator*(double b)
    {
        return Vec3d(
            Vec3d::x * b,
            Vec3d::y * b,
            Vec3d::z * b);
    }
    Vec3d operator*=(double b)
    {
        *this = Vec3d::operator*(b);
        return *this;
    }
    Vec3d operator/(const Vec3d &b)
    {
        return Vec3d(
            Vec3d::x / b.x,
            Vec3d::y / b.y,
            Vec3d::z / b.z);
    }
    Vec3d operator/=(const Vec3d &b)
    {
        *this = Vec3d::operator/(b);
        return *this;
    }
    Vec3d operator/(double b)
    {
        return Vec3d(
            Vec3d::x * b,
            Vec3d::y * b,
            Vec3d::z * b);
    }
    Vec3d operator/=(double b)
    {
        *this = Vec3d::operator/(b);
        return *this;
    }
    Vec3d(double x, double y, double z)
    {
        Vec3d::x = x;
        Vec3d::y = y;
        Vec3d::z = z;
    }
    Vec3d(double x)
    {
        Vec3d::x = x;
        Vec3d::y = x;
        Vec3d::z = x;
    }
    Vec3d()
    {
        //
    }
    ~Vec3d()
    {
        //
    }
};
#define EPSILON 0.000001f
// Check if a line segment intersects with a triangle.
// Möller-Trumbore算法
bool lineSegIntersectTri(Vec3d line[2], Vec3d tri[3], Vec3d *point, double *outDist)
{
    // Calculate the two edges of the triangle and the direction vector of the line segment.
    Vec3d e0 = tri[1] - tri[0];
    Vec3d e1 = tri[2] - tri[0];
    Vec3d dir = line[1] - line[0];
    Vec3d dir_norm = dir.normalize();
    Vec3d h = dir_norm.cross(e1);
    const double a = e0.dot(h);

    if (a > -EPSILON && a < EPSILON)
    {
        return false;
    }

    Vec3d s = line[0] - tri[0];
    const double f = 1.0f / a;
    const double u = f * s.dot(h);

    if (u < 0.0f || u > 1.0f)
    {
        return false;
    }

    Vec3d q = s.cross(e0);
    const double v = f * dir_norm.dot(q);

    if (v < 0.0f || u + v > 1.0f)
    {
        return false;
    }

    const double t = f * e1.dot(q);

    // if (t > EPSILON && t < sqrtf(dir.dot(dir))) // segment intersection
    // {
    //     if (point)
    //     {
    //         *point = line[0] + dir_norm * t;
    //     }

    //     return true;
    // }
    const double segLength = dir.length();
    if (t > EPSILON && t <= segLength)
    { // 线段相交
        if (point)
            *point = line[0] + dir_norm * t;
        if (outDist)
            *outDist = t; // 返回交点距离线段起点的长度
        return true;
    }
    return false;
}

    // This function is a callback used to determine the height at a given (x, y) location.
    // template <class DATATYPE, class ELEMTYPE, int NUMDIMS, class ELEMTYPEREAL, int TMAXNODES, int TMINNODES>

    // RTREE_TEMPLATE
    // inline ELEMTYPE RTREE_QUAL::HeightCallback(const DATATYPE &id, const ELEMTYPE x, const ELEMTYPE y)
    // {
    //     // Create a triangle from the given id, which contains vertex positions.
    //     Vec3d tri[3] =
    //         {
    //             {get<0>(id), get<1>(id), get<2>(id)},
    //             {get<3>(id), get<4>(id), get<5>(id)},
    //             {get<6>(id), get<7>(id), get<8>(id)},
    //         };
    //     // Define a line segment from (x, y, MAXH) to (x, y, MINH).
    //     Vec3d line[2] =
    //         {
    //             {x, y, MAXH},
    //             {x, y, MINH},
    //         };
    //     Vec3d *point = new Vec3d();
    //     // Check if the line segment intersects with the triangle.
    //     if (lineSegIntersectTri(line, tri, point))
    //     {
    //         std::cout<<"222222"<<" ";
    //         return point->z;
    //     }
    //     // If there is no intersection, return a value indicating to continue searching.
    //     return MINH - 1;
    // }

    // // This function is a callback used to check for intersection between a line segment and a triangle.

    // template <class DATATYPE, class ELEMTYPE, int NUMDIMS, class ELEMTYPEREAL, int TMAXNODES, int TMINNODES>
    // inline bool IntersectCallback(const DATATYPE &id, const ELEMTYPE a[3], const ELEMTYPE b[3], ELEMTYPE &outDist)
    // {
    //     outDist = 0;
    //     Vec3d tri[3] =
    //         {
    //             {get<0>(id), get<1>(id), get<2>(id)},
    //             {get<3>(id), get<4>(id), get<5>(id)},
    //             {get<6>(id), get<7>(id), get<8>(id)},
    //         };
    //     Vec3d line[2] =
    //         {
    //             {a[0], a[1], a[2]},
    //             {b[0], b[1], b[2]},
    //         };
    //     // If an intersection is found, return false to indicate the intersection.
    //     if (lineSegIntersectTri(line, tri, NULL))
    //     {
    //         return false;
    //     }
    //     // If no intersection is found, return true to indicate continue searching.
    //     return true;
    // }

    double HeightCallback(DATATYPE id, const double x, const double y)
    {
        // Create a triangle from the given id, which contains vertex positions.
        Vec3d tri[3] =
            {
                {get<0>(id), get<1>(id), get<2>(id)},
                {get<3>(id), get<4>(id), get<5>(id)},
                {get<6>(id), get<7>(id), get<8>(id)},
            };
        Vec3d line[2] =
            {
                {x, y, MAXH},
                {x, y, MINH},
            };
        //Vec3d *point = new Vec3d();
        Vec3d point; // 这是正确的代码，在栈上创建对象
        // Check if the line segment intersects with the triangle.
        // if (lineSegIntersectTri(line, tri, point,nullptr))
        // {
        //     return point->z;
        // }
        if (lineSegIntersectTri(line, tri, &point, nullptr)) // 这是正确的代码
        {
            return point.z; // // 访问成员变量的方式从 -> 变为 .
        }
        // If there is no intersection, return a value indicating to continue searching.
        return MINH - 1;
    }

    // This function is a callback used to check for intersection between a line segment and a triangle.

    bool IntersectCallback(DATATYPE id, const double a[3], const double b[3], double &outDist)
    {

        Vec3d tri[3] =
            {
                {get<0>(id), get<1>(id), get<2>(id)},
                {get<3>(id), get<4>(id), get<5>(id)},
                {get<6>(id), get<7>(id), get<8>(id)},
            };
        Vec3d line[2] =
            {
                {a[0], a[1], a[2]},
                {b[0], b[1], b[2]},
            };
        Vec3d intersectionPoint;
        double dis=-0.1;
        // If an intersection is found, return false to indicate the intersection.
        if (lineSegIntersectTri(line, tri, &intersectionPoint, &dis))
        {
            outDist = dis; // 将距离赋值给 outDist
            return false;
        }
        // If no intersection is found, return true to indicate continue searching.
        return true;
    }

    // 运行Q-view.cpp文件时需要打开此代码
    //  auto heightCallback = [](const DATATYPE &id, double x, double y)
    //  {
    //      return HeightCallback(id, x, y);
    //  };
    //  auto intersectCallback = [](const auto &id, const auto &a, const auto &b, auto &outDist)
    //  {
    //      return IntersectCallback(id, a, b, outDist);
    //  };

#endif
    // QUERY3DRTREE_H
