#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "query3dRtree.h"
#include <pybind11/numpy.h>
#include <pybind11/functional.h>

namespace py = pybind11;
// 假设 ELEMTYPE 和 DATATYPE 是已定义的类型，这里简单定义一下
//using DATATYPE = std::array<double, 9>; // 通用的字节数组
using ELEMTYPE = double;                // 坐标类型
constexpr int NUMDIMS = 3;              // 空间维度
using ELEMTYPEREAL = double;            // 实数类型
constexpr int TMAXNODES = 8;            // 最大节点数
constexpr int TMINNODES = 4;            // 最小节点数

class RTreeWrapper : public RTree<DATATYPE, ELEMTYPE, NUMDIMS, ELEMTYPEREAL, TMAXNODES, TMINNODES>
{
public:
     // RTreeType 是 RTree 的特定实例，配置了数据类型、维度等。
     //using RTreeType = RTree<DATATYPE, ELEMTYPE, NUMDIMS, ELEMTYPEREAL, TMAXNODES, TMINNODES>;
     using RTreeType = RTree<DATATYPE, double,3, double,8, 4>;
     bool LoadFromFile(const std::string &file_name)
     {
          std::cout<<"loadFilefunction"<<std::endl;
          return RTreeType::Load(file_name.c_str());
     }
     // 预定义的回调函数
     static double DefaultHeightCallback(const DATATYPE &id, ELEMTYPE x, ELEMTYPE y)
     {
         
          return HeightCallback(id, x, y);
     }

     static bool DefaultIntersectCallback(const DATATYPE &id, const ELEMTYPE *a, const ELEMTYPE *b, double &outDist)
     {    
          return IntersectCallback(id, a, b, outDist);
     }

     std::tuple<bool,double> Intersect3d(const std::array<ELEMTYPE, NUMDIMS> &s_start,
                      const std::array<ELEMTYPE, NUMDIMS> &s_end) const
     {
          // 直接调用内部实现，不再需要回调
          //return RTreeType::Intersect3d(s_start.data(), s_end.data(), DefaultIntersectCallback,outDist);
          double outDist = 0;
          bool result = RTreeType::Intersect3d(s_start.data(), s_end.data(), DefaultIntersectCallback, outDist);
          return std::make_tuple(result, outDist);
     }

     ELEMTYPE GetHeight3d(ELEMTYPE x, ELEMTYPE y) const
     {
          // 直接调用内部实现，不再需要回调
          return RTreeType::Getheight3d(x, y, DefaultHeightCallback);
     }


};

PYBIND11_MODULE(my_rtree, m)
{
     // 绑定 RTree
     py::class_<RTreeWrapper>(m, "RTree")
         .def(py::init<>(), "Initialize an RTree object.")
         .def("load_from_file", &RTreeWrapper::LoadFromFile, py::arg("file_name"))
         .def("intersect3d", &RTreeWrapper::Intersect3d, py::arg("s_start"), py::arg("s_end"))
         .def("get_height_3d", &RTreeWrapper::GetHeight3d, py::arg("x"), py::arg("y"));

         
}

