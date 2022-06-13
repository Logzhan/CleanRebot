/*
 * Copyright 2016 The Cartographer Authors
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "cartographer/mapping/2d/probability_grid_range_data_inserter_2d.h"
#include <cstdlib>
#include "Eigen/Core"
#include "Eigen/Geometry"
#include "cartographer/mapping/2d/xy_index.h"
#include "cartographer/mapping/internal/2d/ray_to_pixel_mask.h"
#include "cartographer/mapping/probability_values.h"
#include "glog/logging.h"
extern int cnew[4000000];
// static int x_last=0;
// static int y_last=0;
// static int xm_last=0;
// static int ym_last=0;

namespace cartographer {
namespace mapping {
namespace {

// Factor for subpixel accuracy of start and end point for ray casts.
constexpr int kSubpixelScale = 1000;

// 根据点云的bounding box, 看是否需要对地图进行扩张
void GrowAsNeeded(const sensor::RangeData& range_data,
                  ProbabilityGrid* const probability_grid) {
  // 找到点云的bounding_box
  Eigen::AlignedBox2f bounding_box(range_data.origin.head<2>());
  // Padding around bounding box to avoid numerical issues at cell boundaries.
  constexpr float kPadding = 1e-6f;
  for (const sensor::RangefinderPoint& hit : range_data.returns) {
    bounding_box.extend(hit.position.head<2>());
  }
  for (const sensor::RangefinderPoint& miss : range_data.misses) {
    bounding_box.extend(miss.position.head<2>());
  }
  // 是否对地图进行扩张
  probability_grid->GrowLimits(bounding_box.min() -
                               kPadding * Eigen::Vector2f::Ones());
  probability_grid->GrowLimits(bounding_box.max() +
                               kPadding * Eigen::Vector2f::Ones());
}

/**
 * @brief 根据雷达点对栅格地图进行更新
 * 
 * @param[in] range_data 
 * @param[in] hit_table 更新占用栅格时的查找表
 * @param[in] miss_table 更新空闲栅格时的查找表
 * @param[in] insert_free_space 
 * @param[in] probability_grid 栅格地图
 */

//...

Eigen::Affine3d ToEigen(const ::cartographer::transform::Rigid3d& rigid3) {
  return Eigen::Translation3d(rigid3.translation()) * rigid3.rotation();
}
void CastRays(const sensor::RangeData& range_data,
              const std::vector<uint16>& hit_table,
              const std::vector<uint16>& miss_table,
              const bool insert_free_space, 
              ProbabilityGrid* probability_grid,
              const transform::Rigid3d local_pose_) {
  // 根据雷达数据调整地图范围
  GrowAsNeeded(range_data, probability_grid);

  const MapLimits& limits = probability_grid->limits();
  const double superscaled_resolution = limits.resolution() / kSubpixelScale;
  const MapLimits superscaled_limits(
      superscaled_resolution, limits.max(),
      CellLimits(limits.cell_limits().num_x_cells * kSubpixelScale,
                 limits.cell_limits().num_y_cells * kSubpixelScale));
  // 雷达原点在地图中的像素坐标, 作为画线的起始坐标
  const Eigen::Array2i begin =
      superscaled_limits.GetCellIndex(range_data.origin.head<2>());
  // Compute and add the end points.
  std::vector<Eigen::Array2i> ends;
  ends.reserve(range_data.returns.size());
  //。。。for (const sensor::RangefinderPoint& hit : range_data.returns) {
  for (int i=0;i<range_data.returns.size();i++) {
    // 计算hit点在地图中的像素坐标, 作为画线的终止点坐标
    //。。。ends.push_back(superscaled_limits.GetCellIndex(hit.position.head<2>()));
    ends.push_back(superscaled_limits.GetCellIndex((range_data.returns.points())[i].position.head<2>()));
    // 更新hit点的栅格值
    //新写函数HitpplyLookupTable，增加入口参数range_data.returns
    //。。。probability_grid->ApplyLookupTable(ends.back()/kSubpixelScale, hit_table);
    probability_grid->ApplyLookupTableHit(ends.back()/kSubpixelScale, hit_table,uint8(range_data.returns.intensities()[i]));
  }
  
  // 如果不插入free空间就可以结束了
  if (!insert_free_space) {
    return;
  }

  const Eigen::Matrix4d homo =ToEigen(local_pose_).matrix();           
  int o_x=1000;
  int o_y=1000;
  int xm;
  int ym;
  if (begin.y()/1000>150){
   xm=-begin.y()/1000+20*homo(0,3)+(o_x+200);
   ym=begin.x()/1000-20*homo(1,3)+(o_y-200);

  // LOG(INFO)<<"wyk--guiji jububiao:"<<begin.x()/1000<<"         "<<begin.y()/1000;
  // LOG(INFO)<<"wyk--guiji homo:"<<homo(0,3)<<"      "<<homo(1,3);
  // LOG(INFO)<<"wyk--guiji quanjuzuobiao:"<<xm<<"         "<<ym;
  }else{
     xm=-begin.y()/1000+20*homo(0,3)+(o_x+100);
     ym=begin.x()/1000-20*homo(1,3)+(o_y-100);

  // LOG(INFO)<<"wyk--guiji jububiao:"<<begin.x()/1000<<"         "<<begin.y()/1000;
  // LOG(INFO)<<"wyk--guiji homo:"<<homo(0,3)<<"      "<<homo(1,3);
  // LOG(INFO)<<"wyk--guiji quanjuzuobiao:"<<xm<<"         "<<ym;
  }
  


  
 
  //LOG(INFO)<<"wyk--guiji zuobiao_last:"<<x_last<<"         "<<y_last;


    // if(x_last==0 && y_last==0){
    //       x_last=xm;
    //       y_last=ym;
    //     }else{
    // if(1){
    //         x_last=xm;
    //         y_last=ym;
    //         LOG(INFO)<<"wyk--guiji quanjuzuobiao:"<<xm<<"         "<<ym;
    //         LOG(INFO)<<"wyk--guiji zuobiao_last:"<<x_last<<"         "<<y_last;
    //         }
    //     }



  std::vector<Eigen::Array2i> end;
  end.reserve(range_data.returns.size());
  // Now add the misses.
  for (int i=0;i<range_data.returns.size();i++) {
    //返回一个vector miss点坐标的集合 如果在视线范围内看的到就进行光线跟踪，否则返回空vector 不做光线跟踪
    end.push_back(superscaled_limits.GetCellIndex((range_data.returns.points())[i].position.head<2>()));
    //如果需要使非摄像头视域的小物体保留下来，则需要当前时刻的整张地图信息以及当前时刻该子图的位姿
    // LOG(INFO)<<"wyk--guiji begin:"<<begin.x()/1000<<"         "<<begin.y()/1000;
    std::vector<Eigen::Array2i> ray =
    RayToPixelMaskVisualNew(begin, end.back(), kSubpixelScale,uint8(range_data.returns.intensities()[i]),local_pose_,0,hit_table,probability_grid);
    if (ray.size()==0)
        continue;
    //if(ray.size()==0)
    //LOG(INFO)<<"test!!!!!!!!!!!!!!!!!is:"<<ray.size();
    for (const Eigen::Array2i& cell_index : ray) {
      // 从起点到end点之前, 更新miss点的栅格值
      probability_grid->ApplyLookupTable(cell_index, miss_table);
    }
  }

  // // Now add the misses.
  // for (const Eigen::Array2i& end : ends) {
  //   //返回一个vector miss点坐标的集合 如果在视线范围内看的到就进行光线跟踪，否则返回空vector 不做光线跟踪
  //   std::vector<Eigen::Array2i> ray =
  //       RayToPixelMask(begin, end, kSubpixelScale);
  //   for (const Eigen::Array2i& cell_index : ray) {
  //     // 从起点到end点之前, 更新miss点的栅格值
  //     probability_grid->ApplyLookupTable(cell_index, miss_table);
  //   }
  // }

  //Finally, compute and add empty rays based on misses in the range data.
  for (const sensor::RangefinderPoint& missing_echo : range_data.misses) {
    std::vector<Eigen::Array2i> ray = RayToPixelMask(
        begin, superscaled_limits.GetCellIndex(missing_echo.position.head<2>()),
        kSubpixelScale);
    for (const Eigen::Array2i& cell_index : ray) {
      // 从起点到misses点之前, 更新miss点的栅格值
      probability_grid->ApplyLookupTable(cell_index, miss_table);
    }
  }
}



void CastRays(const sensor::RangeData& range_data,
              const std::vector<uint16>& hit_table,
              const std::vector<uint16>& miss_table,
              const bool insert_free_space, ProbabilityGrid* probability_grid) {
  // 根据雷达数据调整地图范围
  GrowAsNeeded(range_data, probability_grid);

  const MapLimits& limits = probability_grid->limits();
  const double superscaled_resolution = limits.resolution() / kSubpixelScale;
  const MapLimits superscaled_limits(
      superscaled_resolution, limits.max(),
      CellLimits(limits.cell_limits().num_x_cells * kSubpixelScale,
                 limits.cell_limits().num_y_cells * kSubpixelScale));
  // 雷达原点在地图中的像素坐标, 作为画线的起始坐标
  const Eigen::Array2i begin =
      superscaled_limits.GetCellIndex(range_data.origin.head<2>());
  // Compute and add the end points.
  std::vector<Eigen::Array2i> ends;
  ends.reserve(range_data.returns.size());
  //。。。for (const sensor::RangefinderPoint& hit : range_data.returns) {
  for (int i=0;i<range_data.returns.size();i++) {
    // 计算hit点在地图中的像素坐标, 作为画线的终止点坐标
    //。。。ends.push_back(superscaled_limits.GetCellIndex(hit.position.head<2>()));
    ends.push_back(superscaled_limits.GetCellIndex((range_data.returns.points())[i].position.head<2>()));
    // 更新hit点的栅格值
    //新写函数HitpplyLookupTable，增加入口参数range_data.returns
    //。。。probability_grid->ApplyLookupTable(ends.back()/kSubpixelScale, hit_table);
    probability_grid->ApplyLookupTableHit(ends.back()/kSubpixelScale, hit_table,uint8(range_data.returns.intensities()[i]));
  }
  
  // 如果不插入free空间就可以结束了
  if (!insert_free_space) {
    return;
  }
    std::vector<Eigen::Array2i> end;
    end.reserve(range_data.returns.size());
    for (int i=0;i<range_data.returns.size();i++) {
    //返回一个vector miss点坐标的集合 如果在视线范围内看的到就进行光线跟踪，否则返回空vector 不做光线跟踪
    end.push_back(superscaled_limits.GetCellIndex((range_data.returns.points())[i].position.head<2>()));
    std::vector<Eigen::Array2i> ray =
        RayToPixelMaskVisual(begin, end.back(), kSubpixelScale,uint8(range_data.returns.intensities()[i]));
    if (ray.size()==0)
        continue;
    for (const Eigen::Array2i& cell_index : ray) {
      // 从起点到end点之前, 更新miss点的栅格值
      probability_grid->ApplyLookupTable(cell_index, miss_table);
    }
  }

  // Finally, compute and add empty rays based on misses in the range data.
  for (const sensor::RangefinderPoint& missing_echo : range_data.misses) {
    std::vector<Eigen::Array2i> ray = RayToPixelMask(
        begin, superscaled_limits.GetCellIndex(missing_echo.position.head<2>()),
        kSubpixelScale);
    for (const Eigen::Array2i& cell_index : ray) {
      // 从起点到misses点之前, 更新miss点的栅格值
      probability_grid->ApplyLookupTable(cell_index, miss_table);
    }
  }
}
}  // namespace

proto::ProbabilityGridRangeDataInserterOptions2D
CreateProbabilityGridRangeDataInserterOptions2D(
    common::LuaParameterDictionary* parameter_dictionary) {
  proto::ProbabilityGridRangeDataInserterOptions2D options;
  options.set_hit_probability(
      parameter_dictionary->GetDouble("hit_probability"));
  options.set_miss_probability(
      parameter_dictionary->GetDouble("miss_probability"));
  options.set_insert_free_space(
      parameter_dictionary->HasKey("insert_free_space")
          ? parameter_dictionary->GetBool("insert_free_space")
          : true);
  CHECK_GT(options.hit_probability(), 0.5);
  CHECK_LT(options.miss_probability(), 0.5);
  return options;
}

// 写入器的构造, 新建了2个查找表
ProbabilityGridRangeDataInserter2D::ProbabilityGridRangeDataInserter2D(
    const proto::ProbabilityGridRangeDataInserterOptions2D& options)
    : options_(options),
      // 生成更新占用栅格时的查找表 // param: hit_probability
      hit_table_(ComputeLookupTableToApplyCorrespondenceCostOdds(
          Odds(options.hit_probability()))),    // 0.55
      // 生成更新空闲栅格时的查找表 // param: miss_probability
      miss_table_(ComputeLookupTableToApplyCorrespondenceCostOdds(
          Odds(options.miss_probability()))) {} // 0.49

/**
 * @brief 将点云写入栅格地图
 * 
 * @param[in] range_data 要写入地图的点云
 * @param[in] grid 栅格地图
 */
void ProbabilityGridRangeDataInserter2D::Insert(
    const sensor::RangeData& range_data, GridInterface* const grid) const {
  //将Grid2D类型强制转化成其子类
  ProbabilityGrid* const probability_grid = static_cast<ProbabilityGrid*>(grid);
  CHECK(probability_grid != nullptr);
  // By not finishing the update after hits are inserted, we give hits priority
  // (i.e. no hits will be ignored because of a miss in the same cell).
  // param: insert_free_space
  CastRays(range_data, hit_table_, miss_table_, options_.insert_free_space(),
           probability_grid);
  probability_grid->FinishUpdate();
}

void ProbabilityGridRangeDataInserter2D::Insert_new(
    const sensor::RangeData& range_data, GridInterface* const grid,const transform::Rigid3d local_pose_) const {
  //将Grid2D类型强制转化成其子类 Grid2D没有子图位姿 该函数需要改变接口，得到submap的local_pose
  ProbabilityGrid* const probability_grid = static_cast<ProbabilityGrid*>(grid);
  CHECK(probability_grid != nullptr);
  // By not finishing the update after hits are inserted, we give hits priority
  // (i.e. no hits will be ignored because of a miss in the same cell).
  // param: insert_free_space
  CastRays(range_data, hit_table_, miss_table_, options_.insert_free_space(),
           probability_grid,local_pose_);
  probability_grid->FinishUpdate();
     
}


}  // namespace mapping
}  // namespace cartographer

