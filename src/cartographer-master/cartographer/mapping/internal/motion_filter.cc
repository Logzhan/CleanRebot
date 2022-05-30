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
 * Description : cartographer运动滤波器相关功能
                 1、判断当前状态是否和上次状态非常接近，根据判断的结果决定是否插入
				 子图
 * Date        : 2022-05-29 wyk & zhanli
 */

#include "cartographer/mapping/internal/motion_filter.h"

#include "cartographer/transform/transform.h"
#include "glog/logging.h"
#include <iostream>
#include <math.h>
namespace cartographer {
namespace mapping {

proto::MotionFilterOptions CreateMotionFilterOptions(
    common::LuaParameterDictionary* const parameter_dictionary) {
  proto::MotionFilterOptions options;
  options.set_max_time_seconds(
      parameter_dictionary->GetDouble("max_time_seconds"));
  options.set_max_distance_meters(
      parameter_dictionary->GetDouble("max_distance_meters"));
  options.set_max_angle_radians(
      parameter_dictionary->GetDouble("max_angle_radians"));
  return options;
}

MotionFilter::MotionFilter(const proto::MotionFilterOptions& options)
    : options_(options) {}


bool MotionFilter::IsSimilar(const common::Time time,
                             const transform::Rigid3d& pose) {
							 
  LOG_IF_EVERY_N(INFO, num_total_ >= 500, 500)
      << "Motion filter reduced the number of nodes to "
      << 100. * num_different_ / num_total_ << "%.";
  ++num_total_;

  // 如果角度变化较大，此时摄像头结合雷达生成的语义标准，此时插入到Submap中
  // 容易产生错误的地图
  if (transform::GetAngle(pose.inverse() * last_pose_) > 0.02) {
   last_time_ = time;
   last_pose_ = pose;
    return true;
  }
  
  if (num_total_ > 1 &&
      // 原始版本是common::FromSeconds(options_.max_time_seconds())
      time - last_time_ <= common::FromSeconds(0.8) &&
      (pose.translation() - last_pose_.translation()).norm() <=
          options_.max_distance_meters() &&
	  // 计算上当前角度和上一次角度的差值(弧度radians)
      transform::GetAngle(pose.inverse() * last_pose_) <=
          options_.max_angle_radians()) {
    return true;
  }
  last_time_ = time;
  last_pose_ = pose;
  // 统计不同雷达帧的数量，只有返回false的时候才更新num_different_
  ++num_different_;
  return false;
}

bool MotionFilter::IsSimilarNew(const common::Time time,
                             const transform::Rigid3d& pose,
                             const sensor::RangeData& range_data) {
  LOG_IF_EVERY_N(INFO, num_total_ >= 500, 500)
      << "Motion filter reduced the number of nodes to "
      << 100. * num_different_ / num_total_ << "%.";
  ++num_total_;
  int flag=0;
  for(int i=0;i<range_data.returns.intensities().size();i++){
      if(range_data.returns.intensities()[i]==115 || range_data.returns.intensities()[i]==117 || range_data.returns.intensities()[i]==200 ||
      range_data.returns.intensities()[i]==201 || range_data.returns.intensities()[i]==121  )
        flag=1;
  }
  if (num_total_ > 1 &&
      time - last_time_ <= common::FromSeconds(options_.max_time_seconds()) &&
      (pose.translation() - last_pose_.translation()).norm() <=
          options_.max_distance_meters() &&
      transform::GetAngle(pose.inverse() * last_pose_) <=
          options_.max_angle_radians() &&  flag==0) {
    return true;
  }
 if ( time - last_time_ >= common::FromSeconds(0.1) &&  flag==1) {
        last_time_ = time;
        last_pose_ = pose;
        ++num_different_;
        return false;
  }else{
      return true;
  }
  last_time_ = time;
  last_pose_ = pose;
  ++num_different_;
  return false;
}

}  // namespace mapping
}  // namespace cartographer
