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

#ifndef CARTOGRAPHER_MAPPING_INTERNAL_2D_RAY_TO_PIXEL_MASK_H_
#define CARTOGRAPHER_MAPPING_INTERNAL_2D_RAY_TO_PIXEL_MASK_H_

#include <vector>

#include "cartographer/transform/transform.h"
#include "cartographer/mapping/2d/probability_grid.h"
namespace cartographer {
namespace mapping {

// Compute all pixels that contain some part of the line segment connecting
// 'scaled_begin' and 'scaled_end'. 'scaled_begin' and 'scaled_end' are scaled
// by 'subpixel_scale'. 'scaled_begin' and 'scaled_end' are expected to be
// greater than zero. Return values are in pixels and not scaled.
std::vector<Eigen::Array2i> RayToPixelMask(const Eigen::Array2i& scaled_begin,
                                           const Eigen::Array2i& scaled_end,
                                           int subpixel_scale);
std::vector<Eigen::Array2i> RayToPixelMaskVisualNew(const Eigen::Array2i& scaled_begin,
                                           const Eigen::Array2i& scaled_end,
                                           int subpixel_scale,
                                           const uint8 intensities,
                                           const transform::Rigid3d local_pose_,
                                           int flag,
                                           const std::vector<uint16>& hit_table,
                                           cartographer::mapping::ProbabilityGrid* probability_grid);
std::vector<Eigen::Array2i> RayToPixelMaskVisual(const Eigen::Array2i& scaled_begin,
                                           const Eigen::Array2i& scaled_end,
                                           int subpixel_scale,
                                           const uint8 intensities);
std::vector<Eigen::Array2i> ProcessPixelMask(std::vector<Eigen::Array2i> pixel_mask,
                                           const Eigen::Array2i& scaled_begin,
                                           const Eigen::Array2i& scaled_end,
                                           int subpixel_scale,
                                           int flag,
                                           const uint8 intensities,
                                           const transform::Rigid3d local_pose_);
std::vector<Eigen::Array2i> ProcessPixelMaskHit(std::vector<Eigen::Array2i> pixel_mask,
                                           const Eigen::Array2i& scaled_begin,
                                           const Eigen::Array2i& scaled_end,
                                           int subpixel_scale,
                                           int flag,
                                           const uint8 intensities,
                                           const transform::Rigid3d local_pose_,
                                            const std::vector<uint16>& hit_table,
                                         cartographer::mapping::ProbabilityGrid* probability_grid);
std::vector<Eigen::Array2i> LocalToMap(const Eigen::Array2i& pixel,
              const Eigen::Array2i& begin,
              const Eigen::Array2i& end,
              int subpixel_scale,
              int flag,
              const transform::Rigid3d local_pose_);


}  // namespace mapping
}  // namespace cartographer

#endif  // CARTOGRAPHER_MAPPING_INTERNAL_2D_RAY_TO_PIXEL_MASK_H_
