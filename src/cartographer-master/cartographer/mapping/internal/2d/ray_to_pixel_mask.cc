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

#include "cartographer/mapping/internal/2d/ray_to_pixel_mask.h"

#include "Eigen/Dense"
// 全局变量：从语义地图话题获取到的语义地图信息
int cnew[4000000]={0};
namespace cartographer {
namespace mapping {
namespace {

bool isEqual(const Eigen::Array2i& lhs, const Eigen::Array2i& rhs) {
  return ((lhs - rhs).matrix().lpNorm<1>() == 0);
}
}  // namespace

// Compute all pixels that contain some part of the line segment connecting
// 'scaled_begin' and 'scaled_end'. 'scaled_begin' and 'scaled_end' are scaled
// by 'subpixel_scale'. 'scaled_begin' and 'scaled_end' are expected to be
// greater than zero. Return values are in pixels and not scaled.
std::vector<Eigen::Array2i> RayToPixelMask(const Eigen::Array2i& scaled_begin,
                                           const Eigen::Array2i& scaled_end,
                                           int subpixel_scale) {
  // For simplicity, we order 'scaled_begin' and 'scaled_end' by their x
  // coordinate.
  if (scaled_begin.x() > scaled_end.x()) {
    return RayToPixelMask(scaled_end, scaled_begin, subpixel_scale);
  }

  CHECK_GE(scaled_begin.x(), 0);
  CHECK_GE(scaled_begin.y(), 0);
  CHECK_GE(scaled_end.y(), 0);
  std::vector<Eigen::Array2i> pixel_mask;
  // Special case: We have to draw a vertical line in full pixels, as
  // 'scaled_begin' and 'scaled_end' have the same full pixel x coordinate.
  if (scaled_begin.x() / subpixel_scale == scaled_end.x() / subpixel_scale) {
    Eigen::Array2i current(
        scaled_begin.x() / subpixel_scale,
        std::min(scaled_begin.y(), scaled_end.y()) / subpixel_scale);
    pixel_mask.push_back(current);
    const int end_y =
        std::max(scaled_begin.y(), scaled_end.y()) / subpixel_scale;
    for (; current.y() <= end_y; ++current.y()) {
      if (!isEqual(pixel_mask.back(), current)) pixel_mask.push_back(current);
    }
    return pixel_mask;
  }

  const int64 dx = scaled_end.x() - scaled_begin.x();
  const int64 dy = scaled_end.y() - scaled_begin.y();
  const int64 denominator = 2 * subpixel_scale * dx;

  // The current full pixel coordinates. We scaled_begin at 'scaled_begin'.
  Eigen::Array2i current = scaled_begin / subpixel_scale;
  pixel_mask.push_back(current);

  // To represent subpixel centers, we use a factor of 2 * 'subpixel_scale' in
  // the denominator.
  // +-+-+-+ -- 1 = (2 * subpixel_scale) / (2 * subpixel_scale)
  // | | | |
  // +-+-+-+
  // | | | |
  // +-+-+-+ -- top edge of first subpixel = 2 / (2 * subpixel_scale)
  // | | | | -- center of first subpixel = 1 / (2 * subpixel_scale)
  // +-+-+-+ -- 0 = 0 / (2 * subpixel_scale)

  // The center of the subpixel part of 'scaled_begin.y()' assuming the
  // 'denominator', i.e., sub_y / denominator is in (0, 1).
  int64 sub_y = (2 * (scaled_begin.y() % subpixel_scale) + 1) * dx;

  // The distance from the from 'scaled_begin' to the right pixel border, to be
  // divided by 2 * 'subpixel_scale'.
  const int first_pixel =
      2 * subpixel_scale - 2 * (scaled_begin.x() % subpixel_scale) - 1;
  // The same from the left pixel border to 'scaled_end'.
  const int last_pixel = 2 * (scaled_end.x() % subpixel_scale) + 1;

  // The full pixel x coordinate of 'scaled_end'.
  const int end_x = std::max(scaled_begin.x(), scaled_end.x()) / subpixel_scale;

  // Move from 'scaled_begin' to the next pixel border to the right.
  sub_y += dy * first_pixel;
  if (dy > 0) {
    while (true) {
      if (!isEqual(pixel_mask.back(), current)) pixel_mask.push_back(current);
      while (sub_y > denominator) {
        sub_y -= denominator;
        ++current.y();
        if (!isEqual(pixel_mask.back(), current)) pixel_mask.push_back(current);
      }
      ++current.x();
      if (sub_y == denominator) {
        sub_y -= denominator;
        ++current.y();
      }
      if (current.x() == end_x) {
        break;
      }
      // Move from one pixel border to the next.
      sub_y += dy * 2 * subpixel_scale;
    }
    // Move from the pixel border on the right to 'scaled_end'.
    sub_y += dy * last_pixel;
    if (!isEqual(pixel_mask.back(), current)) pixel_mask.push_back(current);
    while (sub_y > denominator) {
      sub_y -= denominator;
      ++current.y();
      if (!isEqual(pixel_mask.back(), current)) pixel_mask.push_back(current);
    }
    CHECK_NE(sub_y, denominator);
    CHECK_EQ(current.y(), scaled_end.y() / subpixel_scale);
    return pixel_mask;
  }

  // Same for lines non-ascending in y coordinates.
  while (true) {
    if (!isEqual(pixel_mask.back(), current)) pixel_mask.push_back(current);
    while (sub_y < 0) {
      sub_y += denominator;
      --current.y();
      if (!isEqual(pixel_mask.back(), current)) pixel_mask.push_back(current);
    }
    ++current.x();
    if (sub_y == 0) {
      sub_y += denominator;
      --current.y();
    }
    if (current.x() == end_x) {
      break;
    }
    sub_y += dy * 2 * subpixel_scale;
  }
  sub_y += dy * last_pixel;
  if (!isEqual(pixel_mask.back(), current)) pixel_mask.push_back(current);
  while (sub_y < 0) {
    sub_y += denominator;
    --current.y();
    if (!isEqual(pixel_mask.back(), current)) pixel_mask.push_back(current);
  }
  CHECK_NE(sub_y, 0);
  CHECK_EQ(current.y(), scaled_end.y() / subpixel_scale);
  return pixel_mask;
}

Eigen::Affine3d ToEigen(const ::cartographer::transform::Rigid3d& rigid3) {
  return Eigen::Translation3d(rigid3.translation()) * rigid3.rotation();
}

/**----------------------------------------------------------------------
* Function    : LocalToMap
* Description : 根据local_pose将miss点转换成全局坐标
* Date        : 2022/05/25 wangyukun 
*---------------------------------------------------------------------**/
std::vector<Eigen::Array2i> LocalToMap(const Eigen::Array2i& pixel,
              const Eigen::Array2i& scaled_begin,
              const Eigen::Array2i& scaled_end,
              int subpixel_scale,
              int flag,
              const transform::Rigid3d local_pose_){
			  
	int o_x=1000;
	int o_y=1000;
	int xm;
	int ym;
	int xm_begin;
	int ym_begin;
	int xm_end;
	int ym_end;
  
	const Eigen::Matrix4d homo = ToEigen(local_pose_).matrix();
    if (flag == 0){
        if (scaled_begin.y()/1000 > 150){
			xm=-(pixel).y()+20*homo(0,3)+(o_x+200);
          	ym=(pixel).x()-20*homo(1,3)+(o_y-200);
			
          	xm_begin=-scaled_begin.y()/1000+20*homo(0,3)+(o_x+200);
          	ym_begin=scaled_begin.x()/1000-20*homo(1,3)+(o_y-200);
          	xm_end=-scaled_begin.y()/1000+20*homo(0,3)+(o_x+200);
          	ym_end=scaled_begin.x()/1000-20*homo(1,3)+(o_y-200);
		}else{
          	xm=-(pixel).y()+20*homo(0,3)+(o_x+100);
          	ym=(pixel).x()-20*homo(1,3)+(o_y-100); 
          	xm_begin=-scaled_begin.y()/1000+20*homo(0,3)+(o_x+100);
          	ym_begin=scaled_begin.x()/1000-20*homo(1,3)+(o_y-100);
          	xm_end=-scaled_end.y()/1000+20*homo(0,3)+(o_x+100);
          	ym_end=scaled_end.x()/1000-20*homo(1,3)+(o_y-100);
  		}
	}else{
		if(scaled_end.y()/1000>150){
			xm=-(pixel).y()+20*homo(0,3)+(o_x+200);
			ym=(pixel).x()-20*homo(1,3)+(o_y-200);
			xm_begin=-scaled_end.y()/1000+20*homo(0,3)+(o_x+200);
			ym_begin=scaled_end.x()/1000-20*homo(1,3)+(o_y-200);
			xm_end=-scaled_begin.y()/1000+20*homo(0,3)+(o_x+200);
			ym_end=scaled_begin.x()/1000-20*homo(1,3)+(o_y-200);
		}else{
			xm=-(pixel).y()+20*homo(0,3)+(o_x+100);
			ym=(pixel).x()-20*homo(1,3)+(o_y-100); 
			xm_begin=-scaled_end.y()/1000+20*homo(0,3)+(o_x+100);
			ym_begin=scaled_end.x()/1000-20*homo(1,3)+(o_y-100);
			xm_end=-scaled_begin.y()/1000+20*homo(0,3)+(o_x+100);
			ym_end=scaled_begin.x()/1000-20*homo(1,3)+(o_y-100);
			//  LOG(INFO)<<"wyk--begin is:"<<xm_begin<<" ni      "<<ym_begin;
			//  LOG(INFO)<<"wyk--end   is:"<<xm_end<<" end       "<<ym_end;
			//  LOG(INFO)<<"wyk--xiaowuti is:"<<xm<<"        "<<ym;
		}
      }
	  
	std::vector<Eigen::Array2i> map_position;
	map_position.push_back(Eigen::Array2i(xm_begin, ym_begin));
	map_position.push_back(Eigen::Array2i(xm,ym));
	map_position.push_back(Eigen::Array2i(xm_end,ym_end));
	return map_position;
}

/**
 * @brief 非摄像头视域或摄像头视域超过两米范围的小物体不更新，即该范围内的删去miss点
 * 
 */
std::vector<Eigen::Array2i> ProcessPixelMask(std::vector<Eigen::Array2i> pixel_mask,
                                           const Eigen::Array2i& scaled_begin,
                                           const Eigen::Array2i& scaled_end,
                                           int subpixel_scale,
                                           int flag,
                                           const uint8 intensities,
                                           const transform::Rigid3d local_pose_){
     
      int xm;
      int ym;
      int xm_begin=LocalToMap(scaled_begin,scaled_begin,scaled_end,subpixel_scale,flag,local_pose_)[0].x();
      int ym_begin=LocalToMap(scaled_begin,scaled_begin,scaled_end,subpixel_scale,flag,local_pose_)[0].y();
      int xm_end;
      int ym_end;
 for(std::vector<Eigen::Array2i>::iterator pixel=pixel_mask.begin();pixel!=pixel_mask.end(); ){
      xm=LocalToMap(*pixel,scaled_begin,scaled_end,subpixel_scale,flag,local_pose_)[1].x();
      ym=LocalToMap(*pixel,scaled_begin,scaled_end,subpixel_scale,flag,local_pose_)[1].y();
        if(0<ym*2000+xm && ym*2000+xm<4000000){
          //标志是否删去miss 0代表没有删去
          int flag_miss=0;
          //以小物体为中心2*2搜索
          for (int i=-1;i<=1;i++){
            for (int j=-1;j<=1;j++){
              if (cnew[(ym+j)*2000+xm+i]>100){
              //如果在非摄像头视域
              if (intensities == 0) {     
                pixel=pixel_mask.erase(pixel);
                //probability_grid->ApplyLookupTable(*pixel, hit_table);
                flag_miss=1;
                break;
              //在摄像头视域
              }else{
                //超过两米 小物体不更新
                if ((xm-xm_begin)*(xm-xm_begin)+(ym-ym_begin)*(ym-ym_begin)>1600){
                  pixel=pixel_mask.erase(pixel);
                  flag_miss=1;
                  break;
                }
              } 
            }
          }
           if (flag_miss==1)
              break;
        }
        if (flag_miss == 0)
              ++pixel;  
        } 
        else{
          LOG(INFO)<<"frame erro:"<<xm<<"      "<<ym;
        
  } 
} 
  return pixel_mask;
}




std::vector<Eigen::Array2i> RayToPixelMaskVisualNew(const Eigen::Array2i& scaled_begin,
                                           const Eigen::Array2i& scaled_end,
                                           int subpixel_scale,
                                           const uint8 intensities,
                                           const transform::Rigid3d local_pose_,
                                           int flag,
                                           const std::vector<uint16>& hit_table,
                                         cartographer::mapping::ProbabilityGrid* probability_grid) {
	int o_x=1000;
	int o_y=1000;
	int xm;
	int ym;
	int xm_begin;
	int ym_begin;
	int xm_end;
	int ym_end;
  
  // For simplicity, we order 'scaled_begin' and 'scaled_end' by their x
  // coordinate.
  //测试地图数据是否收到
  // int flag = 0;
  // for (int y = 2000 - 1; y >= 0;--y) {
  //       for (int x = 0; x < 2000; ++x) {
  //            if (cnew[y*2000+x]>100){
  //            std::cout<<x<<"x   "<<y<<"y    "<<std::endl;
  //            flag=1;}
  //       }
  // }
  
  //测试坐标转换是否正确
  // if(intensities==200){
    //const Eigen::Matrix4d homo =ToEigen(local_pose_).matrix();
    // LOG(INFO)<<"wyk--guiji scal_begin:"<<scaled_begin.x()/1000<<"         "<<scaled_begin.y()/1000;
    //  LOG(INFO)<<"wyk--guiji end:"<<scaled_end.x()/1000<<"         "<<scaled_end.y()/1000;
//  if (scaled_begin.y()/1000>150){
//    int o_x=1000;
//    int o_y=1000;
//    int xm=-scaled_begin.y()/1000+20*homo(0,3)+(o_x+200);
//    int ym=scaled_begin.x()/1000-20*homo(1,3)+(o_y-200);

//   LOG(INFO)<<"wyk--guiji jububiao:"<<scaled_begin.x()/1000<<"         "<<scaled_begin.y()/1000;
//   LOG(INFO)<<"wyk--guiji homo:"<<homo(0,3)<<"      "<<homo(1,3);
//   LOG(INFO)<<"wyk--guiji quanjuzuobiao:"<<xm<<"         "<<ym;
//   }else{
//      int o_x=1000;
//      int o_y=1000;
//      int xm=-scaled_begin.y()/1000+20*homo(0,3)+(o_x+100);
//      int ym=scaled_begin.x()/1000-20*homo(1,3)+(o_y-100);
//   LOG(INFO)<<"wyk--guiji jububiao:"<<scaled_begin.x()/1000<<"         "<<scaled_begin.y()/1000;
//   LOG(INFO)<<"wyk--guiji homo:"<<homo(0,3)<<"      "<<homo(1,3);
//   LOG(INFO)<<"wyk--guiji quanjuzuobiao:"<<xm<<"         "<<ym;
//   }
  

  // }
  //小物体测试
  // for (int y = 1200; y >=900 ;--y) {
  //       for (int x = 1000; x <1200 ; ++x) {
  //            if(cnew[y*2000+x]==200){
  //              std::cout<<x<<"c1    ";
  //             std::cout<<y<<"c1    "<<std::endl;
  //            }

  //       }
  // }

    if (scaled_begin.x() > scaled_end.x()) {
    	//flag = 1代表end和begin交换
    	flag = 1;
    	return RayToPixelMaskVisualNew(scaled_end, scaled_begin, subpixel_scale,intensities,local_pose_,flag,hit_table,probability_grid);
 	}
  
    // if(intensities >100){
    //    std::cout<<int(intensities)<<std::endl;
    //   std::cout<<scaled_end.x()<<"end    ";
    //   std::cout<<scaled_end.y()<<"end    "<<std::endl;
    //   std::cout<<scaled_begin.x()<<"begin   ";
    //   std::cout<<scaled_begin.y()<<"begin  "<<std::endl;
    // }
    

  CHECK_GE(scaled_begin.x(), 0);
  CHECK_GE(scaled_begin.y(), 0);
  CHECK_GE(scaled_end.y(), 0);
  std::vector<Eigen::Array2i> pixel_mask;
  // Special case: We have to draw a vertical line in full pixels, as
  // 'scaled_begin' and 'scaled_end' have the same full pixel x coordinate.
  if (scaled_begin.x() / subpixel_scale == scaled_end.x() / subpixel_scale) {
    Eigen::Array2i current(
        scaled_begin.x() / subpixel_scale,
        std::min(scaled_begin.y(), scaled_end.y()) / subpixel_scale);
    pixel_mask.push_back(current);
    const int end_y =
        std::max(scaled_begin.y(), scaled_end.y()) / subpixel_scale;
    for (; current.y() <= end_y; ++current.y()) {
      if (!isEqual(pixel_mask.back(), current)) pixel_mask.push_back(current);
    }

     const Eigen::Matrix4d homo =ToEigen(local_pose_).matrix();
  // std::cout<<homo<<"pose    "<<std::endl;
 

    return ProcessPixelMask(pixel_mask,scaled_begin,scaled_end,subpixel_scale,flag,intensities,local_pose_);
  }

  const int64 dx = scaled_end.x() - scaled_begin.x();
  const int64 dy = scaled_end.y() - scaled_begin.y();
  const int64 denominator = 2 * subpixel_scale * dx;

  // The current full pixel coordinates. We scaled_begin at 'scaled_begin'.
  Eigen::Array2i current = scaled_begin / subpixel_scale;
  pixel_mask.push_back(current);

  // To represent subpixel centers, we use a factor of 2 * 'subpixel_scale' in
  // the denominator.
  // +-+-+-+ -- 1 = (2 * subpixel_scale) / (2 * subpixel_scale)
  // | | | |
  // +-+-+-+
  // | | | |
  // +-+-+-+ -- top edge of first subpixel = 2 / (2 * subpixel_scale)
  // | | | | -- center of first subpixel = 1 / (2 * subpixel_scale)
  // +-+-+-+ -- 0 = 0 / (2 * subpixel_scale)

  // The center of the subpixel part of 'scaled_begin.y()' assuming the
  // 'denominator', i.e., sub_y / denominator is in (0, 1).
  int64 sub_y = (2 * (scaled_begin.y() % subpixel_scale) + 1) * dx;

  // The distance from the from 'scaled_begin' to the right pixel border, to be
  // divided by 2 * 'subpixel_scale'.
  const int first_pixel =
      2 * subpixel_scale - 2 * (scaled_begin.x() % subpixel_scale) - 1;
  // The same from the left pixel border to 'scaled_end'.
  const int last_pixel = 2 * (scaled_end.x() % subpixel_scale) + 1;

  // The full pixel x coordinate of 'scaled_end'.
  const int end_x = std::max(scaled_begin.x(), scaled_end.x()) / subpixel_scale;

  // Move from 'scaled_begin' to the next pixel border to the right.
  sub_y += dy * first_pixel;
  if (dy > 0) {
    while (true) {
      if (!isEqual(pixel_mask.back(), current)) pixel_mask.push_back(current);
      while (sub_y > denominator) {
        sub_y -= denominator;
        ++current.y();
        if (!isEqual(pixel_mask.back(), current)) pixel_mask.push_back(current);
      }
      ++current.x();
      if (sub_y == denominator) {
        sub_y -= denominator;
        ++current.y();
      }
      if (current.x() == end_x) {
        break;
      }
      // Move from one pixel border to the next.
      sub_y += dy * 2 * subpixel_scale;
    }
    // Move from the pixel border on the right to 'scaled_end'.
    sub_y += dy * last_pixel;
    if (!isEqual(pixel_mask.back(), current)) pixel_mask.push_back(current);
    while (sub_y > denominator) {
      sub_y -= denominator;
      ++current.y();
      if (!isEqual(pixel_mask.back(), current)) pixel_mask.push_back(current);
    }
    CHECK_NE(sub_y, denominator);
    CHECK_EQ(current.y(), scaled_end.y() / subpixel_scale);
  
    return ProcessPixelMask(pixel_mask,scaled_begin,scaled_end,subpixel_scale,flag,intensities,local_pose_);
  }

  // Same for lines non-ascending in y coordinates.
  while (true) {
    if (!isEqual(pixel_mask.back(), current)) pixel_mask.push_back(current);
    while (sub_y < 0) {
      sub_y += denominator;
      --current.y();
      if (!isEqual(pixel_mask.back(), current)) pixel_mask.push_back(current);
    }
    ++current.x();
    if (sub_y == 0) {
      sub_y += denominator;
      --current.y();
    }
    if (current.x() == end_x) {
      break;
    }
    sub_y += dy * 2 * subpixel_scale;
  }
  sub_y += dy * last_pixel;
  if (!isEqual(pixel_mask.back(), current)) pixel_mask.push_back(current);
  while (sub_y < 0) {
    sub_y += denominator;
    --current.y();
    if (!isEqual(pixel_mask.back(), current)) pixel_mask.push_back(current);
  }
  CHECK_NE(sub_y, 0);
  CHECK_EQ(current.y(), scaled_end.y() / subpixel_scale);
  
  return ProcessPixelMask(pixel_mask,scaled_begin,scaled_end,subpixel_scale,flag,intensities,local_pose_);

 }



std::vector<Eigen::Array2i> RayToPixelMaskVisual(const Eigen::Array2i& scaled_begin,
                                           const Eigen::Array2i& scaled_end,
                                           int subpixel_scale,
                                           const uint8 intensities) {
  // For simplicity, we order 'scaled_begin' and 'scaled_end' by their x
  // coordinate.
  if (intensities != 0){
    if (scaled_begin.x() > scaled_end.x()) {
    return RayToPixelMask(scaled_end, scaled_begin, subpixel_scale);
  }

  CHECK_GE(scaled_begin.x(), 0);
  CHECK_GE(scaled_begin.y(), 0);
  CHECK_GE(scaled_end.y(), 0);
  std::vector<Eigen::Array2i> pixel_mask;
  // Special case: We have to draw a vertical line in full pixels, as
  // 'scaled_begin' and 'scaled_end' have the same full pixel x coordinate.
  if (scaled_begin.x() / subpixel_scale == scaled_end.x() / subpixel_scale) {
    Eigen::Array2i current(
        scaled_begin.x() / subpixel_scale,
        std::min(scaled_begin.y(), scaled_end.y()) / subpixel_scale);
    pixel_mask.push_back(current);
    const int end_y =
        std::max(scaled_begin.y(), scaled_end.y()) / subpixel_scale;
    for (; current.y() <= end_y; ++current.y()) {
      if (!isEqual(pixel_mask.back(), current)) pixel_mask.push_back(current);
    }
    return pixel_mask;
  }

  const int64 dx = scaled_end.x() - scaled_begin.x();
  const int64 dy = scaled_end.y() - scaled_begin.y();
  const int64 denominator = 2 * subpixel_scale * dx;

  // The current full pixel coordinates. We scaled_begin at 'scaled_begin'.
  Eigen::Array2i current = scaled_begin / subpixel_scale;
  pixel_mask.push_back(current);

  // To represent subpixel centers, we use a factor of 2 * 'subpixel_scale' in
  // the denominator.
  // +-+-+-+ -- 1 = (2 * subpixel_scale) / (2 * subpixel_scale)
  // | | | |
  // +-+-+-+
  // | | | |
  // +-+-+-+ -- top edge of first subpixel = 2 / (2 * subpixel_scale)
  // | | | | -- center of first subpixel = 1 / (2 * subpixel_scale)
  // +-+-+-+ -- 0 = 0 / (2 * subpixel_scale)

  // The center of the subpixel part of 'scaled_begin.y()' assuming the
  // 'denominator', i.e., sub_y / denominator is in (0, 1).
  int64 sub_y = (2 * (scaled_begin.y() % subpixel_scale) + 1) * dx;

  // The distance from the from 'scaled_begin' to the right pixel border, to be
  // divided by 2 * 'subpixel_scale'.
  const int first_pixel =
      2 * subpixel_scale - 2 * (scaled_begin.x() % subpixel_scale) - 1;
  // The same from the left pixel border to 'scaled_end'.
  const int last_pixel = 2 * (scaled_end.x() % subpixel_scale) + 1;

  // The full pixel x coordinate of 'scaled_end'.
  const int end_x = std::max(scaled_begin.x(), scaled_end.x()) / subpixel_scale;

  // Move from 'scaled_begin' to the next pixel border to the right.
  sub_y += dy * first_pixel;
  if (dy > 0) {
    while (true) {
      if (!isEqual(pixel_mask.back(), current)) pixel_mask.push_back(current);
      while (sub_y > denominator) {
        sub_y -= denominator;
        ++current.y();
        if (!isEqual(pixel_mask.back(), current)) pixel_mask.push_back(current);
      }
      ++current.x();
      if (sub_y == denominator) {
        sub_y -= denominator;
        ++current.y();
      }
      if (current.x() == end_x) {
        break;
      }
      // Move from one pixel border to the next.
      sub_y += dy * 2 * subpixel_scale;
    }
    // Move from the pixel border on the right to 'scaled_end'.
    sub_y += dy * last_pixel;
    if (!isEqual(pixel_mask.back(), current)) pixel_mask.push_back(current);
    while (sub_y > denominator) {
      sub_y -= denominator;
      ++current.y();
      if (!isEqual(pixel_mask.back(), current)) pixel_mask.push_back(current);
    }
    CHECK_NE(sub_y, denominator);
    CHECK_EQ(current.y(), scaled_end.y() / subpixel_scale);
    return pixel_mask;
  }

  // Same for lines non-ascending in y coordinates.
  while (true) {
    if (!isEqual(pixel_mask.back(), current)) pixel_mask.push_back(current);
    while (sub_y < 0) {
      sub_y += denominator;
      --current.y();
      if (!isEqual(pixel_mask.back(), current)) pixel_mask.push_back(current);
    }
    ++current.x();
    if (sub_y == 0) {
      sub_y += denominator;
      --current.y();
    }
    if (current.x() == end_x) {
      break;
    }
    sub_y += dy * 2 * subpixel_scale;
  }
  sub_y += dy * last_pixel;
  if (!isEqual(pixel_mask.back(), current)) pixel_mask.push_back(current);
  while (sub_y < 0) {
    sub_y += denominator;
    --current.y();
    if (!isEqual(pixel_mask.back(), current)) pixel_mask.push_back(current);
  }
  CHECK_NE(sub_y, 0);
  CHECK_EQ(current.y(), scaled_end.y() / subpixel_scale);
  return pixel_mask;
  }else {
    return  std::vector<Eigen::Array2i> ();
  }
  
}



}  // namespace mapping
}  // namespace cartographer
