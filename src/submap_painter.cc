/*
 * Copyright 2017 The Cartographer Authors
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

#include "cartographer/io/submap_painter.h"
#include "cartographer/io/user_sem_map.h"
#include <iostream>
#include "cartographer/mapping/2d/submap_2d.h"
#include "cartographer/mapping/3d/submap_3d.h"
#include "cairo/cairo.h"
#include <map>
#include <utility>

int fusion_sem_map[4000000]={0};
int origin_sem_map[4000000]={0};

//如何初始化（0，0）？
// semmap : record different obj map
std::vector < std::map<unsigned char,int> > semmap(4000000); 
//存放计算过的子图
std::map<int,int> semmap_id;
UserSemMap user_sem_map(100,100,0.5);

namespace cartographer {
namespace io {
namespace {

Eigen::Affine3d ToEigen(const ::cartographer::transform::Rigid3d& rigid3) {
  return Eigen::Translation3d(rigid3.translation()) * rigid3.rotation();
}

void CairoPaintSubmapSlices(
    const double scale,
    const std::map<::cartographer::mapping::SubmapId, SubmapSlice>& submaps,
    cairo_t* cr, std::function<void(const ::cartographer::mapping::SubmapId&,const SubmapSlice&)> draw_callback) {
  cairo_scale(cr, scale, scale);

  for (auto& pair : submaps) {
    const auto& submap_slice = pair.second;
    const auto& submap_id = pair.first;
    if (submap_slice.surface == nullptr) {
      return;
    }
    const Eigen::Matrix4d homo =
        ToEigen(submap_slice.pose * submap_slice.slice_pose).matrix();

    cairo_save(cr);
    cairo_matrix_t matrix;
    cairo_matrix_init(&matrix, homo(1, 0), homo(0, 0), -homo(1, 1), -homo(0, 1),
                      homo(0, 3), -homo(1, 3));
    cairo_transform(cr, &matrix);

    const double submap_resolution = submap_slice.resolution;
    cairo_scale(cr, submap_resolution, submap_resolution);

    // Invokes caller's callback to utilize slice data in global cooridnate
    // frame. e.g. finds bounding box, paints slices.
    draw_callback(submap_id,submap_slice);
    cairo_restore(cr);
  }
}

bool Has2DGrid(const mapping::proto::Submap& submap) {
  return submap.has_submap_2d() && submap.submap_2d().has_grid();
}

bool Has3DGrids(const mapping::proto::Submap& submap) {
  return submap.has_submap_3d() &&
         submap.submap_3d().has_low_resolution_hybrid_grid() &&
         submap.submap_3d().has_high_resolution_hybrid_grid();
}

}  // namespace

//按value值比较
bool cmp_value(const std::pair<int, int> left,const std::pair<int,int> right){
	return left.second < right.second;
}

int* Slices(
    const std::map<::cartographer::mapping::SubmapId, SubmapSlice>& submaps,
    const double resolution) {
  Eigen::AlignedBox2f bounding_box;
  {
    auto surface = MakeUniqueCairoSurfacePtr(
        cairo_image_surface_create(kCairoFormat, 1, 1));
    auto cr = MakeUniqueCairoPtr(cairo_create(surface.get()));
    const auto update_bounding_box = [&bounding_box, &cr](double x, double y) {
      cairo_user_to_device(cr.get(), &x, &y);
      bounding_box.extend(Eigen::Vector2f(x, y));
    };
    
    
    CairoPaintSubmapSlices(
        1. / resolution, submaps, cr.get(),
        [&update_bounding_box](const ::cartographer::mapping::SubmapId& submap_id,const SubmapSlice& submap_slice) {
          update_bounding_box(0, 0);
          update_bounding_box(submap_slice.width, 0);
          update_bounding_box(0, submap_slice.height);
          update_bounding_box(submap_slice.width, submap_slice.height);
        });
  }

  const int kPaddingPixel = 5;
  const Eigen::Array2i size(
      std::ceil(bounding_box.sizes().x()) + 2 * kPaddingPixel,
      std::ceil(bounding_box.sizes().y()) + 2 * kPaddingPixel);
  const Eigen::Array2f origin(-bounding_box.min().x() + kPaddingPixel,
                              -bounding_box.min().y() + kPaddingPixel);

  auto surface = MakeUniqueCairoSurfacePtr(
      cairo_image_surface_create(kCairoFormat, size.x(), size.y()));
    int widthmap=cairo_image_surface_get_width(surface.get());
    int heightmap=cairo_image_surface_get_height(surface.get());
   //int a[widthmap*heightmap]={0};
  {
    auto cr = MakeUniqueCairoPtr(cairo_create(surface.get()));
    cairo_set_source_rgba(cr.get(), 0.5,0.0, 0., 1.0);
    cairo_paint(cr.get());
   
    //...cairo_translate(cr.get(), origin.x(), origin.y());
    cairo_translate(cr.get(), 1000,1000);
   // std::cout<< origin.x()<<"x    "<<origin.y()<<"y   "<<std::endl;
    int o_x=1000; 
    int o_y=1000; 
    CairoPaintSubmapSlices(1. / resolution, submaps, cr.get(),
		[&cr,&surface,&widthmap,&heightmap,&o_x,&o_y](const ::cartographer::mapping::SubmapId& submap_id,const SubmapSlice& submap_slice) {

		int height = submap_slice.height;  //确保width是像素点尺寸长
		int width = submap_slice.width;

		const uint32_t* pixel_data = reinterpret_cast<uint32_t*>(cairo_image_surface_get_data(submap_slice.surface.get()));
		const Eigen::Matrix4d homo = ToEigen(submap_slice.pose * submap_slice.slice_pose).matrix();


		int id = submap_id.submap_index*100 + submap_slice.version;

		auto iter_id=semmap_id.find(id);
		//如果semmap_id中找不到该子图的id，则将id加到该semmap中，并进行语义的计算
		if (iter_id == semmap_id.end()){
			semmap_id.insert(std::pair<int,int>(id,1));

			for (int y = height - 1; y >= 0; --y) {
				for (int x = 0; x < width; ++x) {
					const uint32_t packed = pixel_data[y * width + x];
					const unsigned char sem =(uint8_t)(packed & 0xFFu);

					int xm=homo(1,0)*x - homo(1,1)*y + 20*homo(0,3) + o_x - 1;
					int ym=homo(0,0)*x - homo(0,1)*y - 20*homo(1,3) + o_y + 1;

					std::map <unsigned char,int>::iterator iter;

					//记得给semmap初始化
					//vector.at()会抛出越界访问,vector[]不会，注意vector的越界调用，编译可能不会报错，读取到的数据是随机的
					iter=semmap.at(ym * 2000 + xm).find(sem);

					if(iter != semmap.at(ym * 2000 + xm).end()){
						if(iter->second <= 5){
							semmap.at(ym * 2000 + xm)[sem]=semmap.at(ym * 2000 + xm)[sem]+1;
							//在摄像头视域时，for循环给其他语义减一
							if (sem!=0){
								iter=semmap.at(ym * 2000 + xm).begin();
								for(;iter !=  semmap.at(ym * 2000 + xm).end();iter++){
									if(iter->first!=sem && iter->second>0 ){
										iter->second=iter->second-1;
									}
								}
							}                             
						}  
					}else{
						semmap.at(ym * 2000 + xm).insert(std::pair<unsigned char,int>(sem,1));
					}
					//策略2 保存0的键值对
					int value_0;
					auto iter0=semmap.at(ym * 2000 + xm).find(0);

					//如果存在0的键值对  因为不知道怎么给全局变量赋值0：0
					if (iter0 != semmap.at(ym*2000+xm).end() ){
					//如果0的value值大于等于1
						value_0=iter0->second;
						//删除0的键值对
						semmap.at(ym*2000+xm).erase(iter0);
						//返回除0外value值最大的迭代器
						auto max_probability= max_element(semmap.at(ym * 2000 + xm).begin(),semmap.at(ym * 2000 + xm).end(),cmp_value);
						if ((max_probability->second)>=3 && max_probability->first!=1 && xm<2000 && ym<2000 && xm>=0 && ym>=0){
							//std::cout<<max_probability->second<<"sec    ";
							origin_sem_map[ym * 2000 + xm]=max_probability->first;
						}else{
							origin_sem_map[ym * 2000 + xm]=0;                                                     
						}
						//恢复0键值对
						semmap.at(ym * 2000 + xm).insert(std::pair<unsigned char,int>(0,value_0));                                 
					}
				}
			}
		}
	});
   	cairo_surface_flush(surface.get());
  	}
	return origin_sem_map;
}

int* GetFusionSemMap(){
	return fusion_sem_map;
}

PaintSubmapSlicesResult PaintSubmapSlices(
    const std::map<::cartographer::mapping::SubmapId, SubmapSlice>& submaps,
    const double resolution) {
  Eigen::AlignedBox2f bounding_box;
  {
    auto surface = MakeUniqueCairoSurfacePtr(
        cairo_image_surface_create(kCairoFormat, 1, 1));
    auto cr = MakeUniqueCairoPtr(cairo_create(surface.get()));
    const auto update_bounding_box = [&bounding_box, &cr](double x, double y) {
      cairo_user_to_device(cr.get(), &x, &y);
      bounding_box.extend(Eigen::Vector2f(x, y));
    };

    CairoPaintSubmapSlices(
        1. / resolution, submaps, cr.get(),
        [&update_bounding_box](const ::cartographer::mapping::SubmapId& submap_id,const SubmapSlice& submap_slice) {
          update_bounding_box(0, 0);
          update_bounding_box(submap_slice.width, 0);
          update_bounding_box(0, submap_slice.height);
          update_bounding_box(submap_slice.width, submap_slice.height);
        });
   }

	const int kPaddingPixel = 5;
	const Eigen::Array2i size(
		std::ceil(bounding_box.sizes().x()) + 2 * kPaddingPixel,
		std::ceil(bounding_box.sizes().y()) + 2 * kPaddingPixel);
	const Eigen::Array2f origin(-bounding_box.min().x() + kPaddingPixel,
								-bounding_box.min().y() + kPaddingPixel);

  	auto surface = MakeUniqueCairoSurfacePtr(
      	cairo_image_surface_create(kCairoFormat, size.x(), size.y()));
  	{
		auto cr = MakeUniqueCairoPtr(cairo_create(surface.get()));
		cairo_set_source_rgba(cr.get(),  0.5, 0.0, 0.0, 1.0);
		cairo_paint(cr.get());
		cairo_translate(cr.get(), origin.x(), origin.y());
		CairoPaintSubmapSlices(1. / resolution, submaps, cr.get(),
							[&cr,&surface](const ::cartographer::mapping::SubmapId& submap_id,const SubmapSlice& submap_slice) {
								cairo_set_source_surface(
									cr.get(), submap_slice.surface.get(), 0., 0.);
								cairo_paint(cr.get());
							});
		// zhanli:@newAdd: update user fusion sem map
		user_sem_map.UpdateFusionSemMap(reinterpret_cast<uint32_t*>(cairo_image_surface_get_data(surface.get())), 
			origin,
			cairo_image_surface_get_height(surface.get()),
			cairo_image_surface_get_width(surface.get()),
			fusion_sem_map,origin_sem_map,0);

		cairo_surface_flush(surface.get());
  	}
  	return PaintSubmapSlicesResult(std::move(surface), origin);
}

void FillSubmapSlice(
    const ::cartographer::transform::Rigid3d& global_submap_pose,
    const ::cartographer::mapping::proto::Submap& proto,
    SubmapSlice* const submap_slice,
    mapping::ValueConversionTables* conversion_tables) {
  ::cartographer::mapping::proto::SubmapQuery::Response response;
  ::cartographer::transform::Rigid3d local_pose;
  if (proto.has_submap_3d()) {
    mapping::Submap3D submap(proto.submap_3d());
    local_pose = submap.local_pose();
    submap.ToResponseProto(global_submap_pose, &response);
  } else {
    ::cartographer::mapping::Submap2D submap(proto.submap_2d(),
                                             conversion_tables);
    local_pose = submap.local_pose();
    submap.ToResponseProto(global_submap_pose, &response);
  }
  submap_slice->pose = global_submap_pose;

  auto& texture_proto = response.textures(0);
  const SubmapTexture::Pixels pixels = UnpackTextureData(
      texture_proto.cells(), texture_proto.width(), texture_proto.height());
  submap_slice->width = texture_proto.width();
  submap_slice->height = texture_proto.height();
  submap_slice->resolution = texture_proto.resolution();
  submap_slice->slice_pose =
      ::cartographer::transform::ToRigid3(texture_proto.slice_pose());
  submap_slice->surface =
      DrawTexture(pixels.intensity, pixels.alpha, pixels.color,texture_proto.width(),
                  texture_proto.height(), &submap_slice->cairo_data);
}

void DeserializeAndFillSubmapSlices(
    ProtoStreamDeserializer* deserializer,
    std::map<mapping::SubmapId, SubmapSlice>* submap_slices,
    mapping::ValueConversionTables* conversion_tables) {
  std::map<mapping::SubmapId, transform::Rigid3d> submap_poses;
  for (const auto& trajectory : deserializer->pose_graph().trajectory()) {
    for (const auto& submap : trajectory.submap()) {
      submap_poses[mapping::SubmapId(trajectory.trajectory_id(),
                                     submap.submap_index())] =
          transform::ToRigid3(submap.pose());
    }
  }
  mapping::proto::SerializedData proto;
  while (deserializer->ReadNextSerializedData(&proto)) {
    if (proto.has_submap() &&
        (Has2DGrid(proto.submap()) || Has3DGrids(proto.submap()))) {
      const auto& submap = proto.submap();
      const mapping::SubmapId id{submap.submap_id().trajectory_id(),
                                 submap.submap_id().submap_index()};
      FillSubmapSlice(submap_poses.at(id), submap, &(*submap_slices)[id],
                      conversion_tables);
    }
  }
}

SubmapTexture::Pixels UnpackTextureData(const std::string& compressed_cells,
                                        const int width, const int height) {
	SubmapTexture::Pixels pixels;
	std::string cells;
	// 将压缩后的地图栅格数据 解压成 字符串
	::cartographer::common::FastGunzipString(compressed_cells, &cells);

	const int num_pixels = width * height;
	CHECK_EQ(cells.size(), 2 * num_pixels);
	pixels.intensity.reserve(num_pixels);
	pixels.alpha.reserve(num_pixels);
	// wangyukun @newAdd : 增加color
	pixels.color.reserve(num_pixels);

	// 填充数据
	for (int i = 0; i < height; ++i) {
		for (int j = 0; j < width; ++j) {
			const int8  delta = cells[(i * width + j) * 2];
			const uint8 value = delta > 0 ? delta : 0;
			const uint8 alpha = delta > 0 ? 0 : -delta;
			pixels.intensity.push_back(value);
			pixels.alpha.push_back(alpha);
			pixels.color.push_back(cells[(i * width + j) * 2 + 1]);
		}
	}
  	return pixels;
}

/**
 * @brief 指向新创建的图像的指针
 * 
 * @param[in] intensity 地图栅格数据
 * @param[in] alpha 地图栅格的透明度
 * @param[in] width 地图的宽
 * @param[in] height 地图的高
 * @param[out] cairo_data 4字节的值, 左边3个字节分别存储了alpha_value intensity_value 与 observed
 * @return UniqueCairoSurfacePtr 指向新创建的图像的指针
 */
//....
UniqueCairoSurfacePtr DrawTexture(const std::vector<char>& intensity,
                                  const std::vector<char>& alpha,
                                  //。。。增加color传参
                                   const std::vector<char>& color,
                                  const int width, const int height,
                                  std::vector<uint32_t>* const cairo_data) {
  CHECK(cairo_data->empty());

  // Properly dealing with a non-common stride would make this code much more
  // complicated. Let's check that it is not needed.
  const int expected_stride = 4 * width;
  CHECK_EQ(expected_stride, cairo_format_stride_for_width(kCairoFormat, width));

  // 对cairo_data进行填充
  for (size_t i = 0; i < intensity.size(); ++i) {
    // We use the red channel to track intensity information. The green
    // channel we use to track if a cell was ever observed.
    // 使用红色通道来跟踪强度信息 绿色通道来追踪栅格是否被观察到
    //。。。填补最后一个字节 color
    const uint8_t color_value = color.at(i);
    const uint8_t intensity_value = intensity.at(i);
    const uint8_t alpha_value = alpha.at(i);
    const uint8_t observed = (intensity_value == 0 && alpha_value == 0) ? 0 : 255;
    // tag: 这里需要确认一下
    cairo_data->push_back((alpha_value << 24) |     // 第一字节 存储透明度
                          (intensity_value << 16) | // 第二字节 存储栅格值
                          (observed << 8) |         // 第三字节 存储是否被更新过
                           color_value);            // 第四字节 始终为0
   
  }

  // c++11: reinterpret_cast 用于进行各种不同类型的指针之间、不同类型的引用之间以及指针和能容纳指针的整数类型之间的转换

  // MakeUniqueCairoSurfacePtr 生成一个指向cairo_surface_t数据的指针
  auto surface = MakeUniqueCairoSurfacePtr(
    // cairo_image_surface_create_for_data: 根据提供的像素数据创建surface, 返回指向新创建的surface的指针
    cairo_image_surface_create_for_data(
      reinterpret_cast<unsigned char*>(cairo_data->data()), kCairoFormat, width,
      height, expected_stride) );
        
  CHECK_EQ(cairo_surface_status(surface.get()), CAIRO_STATUS_SUCCESS)
      << cairo_status_to_string(cairo_surface_status(surface.get()));
  return surface;
}

UniqueCairoSurfacePtr DrawTexture(const std::vector<char>& intensity,
                                  const std::vector<char>& alpha,
                                  const int width, const int height,
                                  std::vector<uint32_t>* const cairo_data) {
  CHECK(cairo_data->empty());

  // Properly dealing with a non-common stride would make this code much more
  // complicated. Let's check that it is not needed.
  const int expected_stride = 4 * width;
  CHECK_EQ(expected_stride, cairo_format_stride_for_width(kCairoFormat, width));
  for (size_t i = 0; i < intensity.size(); ++i) {
    // We use the red channel to track intensity information. The green
    // channel we use to track if a cell was ever observed.
    const uint8_t intensity_value = intensity.at(i);
    const uint8_t alpha_value = alpha.at(i);
    const uint8_t observed =
        (intensity_value == 0 && alpha_value == 0) ? 0 : 255;
    cairo_data->push_back((alpha_value << 24) | (intensity_value << 16) |
                          (observed << 8) | 0);
    
  }

  auto surface = MakeUniqueCairoSurfacePtr(cairo_image_surface_create_for_data(
      reinterpret_cast<unsigned char*>(cairo_data->data()), kCairoFormat, width,
      height, expected_stride));
  CHECK_EQ(cairo_surface_status(surface.get()), CAIRO_STATUS_SUCCESS)
      << cairo_status_to_string(cairo_surface_status(surface.get()));
  return surface;
}

}  // namespace io
}  // namespace cartographer


