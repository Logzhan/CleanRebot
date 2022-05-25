#include <iostream>
#include "cartographer/io/user_sem_map.h"


void test_function(){
	std::cout << "test user_seem_map" << std::endl;
}

// UserSemMap construct function
UserSemMap::UserSemMap(size_t width, size_t height, float resolution){

	std::cout << "SemMap::" << "construct user sem map" << std::endl;
	// resolution can not to small, else using the default value
	if(resolution < 0.01){
		resolution = 0.5;
	}
	map_size_h = height / resolution;
	map_size_w = width  / resolution;
	map_resoution = resolution;

	cout << "SemMap:: map size width  = " << map_size_w << endl;
	cout << "SemMap:: map size height = " << map_size_h << endl;
	cout << "SemMap:: map size resoution = " << map_resoution << endl;
}

/**----------------------------------------------------------------------
* Function    : UpdateFusionSemMap
* Description : 更新融合语义地图，结合原始语义地图和cartographer原始地图，生成
                最终呈现的融合语义地图
* Date        : 2022/05/25 wangyukun & zhanli
*---------------------------------------------------------------------**/
void UserSemMap::UpdateFusionSemMap(const uint32_t* pixel_data, Eigen::Array2f origin,
	int height_map, int width_map, int* fusion_sem_map, int* origin_sem_map,
	int sem_map_width){

	for (int y = height_map - 1; y >= 0; --y) {
        for (int x = 0; x < width_map; ++x) {
            const uint32_t packed = pixel_data[y * width_map + x];
            const unsigned char probality = packed >> 16;
            const unsigned char observed = packed >> 8;
            const unsigned char value=std::lround((1. - probality / 255.) * 100.);

            int xm = x + (1000-origin.x());
            int ym = y + (1000-origin.y());
            
			if(xm > 2000 || ym > 2000 || xm < 0 || ym < 0)continue;
			
			fusion_sem_map[ym * 2000 + xm] = -1;    

            if (observed != 0){	
			  	fusion_sem_map[ym * 2000 + xm] = value;
              	if(value > 65 && value <= 100 && origin_sem_map[ym * 2000 + xm] > 100){
                	fusion_sem_map[ym * 2000 + xm] = origin_sem_map[ym * 2000 + xm];
              	}
            } 
       }
   	}
}


UserSemMap::~UserSemMap(){

}