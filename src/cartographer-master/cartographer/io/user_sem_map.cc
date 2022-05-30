#include <iostream>
#include "cartographer/io/user_sem_map.h"

bool CompareValue(const std::pair<int, int> left,const std::pair<int,int> right){
	return left.second < right.second;
}

// UserSemMap construct function
UserSemMap::UserSemMap(size_t width, size_t height, float resolution){

	std::cout << "SemMap::" << "construct user sem map" << std::endl;
	// resolution can not to small, else using the default value
	if(resolution < 0.001){
		resolution = 0.05;
	}
	map_size_h = height / resolution;
	map_size_w = width  / resolution;
	map_resoution = resolution;
	map_scale = std::lround((1. /  resolution));

	cout << "SemMap:: map size width  = " << map_size_w << endl;
	cout << "SemMap:: map size height = " << map_size_h << endl;
	cout << "SemMap:: map size resoution = " << map_resoution << endl;
	cout << "SemMap:: map scale = " << map_scale << endl;
	
	// create the inital size
	fusion_sem_map = new int[map_size_h * map_size_w];
	origin_sem_map = new int[map_size_h * map_size_w];
	type_score_map.resize(map_size_h * map_size_w);
}

/**----------------------------------------------------------------------
* Function    : UpdateFusionSemMap
* Description : 更新融合语义地图，结合原始语义地图和cartographer原始地图，生成
                最终呈现的融合语义地图
* Date        : 2022/05/25 wangyukun & zhanli
*---------------------------------------------------------------------**/
void UserSemMap::UpdateFusionSemMap(const uint32_t* pixel_data, Eigen::Array2f origin,
	int height_map, int width_map){

	for (int y = height_map - 1; y >= 0; --y) {
        for (int x = 0; x < width_map; ++x) {
            const uint32_t packed = pixel_data[y * width_map + x];
            const unsigned char probality = packed >> 16;
            const unsigned char observed = packed >> 8;
            const unsigned char value=std::lround((1. - probality / 255.) * 100.);

            int xm = x + (GetSemMapCenterPosX() - origin.x());
            int ym = y + (GetSemMapCenterPosY() - origin.y());
            
			if(!CheckSemMapLocIsVaild(xm, ym))continue;
			
			fusion_sem_map[ym * map_size_w + xm] = -1;    

            if (observed != 0){	
			  	fusion_sem_map[ym * map_size_w + xm] = value;
              	if(value > 65 && value <= 100 && origin_sem_map[ym * map_size_w + xm] > 100){
                	fusion_sem_map[ym * map_size_w + xm] = origin_sem_map[ym * map_size_w + xm];
              	}
            } 
       }
   	}
}

void UserSemMap::UpdateOriginSemMap(const uint32_t* pixel_data, 
				const Eigen::Matrix4d homo,int submap_slice_height,
				int submap_slice_width,int submap_index, 
				int submap_version)
{
	// check this submap weather is had been processed 
	if(isSumMapIdIsProcessed(submap_index, submap_version)){
    	return;
    }
	for (int y = submap_slice_height - 1; y >= 0; --y) {
    	for (int x = 0; x < submap_slice_width; ++x) {
			// get sem data from data
			uint8_t sem = (uint8_t)(pixel_data[y * submap_slice_width + x] & 0xFFu);
			// if sem == 0 don not proc anything
			if(sem == 0){
				continue;
			}

			int o_x = GetSemMapCenterPosX();
			int o_y = GetSemMapCenterPosY();

			// cal the location
			int xm = homo(1,0)*x - homo(1,1)*y + map_scale * homo(0,3) + o_x - 1;
			int ym = homo(0,0)*x - homo(0,1)*y - map_scale * homo(1,3) + o_y + 1;

			// check location is vaild
			if(!CheckSemMapLocIsVaild(xm, ym)){
				continue;
			}
			// Update sem map filter
			UpdateSemMapFilter(type_score_map.at(ym * map_size_w + xm), sem);
			// Select the best sem to origin sem map
			SetOriginSemMapBestMatch(type_score_map.at(ym * map_size_w + xm), xm, ym);
		}
    }
	AddProcessedSubMapIdToMap(submap_index, submap_version);
}

void UserSemMap::SetOriginSemMapBestMatch(std::map <unsigned char,int>& sem_list,
										int pos_x, int pos_y)
{
	auto max_probability = max_element(sem_list.begin(), sem_list.end(), CompareValue);

	if ((max_probability->second) >= 3 && max_probability->first != 1){
		origin_sem_map[pos_y * map_size_w + pos_x] = max_probability->first;
	}else{
		origin_sem_map[pos_y * map_size_w + pos_x] = 0;                                                     
	}
}

void UserSemMap::UpdateSemMapFilter(std::map <unsigned char,int>& sem_list,uint8_t sem){
	auto iter = sem_list.find(sem);

	if(iter == sem_list.end()){
		sem_list.insert(std::pair<unsigned char,int>(sem,1));
	}

	if(iter->second <= 5){
		sem_list[sem]++;
		for(iter = sem_list.begin();iter != sem_list.end();iter++){
			if(iter->first != sem && iter->second > 0){
				iter->second--;
			}
		}
	}
}

int UserSemMap::CalSumMapId(int submap_index, int submap_version){
	return submap_index * 100 + submap_version;
}

void UserSemMap::AddProcessedSubMapIdToMap(int submap_index, int submap_version){
	int id = CalSumMapId(submap_index, submap_version);
	processed_submap_id.insert(std::pair<int,int>(id, 1));
}

bool UserSemMap::isSumMapIdIsProcessed(int submap_index, int submap_version){
	int id = CalSumMapId(submap_index, submap_version);
	if(processed_submap_id.find(id) != processed_submap_id.end()){
      return true;
    }
	return false;
}

int UserSemMap::GetSemMapCenterPosX(void){
	return (int)(map_size_w / 2);
}

int UserSemMap::GetSemMapCenterPosY(void){
	return (int)(map_size_h / 2);
}

bool UserSemMap::CheckSemMapLocIsVaild(int pos_x, int pos_y){
	if(pos_x >= map_size_w || pos_y >= map_size_h || pos_x < 0 || pos_y < 0){
		return false;
	}
	return true;
}

UserSemMap::~UserSemMap(){
	delete[] fusion_sem_map;
	delete[] origin_sem_map;
}