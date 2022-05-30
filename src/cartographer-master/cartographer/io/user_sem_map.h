#ifndef _USER_SEM_MAP_H_
#define _USER_SEM_MAP_H_

#include <vector>
#include <map>
#include "Eigen/Core"

using namespace std;
void test_function(void);

class UserSemMap{
public:
	int* fusion_sem_map;
	int* origin_sem_map;
	std::vector<std::map<unsigned char,int>> type_score_map; 
	std::map<int,int> processed_submap_id;
	
	// construct user sem map, must set width 
	UserSemMap(size_t width, size_t height, float resolution);

	int SetOriginMapData(size_t x, size_t y, int valule);

	int GetOriginMapData(size_t x, size_t y, int valule);

	void UpdateFusionSemMap(const uint32_t* pixel_data, Eigen::Array2f origin,
		int height_map, int width_map);

	void UpdateOriginSemMap(const uint32_t* pixel_data, 
		const Eigen::Matrix4d homo,int submap_slice_height,
		int submap_slice_width,int submap_index, 
		int submap_version);

	void AddProcessedSubMapIdToMap(int submap_index, int submap_version);

	void UpdateSemMapFilter(std::map <unsigned char,int>& sem_list,uint8_t sem);

	void SetOriginSemMapBestMatch(std::map <unsigned char,int>& sem_list,int pos_x, int pos_y);

	int GetSemMapCenterPosX(void);

	int GetSemMapCenterPosY(void);

	bool CheckSemMapLocIsVaild(int pos_x, int pos_y);

	bool isSumMapIdIsProcessed(int submap_index, int submap_version);

	int UpdateSemMapData();

	~UserSemMap();
private:
	size_t map_size_w;
	size_t map_size_h;
	float  map_resoution;
	size_t map_scale;
	int CalSumMapId(int submap_index, int submap_version);
};

#endif