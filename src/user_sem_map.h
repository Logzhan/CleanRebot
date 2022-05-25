#ifndef _USER_SEM_MAP_H_
#define _USER_SEM_MAP_H_

#include <vector>
#include <map>
#include "Eigen/Core"

using namespace std;
void test_function(void);

class UserSemMap{
public:
	struct semMapInfo{
		int origin_map;
		int carto_map;
		int fusion_map;
		map <unsigned char,int> type_score_map;
	};

	vector <semMapInfo> sem_map; 
	// construct user sem map, must set width 
	UserSemMap(size_t width, size_t height, float resolution);
	// return the sem map pos data resiult
	//int GetSemMapDataWithPos(int pos_x, int pos_y);

	int SetOriginMapData(size_t x, size_t y, int valule);
	int GetOriginMapData(size_t x, size_t y, int valule);

	void UpdateFusionSemMap(const uint32_t* pixel_data, Eigen::Array2f origin,
		int height_map, int width_map, int* fusion_sem_map, int* origin_sem_map,
		int sem_map_width);

	int UpdateSemMapData();

	~UserSemMap();
private:
	size_t map_size_w;
	size_t map_size_h;
	float  map_resoution;
};

#endif