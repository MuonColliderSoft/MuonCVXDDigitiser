#ifndef FindUnionAlgorithm_h
#define FindUnionAlgorithm_h 1

#include "AbstractSensor.h"
#include <unordered_map>

using std::vector;
using std::unordered_map;

using ClusterOfPixel = vector<LinearPosition>;

using ClusterOfCoordinate = vector<GridCoordinate>;

class FindUnionAlgorithm
{
public:
    FindUnionAlgorithm(int n_row, int n_col);
    virtual ~FindUnionAlgorithm() {}
    int find(int x, int y);
    void merge(int x1, int y1, int x2, int y2);
    void init();
    void close();
    void invalidate(int x, int y);
    vector<ClusterOfPixel> get_clusters();
    vector<ClusterOfCoordinate> list_clusters();

private:
    int rows;
    int columns;
    unordered_map<int, ClusterOfPixel> c_table;

    GridPosition locate;
    vector<int> data;
};


#endif //FindUnionAlgorithm_h