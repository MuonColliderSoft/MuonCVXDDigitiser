#include "FindUnionAlgorithm.h"

FindUnionAlgorithm::FindUnionAlgorithm(int n_row, int n_col) :
    rows(n_row),
    columns(n_col),
    c_table(),
    locate(n_row, n_col),
    data(rows * columns, 0)
{}

void FindUnionAlgorithm::init()
{
    for (size_t k = 0; k < data.size(); k++) data[k] = k;
    c_table.clear();
}

void FindUnionAlgorithm::close()
{
    /*
     * Remove internal references in the grid, each cell contains the cluster ID
     */
    for (int h = 0; h < rows; h++)
    {
        for (int k = 0; k < columns; k++)
        {
            if (data[locate(h, k)] >= 0) find(h, k);
        }
    }

    /*
     * Prepare the set of clusters
     */
    for (size_t k = 0; k < data.size(); k++)
    {
        auto c_item = c_table.find(data[k]);
        if (c_item != c_table.end())
        {
            (c_item->second).push_back(k);
        }
        else if (data[k] >= 0)
        {
            ClusterOfPixel p_list { k };
            c_table.emplace(data[k], p_list);
        }
    }
}

int FindUnionAlgorithm::find(int x, int y)
{
    int res = locate(x, y);
    while (data[res] != res)
    {
        if (res < 0) return -1;
        res = data[res];
    }

    int pos = locate(x, y);
    while (data[pos] != pos)
    {
        int curr = data[pos];
        data[pos] = res;
        pos = curr;
    }
    return res;
}

void FindUnionAlgorithm::merge(int x1, int y1, int x2, int y2)
{
    int pset1 = find(x1, y1);
    int pset2 = find(x2, y2);

    if (pset1 < 0 || pset2 < 0) return;

    if (pset1 >= pset2)
    {
        data[pset1] = pset2;
    }
    else
    {
        data[pset2] = pset1;
    }
}

void FindUnionAlgorithm::invalidate(int x, int y)
{
    data[locate(x, y)] = -1;
}

vector<ClusterOfPixel> FindUnionAlgorithm::get_clusters()
{
    vector<ClusterOfPixel> result;
    for (auto c_item : c_table)
    {
        result.push_back(std::move(c_item.second));
    }
    return result;
}

vector<ClusterOfCoordinate> FindUnionAlgorithm::list_clusters()
{
    vector<ClusterOfCoordinate> result;
    for (auto c_item : c_table)
    {
        ClusterOfCoordinate tmp_item;
        for (LinearPosition curr_pos : c_item.second)
        {
            tmp_item.push_back(locate(curr_pos));
        }
        result.push_back(std::move(tmp_item));
    }
    return result;
}