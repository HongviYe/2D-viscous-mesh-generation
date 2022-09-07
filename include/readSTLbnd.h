#ifndef _READ_STL_BND_H_
#define _READ_STL_BND_H_

#include <igl/list_to_matrix.h>
#include <igl/readSTL.h>
#include <igl/readWRL.h>
#include <igl/boundary_loop.h>


#include <Eigen/Core>
#include <string>
#include <vector>
#include <array>
#include <algorithm>
#include <fstream>
#include <cstdio>

namespace RIGIDT {
    template <typename DerivedV, typename DerivedF>
    bool readSTLbnd(
        const std::string str,
        Eigen::PlainObjectBase<DerivedV>& V,
        Eigen::PlainObjectBase<DerivedF>& F)
    {



        std::vector<std::vector<int> >index;
        std::vector<std::vector<double > >points;

        igl::readWRL(str, points, index);


        std::vector<std::vector<int>> all_bnds;
        Eigen::MatrixXi indexMatrix;
        igl::list_to_matrix(index, indexMatrix);
        igl::boundary_loop(indexMatrix, all_bnds);
        std::vector<std::array<double, 2>> vertexs;
        std::vector<std::array<int, 2>> edges;
        std::map<int, int> index_map; int pnum = 0;
        for (int i = 0; i < all_bnds.size(); i++) {
            std::reverse(all_bnds[i].begin(), all_bnds[i].end());
            for (int j = 0; j < all_bnds[i].size(); j++) {
                if (index_map.find(all_bnds[i][j]) == index_map.end()) {
                    index_map[all_bnds[i][j]] = pnum++;
                    vertexs.push_back({ points[all_bnds[i][j]][0],points[all_bnds[i][j]][1]});
                }
            }
            for (int j = 0; j < all_bnds[i].size(); j++) {
                edges.push_back({ index_map[all_bnds[i][j]],index_map[all_bnds[i][(j+1)% all_bnds[i].size()]] });
            }
        }





       
        bool V_rect = igl::list_to_matrix(vertexs, V);
        if (!V_rect)
        {
        //    printf(format, "V", igl::min_size(vertexs), igl::max_size(vertexs));
            return false;
        }

        bool F_rect = igl::list_to_matrix(edges, F);
        if (!F_rect)
        {
        //   printf(format, "F", igl::min_size(edges), igl::max_size(edges));
            return false;
        }

        

        return true;

    }
}
#endif