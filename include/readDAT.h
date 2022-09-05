#ifndef _READ_DAT_H_
#define _READ_DAT_H_

#include <igl/list_to_matrix.h>

#include <Eigen/Core>
#include <string>
#include <vector>
#include <array>
#include <algorithm>
#include <fstream>
#include <cstdio>

namespace expan {
    template <typename DerivedV, typename DerivedF>
    bool readNACADAT(
        const std::string str,
        Eigen::PlainObjectBase<DerivedV>& V,
        Eigen::PlainObjectBase<DerivedF>& F)
    {
        std::vector < std::array < double,2 >> vertexs;
        std::vector<std::array<int,2>> edges;

        std::ifstream fin(str);
        std::string str1;
        std::getline(fin,str1);
        int number_of_point_count = 0;
        while (true) {
            std::array<double, 2> vs;
            if (!(fin >> vs[0] >> vs[1])) {
                edges.pop_back();
                edges.push_back({number_of_point_count-1, 0});
                break;
            }
            else {
                vertexs.push_back(vs);
                edges.push_back({number_of_point_count,number_of_point_count+1});
                number_of_point_count++;
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