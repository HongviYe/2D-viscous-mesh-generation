#ifndef _NEAREST_SEARCH_H_
#define _NEAREST_SEARCH_H_

#include <Eigen/core>
/*
* * this file contain the easiest implementation of nearest search
* We should use AABB tree or kd tree to do that if we have free time
*/

static Eigen::VectorXd nearestDistance(const Eigen::MatrixXd& V, const Eigen::MatrixXd& V__compare, const Eigen::MatrixXi& F) {
	Eigen::VectorXd ans(F.rows());
	Eigen::VectorXd pnt(V.rows());
	for (int i = 0; i < V.rows(); i++) {
		pnt(i) = std::numeric_limits<double>::max();
		for (int j = 0; j < V__compare.rows(); j++) {
			pnt(i) = std::min(pnt(i),sqrt(pow(V(j,0)- V__compare(i,0),2)+ pow(V(j, 1) - V__compare(i, 1), 2)));
		}
	}
	for (int i = 0; i < F.rows(); i++) {
		ans(i) = pnt(F(i, 0));
	}
	for (int i = 0; i < F.rows(); i++) {
		ans(i) += pnt(F(i, 1));
	}
	for (int i = 0; i < F.rows(); i++) {
		ans(i) /=2;
	}

	return ans;
}

#endif //! _NEAREST_SEARCH_H_