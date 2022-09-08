#ifndef _WRITE_VTK_H_
#define _WRITE_VTK_H_
#include <Eigen/Core>
#include <string>
#include <fstream>
namespace RIGIDT {
	void header(std::ostream& fout) {		
		fout << "# vtk DataFile Version 2.0" << std::endl;
		fout << "boundary layer mesh" << std::endl;
		fout << "ASCII" << std::endl;
		fout << "DATASET UNSTRUCTURED_GRID" << std::endl;
	}
	template <typename DerivedV, typename DerivedF>
	bool writeVTK(
		const std::string str,
		const Eigen::MatrixBase<DerivedV>& V,
		const Eigen::MatrixBase<DerivedF>& F)
	{
		using namespace std;
		using namespace Eigen;
		std::ofstream fout(str);
		if (!fout.is_open())
			return false;
		header(fout);
		fout << "POINTS " << V.size() << " double" << std::endl;
		fout <<V.format(IOFormat(FullPrecision, DontAlignCols, " ", "\n", "", "", "", "")); 
		fout << "CELLS " << F.rows() << " " << F.rows() * (F.cols() + 1) << std::endl;
		fout <<(F.array()).format(IOFormat(FullPrecision, DontAlignCols, " ", "\n", "3 ", "", "", ""));
		fout << "CELL_TYPES " << F.rows() << std::endl;
		for (int i = 0; i < F.rows(); i++) {
			fout << 5 << std::endl;
		}
		fout.close();
		return true;
	}
}
#endif //! _WRITE_VTK_H_