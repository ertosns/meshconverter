//#define DEBUG 1

#include <cmath>
#include <random>
#include <memory>
#include <vector>
#include "converters/convert.h"
#include "io/AsciiStl.h"
#include <cassert>
#include <exception>



int main(int argc, char **args)
{
	std::string nef_in_obj("../resources/N_E_F_E_R_T_I_T_I_.obj");
	std::string nef_in_stl("../resources/N_E_F_E_R_T_I_T_I_.stl");

	std::string car_in_stl("../resources/car1.obj");
	std::string cube_in_stl("../resources/cube.stl");
	std::string cube_in_obj("../resources/cube.obj");
	std::string shapr_cube_in_obj("../resources/shapr_cube.obj");


	std::string cup_in_obj("../resources/Cup.obj");

	std::string cup_out_stl("../resources/Cup.stl");
	
	std::string out_obj("../resources/out.obj");
	std::string out_stl("../resources/out.stl");
	
	auto mesh = make_shared<Mesh<double>>();
	
	Converter<double> oneway;

	try 
	{
		oneway = Converter<double>(MeshType::kObjText, MeshType::kSTLText, mesh);
		//twoway = Converter<float>(MeshType::kObjText,  MeshType::kSTLText);
	}
	catch (exception &e)
	{
		std::cout << "ERROR: " << e.what() << std::endl;
		std::this_thread::sleep_for(std::chrono::seconds(3));

	}
	catch (...)
	{
		std::cout << "Error happend converting..." << std::endl;
		std::this_thread::sleep_for(std::chrono::seconds(3));

	}
	oneway.Convert(shapr_cube_in_obj, out_stl);
	std::cout << "finished..." << std::endl;
	return 0;
}