/*
 *------------------------------------------------------
 * (c) 2021, mohab metwally All Rights Reserved
 *
 * the content of this file may not be disclosed, copied, or duplicated in
 * any form, in whole, or part, without the prior wirtten permission of
 *  mohab metwally
 *------------------------------------------------------
 *
 * @author: mohab metwally
 * @date: 2021/12/1
*/

#pragma once

#include "..\Command.h"
#include "..\Definitions.h"
#include "..\Mesh.h"
#include "Point.h"
#include "..\Face.h"
#include <memory>
#include <fstream>
#include <ostream>
#include <iostream>

#include <chrono>

using namespace std;


template<class T>
class AsciiSTLReader : public Command<T>
{
public:
	AsciiSTLReader() = delete;
	AsciiSTLReader(shared_ptr<Mesh<T>> t) : Command<T>(t)
	{

	}

	~AsciiSTLReader()
	{

	}

	bool Exec(std::string file_name) 
	{


		std::shared_ptr<Mesh<T>> mesh = GetData();

		std::ifstream stlfile(file_name);
		if (!stlfile.is_open())
		{
			std::cerr << "can't open file : " + file_name;
			return false;
		}

		// check if binary
		std::string tmp, name;
		stlfile >> tmp;
		auto StrToLower = [](const std::string& str) 
		{
			std::string res = "";
			for (auto c : str) 
			{
				res += tolower(c);
			}
			return res;
		};


		std::unordered_map<std::string, int> vMap;
		T x, y, z, nx, ny, nz;

		while (stlfile >> tmp)
		{
			if (tmp == "facet") 
			{
				// read normal
				stlfile >> tmp >> nx >> ny >> nz;
				std::getline(stlfile, tmp); // "outer loop"
				std::vector<int> face_vertices;
				while (face_vertices.size() < 3) 
				{
					// read vertices
					if (stlfile >> tmp) {
						if (tmp == "vertex") {
							std::getline(stlfile, tmp);
							if (vMap.find(tmp) == vMap.end()) 
							{
								// add new vertex
								std::istringstream istr(tmp);
								istr >> x >> y >> z;
								vMap[tmp] = mesh->AddVertex(x, y, z);
							}
							// prepare the face vertices
							face_vertices.push_back(vMap[tmp]);
						}
					}
				}
				if (face_vertices.size() != 3) 
				{
					return false;
				}
				// add face
				mesh->AddFace(face_vertices);
			}
		}
		stlfile.close();
		return true;
	}
};

template<class T>
class AsciiSTLWriter : public Command<T>
{
public:
	AsciiSTLWriter() = delete;
	AsciiSTLWriter(shared_ptr<Mesh<T>> t) : Command<T>(t)
	{

	}
	~AsciiSTLWriter()
	{

	}

	bool Exec(std::string file_name)
	{
		std::cout << "writing ascii stl.." << std::endl;
		std::this_thread::sleep_for(std::chrono::seconds(3));
		auto mesh = GetData();
		std::ofstream out(file_name.c_str());
		if (!out.is_open())
		{
			return false;
		}
		const auto& faces = mesh->GetFaces();

		auto write_face = [&out, &mesh](shared_ptr<Face<T>> face)
		{
			out << "facet normal " << mesh->GetFaceNormal(face->GetIndex()) << std::endl;
			out << "outer loop" << std::endl;
			for (int j = 0; j < 3; ++j) {
				out << "vertex " << *face->GetVertex(j) << std::endl;
			}
			out << "endloop" << std::endl;
			out << "endfacet" << std::endl;
		};
		//TODO (fix) viewing this in meshlab results in missing face(the first face!), unless it's written twice it doesn't appear! i have no idea why this is happening
		write_face(faces[0]);
		for (auto face : faces)
		{
			write_face(face);
		}
		
		out << "endsolid" << std::endl;
		out.close();
		return true;
	}
};