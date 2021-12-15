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

#include <memory>
#include "..\Definitions.h"
#include "..\Command.h"
#include "..\Point.h"
#include "..\Face.h"
#include "..\Mesh.h"

using namespace std;



template<class T>
class BinarySTLReader : public Command<T>
{
public:
	BinarySTLReader() = delete;

	BinarySTLReader(shared_ptr<Mesh<T>> t) : Command<T>(t)
	{

	}
	~BinarySTLReader()
	{

	}

	bool Exec(std::string file_name) 
	{
		std::cout << "reading binary stl" << std::endl;
		std::shared_ptr<Mesh<T>> mesh = GetData();
		std::ifstream inp(file_name.c_str(), std::ios::in | std::ios::binary);
		if (!inp.is_open())
		{
			std::cerr << "can't ope the file" << std::endl;
			return false;
		}
		char buffer[80] = { 0 };
		inp.read(buffer, 80);

		unsigned int numFaces = 0;
		// making sure it is little endian
		for (int i = 0; i < 4; ++i) {
			inp.read(buffer, 1);
			unsigned int tmp = (unsigned char)(buffer[0]);
			tmp <<= (i * 8);
			numFaces = numFaces | tmp;
		}
		std::cout << "numFaces = " << numFaces << std::endl;
		std::string key;
		std::vector<T> tmp_data(3);
		std::vector<int> face_indeces(3);
		std::unordered_map<std::string, int> vertexIndex;

		auto read = [&](bool add) {
			inp.read(buffer, 12);

			tmp_data[0] = *((float*)&(buffer[0]));
			tmp_data[1] = *((float*)&(buffer[4]));
			tmp_data[2] = *((float*)&(buffer[8]));

			if (add) {
				key = std::string(buffer, buffer + 12);
				if (vertexIndex.find(key) == vertexIndex.end()) 
				{
					vertexIndex[key] = mesh->AddVertex(tmp_data);
				}
				return vertexIndex[key];
			}
			return -1;
		};

		for (int i = 0; i < numFaces; ++i)
		{
			// normal
			read(false);
			// point 1
			face_indeces[0] = read(true);
			// point 2
			face_indeces[1] = read(true);
			// point 3
			face_indeces[2] = read(true);

			inp.read(buffer, 2);
			// add face
			mesh->AddFace(face_indeces);
		}
		mesh->CalculateFaceNormals();
		std::cout << "calculate face normals" << std::endl;

		inp.close();
		return true;
	}
};

template<class T>
class BinarySTLWriter : public Command<T>
{
public:
	BinarySTLWriter() = delete;
	BinarySTLWriter(shared_ptr<Mesh<T>> t) : Command<T>(t)
	{

	}
	~BinarySTLWriter()
	{

	}
	bool Exec(std::string file_name) 
	{ 
		std::cout << "writing binary stl" << std::endl;
		std::this_thread::sleep_for(std::chrono::seconds(3));
		
		std::shared_ptr<Mesh<T>> mesh = GetData();
		std::ofstream out(file_name.c_str(), std::ios::out | std::ios::binary);
		if (!out.is_open())
		{
			std::cout << "can't open file for writing: " << file_name << std::endl;
			return false;
		}
		char header[80] = { 0 };
		out.write(header, 80);

		vector<shared_ptr<Face<T>>> faces = mesh->GetFaces();
		unsigned int numFaces = faces.size();
		// making sure it is little endian
		for (int i = 0; i < 4; ++i) {
			out << (unsigned char)(numFaces & 0xff);
			numFaces >>= 8;
		}

		for (int i =0; i < faces.size(); i++)
		{
			shared_ptr<Vertex<T>> p1 = faces[i]->Getv1();
			shared_ptr<Vertex<T>> p2 = faces[i]->Getv2();
			shared_ptr<Vertex<T>> p3 = faces[i]->Getv3();
			Vertex<T> n = mesh->GetFaceNormal(faces[i]->GetIndex());
			auto write = [&](const auto& x, const auto& y, const auto& z) {
				float tmpX = x;
				float tmpY = y;
				float tmpZ = z;
				out.write(reinterpret_cast<const char*>(&tmpX), sizeof(float));
				out.write(reinterpret_cast<const char*>(&tmpY), sizeof(float));
				out.write(reinterpret_cast<const char*>(&tmpZ), sizeof(float));
			};
			write(n.x(), n.y(), n.z());
			write(p1->x(), p1->y(), p1->z());
			write(p2->x(), p2->y(), p2->z());
			write(p3->x(), p3->y(), p3->z());
			out.write(header, 2);
		}
		out.close();
		std::cout << "exiting Write Exec" << std::endl;
		return true;
	}
};
