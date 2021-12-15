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

#include <string>
#include <cassert>
#include <numeric>
#include <algorithm>
#include "..\Command.h"
#include "..\Definitions.h"
#include "..\Mesh.h"
#include "obj_io.h"
#include "VertexParameters.h"

template<class T>
class ObjReader : public Command<T>
{
public:
	ObjReader(shared_ptr<Mesh<T>> t) : Command<T>(t)
	{

	}
	~ObjReader()
	{

	}
	
	bool Exec(std::string file_name)
	{
		std::cout << "Obj Read" << std::endl;
		//TODO use member function instead of lambdas
		auto mesh = GetData();
		// We cannot assume the order in which callbacks are invoked,
  		// so we need to keep track of which vertex to add properties to.
  		// The first encountered position gets added to the first vertex, etc.
  		auto pos_count = std::uint32_t{0};
  		auto tex_count = std::uint32_t{0};
  		auto nml_count = std::uint32_t{0};
  		// V
  		// position
  		// wrap add position addition
  		auto add_position = thinks::MakeObjAddFunc<thinks::ObjPosition<T,3>>([&mesh, &pos_count](const auto &pos)
  		{
  			Vertex<T> v(pos.values[0], pos.values[1], pos.values[2]);
  			mesh->AddVertex(v);
			pos_count++;

  		});
  		// F
  		// Face
  		// wrap the Face addition
		using ObjFaceType = thinks::ObjPolygonFace<thinks::ObjIndex<std::uint32_t>>;
  		auto add_face = thinks::MakeObjAddFunc<ObjFaceType>([&mesh](const auto&face)
  		{
  			//auto poly_len = face.size();
  			std::vector<shared_ptr<Vertex<T>>> polygon;
  			for (auto index : face.values)
  			{
				auto v = mesh->GetVertex(index.value);
				polygon.push_back(v);
  			}
#ifdef DEBUG
			for (auto v : polygon)
			{
				std::cout << *v << std::endl;
			}
			std::cout << "------------" << std::endl;
#endif
  			mesh->AddFace(polygon);
  		});

		vector<VertexParameters<T>> params(pos_count);

		// Vt
		// Texture coordinates
		// wrap the vertex texture addition
		auto add_tex_coord = thinks::MakeObjAddFunc<thinks::ObjTexCoord<T,2>>([&mesh, &tex_count, &params](const auto& tex)
		{
			if (params.size() <= tex_count)
			{
				params.push_back(VertexParameters<T>());
			}
			//vector<VertexParameters<T>> params(pos_count);

			Vertex2<T> tex_coord({ tex.values[0], tex.values[1] });
			params[tex_count++].AddTexture(tex_coord);
		});

		//Vn
		// Normals
		// wrap the vertex normal addition
		auto add_normal = thinks::MakeObjAddFunc<thinks::ObjNormal<T>>([&mesh, &nml_count, &params](const auto&nml){
			if (params.size() <= nml_count)
			{
				params.push_back(VertexParameters<T>());
			}
			Vertex<T> nml_coord(nml.values[0], nml.values[1], nml.values[2]);
			params[nml_count++].SetNormal(nml_coord);
		});
		auto ifs = std::ifstream(file_name.c_str());
		if (!ifs) {
			return false;
		}
		const auto result = thinks::ReadObj(ifs, add_position, add_face, add_tex_coord, add_normal);
		ifs.close();
		mesh->SetVertexParamters(params);

		if (!(result.position_count == mesh->GetNumberOfVertices())) {
			return false;
		}
		return true;
	}
	
};

template<class T>
class ObjWriter : public Command<T>
{
public:
	ObjWriter(shared_ptr<Mesh<T>> t) : Command<T>(t)
	{

	}
	~ObjWriter()
	{

	}
	bool Exec(std::string file_name)
	{
		std::cout << "Obj Write" << std::endl;
		//TODO use member function instead of lambdas
		shared_ptr<Mesh<T>> mesh = GetData();

		int vertex_count = mesh->GetNumberOfVertices();
		int face_count = mesh->GetNumberOfFaces();
		// V
		// Position 
		// Vertex writer wrapper
		int vertex_counter = 0;
		auto pos_mapper = [&vertex_counter, &mesh]() {
			using ObjPositionType = thinks::ObjPosition<T,3>;
			if (vertex_counter >= mesh->GetNumberOfVertices())
			{
				return thinks::ObjEnd<ObjPositionType>();
			}
			Vertex<T> v = *mesh->GetVertex(vertex_counter++);
			return thinks::ObjMap(ObjPositionType(v[0], v[1], v[2]));
		};
		// F
		// Face
		// Face Writer wrapper
		int face_counter = 0;
		auto face_mapper = [&face_counter, &mesh]() {
			using ObjIndexType = thinks::ObjIndex<std::uint32_t>;
			using ObjFaceType = thinks::ObjPolygonFace<ObjIndexType>;
			if (face_counter >= mesh->GetNumberOfFaces())
			{
				return thinks::ObjEnd<ObjFaceType>();
			}
			std::shared_ptr<Face<T>> f = mesh->GetFace(face_counter++);
			//std::vector<ObjIndexType> face_indices;
			auto obj_face = ObjFaceType{};
			for (auto v : f->GetVertices())
			{
				int v_idx = mesh->GetIndex(v);
				obj_face.values.push_back(ObjIndexType(v_idx));
			}
			return thinks::ObjMap(obj_face);

		};


		Vertex2<T> dummy_tex({ 0,0 });
		// Vt
		// Vertex Texture
		// Vertex Texture wrapper
		int vertex_tex_counter = 0;
		auto tex_mapper = [&vertex_tex_counter, &mesh]() {
			using ObjTexCoordType = thinks::ObjTexCoord<T,2>;
			if (vertex_tex_counter >= mesh->GetVertexParametersSize())
			{
				return thinks::ObjEnd<ObjTexCoordType>();
			}
			VertexParameters<T> tex = mesh->GetVertexParamter(vertex_tex_counter++);
			auto tex_vertex = tex.GetTexture();
			//auto v_ptr = mesh->GetVertex(vertex_tex_counter++);
			//const auto tex = v_ptr->GetTexture();
			return thinks::ObjMap(ObjTexCoordType(tex_vertex[0], tex_vertex[1]));
		};

		// Vn
		// Vertex Normal
		// Vertex Normal wrapper
		int vertex_nml_counter = 0;
		auto nml_mapper = [&vertex_nml_counter, &mesh]() {
			using ObjNormalType = thinks::ObjNormal<T>;
			if (vertex_nml_counter >= mesh->GetVertexParametersSize())
			{
				return thinks::ObjEnd<ObjNormalType>();
			}
			VertexParameters<T> nml = mesh->GetVertexParamter(vertex_nml_counter++);
			Vertex<T> nml_vertex = nml.GetNormal();
			//auto v_ptr = mesh->GetVertex(vertex_nml_counter++);
			//const auto nml = v_ptr->GetNormal();
			return thinks::ObjMap(ObjNormalType(nml_vertex[0], nml_vertex[1], nml_vertex[2]));
		};

 		// Open the OBJ file and pass in the mappers, which will be called
  		// internally to write the contents of the mesh to the file.
		
		std::cout << "writing Obj to file_name: " << file_name << std::endl;
  		auto ofs = std::ofstream(file_name.c_str());
		if (!ofs)
		{
			return false;
		}
  		const auto result = thinks::WriteObj(ofs, pos_mapper, face_mapper, tex_mapper, nml_mapper);
 	 	ofs.close();
		std::cout << "mesh writen to Obj file"  << std::endl;
		return true;
	}
};

