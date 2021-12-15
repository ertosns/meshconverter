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

/* \brief class convert 3d mesh from type A to type B
*/

#include <string>
#include <memory>
#include <thread>
#include <chrono>

#include "../Mesh.h"
#include "../Definitions.h"
#include "../Command.h"
#include "../io/Obj.h"
#include "../io/AsciiStl.h"
#include "../io/BinaryStl.h"

enum MeshType
{
	kSTLText,
	kSTLBinary,
	kObjText,
	kObjBinary
};

template<class T>
class CommandBuilder
{
	typedef std::shared_ptr<Command<T>> CommandPtr;
public:
	CommandBuilder() : command_(nullptr)
	{

	}
	CommandBuilder(CommandBuilder& builder)
	{
		type_ = builder.GetType();
		command_ = builder.GetCommand();
	}
	CommandBuilder(MeshType mesh_type) : type_(mesh_type), command_(nullptr)
	{
		
	}
	void operator=(const CommandBuilder& builder)
	{
		type_ = builder.type_;
		command_ = builder.command_;
	}
	void Exec(std::string file_path)
	{
		std::cout << "command builder exec" << std::endl;
		command_->Exec(file_path);
	}
	MeshType GetType()
	{
		return type_;
	}
	CommandPtr GetCommand()
	{
		return command_;
	}

	shared_ptr<Mesh<T>> GetData()
	{
		return command_->GetData();
	}
protected:
	MeshType type_;
	CommandPtr command_;
};

template<class T>
class ReaderBuilder : public CommandBuilder<T>
{
public:
	ReaderBuilder() : CommandBuilder<T>()
	{

	}
	ReaderBuilder(ReaderBuilder &builder) : CommandBuilder<T>(builder)
	{

	}
	ReaderBuilder(std::shared_ptr<Mesh<T>> mesh, MeshType mesh_type) : CommandBuilder<T>(mesh_type)
	{
		switch (type_)
		{
		case kSTLText: 
		{
			command_ = std::make_shared<AsciiSTLReader<T>>(mesh);
			break;
		}
		case kSTLBinary:
		{
			command_ = std::make_shared<BinarySTLReader<T>>(mesh);
			break;
		}
		case kObjText:
		{
			command_ = std::make_shared<ObjReader<T>>(mesh);
			break;
		}
		default:
			throw "this type isn't supported yet";
		}
	}
};

template<class T>
class WriterBuilder : public CommandBuilder<T>
{
public:
	WriterBuilder() : CommandBuilder<T>()
	{

	}
	WriterBuilder(WriterBuilder& builder) : CommandBuilder<T>(builder)
	{

	}
	WriterBuilder(std::shared_ptr<Mesh<T>> mesh, MeshType mesh_type) : CommandBuilder<T>(mesh_type)
	{
		switch (type_)
		{
		case kSTLText: 
		{
			command_ = std::make_shared<AsciiSTLWriter<T>>(mesh);
			break;
		}
		case kSTLBinary:
		{
			command_ = std::make_shared<BinarySTLWriter<T>>(mesh);

			break;
		}
		case kObjText:
		{
			command_ = std::make_shared<ObjWriter<T>>(mesh);
			break;
		}

		default:
			throw "this type isn't supported yet";
		}
	}
};

template<class T>
class Converter
{
public:
	/*
	* Converter constructor from 3d
	*/
	Converter()
	{

	}
	Converter(Converter<T>& converter) : reader_(converter.GetReader()), writer_(converter.GetWriter())
	{

	}
	Converter(MeshType in_type, MeshType out_type, shared_ptr<Mesh<T>> mesh)
		: mesh_(mesh), reader_(mesh, in_type), writer_(mesh, out_type)
	{
	}
	~Converter()
	{
		//
	}

	void Convert(std::string in, std::string out)
	{
		//
		std::cout << "reading..." << std::endl;
		//
		reader_.Exec(in);
		//

		//
		ProcessMesh();




		//
		std::cout << "writing..." << std::endl;
		writer_.Exec(out);


	}

	void operator=(const Converter<T>& converter)
	{
		reader_ = converter.reader_;
		writer_ = converter.writer_;
	}

	ReaderBuilder<T> GetReader()
	{
		return reader_;
	}
	WriterBuilder<T> GetWriter()
	{
		return writer_;
	}

private:
	ReaderBuilder<T> reader_;
	WriterBuilder<T> writer_;
	shared_ptr<Mesh<T>> mesh_;
	/*
	MeshType GetMeshExtType(std::string file_name)
	{
		//parse the files, and detect their types.
	}
	*/
	void ProcessMesh() 
	{
		auto mesh = reader_.GetData();
		mesh->Triangulate();
		//reader_.GetData()->CalculateFaceNormals();

		std::cout << "vertices: " << mesh->GetNumberOfVertices() << std::endl;
		std::cout << "faces: " << mesh->GetNumberOfFaces() << std::endl;
		std::cout << "volume: " << mesh->CalculateWatertightVolume() << std::endl;
		std::cout << "area: " << mesh->GetArea() << std::endl;
	}
};