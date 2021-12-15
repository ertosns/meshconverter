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

#include <fstream>
#include <memory>
#include <vector>
#include <string>
#include "Mesh.h"

using namespace std;

template<class T>
class Command
{
public:
	Command(shared_ptr<Mesh<T>> t) : data_(t)
	{

	}
	~Command()
	{

	}

	Command(const Command<T>& command) : data_(command.data_)
	{
	
	}

	Command(const std::shared_ptr<Command<T>>& command) : data_(command->data_)
	{

	}

	virtual bool Exec(std::string) = 0;

	shared_ptr<Mesh<T>> GetData()
	{
		return data_;
	}
protected:
	
private:
	shared_ptr<Mesh<T>> data_; /*representation of data*/
};




