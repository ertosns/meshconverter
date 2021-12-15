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
#include "Transformation.h"

template<class T>
class Element 
{
public:

	/**
	 * Empty Constructor
	 */
	Element(){}

	/**
	 * Empty Destructor
	 */
  virtual ~Element() {}

    
  void SetIndex(int value)
  {
    index_ = value;
  }

  
  int GetIndex() const
  {
    return index_;
  }

  virtual void Transform(Transformation<T>&) = 0;

private:
	// index of the Element
	int index_;


};

