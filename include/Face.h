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
#include <vector>
#include "Point.h"
#include "Edge.h"
#include "Element.h"
#include "Transformation.h"

using namespace std;

template<class T>

class Face : public Element<T>
{

public:

  /**
   * Empty Constructor
   */
  Face() : Element<T>()
  {
  }

  /**
   * Empty Destructor
   */
  virtual ~Face()
  {
  }


  /**
  * This is a constructor using array of vertices
   * @param  vertices are the vertices making this face
   */
  Face(const std::vector<std::shared_ptr<Vertex<T>>>& vertices)
    : vertices_(vertices), Element<T>()
  {
  }

  /*
  * Get first vertex
  * @return the first vertex in this face
  */
  std::shared_ptr<Vertex<T>> Getv1()  const
  {
    return vertices_[0];
  }

  /*
  * Get second vertex
  * @return the first vertex in this face
  */
  std::shared_ptr<Vertex<T>> Getv2()  const
  {
    return vertices_[1];
  }

  /*
  * Get third vertex
  * @return the first vertex in this face
  */
  std::shared_ptr<Vertex<T>> Getv3()  const
  {
    return vertices_[2];
  }

  /*
   * Get vertex index
   * @param index is the index of the veretex to be returned
   * @return the vertex at index index in this face
   */
  std::shared_ptr<Vertex<T>> GetVertex(int index) const
  {
    if (index < 0 || index >= vertices_.size()) 
    {
      throw "Out of range!";
    }
    return vertices_[index];
  }

  /**
   * Set the value of vertices_
   * vertices needs to be orderd
   * @param value the new value of vertices_
   */
  void SetVertices(const std::vector<std::shared_ptr<Vertex<T>>>& value)
  {
    vertices_ = value;
  }

  /**
   * Get the value of vertices_
   * vertices needs to be orderd
   * @return the value of vertices_
   */
  std::vector<std::shared_ptr<Vertex<T>>>& GetVertices()
  {
    return vertices_;
  }

  /*
  * This function  checks if this object is connected to the given vertex
  * @return true if the given vertex is directly connected to this object
  * @param shared_ptr that contains the base class to vertex
  */
  virtual bool IsConnectedToVertex(const std::shared_ptr<Vertex<T> >& vertex) const 
  {
    for (auto& v : vertices_) 
    {
      if (vertex == v) {
        return true;
      }
    }
    return false;
  }

  /**
  * This function transforms the face using the given transformation matrix
  * @param  transformation_matrix
  */
  virtual void Transform(Transformation<T>& transformation) 
  {
    for (auto v : vertices_) {
      v->Transform(transformation);
    }
  }

  /*
  * Compare if the two faces have the same vertex indecs
  */
  bool operator==(const Face& face) 
  {
    // number of vertices
    if (vertices_.size() != face.vertices_.size()) 
    {
      return false;
    }
    const int size = vertices_.size();
    // get the starting point
    int fInd = -1;
    for (int i = 0; i < size; ++i) {
      if (face.vertices_[i]->GetIndex() == vertices_[0]->GetIndex()) 
      {
        fInd = i;
        break;
      }
    }
    // no match found
    if (fInd == -1) 
    {
      return false;
    }
    // check all indeces
    for (int i = 0; i < size; ++i) {
      if (face.vertices_[(fInd + i) % size]->GetIndex() != vertices_[i]->GetIndex()) 
      {
        return false;
      }
    }
    // all indeces are the same
    // not necessarily starting from the same index inside the face
    // but they have the same indeces in the same order
    return true;
  }

  /*
  * This function is mainly used to reverse the order of the vertices in this face
  * It will be mainly used for flipping the normal of this face
  */
  void ReverseOrder() 
  {
    std::reverse(vertices_.begin(), vertices_.end());
  }

  /*
  * This function returns the number of vertices making this face
  */
  int GetNumberOfVertices() const 
  {
    return vertices_.size();
  }


private:

  std::vector<std::shared_ptr<Vertex<T>>> vertices_;

};
