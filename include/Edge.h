#pragma once

#include <unordered_set>
#include <memory>

#include "Element.h"
#include "Point.h"
#include "Transformation.h"

using namespace std;
template<class T>
class Edge : public Element<T>
{
 
public:


  /**
   * Empty Constructor
   */
  Edge() : Element<T>()
  {
  }

  /**
   * Empty Destructor
   */
  virtual ~Edge()
  {
  }


  /**
  * Constructor using two vertices
   * @param  v1 The first vertex
   * @param  v2 The second vertex
   */
  Edge(const shared_ptr<Vertex<T>>& v1, const shared_ptr<Vertex<T>>& v2) :  v1_(v1), v2_(v2), Element<T>()
  {
  }

  /**
  * Set the value of the first vertex
  * @param value the new value of v1_
  */
  void SetV1(shared_ptr<Vertex<T>>& value)
  {
      v1_ = value;
  }
  /**
   * Get the value of the first vertex
   * @return the value of v1_
   */
  shared_ptr<Vertex<T>> GetV1() const
  {
    return v1_;
  }

  /**
   * Set the value of the second vertex
   * @param value the new value of v2_
   */
  void SetV2(shared_ptr<Vertex<T>>& value)
  {
    v2_ = value;
  }

  /**
   * Get the value of the second veretx
   * @return the value of v2_
   */
  shared_ptr<Vertex<T>> GetV2() const
  {
    return v2_;
  }

  /*
  * This function  checks if this object is connected to the given vertex
  * @return true if the given vertex is directly connected to this object
  * @param shared_ptr that contains the base class to vertex
  */
  virtual bool IsConnectedToVertex(const std::shared_ptr<Vertex<T>>& vertex) const 
  {
    if (GetV1() == vertex || GetV2() == vertex) 
    {
      return true;
    }
    return false;
  }

  /* 
  * This function taked a vertex that has to be one of the 2 ends of the edge and returns the other vertex
  * @param vertex is the first vertex in this edge
  * @returns the other end of the current edge
  */
  shared_ptr<Vertex<T>> GetOtherVertex(const shared_ptr<Vertex<T>>& vertex) 
  {
    if (v1_ == vertex) 
    {
      return v2_;
    }
    // vertex has to be v1_
    return v1_;
  }

  /**
  * This function transforms the edge using the given transformation matrix
  * @param  transformation
  */
  virtual void Transform(Transformation<T>& transformation) 
  {
    v1_->Transform(transformation);
    v2_->Transform(transformation);
  }

  /*
  * Compare if the two edges have the same two vertex indecs
  */
  bool operator==(const Edge<T>& edge) {
    return (v1_->GetIndex() == edge->v1_->GetIndex() && v2_->GetIndex() == edge->v2_->GetIndex()) 
      || (v2_->GetIndex() == edge->v1_->GetIndex() && v1_->GetIndex()== edge->v2_->GetIndex());
  }

private:
  shared_ptr<Vertex<T>> v1_;
  shared_ptr<Vertex<T>> v2_;
};
