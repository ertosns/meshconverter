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

#include <algorithm>
#include <cmath>
#include <vector>
#include "Transformation.h"
#include "Element.h"
#include <iostream>
#include <ostream>
#include <fstream>

template<class T, unsigned int N>
class Point : public Element<T>
{
public:
  Point() 
  {
  }

  /**
   * This constructor initializes the point using vector of data
   * @param  values
   */
  Point(const std::vector<T>& values) : Element<T>()
  {
    for (unsigned int i = 0; i < N; ++i)
    {
      data_[i] = values[i];
    }
  }


  /**
  *  This is the copy constructor
   * @param  point
   */
  Point(const Point& point) : Element<T>()
  {
    for (unsigned int i = 0; i < N; ++i) 
    {
      data_[i] = point.data_[i];
    }
  }


  /**
   * Empty Destructor
   */
  virtual ~Point()
  {
  }



  /**
   * Compare two points
   * @return bool
   * @param  point
   */
  bool operator==(const Point& point) const
  {
    for (unsigned int i = 0; i < N; ++i) 
    {
      if (data_[i] != point.data_[i]) {
        return false;
      }
    }

    return true;
  }

  /**
  * Compare two points for inequality
  * @return bool
  * @param  point
  */
  bool operator!=(const Point& point) const
  {
	  for (unsigned int i = 0; i < N; ++i) 
    {
		  if (data_[i] == point.data_[i]) {
			  return false;
		  }
	  }

	  return true;
  }

  /**
   * Assignment operator
   * @return the current point
   * @param  point
   */
  template<class type>
  void operator=(const type& point)
  {
    for (unsigned int i = 0; i < N; ++i)
    {
      data_[i] = point.data_[i];
    }
  }

  /**
   * Add two points
   * @return Point
   * @param  point
   */
  //template<class type>
  Point operator+(const Point& point) const
  {
    Point res;
    for (unsigned int i = 0; i < N; ++i) 
    {
      res.data_[i] = data_[i] + point.data_[i];
    }

    return res;
  }

  /*
  * Get distance to the point 
  * @param pt is the point to get distance to 
  * @return the distance
  */
  T DistanceTo(const Point& pt) const 
  {
    return (pt - *this).GetL2Norm();
  }

  /**
   * Add a point to this point
   * @param  point
   */
  template<class type>
  void operator+=(const type& point)
  {
    for (unsigned int i = 0; i < N; ++i) 
    {
      data_[i] += point.data_[i];
    }
  }

  /**
   * Subtract two points
   * @return Point
   * @param  point
   */
  Point operator-(const Point& point) const
  {
    Point res;
    for (unsigned int i = 0; i < N; ++i) {
      res.data_[i] = data_[i] - point.data_[i];
    }

    return res;
  }

  /**
   * Subtract a point from this point
   * @param  point
   */
  template<class type>
  void operator-=(const type& point)
  {
    for (unsigned int i = 0; i < N; ++i) {
      data_[i] -= point.data_[i];
    }
  }

  /**
   * Negate point
   * @return Point multiplied by -1
   */
  Point operator-()
  {
    Point res;
    for (unsigned int i = 0; i < N; ++i) {
      res.data_[i] = -data_[i];
    }

    return res;
  }


  /**
   * Multiply point with factor
   * @return Point
   * @param  factor
   */
  template<class type>
  Point operator*(type factor) const
  {
    Point res;
    for (unsigned int i = 0; i < N; ++i) 
    {
      res.data_[i] = data_[i] * factor;
    }

    return res;
  }

  /**
   * Multiply this point with a factor
   * @param  factor
   */
  template<class type>
  void operator*=(type factor)
  {
    for (unsigned int i = 0; i < N; ++i) 
    {
      data_[i] *= factor;
    }
  }


  /**
   * Divide point by a factor
   * @return Point
   * @param  denominator
   */
  template<class type>
  Point operator/(type denominator) const
  {
    if (denominator == 0) {
      throw "Can not divide by Zero!";
    }

    Point res;
    for (unsigned int i = 0; i < N; ++i) {
      res.data_[i] = data_[i] / denominator;
    }

    return res;
  }

  /**
   * Divide this point by a factor
   * @param  denominator
   */
  template<class type>
  void operator/=(type denominator)
  {
    if (denominator == 0) {
      throw "Can not divide by Zero!";
    }

    for (unsigned int i = 0; i < N; ++i) {
      data_[i] /= denominator;
    }
  }

  Point Cross(const Point& point_b) const
  {
      Point<T, 3> res;

      res[0] = Point<T, 3>::GetData(1) * point_b.GetData(2) - Point<T, 3>::GetData(2) * point_b.GetData(1);
      res[1] = Point<T, 3>::GetData(2) * point_b.GetData(0) - Point<T, 3>::GetData(0) * point_b.GetData(2);
      res[2] = Point<T, 3>::GetData(0) * point_b.GetData(1) - Point<T, 3>::GetData(1) * point_b.GetData(0);

      return res;
  }


  /**
  *  This function performs the dot product between the current point and another point of the same size and type
   * @return T
   * @param  point_b
   */
  T Dot(const Point& point_b) const
  {
    T res = 0;
    for (unsigned int i = 0; i < N; ++i) 
    {
      res += data_[i] * point_b.data_[i];
    }

    return res;
  }


  /**
  *  This function gets the L2 norm of the current point
   * @return double
   */
  T GetL2Norm() const
  {
    return std::sqrt(Dot(*this));
  }


  /**
  *  This function normalizes the current point
   */
  void Normalize()
  {
    T normal = GetL2Norm();
    if (normal > 0) {
      for (unsigned int i = 0; i < N; ++i) 
      {
        data_[i] = data_[i] / normal;
      }
    }
  }

  /*
  * Get normalized version of this points
  */
  Point GetNormalized()
  {
    Point res(*this);
    res.Normalize();
    return res;
  }

  /**
  * This operator is indexing the point coefficients
   * @param  index
   */
  T& operator[](int index)
  {
    if (index < 0 || index >= N) 
    {
      throw "Subscript out of range!";
    }
    return data_[index];
  }

  /*
  * This function Tranforms the point using transformation matrix
  */
  virtual void Transform(Transformation<T>& transformation) 
  {
    Point<T, 4> pt({0, 0, 0, 1});
    auto size = std::min(3u, N);
    for (unsigned int i = 0; i < size; ++i) {
      pt[i] = data_[i];
    }
    Point<T, 4> tmp = transformation * pt;
    for (unsigned int i = 0; i < size; ++i) {
      data_[i] = tmp[i];
    }
  }


  /*
  * Returns true if all elements are zeros and false othjerwise
  */
  bool IsZero() const 
  {
    for (int i = 0; i < N; ++i) 
    {
      if (std::abs(data_[i]) > T(0.0001)) 
      {
        return false;
      }
    }
    return true;
  }
  /*
  * This function returns the value at the index specified
  * @param index is the index spcified inside the point
  */
  T GetData(int index) const 
  {
    return data_[index];
  }

  /*
  * Get teh absolute of a point
  */
  void MakeAbs() 
  {
    for (int i = 0; i < N; ++i) 
    {
      if (data_[i] < 0) {
        data_[i] = -data_[i];
      }
    }
  }

  /*
  * This function prints the point
  */
  void Print() const 
  {
    std::cout << "[" << data_[0];
    for (int i = 1; i < N; ++i) 
    {
      std::cout << ", " << data_[i];
    }
    std::cout << "]" << std::endl;
  }

  //friend std::ostream operator<<(std::ostream& os, Point<T,N>& v);

 size_t Size() const
 {
  return N;
 }
protected:
  T data_[N];
};

template<class T> using Vertex2 = Point<T,2>; // for parameters
//template<class T> using Vertex = Point<T,3>; // for mesh coordinates
template<class T> using Vertex4 = Point<T,4>; // for transformation vector

template<class T>
class Vertex : public Point<T, 3>
{
public:
    Vertex() : Point<T,3>()
    {

    }


    Vertex(const Point<T,3>& p) : Point<T,3>()
    {
        data_[0] = p.GetData(0);
        data_[1] = p.GetData(1);
        data_[2] = p.GetData(2);
    }


    Vertex(T x, T y, T z) : Point<T,3>()
    {

        data_[0] = x;
        data_[1] = y;
        data_[2] = z;
    }

    Vertex(const Vertex<T>& p) : Point<T, 3>()
    {

        data_[0] = p.x();
        data_[1] = p.y();
        data_[2] = p.z();
    }


    Vertex(const std::vector<T>& values) : Point<T, 3>(values)
    {

    }
    //TODO validate
    T x() const
    {
        return data_[0];
    }

    T y() const
    {
        return data_[1];
    }

    T z() const
    {
        return data_[2];
    }

    template<class T>
    friend std::ostream& operator<<(std::ostream& os, const Vertex<T>& v);
};


template<class T>
std::ostream& operator<<(std::ostream& os, const Vertex<T>& v)
{
    size_t len = v.Size();
    os << v.x() << kVertexSeparator << v.y() << kVertexSeparator << v.z();
    return os;
};