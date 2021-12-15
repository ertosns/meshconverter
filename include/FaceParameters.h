#pragma once

#include "Point.h"

template<class T>
class FaceParameters
{
public:
  /**
   * Empty Constructor
   */
   // Constructors/Destructors
   //  
  FaceParameters()
  {
  }

  /**
   * Empty Destructor
   */
  virtual ~FaceParameters()
  {
  }

  /**
   * Set the value of area_
   * @param value the new value of area_
   */
  void SetArea(const T& value)
  {
    area_ = value;
  }

  /**
   * Get the value of area_
   * @return the value of area_
   */
  T GetArea() const
  {
    return area_;
  }

  /**
   * Set the value of normal_
   * @param value the new value of normal_
   */
  void SetNormal(const Vertex<T>& value)
  {
    normal_ = value;
  }

  /**
   * Get the value of normal_
   * @return the value of normal_
   */
  Vertex<T> GetNormal() const
  {
    return normal_;
  }

  /**
   * Flip the normal 
   */
  void FlipNormal()
  {
    normal_ = -normal_;
  }


  /**
   * Set the value of mean_
   * @param value the new value of mean_
   */
  void SetMean(const Vertex<T>& value)
  {
    mean_ = value;
  }

  /**
   * Get the value of mean_
   * @return the value of mean_
   */
  Vertex<T> GetMean() const
  {
    return mean_;
  }

  /**
   * Set the value of centroid_
   * @param value the new value of centroid_
   */
  void SetCentroid(const Vertex<T>& value)
  {
    centroid_ = value;
  }

  /**
   * Get the value of centroid_
   * @return the value of centroid_
   */
  Vertex<T> GetCentroid() const
  {
    return centroid_;
  }

private:
  /*
  * Face area
  */
  T area_;
  /*
  * face normal
  */
  Vertex<T> normal_;
  /*
  * face mean
  */
  Vertex<T> mean_;
  /*
  * face centroid
  */
  Vertex<T> centroid_;

};