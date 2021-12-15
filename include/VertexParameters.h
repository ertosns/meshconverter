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

#include "Point.h"

enum VertexParamType
{
  kNormal,
  kTexture,
};


template<class T>
class VertexParameters
{
    //TODO (fix) use variant instead
    // also note vertex has nml, and tex decide which to hold those params
  //typedef std::variant<Vertex3<T>, Vertex2<T>, Curvature> ParameterType;

public:

    VertexParameters() : tex_(Vertex2<T>({ 0,0 })), nml_(Vertex<T>{0, 0, 0})
  {
  }

  virtual ~VertexParameters()
  {
  }


  /**
   * Set the value of normal_
   * @param value the new value of normal_
   */
  void SetNormal(const Vertex<T>& value)
  {
    nml_ = value;
    type_ = VertexParamType::kNormal;
  }

  /**
   * Get the value of normal_
   * @return the value of normal_
   */
  Vertex<T> GetNormal() const
  {
    return nml_;
  }

  /**
   * reverse normal
   */
  void FlipNormal()
  {
    nml_ = -nml_;
  }

  void AddTexture(const Vertex2<T> & value)
  {
    tex_ = value;
    type_ = VertexParamType::kTexture;
  }

  Vertex2<T> GetTexture()
  {
      return tex_;
  }

private:
  VertexParamType type_;
  Vertex2<T> tex_;
  Vertex<T> nml_ ;
};
