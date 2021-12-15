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

#include "VertexParameters.h"
#include "FaceParameters.h"
#include "Point.h"
#include "Face.h"

template<class T>
class MeshParameters
{

public:

  MeshParameters()
  {
  }

  /**
   * Empty Destructor
   */
  virtual ~MeshParameters()
  {
  }

  /**
  * Set the value of vertex_paramters_
  * @param value the new value of vertex_paramters_
  */
  void SetVertexParamters(const std::vector<VertexParameters<T> >& value)
  {
    vertex_parameters_ = value;
  }

  void SetVertexParamter(const VertexParameters<T> & value)
  {
      vertex_parameters_.push_back(value);
  }

  /**
   * Get the value of vertex_paramters_
   * @return the value of vertex_paramters_
   */
  std::vector<VertexParameters<T> >& GetVertexParamters()
  {
    return vertex_parameters_;
  }

  VertexParameters<T>& GetVertexParamter(int index)
  {
      return vertex_parameters_[index];
  }
  /**
   * Get the value of vertex_paramters_
   * @return the value of vertex_paramters_
   */
  Vertex<T> GetVertexNormal(unsigned int index) const
  {
    if (index < 0 || index >= vertex_parameters_.size()) 
    {
      throw "vertex normals are not estimated!";
    }

    return vertex_parameters_[index].GetNormal();
  }



  /**
   * Set the value of face_parameters_
   * @param value the new value of face_parameters_
   */
  void SetFaceParameters(const std::vector<FaceParameters<T> >& value)
  {
    face_parameters_ = value;
  }

  /**
   * Get the value of face_parameters_
   * @return the value of face_parameters_
   */
  std::vector<FaceParameters<T> >& GetFaceParameters()
  {
    return face_parameters_;
  }


  /**
  * This function estimates face normals for all mesh faces
  */
  void CalculateFaceNormals(std::vector<shared_ptr<Face<T>>>& faces)
  {

    if (face_parameters_.size() != faces.size()) 
    {
        const auto size = faces.size();
      face_parameters_.resize(size);
      for (auto i = 0; i < size; ++i) {
        T area = 0;
        auto& verts = faces[i]->GetVertices();
        Vertex<T> vec_prev = *verts[1] - *verts[0];
        for (int j = 2; j < verts.size(); ++j) {
          Vertex<T> vec_n = *verts[j] - *verts[0];
          Vertex<T> normal = vec_prev.Cross(vec_n);
          area += normal.GetL2Norm();
          if (j == 2) {
            face_parameters_[i].SetNormal(normal);
          }
        }
       
        face_parameters_[i].SetArea(area / 2);
      }
    }
  }



  /*
  * This function returns the normal of the face given by the face_index
  * @param face_index is the index of the face that we need to get its normal
  * @returns Vertex representing normal of the face
  */
  Vertex<T> GetFaceNormal(int face_index) const
  {
    if (face_index < 0 || face_index >= face_parameters_.size()) 
    {
      throw "Face normals are not estimated!";
    }

    return face_parameters_[face_index].GetNormal();
  }

  /*
  * This function returns the area of the face given by the face_index
  * @param face_index is the index of the face that we need to get its normal
  * @returns T the area of the face
  */
  T GetFaceArea(int face_index) const
  {
    if (face_index < 0 || face_index >= face_parameters_.size()) 
    {
      throw "Face normals are not estimated!";
    }

    return face_parameters_[face_index].GetArea();
  }

  /**
  * This function estiamtes vertex normals as the average of face normals
  */
  void CalculateVertexNormals(const std::vector<shared_ptr<Vertex<T>>>& vertices, const std::vector<shared_ptr<Face<T>>>& faces)
  {
    if (vertices.size() != vertex_parameters_.size()) 
    {
      // make sure face normals are caculated
      CalculateFaceNormals(faces);

      auto nFaces = faces.size();
      auto nVertices = vertices.size();
      vertex_parameters_.clear();
      // init arrays
      vertex_parameters_.resize(nVertices);
      std::vector<short> counts(nVertices, 0);
      std::vector<Vertex<T> > normals(nVertices, Vertex<T>(0, 0, 0));
      // sum connected faces normals
      for (auto i = 0; i < nFaces; ++i) {
        for (auto v : faces[i]->GetVertices()) {
          normals[v->GetIndex()] += face_parameters_[i].GetNormal();
          counts[v->GetIndex()]++;
        }
      }
      //average
      for (auto i = 0; i < nVertices; ++i) {
        normals[i] /= counts[i];
        vertex_parameters_[i].SetNormal(normals[i]);
      }
    }
  }
  
  /*
  * This function clears all interinsic mesh paramaters
  */
  void Clear() 
  {
    this->vertex_parameters_.clear();
    this->face_parameters_.clear();
  }

  /*
  * init with the given number of vertices
  */
  void init(int num_of_verts) 
  {
    vertex_parameters_.resize(num_of_verts);
  }

  size_t GetVertexParametersSize() const
  {
      return vertex_parameters_.size();
  }


  size_t GetFaceParametersSize() const
  {
      return face_parameters_.size();
  }

private: 
  /*
  * vector of vertex interinsic paramters
  */
  std::vector<VertexParameters<T> > vertex_parameters_;
  /*
  * vector of vertex interinsic paramters
  */
  std::vector<FaceParameters<T> > face_parameters_;
};