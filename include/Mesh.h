#pragma once

#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "Point.h"
#include "Face.h"
#include "Edge.h"
#include "MeshParameters.h"
#include "Element.h"

/*
* This enumeration is used to set the unit use in the Mesh MM or CM
*/
enum Unit 
{
  kMM = 0,
  kCM
};

enum class BooleanOperation 
{
  kIntersection,
  kUnion,
  kDifference
};

/**
  * class Mesh
  * This classreprents a the mesh inside the library. It has three arrays: vertices, edges, and faces
  */
template<class T>
class Mesh
{
  typedef std::shared_ptr<Element<T> > ObjectPtr;
  /*
  * type definition to shared pointer of the type Vertex
  */
  typedef std::shared_ptr<Vertex<T> > VertexPtr;
  /*
  * type definition to shared pointer of the type Edge
  */
  typedef std::shared_ptr<Edge<T> > EdgePtr;
  /*
  * type definition to shared pointer of the type Face
  */
  typedef std::shared_ptr<Face<T> > FacePtr;
public:

  Mesh()
  {
  }


  /**
  * Copy constructor for the mesh
   * @param  mesh
   */
  Mesh(const Mesh& mesh)
  {
    const int nV = mesh.GetNumberOfVertices();
    // add vertices
    for (int i = 0; i < nV; ++i) {
      auto v = mesh.GetVertex(i);
      AddVertex(v->x(), v->y(), v->z());
    }
    // add edges
    const int nE = mesh.GetNumberOfEdges();
    for (int i = 0; i < nE; ++i) {
      EdgePtr edge = mesh.GetEdge(i);

      AddEdge(
        vertices_[mesh.GetIndex(edge->GetV1())],
        vertices_[mesh.GetIndex(edge->GetV2())]);
    }

    // Add faces
    const int nF = mesh.GetNumberOfFaces();
    for (int i = 0; i < nF; ++i) {
      FacePtr face = mesh.GetFace(i);
      const std::vector<VertexPtr>& old_verts = face->GetVertices();
      std::vector<int> old_verts_indeces;
      for (const auto& v : old_verts) {
        old_verts_indeces.push_back(v->GetIndex());
      }
      AddFace(old_verts_indeces);
    }
  }

  /**
   * Empty Destructor
   */
  virtual ~Mesh()
  {
  }

  /**
   * @return number of faces
   */
  unsigned int GetNumberOfFaces() const
  {
    return faces_.size();
  }

  /**
   * @return number of vertices
   */
  unsigned int GetNumberOfVertices() const
  {
    return vertices_.size();
  }

  /**
   * @return number of edges
   */
  unsigned int GetNumberOfEdges() const
  {
    return edges_.size();
  }


  /**
  * Compare two meshes to have same vertices and faces 
  * by comparing the values not pointers the same objects
  * @return bool
  * @param  mesh_b
  */
  bool operator==(const Mesh& mesh_b) const
  {
    const auto nVertices = GetNumberOfVertices();
    const auto nFaces = GetNumberOfFaces();
    if (nVertices != mesh_b.GetNumberOfVertices()) {
      return false;
    }
    if (GetNumberOfEdges() != mesh_b.GetNumberOfEdges()) {
      return false;
    }
    if (nFaces != mesh_b.GetNumberOfFaces()) {
      return false;
    }
    // check vertices
    for (int i = 0;i <  GetNumberOfVertices(); ++i) {
      if (!(*GetVertex(i) == *mesh_b.GetVertex(i))) {
        return false;
      }
    }
    // check faces
    for (auto i = 0; i < nFaces; ++i) {
      if (!(*GetFace(i) == *mesh_b.GetFace(i))) {
        return false;
      }
    }
    return true;
  }



  /*
  * Get edges connected to a specific face
  */
  void GetFaceEdges(const FacePtr& face, std::unordered_set<EdgePtr>& edges) const {
    auto& verts = face->GetVertices();
    const auto n = verts.size();
    for (auto i = 0; i < n; ++i) {
      edges.insert(FindEdge(verts[i], verts[(i + 1) % n]));
    }
  }
  /**
 * This function adds a vertex to the mesh
  * @param  x x-coordinate of the new Vertex
  * @param  y y-coordinate of the new Vertex
  * @param  z z-coordinate of the new Vertex
  * @return the index of this vertex
  */
  int AddVertex(T x, T y, T z)
  {
      Vertex<T> vertex({ x, y, x });
      return AddVertex(vertex);
  }
 
  /**
 * This function adds a vertex to the mesh
  * @param  vertex object to be added of type Vertex
  * @return the index of this vertex
  */
  int AddVertex(Vertex<T> vertex)
  {
      int vertexIndex = vertices_.size();
      auto v = std::make_shared<Vertex<T> >(vertex);
      v->SetIndex(vertexIndex);
      vertices_.push_back(v);
      //index_[v] = vertexIndex;
      return vertexIndex;
  }

  /**
   * @return the edge index (if a new edge is added it will return its index,
   * otherwise the edge that was found will be returned). In case the edge cannot be created nor found return -1;
   * @param  v1 the first vertex index in the edge
   * @param  v2 the second vertex index in the edge
   */
  int AddEdge(int v1_index, int v2_index)
  {
    if (v1_index < 0 || v1_index >= vertices_.size() || v2_index < 0 || v2_index >= vertices_.size()) {
      return -1;
    }
    VertexPtr vp1 = vertices_[v1_index];
    VertexPtr vp2 = vertices_[v2_index];
    return AddEdge(vp1, vp2);
  }


  /**
   * The function takes the vertex indeces in the face order and creates/finds the edges
   * and add them when creating the face
   * @param  vertices
   * @return the function return the face index or -1 in case of the face could not be created
   */
  int AddFace(const std::vector<int>& vertices)
  {
    const auto numVertices = vertices.size();
    std::vector<VertexPtr> vertPtrs(numVertices);
    for (auto i = 0; i < numVertices; ++i) {
      vertPtrs[i] = GetVertex(vertices[i]);
    }
    std::cout << "added face of size: " << vertices.size() << std::endl;
    return AddFace(vertPtrs);
  }

  /**
   * // The function takes the vertex pointers in the face order and creates/finds the edges
   * and add them when creating the face
   * @param  vertices
   * @return the function return the face index or -1 in case of the face could not be created
   */
  int AddFace(const std::vector<VertexPtr>& vertices)
  {

    const auto numVertices = vertices.size();
    std::vector<EdgePtr> edges;
    for (auto i = 0; i < numVertices; ++i) {
      int j = (i + 1) % numVertices;
      // find/create edge
      auto edgeIndex = AddEdge(GetIndex(vertices[i]), GetIndex(vertices[j]));
      if (edgeIndex > -1) {
        edges.push_back(edges_[edgeIndex]);
      }
      else {
        return -1;
      }
    }
    return AddFace(vertices, edges);
  }

  /**
  * This function returns a vertex pointer by index
   * @param  index the index of the vertex, no check is made for validitiy of the index
   */
  VertexPtr GetVertex(int index) const
  {
    return vertices_[index];
  }

  /**
 * This function set a vertex by index
  * @param  index the index of the vertex, no check is made for validitiy of the index
  * @param  vertex : the pointer to the new vertex
  */
  void SetVertex(int index, VertexPtr vertex)
  {
	  vertices_[index]->SetX(vertex->x());
    vertices_[index]->SetY(vertex->y());
    vertices_[index]->SetZ(vertex->z());
  }

  /**
* This function returns an edge pointer by index
 * @param  index the index of the edge, no check is made for validitiy of the index
 */
  EdgePtr GetEdge(int index) const
  {
    return edges_[index];
  }


  /**
* This function returns a face pointer by index
 * @param  index the index of the face, no check is made for validitiy of the index
 */
  FacePtr GetFace(int index) const
  {
    return faces_[index];
  }

  /*
  * Get the index of the passed shared_ptr in the corresponding array (vertices_, edges_, faces_)
  * @param obj_ptr is a base pointer for all object types
  * @return the index of the passed pointer if exists, otherwise returns -1
  */
 
  inline int GetIndex(const ObjectPtr& obj_ptr) const
  {  
    if (obj_ptr != NULL) {
      return obj_ptr->GetIndex();
    }
    return -1;
  }

  /**
  * Delete vertex by the given index
   * @param  index
   * @param update is a flag to update the vertices array
   */
  void DeleteVertex(int index, bool update = true)
  {
    auto& vertex = vertices_[index];
    DeleteVertex(vertex, update);
  }

  /*
  * Delete vertex by the given pointer
  * @param vertex is the vertex to be deleted
  * @param update is a flag to update vertices_ array
  */
  void DeleteVertex(VertexPtr vertex, bool update = true)
  {
    if (vertex != NULL) {
      const auto& entry = vertex_edge_.find(vertex);
      if (entry != vertex_edge_.end()) {
        auto edges = entry->second;
        DeleteEdges(edges, false);
      }
      vertex_edge_.erase(vertex);
      int index = GetIndex(vertex);
      if (index > -1 && vertices_[index]) {
        vertices_[index].reset();
      }
      vertex->SetIndex(-1);
    }
    deleted_vertex_ = true;
    if (update) {
      OptimizeArray(vertices_, deleted_vertex_);
    }
  }

  /*
  * This function deletes set of vertices sent by indeces
  * @param vertex_indeces is a vector of indeces to be deleted
  */
  void DeleteVertices(const std::vector<int>& vertex_indeces) {
    for (auto v : vertex_indeces) {
      DeleteVertex(v, false);
    }
    OptimizeArray(vertices_, deleted_vertex_);
  }

  /*
  * This function deletes set of vertices sent by pointers
  * @param vertices is a vector of vertices to be deleted
  */
  void DeleteVertices(const std::vector<VertexPtr>& vertices) {
    for (auto v : vertices) {
      DeleteVertex(v, false);
    }
    OptimizeArray(vertices_, deleted_vertex_);
  }

  /**
  * Delete edge by the given index
  * @param  index
  * @param update is a falg to update edges array
  * @param delete_vertex_connection do we need to delete the connection saved from vertex to edge
  */
  void DeleteEdge(int index, bool update = true, bool delete_vertex_connection = true)
  {
    auto& edge = edges_[index];
    DeleteEdge(edge, update, delete_vertex_connection);
  }

  /**
  * Delete edge by the given Edge pointer
  * @param  edge is the edge to be deleted
  * @param update is a flag to update the edges array
  * @param delete_vertex_connection do we need to delete the connection saved from vertex to edge
  */
  void DeleteEdge(EdgePtr edge, bool update = true, bool delete_vertex_connection = true)
  {
    if (edge != NULL) {
      if (delete_vertex_connection) {
        DeleteVertexEdgeConnection(edge->GetV1(), edge);
        DeleteVertexEdgeConnection(edge->GetV2(), edge);
      }

      auto entry = edge_face_.find(edge);
      if (entry != edge_face_.end()) {
        auto& faces = entry->second;
        DeleteFaces(faces, false);
      }
      edge_face_.erase(edge);
      int index = GetIndex(edge);
      if (index > -1 && edges_[index]) {
        edges_[index].reset();
      }
      edge->SetIndex(-1);
    }
    
    deleted_edge_ = true;
    if (update) {
      OptimizeArray(edges_, deleted_edge_);
    }
  }

  /*
  * This function deletes the connection between the edge and the vertex
  * the function assumes a connection already exists
  * @param v1 is a pointer to the vertex
  * @param edge is a pointer to the edge to be deleted
  */
  void DeleteVertexEdgeConnection(const std::shared_ptr<Vertex<T>>& v1, const std::shared_ptr<Edge<T>>& edge)
  {
    auto end = vertex_edge_[v1].end();
    for (typename std::vector<EdgePtr>::iterator it = vertex_edge_[v1].begin(); it != end; ++it) {
      if (*it == edge) {
        vertex_edge_[v1].erase(it);
        return;
      }
    }
  }

  /*
  * This function deletes set of edges sent by indeces
  * @param edge_indeces is a vector of indeces to be deleted
  * @param delete_vertex_connection do we need to delete the connection saved from vertex to edge
  */
  void DeleteEdges(const std::vector<int>& edge_indeces, bool delete_vertex_connection = true) {
    for (auto e : edge_indeces) {
      DeleteEdge(e, false, delete_vertex_connection);
    }
    OptimizeArray(edges_, deleted_edge_);
  }

  /*
  * This function deletes set of edges sent by pointers
  * @param edges is a vector of edges to be delete
  * @param delete_vertex_connection do we need to delete the connection saved from vertex to edge
  */
  void DeleteEdges(const std::vector<EdgePtr>& edges, bool delete_vertex_connection = true) {
    for (auto e : edges) {
      DeleteEdge(e, false, delete_vertex_connection);
    }
    OptimizeArray(edges_, deleted_edge_);
  }

  /**
  * Delete face by index
  * @param  index of the face to be deleted
  * @param delete_edge_connection do we need to delete the connection saved from edge to face
  */
  void DeleteFace(int index, bool update = true, bool delete_edge_connection = true)
  {
    auto& face = faces_[index];
    DeleteFace(face, update, delete_edge_connection);
  }

  /**
  * Delete face by pointer
  * @param face is the face to be deleted
  * @param update is a flag to update faces_ array
  * @param delete_edge_connection do we need to delete the connection saved from edge to face
  */
  void DeleteFace(FacePtr face_ptr, bool update = true, bool delete_edge_connection = true)
  {
    if (face_ptr != NULL) {
      if (delete_edge_connection) {
        auto& verts = face_ptr->GetVertices();
        const auto n = verts.size();
        for (auto i = 0; i < n; ++i) {
          const EdgePtr& edge = FindEdge(verts[i], verts[(i + 1) % n]);
          DeleteEdgeFaceConnection(edge, face_ptr);
        }
      }
      int index = GetIndex(face_ptr);
      if (index > -1 && faces_[index]) {
        faces_[index].reset();
      }
      face_ptr->SetIndex(-1);
    }
    
    deleted_face_ = true;
    if (update) {
      OptimizeArray(faces_, deleted_face_);
    }
  }

  /*
  * This function deletes the connection saved from edge to face
  * 
  * @param edge is the edge we are looking for its connection
  * @param face_ptr is the face to be deleted
  */
  void DeleteEdgeFaceConnection(const std::shared_ptr<Edge<T>>& edge,const std::shared_ptr<Face<T> > & face_ptr)
  {
    auto end = edge_face_[edge].end();
    for (typename std::vector<FacePtr>::iterator it = edge_face_[edge].begin(); it != end; ++it) {
      if (*it = face_ptr) {
        edge_face_[edge].erase(it);
        break;
      }
    }
  }

  /*
  * This function deletes set of faces sent by indeces
  * @param face_indeces is a vector of indeces to be deleted
  * @param delete_edge_connection do we need to delete the connection saved from edge to face
  */
  void DeleteFaces(const std::vector<int>& face_indeces, bool delete_edge_connection = true) {
    for (auto f : face_indeces) {
      DeleteFace(f, false, delete_edge_connection);
    }
    OptimizeArray(faces_, deleted_face_);
  }

  /*
  * This function deletes set of faces sent by pointers
  * @param faces is a vector of faces to be deleted
  * @param delete_edge_connection do we need to delete the connection saved from edge to face
  */
  void DeleteFaces(const std::vector<FacePtr>& faces, bool delete_edge_connection = true) {
    for (auto f : faces) {
      DeleteFace(f, false, delete_edge_connection);
    }
    OptimizeArray(faces_, deleted_face_);
  }

  /**
  * Clear all model arrays (vertices, faces and edges)
  * It should also clear the interinisic paramters (NOT YET IMPLEMENTED)
  */
  void ClearModel()
  {
    this->edges_.clear();
    this->faces_.clear();
    this->vertices_.clear();
    ClearModelExceptFacesAndVertices();
  }


  /**
  * Clear interinsic paramters only
   */
  void ClearModelExceptFacesAndVertices()
  {
    mesh_paramters_.Clear();
  }

  size_t GetVertexParametersSize() const
  {
      return mesh_paramters_.GetVertexParametersSize();
  }



  size_t GetFaceParametersSize() const
  {
      return mesh_paramters_.GetFaceParametersS();
  }


  /**
   * @return std::vector<std::shared_ptr<Face> >
   * @param  index
   */
  void GetVertexFaces(int index, std::vector<FacePtr>& result)
  {
    GetVertexFaces(GetVertex(index), result);
  }

  /**
   * @return std::vector<std::shared_ptr<Face> >
   * @param  vertex
   */
  void GetVertexFaces(VertexPtr vertex, std::vector<FacePtr>& result)
  {
    // get vertex edges
    const auto& entry = vertex_edge_.find(vertex);
    if (entry != vertex_edge_.end()) {
      const auto& edges = entry->second;
      std::unordered_set<FacePtr> faceSet;
      // iterate over the edges
      for (const auto& e : edges) {
        const auto& faces = edge_face_[e];
        // loop over 
        for (const auto& f : faces) {
          if (faceSet.find(f) == faceSet.end()) {
            // if the face does not exist
            faceSet.insert(f);
            result.push_back(f);
          }
        }
      }
    }
  }


  /**
  * This function transforms the mesh using the given transformation matrix
  * by transforming all the vertices in this mesh
  * @param transformation_matrix
  */
  void Transform(Transformation<T>& transformation)
  {
    for (auto v : vertices_) {
      v->Transform(transformation);
    }
  }

  /**
  * This function estimates face normals for all mesh faces
  */
  void CalculateFaceNormals()
  {
      std::cout << " faces size: " << faces_.size() << std::endl;
    mesh_paramters_.CalculateFaceNormals(faces_);
  }

  /*
  * check if a point is inside
  * assuming a watertight mesh
  *
  */
  bool IsPointInsideWatertight(Vertex<T> p)
  {
      for(int f_idx=0; f_idx < GetNumberOfFaces(); f_idx++)
      {
          Vertex<T> segment = *face.Getv1() - p;
          double convergence = segment.Dot(mesh_paramters_.GetFaceNormal(f_idx));
          if (convergence < 0)
          {
              return false;
          }
      }
      return true;
  }

  void Triangulate()
  {
      vector<shared_ptr<Face<T>>> faces;
      for (int f_idx = 0; f_idx < GetNumberOfFaces(); f_idx++)
      {
          auto face = GetFace(f_idx);
          vector<shared_ptr<Vertex<T>>> vertices = face->GetVertices();
          auto v1 = vertices[0];
          std::cout << " face of size: " << vertices.size() << std::endl;
          for (int i = 1; i < vertices.size()-1; i++)
          {
              auto v2 = vertices[i];
              auto v3 = vertices[i + 1];
              vector<shared_ptr<Vertex<T>>> triangle{ v1,v2,v3 };
              
              auto face = make_shared<Face<T>>(triangle);
              face->SetIndex(f_idx);
              
              faces.push_back(face);
          }
      }
      SetFaces(faces);
      CalculateFaceNormals();
  }
  
  /*
   * Calculate Mesh Volume
   */

  double CalculateWatertightVolume()
  {
      double volume = 0;
      //bool triang_nml_side = true;
      auto AddTetrehedron = [](
          shared_ptr<Vertex<T>> v0,
          shared_ptr<Vertex<T>> v1,
          shared_ptr<Vertex<T>> v2)
      {
          auto v10 = (*v1 - *v0);
          auto v20 = (*v2 - *v0);
          auto nml = v10.Cross(v20);

          
          //triang_nml_side ^= triang_nml_side;
          double vol = v0->Dot(nml);
          return (vol<0)? vol*-1: vol;
      };
      for (int f_idx = 0; f_idx < GetNumberOfFaces(); f_idx++)
      {
          auto face = GetFace(f_idx);
          vector<shared_ptr<Vertex<T>>> vertices = face->GetVertices();
          //triangulation of the polygon mesh
          //std::cout << "triangulation polygon of size:  " << vertices.size() << std::endl;
          for (int i = 1; i < vertices.size()-1; i++)
          {
              auto v1 = vertices[0];
              auto v2 = vertices[i];
              auto v3 = vertices[i + 1];
              double vol = AddTetrehedron(vertices[0], vertices[i], vertices[i + 1]);
              //double vol = std::abs(-1 * v3->x() * v2->y() * v1->z() + v2->x() * v3->y() * v1->z() + v3->x() * v1->y() * v2->z() - 1 * v1->x() * v3->y() * v2->z() - 1 * v2->x() * v1->y() * v3->z() + v1->x() * v2->y() * v3->z());
              //std::cout << " calculating volume: " << vol << std::endl;
              volume += vol;
          }
          //triang_nml_side = true;
      }
      return 1 / 6.0 * volume;
  }


  /*
  * This function returns the normal of the face given by the face_index
  * @param face_index is the index of the face that we need to get its normal
  * @returns Vertex representing normal of the face
  */
  Vertex<T> GetFaceNormal(int face_index) const
  {
    if (face_index < 0 || face_index >= GetNumberOfFaces()) {
      throw "Out of faces index!";
    }
    
    return mesh_paramters_.GetFaceNormal(face_index);
  }

  /**
  * This function estiamtes vertex normals as the average of face normals
  */
  void CalculateVertexNormals()
  {
    mesh_paramters_.CalculateVertexNormals(vertices_, faces_);
  }

  /*
  * Find vertex neigbours
  * @param index is the index of the vertex
  * @param vertex_neighbors are the neigbours of the vertex
  */
  void GetVertexNeighbors(const int index, std::vector<int>& vertex_neighbors) {
    vertex_neighbors.clear();
    std::vector<VertexPtr> nbrs;
    GetVertexNeighbors(vertices_[index], nbrs);
    for (auto& nbr : nbrs) {
      vertex_neighbors.push_back(GetIndex(nbr));
    }
  }

  /*
  * Find vertex neigbours
  * @param vertex is a shared pointer to the vertex
  * @param vertex_neighbors are the neigbours of the vertex
  */
  void GetVertexNeighbors(const VertexPtr vertex, std::vector<VertexPtr>& vertex_neighbors) {
    vertex_neighbors.clear();
    auto edges = vertex_edge_.find(vertex);
    if (edges != vertex_edge_.end()) {
      // there is something here
      for (auto& e : edges->second) {
        vertex_neighbors.push_back(e->GetOtherVertex(vertex));
      }
    }
  }

  double GetArea()
  {
      double area=0;
      for (int i = 0; i < GetNumberOfFaces(); i++)
      {
          area += mesh_paramters_.GetFaceArea(i);
      }
      return area;
  }
  /*
  * Get Vertex faces area sum
  */
  T GetVertexFacesAreaSum(int index) {
    return GetVertexFacesAreaSum(vertices_[index]);
  }

  /*
  * Get Vertex faces area sum
  */
  T GetVertexFacesAreaSum(const VertexPtr& vertex) {
    std::vector<FacePtr> faces;
    GetVertexFaces(vertex, faces);
    T area = 0;
    for (auto& face : faces) {
      area += mesh_paramters_.GetFaceArea(GetIndex(face));
    }

    return area;
  }


  /**
   * Set the value of unit_
   * @param value the new value of unit_
   */
  void SetUnit(Unit value)
  {
    unit_ = value;
  }

  /**
   * Get the value of unit_
   * @return the value of unit_
   */
  Unit GetUnit() const
  {
    return unit_;
  }

  /**
   * Set the value of vertices_
   * @param value the new value of vertices_
   */
  void SetVertices(std::vector<VertexPtr> value)
  {
    vertices_ = value;
  }

  /**
   * Get the value of vertices_
   * @return the value of vertices_
   */
  std::vector<VertexPtr> GetVertices() const
  {
    return vertices_;
  }

  /**
   * Set the value of edges_
   * @param value the new value of edges_
   */
  void SetEdges(std::vector<EdgePtr> value)
  {
    edges_ = value;
  }

  /**
   * Get the value of edges_
   * @return the value of edges_
   */
  std::vector<EdgePtr> GetEdges() const
  {
    return edges_;
  }

  /*
  * This function removes the deleted objects from the list (edges_, faces_, vertices_) in this mesh
  */
  template<class type>
  void OptimizeArray(std::vector<type>& list, bool& deleted_flag)
  {
    int index = 0;
    for (int i = 0; i < list.size(); ++i)
    {
      if (list[i] != NULL) {
        if (i != index) {
          list[i]->SetIndex(index);
          list[index] = list[i];
        }
        ++index;
      }
    }
    // now resize
    list.resize(index);

    deleted_flag = false;
  }

  /**
   * Set the value of faces_
   * @param value the new value of faces_
   */
  void SetFaces(std::vector<FacePtr> value)
  {
    faces_ = value;
  }

  /**
   * Get the value of faces_
   * @return the value of faces_
   */
  std::vector<FacePtr> GetFaces() const
  {
    return faces_;
  }

  /*
  * Get the mean of the vertices
  */
  Vertex<T> CalculateVerticesMean() {
    Vertex<T> mean(0, 0, 0);
    if (vertices_.size()) {
      for (auto a : vertices_) {
        mean += *a;
      }
      mean /= vertices_.size();
    }
    return mean;
  }

  /**
  * Set the value of vertex_paramters_
  * @param value the new value of vertex_paramters_
  */
  void SetVertexParamters(const std::vector<VertexParameters<T> >& value)
  {
    mesh_paramters_.SetVertexParamters(value);
  }

  void AddVertexParameters(const VertexParameters<T> param)
  {
      mesh_paramters_.SetVertexParamter(param);
  }

  /**
   * Get the value of vertex_paramters_
   * @return the value of vertex_paramters_
   */
  std::vector<VertexParameters<T> >& GetVertexParamters()
  {
    return mesh_paramters_.GetVertexParamters();
  }

  VertexParameters<T>& GetVertexParamter(int index)
  {
      return mesh_paramters_.GetVertexParamter(index);
  }

  /**
   * Set the value of face_parameters_
   * @param value the new value of face_parameters_
   */
  void SetFaceParameters(const std::vector<FaceParameters<T> >& value)
  {
    mesh_paramters_.SetFaceParameters(value);
  }

  /**
   * Get the value of face_parameters_
   * @return the value of face_parameters_
   */
  std::vector<FaceParameters<T> >& GetFaceParameters()
  {
    return mesh_paramters_.GetFaceParameters();
  }

  /*
  * Are the two vertices connectes
  * @param vertex_ptr1 pointer to the first vertex
  * @param vertex_ptr2 pointer to the second vertex
  * @return true if they are connected and false otherwise
  */
  bool AreConnected(const VertexPtr& vertex_ptr1, const VertexPtr& vertex_ptr2) const{
    return FindEdge(vertex_ptr1, vertex_ptr2) != NULL;
  }

  /*
  * Are the vertex and the edge connected
  * @param vertex_ptr pointer to the vertex
  * @param edge_ptr pointer to the edge
  * @return true if they are connected and false otherwise
  */
  bool AreConnected(const VertexPtr& vertex_ptr, const EdgePtr& edge_ptr) const {
    return edge_ptr->IsConnectedToVertex(vertex_ptr);
  }


  /*
  * Are the vertex and the face connected
  * @param vertex_ptr pointer to the vertex
  * @param face_ptr pointer to the face
  * @return true if they are connected and false otherwise
  */
  bool AreConnected(const VertexPtr& vertex_ptr, const FacePtr& face_ptr) const {
    return face_ptr->IsConnectedToVertex(vertex_ptr);
  }

  /*
  * Are the two edges connected
  * @param edge_ptr1 pointer to the first edge
  * @param edge_ptr2 pointer to the second edge
  * @return true if they are connected and false otherwise
  */
  bool AreConnected(const EdgePtr& edge_ptr1, const EdgePtr& edge_ptr2) const {
    return (AreConnected(edge_ptr1->GetV1(), edge_ptr2) || AreConnected(edge_ptr1->GetV2(), edge_ptr2));
  }

  /*
  * Are the edge and the face connected
  * @param edge_ptr pointer to the edge
  * @param face_ptr pointer to the face
  * @return true if they are connected and false otherwise
  */
  bool AreConnected(const EdgePtr& edge_ptr, const FacePtr& face_ptr) const {
    auto& verts = face_ptr->GetVertices();
    const auto n = verts.size();
    for (auto i = 0; i < n; ++i) {
      if (edge_ptr == FindEdge(verts[i], verts[(i + 1) % n])) {
        return true;
      }
    }
    return false;
  }

  /*
  * Are the two faces (edge) connected
  * @param face_ptr1 pointer to the first face
  * @param face_ptr2 pointer to the second face
  * @return true if they are connected and false otherwise
  */
  bool AreConnected(const FacePtr& face_ptr1, const FacePtr& face_ptr2) const 
  {
    const auto& vertices1 = face_ptr1->GetVertices();
    const auto& vertices2 = face_ptr2->GetVertices();
    const auto n1 = vertices1.size();
    const auto n2 = vertices2.size();
    for (auto i = 0; i < n1; ++i) {
      for (auto j = 0; j < n2; ++j) {
        if (vertices1[i] == vertices2[j]) {
          // if we have two vertices equals from first and second faces
          if (vertices1[(i + 1) % n1] == vertices2[(j + 1) % n2]) {
            // return true if the next vertex in both faces are equal
            return true;
          }
          else if (vertices1[(i + 1) % n1] == vertices2[(j + n2 - 1) % n2]) {
            // also return true if the next vertex in the first face equals the previous vertex in the second face
            return true;
          }
        }
      }
    }
    // no connection found
    return false;
  }

  /*
   * Find if the edge using the given two vertices exists
   */
  EdgePtr FindEdge(const VertexPtr& v1, const VertexPtr& v2) const
  {
    const auto& entry = vertex_edge_.find(v1);
    if (entry != vertex_edge_.end()) {
      const auto& edges = entry->second;
      for (auto& e : edges) {
        if (e->GetOtherVertex(v1) == v2) {
          return e;
        }
      }
    }
    return NULL;
  }

  /**
  * @return the edge index (if a new edge is added it will return its index,
  * otherwise the edge that was found will be returned). In case the edge cannot be created nor found return -1;
  * @param  v1 the first vertex in the edge
  * @param  v2 the second vertex in the edge
  */
  int AddEdge(VertexPtr& v1, VertexPtr& v2)
  {
    EdgePtr edge = FindEdge(v1, v2);// std::static_pointer_cast<Edge<T>>(v1->FindEdge(v2));
    if (edge != NULL) 
    {
      return GetIndex(edge);
    }
    // the edge is not found
    // create a new edge
    int edgeIndex = edges_.size();
    edge = std::make_shared<Edge<T> >(v1, v2);
    // update index
    edge->SetIndex(edgeIndex);

    edges_.push_back(edge);
    // update maps
    vertex_edge_[v1].push_back(edge);
    vertex_edge_[v2].push_back(edge);

    return edgeIndex;
  }

  /*
  * This function adds a new face given two vectors, the vertices and the edges around this face
  * @param vertices are the vertices pointers in this face
  * @param edges are the edge pointers in this face
  * @return the index of the added face
  */
  int AddFace(const std::vector<VertexPtr>& vertices, const std::vector<EdgePtr>& edges)
  {
    // create the face
    int faceIndex = faces_.size();
    auto face = std::make_shared<Face<T> >(vertices);
    face->SetIndex(faceIndex);
    faces_.push_back(face);
    // update maps
     // should be added only once, since the creation of the face should happen only once
    for (auto e : edges) {
      edge_face_[e].push_back(face);
    }

    return faceIndex;
  }


 
private:
  /*
  * The unit used in this mesh (CM, MM)
  */
  Unit unit_ = kMM;
  /*
  * Array of shared pointers to the vertices
  */
  std::vector<VertexPtr> vertices_;
  /*
  * Array of shared pointers to edges
  */
  std::vector<EdgePtr> edges_;
  /*
  * Array of shared pointers to faces
  */
  std::vector<FacePtr> faces_;
  /*
  * Map from vertices to edges
  */
  std::unordered_map<VertexPtr, std::vector<EdgePtr> > vertex_edge_;
  /*
  * Map from vertices to edges
  */
  std::unordered_map<EdgePtr, std::vector<FacePtr> > edge_face_;
  /*
  *   All interinsic mesh parameters
  */
  MeshParameters<T> mesh_paramters_;
  /*
  * Dirty flags to updated the vertiecs_, edges_, faces_ arrays
  */
  bool deleted_face_ = false;
  bool deleted_edge_ = false;
  bool deleted_vertex_ = false;;

};
