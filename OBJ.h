#ifndef _OBJ_H_
#define _OBJ_H_

#include <Eigen/Core>
#include <Eigen/Geometry>

typedef double ScalarType;
typedef Eigen::Matrix<ScalarType, 3, 1> Vector3;

#include <algorithm>
#include <map>
#include <list>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream> 
#include <iomanip>

class OBJ
{
private:
  class Graph
  {
  public:
    Graph(const OBJ& obj)
      : numVertices_(obj.numVertices())
      , adjLists_(numVertices_)
      , visited_(numVertices_, false)
      , connectedComponents_()
    {
      for (const auto& face : obj.faces)
      {
        for (unsigned int vIdx = 1; vIdx < face.size(); ++vIdx)
        {
          const unsigned int v1 = face[vIdx - 1];
          const unsigned int v2 = face[vIdx];
          adjLists_[v1].emplace_back(v2);
          adjLists_[v2].emplace_back(v1);
        }
      }
    }

    std::vector<std::vector<unsigned int>> getConnectedComponents()
    {
      if (connectedComponents_.empty())
        connectedComponents();

      return connectedComponents_;
    };

  private:
    unsigned int numVertices_;
    std::vector<std::list<unsigned int>> adjLists_;
    std::vector<bool> visited_;
    std::vector<std::vector<unsigned int>> connectedComponents_;

    void dfsUtil(const unsigned int v)
    {
      visited_[v] = true;
      connectedComponents_.back().emplace_back(v);
      for (auto i = adjLists_[v].begin(); i != adjLists_[v].end(); ++i)
      {
        if (!visited_[*i])
          dfsUtil(*i);
      }
    }

    void connectedComponents()
    {
      connectedComponents_.clear();
      connectedComponents_.emplace_back();
      visited_ = std::vector<bool>(numVertices_, false);
      for (unsigned int v = 0; v < numVertices_; ++v)
      {
        if (!visited_[v])
        {
          dfsUtil(v);
          connectedComponents_.emplace_back();
        }
      }
      if (connectedComponents_.back().empty())
        connectedComponents_.pop_back();
    }
  };

public:
  OBJ()
    : vertices()
    , normals()
    , faces()
  {}

  OBJ(const std::vector<Vector3>& v, const std::vector<Vector3>& vn = std::vector<Vector3>(), std::vector<std::vector<unsigned int>>& f = std::vector<std::vector<unsigned int>>())
    : vertices(v)
    , normals(vn)
    , faces(f)
  {}

  OBJ(const std::string& file)
    : OBJ(fromFile(file))
  {}

  void toFile(const std::string& file)
  {
    toFile(*this, file);
  }

  void toVTPFile(const std::string& file)
  {
    toVTPFile(*this, file);
  }

  void slice(const std::vector<unsigned int>& vertices)
  {
    slice(*this, vertices);
  }

  Vector3 centroid(const std::vector<unsigned int> component)
  {
    centroid(*this, component);
  }

  std::vector<OBJ> connectedComponents(const unsigned int sortingDimension = 0)
  {
    return connectedComponents(*this, sortingDimension);
  }

  void fuse(const OBJ& obj)
  {
    *this = fuse(*this, obj);
  }

  static OBJ fromFile(const std::string& file)
  {
    std::ifstream infile(file.c_str());
    std::string line;

    OBJ ret;
    while (std::getline(infile, line))
    {
      if (line.empty())
        continue;
      else if (line[0] == '#')
        continue;

      const std::vector<std::string> tokens = tokenize(line, ' ');
      if (tokens.empty())
        continue;

      const char cmd = tokens[0] == "vn" ? 'n' : tokens[0].c_str()[0];
      switch (cmd)
      {
      case 'v':
        ret.vertices.emplace_back(Vector3(as<ScalarType>(tokens[1]), as<ScalarType>(tokens[2]), as<ScalarType>(tokens[3])));
        break;
      case 'n':
        ret.normals.emplace_back(Vector3(as<ScalarType>(tokens[1]), as<ScalarType>(tokens[2]), as<ScalarType>(tokens[3])));
        break;
      case 'f':
      {
        std::vector<unsigned int> face;
        for (unsigned int t = 1; t < tokens.size(); ++t)
        {
          const std::string token = tokens[t];
          const unsigned int vertex = as<unsigned int>(tokenize(token, '/')[0]) - 1;
          face.emplace_back(vertex);
        }
        ret.faces.emplace_back(face);
        break;
      }
      default:
        break;
      }
    }

    return ret;
  }

  static void toFile(const OBJ& obj, const std::string& file)
  {
    std::ofstream outfile(file.c_str());
    outfile << "# Vertices: " << obj.numVertices() << ", Normals: " << obj.numNormals() << " Faces: " << obj.numFaces() << std::endl;

    outfile << std::fixed << std::setprecision(18);
    for (const auto& v : obj.vertices)
      outfile << "v " << v.x() << " " << v.y() << " " << v.z() << std::endl;

    outfile << std::fixed << std::setprecision(18);
    for (const auto& vn : obj.normals)
      outfile << "vn " << vn.x() << " " << vn.y() << " " << vn.z() << std::endl;

    outfile << std::fixed << std::setprecision(0);
    for (const auto& f : obj.faces)
    {
      outfile << "f";
      for (const auto& v : f)
        outfile << " " << (v + 1) << "//" << (v + 1);
      outfile << std::endl;
    }

    outfile.close();
  }

  static void toVTPFile(const OBJ& obj, const std::string& file)
  {
    std::ofstream outfile(file.c_str());
    outfile << "<?xml version=\"1.0\"?>" << std::endl;
    outfile << "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">" << std::endl;
    outfile << "<PolyData>" << std::endl;
    outfile << "<Piece NumberOfPoints = \"" + std::to_string(obj.numVertices()) + "\" NumberOfVerts = \"0\" NumberOfLines = \"0\" NumberOfStrips = \"0\" NumberOfPolys = \"" + std::to_string(obj.numFaces()) + "\">" << std::endl;
    outfile << "<PointData Normals = \"Normals\">" << std::endl;
    outfile << "<DataArray type = \"Float32\" Name = \"Normals\" NumberOfComponents = \"3\" format = \"ascii\">" << std::endl;
    outfile << std::fixed << std::setprecision(18);
    for (const auto& vn : obj.normals)
      outfile << vn.x() << " " << vn.y() << " " << vn.z() << std::endl;
    outfile << "</DataArray>" << std::endl;
    outfile << "</PointData>" << std::endl;
    outfile << "<Points>" << std::endl;
    outfile << "<DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;
    for (const auto& v : obj.vertices)
      outfile << v.x() << " " << v.y() << " " << v.z() << std::endl;
    outfile << "</DataArray>" << std::endl;
    outfile << "</Points>" << std::endl;
    outfile << "<Polys>" << std::endl;
    outfile << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">" << std::endl;
    outfile << std::fixed << std::setprecision(0);
    for (const auto& f : obj.faces)
    {
      for (const auto& v : f)
        outfile << " " << v;
      outfile << std::endl;
    }
    outfile << "</DataArray>" << std::endl;
    outfile << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">" << std::endl;
    unsigned int offset = 0;
    for (unsigned int f = 0; f < obj.numFaces(); f++)
    {
      offset += obj.faces[f].size();
      outfile << offset << " ";
    }
    outfile << std::endl;
    outfile << "</DataArray>" << std::endl;
    outfile << "</Polys>" << std::endl;
    outfile << "</Piece>" << std::endl;
    outfile << "</PolyData>" << std::endl;
    outfile << "</VTKFile>" << std::endl;

    outfile.close();
  }

  static Vector3 centroid(const OBJ& obj, const std::vector<unsigned int> component)
  {
    Vector3 ret;
    for (const auto& c : component)
      ret += obj.vertices[c];
    return ret / component.size();
  }

  static OBJ slice(const OBJ& obj, const std::vector<unsigned int>& vertices)
  {
    OBJ ret;
    std::map<unsigned int, unsigned int> old2newVertex;
    for (unsigned int v = 0; v < vertices.size(); ++v)
    {
      ret.vertices.emplace_back(obj.vertices[vertices[v]]);
      ret.normals.emplace_back(obj.normals[vertices[v]]);
      old2newVertex[vertices[v]] = v;
    }

    for (const auto& f : obj.faces)
    {
      for (const auto& v : f)
      {
        if (old2newVertex.count(v) > 0)
        {
          ret.faces.emplace_back(f);
          std::transform(ret.faces.back().begin(), ret.faces.back().end(), ret.faces.back().begin(),
            [&](const unsigned int lv) -> unsigned int { return old2newVertex[lv]; }
          );
          break;
        }
      }
    }

    return ret;
  }

  static OBJ fuse(const OBJ& obj1, const OBJ& obj2)
  {
    OBJ ret = obj1;
    for (const auto& v : obj2.vertices)
      ret.vertices.emplace_back(v);
    for (const auto& vn : obj2.normals)
      ret.normals.emplace_back(vn);

    for (const auto& f : obj2.faces)
    {
      ret.faces.emplace_back(f);
      std::transform(ret.faces.back().begin(), ret.faces.back().end(), ret.faces.back().begin(),
        [&](const unsigned int lv) -> unsigned int { return obj1.numVertices() + lv; }
      );
    }

    return ret;
  }

  static std::vector<OBJ> connectedComponents(const OBJ& obj, const unsigned int sortingDimension = 0)
  {
    std::vector<std::vector<unsigned int>> ccs = Graph(obj).getConnectedComponents();
    std::sort(ccs.begin(), ccs.end(),
      [obj, sortingDimension](const std::vector<unsigned int>& a, const std::vector<unsigned int>& b) -> bool
      {
        return centroid(obj, a)[sortingDimension] < centroid(obj, b)[sortingDimension];
      }
    );

    std::vector<OBJ> ret;
    for (const auto& cc : ccs)
      ret.emplace_back(slice(obj, cc));
    return ret;
  }

  unsigned int numVertices() const { return vertices.size(); }
  unsigned int numNormals() const { return normals.size(); }
  unsigned int numFaces() const { return faces.size(); }

  std::vector<Vector3> vertices;
  std::vector<Vector3> normals;
  std::vector<std::vector<unsigned int>> faces;

private:
  template<typename T>
  static const T as(const std::string& s)
  {
    std::stringstream ss(s);
    T ret;
    ss >> ret;
    return ret;
  }

  static std::vector<std::string> tokenize(const std::string& s, const char token = ',')
  {
    std::vector<std::string> ret;
    std::stringstream ss(s);
    std::string tmp;
    while (std::getline(ss, tmp, token))
      ret.push_back(tmp);
    return ret;
  }
};

#endif