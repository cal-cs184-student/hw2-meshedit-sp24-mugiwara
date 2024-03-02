#include "student_code.h"
#include "mutablePriorityQueue.h"

using namespace std;

namespace CGL
{

  /**
   * Evaluates one step of the de Casteljau's algorithm using the given points and
   * the scalar parameter t (class member).
   *
   * @param points A vector of points in 2D
   * @return A vector containing intermediate points or the final interpolated vector
   */
  std::vector<Vector2D> BezierCurve::evaluateStep(std::vector<Vector2D> const &points)
  { 
    // TODO Part 1.
      std::vector<Vector2D> nextPoints;
      size_t numPoints = points.size();

      if (numPoints <= 1) {
          return points;
      }

      for (size_t i = 0; i < numPoints - 1; ++i) {
          Vector2D lerpPoint = points[i] * (1 - t) + points[i + 1] * t;
          nextPoints.push_back(lerpPoint);
      }
      return nextPoints;
  }

  /**
   * Evaluates one step of the de Casteljau's algorithm using the given points and
   * the scalar parameter t (function parameter).
   *
   * @param points    A vector of points in 3D
   * @param t         Scalar interpolation parameter
   * @return A vector containing intermediate points or the final interpolated vector
   */
  std::vector<Vector3D> BezierPatch::evaluateStep(std::vector<Vector3D> const &points, double t) const
  {
    // TODO Part 2.
      std::vector<Vector3D> nextPoints;
      size_t numPoints = points.size();

      if (numPoints <= 1) {
          return points;
      }
      // LERP-ing the list of points
      for (int i = 0; i < numPoints - 1; i++) {
          Vector3D lerpPoint = (1.0 - t) * points[i] + t * points[i + 1];
          nextPoints.push_back(lerpPoint);
      }
      return nextPoints;
  }

  /**
   * Fully evaluates de Casteljau's algorithm for a vector of points at scalar parameter t
   *
   * @param points    A vector of points in 3D
   * @param t         Scalar interpolation parameter
   * @return Final interpolated vector
   */
  Vector3D BezierPatch::evaluate1D(std::vector<Vector3D> const &points, double t) const
  {
    // TODO Part 2.
      if (points.size() == 1) {
          return points[0];
      }
      return evaluate1D(evaluateStep(points, t), t);
  }

  /**
   * Evaluates the Bezier patch at parameter (u, v)
   *
   * @param u         Scalar interpolation parameter
   * @param v         Scalar interpolation parameter (along the other axis)
   * @return Final interpolated vector
   */
  Vector3D BezierPatch::evaluate(double u, double v) const 
  {  
    // TODO Part 2.
      std::vector<Vector3D> intermediatePoints;
      int rows = controlPoints.size();

      for (int row = 0; row < rows; ++row) {
          intermediatePoints.push_back(evaluate1D(controlPoints[row], u));
      }
      return evaluate1D(intermediatePoints, v);
  }

  Vector3D Vertex::normal( void ) const
  {
    // TODO Part 3.
    // Returns an approximate unit normal at this vertex, computed by
    // taking the area-weighted average of the normals of neighboring
    // triangles, then normalizing.
      Vector3D weighted_normal(0, 0, 0);
      HalfedgeCIter h = halfedge();

      if (h->face()->isBoundary()) {
          do {
              h = h->twin()->next();
          } while (h != halfedge() && h->face()->isBoundary());
      }

      HalfedgeCIter h_orig = h;

      do {
          if (!h->face()->isBoundary()) {
              Vector3D v0 = position;
              Vector3D v1 = h->next()->vertex()->position;
              Vector3D v2 = h->next()->next()->vertex()->position;

              Vector3D face_normal = cross(v1 - v0, v2 - v0);
              double face_area = face_normal.norm() / 2.0;

              face_normal = face_normal.unit();
              weighted_normal += face_area * face_normal;
          }
          h = h->twin()->next();
      } while (h != h_orig);

      return weighted_normal.unit();
  }

  EdgeIter HalfedgeMesh::flipEdge( EdgeIter e0 )
  {
    // TODO Part 4.
    // This method should flip the given edge and return an iterator to the flipped edge.
      if (!e0->isBoundary()) {
          // Define halfedges, vertices, edges, and faces involved in the flip.
          HalfedgeIter h0 = e0->halfedge(), h1 = h0->next(), h2 = h1->next(),
                  h3 = h0->twin(), h4 = h3->next(), h5 = h4->next(),
                  h6 = h1->twin(), h7 = h2->twin(), h8 = h4->twin(), h9 = h5->twin();
          VertexIter v0 = h0->vertex(), v1 = h3->vertex(), v2 = h2->vertex(), v3 = h5->vertex();
          EdgeIter e1 = h1->edge(), e2 = h2->edge(), e3 = h4->edge(), e4 = h5->edge();
          FaceIter f1 = h0->face(), f2 = h3->face();

          // Update neighbors for all halfedges.
          h0->setNeighbors(h1, h3, v2, e0, f1);
          h1->setNeighbors(h2, h9, v3, e4, f1);
          h2->setNeighbors(h0, h6, v1, e1, f1);
          h3->setNeighbors(h4, h0, v3, e0, f2);
          h4->setNeighbors(h5, h7, v2, e2, f2);
          h5->setNeighbors(h3, h8, v0, e3, f2);
          h6->setNeighbors(h6->next(), h2, v2, e1, h6->face());
          h7->setNeighbors(h7->next(), h4, v0, e2, h7->face());
          h8->setNeighbors(h8->next(), h5, v3, e3, h8->face());
          h9->setNeighbors(h9->next(), h1, v1, e4, h9->face());

          // Update vertices, edges, and faces to reflect the flip.
          v0->halfedge() = h5;
          v1->halfedge() = h2;
          v2->halfedge() = h4;
          v3->halfedge() = h1;
          e0->halfedge() = h0;
          e1->halfedge() = h2;
          e2->halfedge() = h4;
          e3->halfedge() = h5;
          e4->halfedge() = h1;
          f1->halfedge() = h0;
          f2->halfedge() = h3;
      }
      return e0;
  }

  VertexIter HalfedgeMesh::splitEdge( EdgeIter e0 )
  {
    // TODO Part 5.
    // This method should split the given edge and return an iterator to the newly inserted vertex.
    // The halfedge of this vertex should point along the edge that was split, rather than the new edges.
      if(e0->isBoundary()) {
          return e0->halfedge()->vertex();
      }

      // Retrieve initial halfedges and faces.
      HalfedgeIter h0 = e0->halfedge(), h1 = h0->next(), h2 = h1->next(),
              h3 = h0->twin(), h4 = h3->next(), h5 = h4->next();
      FaceIter f0 = h0->face(), f1 = h3->face();

      // Retrieve original vertices.
      VertexIter v0 = h0->vertex(), v1 = h3->vertex(),
              v2 = h2->vertex(), v3 = h5->vertex();

      // Create new elements: vertex, edges, halfedges, and faces.
      VertexIter new_vertex = newVertex();
      EdgeIter e5 = newEdge(), e6 = newEdge(), e7 = newEdge();
      HalfedgeIter h6 = newHalfedge(), h7 = newHalfedge(), h8 = newHalfedge(),
              h9 = newHalfedge(), h10 = newHalfedge(), h11 = newHalfedge();
      FaceIter f2 = newFace(), f3 = newFace();

      // Compute the position of the new vertex and assign its halfedge.
      new_vertex->position = (v0->position + v1->position) / 2;
      new_vertex->halfedge() = h0;

      // Update connectivity of new and existing elements.
      h0->setNeighbors(h1, h3, new_vertex, e0, f0);
      h1->setNeighbors(h8, h1->twin(), v1, h1->edge(), f0);
      h2->setNeighbors(h6, h2->twin(), v2, h2->edge(), f2);
      h3->setNeighbors(h11, h0, v1, e0, f1);
      h4->setNeighbors(h10, h4->twin(), v0, h4->edge(), f3);
      h5->setNeighbors(h3, h5->twin(), v3, h5->edge(), f1);
      h6->setNeighbors(h7, h9, v0, e6, f2);
      h7->setNeighbors(h2, h8, new_vertex, e5, f2);
      h8->setNeighbors(h0, h7, v2, e5, f0);
      h9->setNeighbors(h4, h6, new_vertex, e6, f3);
      h10->setNeighbors(h9, h11, v3, e7, f3);
      h11->setNeighbors(h5, h10, new_vertex, e7, f1);

      // Assign halfedges to new edges.
      e5->halfedge() = h8;
      e6->halfedge() = h6;
      e7->halfedge() = h10;

      // Mark the new and existing edges accordingly.
      e0->isNew = false;
      e5->isNew = true;
      e6->isNew = true;
      e7->isNew = true;

      // Update faces to point to one of their halfedges.
      f0->halfedge() = h0;
      f1->halfedge() = h3;
      f2->halfedge() = h2;
      f3->halfedge() = h4;

      // Return the newly inserted vertex.
      return new_vertex;
  }



  void MeshResampler::upsample( HalfedgeMesh& mesh ) {
      // TODO Part 6.
      // This routine should increase the number of triangles in the mesh using Loop subdivision.
      // One possible solution is to break up the method as listed below.

      // 1. Compute new positions for all the vertices in the input mesh, using the Loop subdivision rule,
      // and store them in Vertex::newPosition. At this point, we also want to mark each vertex as being
      // a vertex of the original mesh.
      for (VertexIter v = mesh.verticesBegin(); v != mesh.verticesEnd(); ++v) {
          v->isNew = false;
          float n = static_cast<float>(v->degree());
          float u = (n == 3.0f) ? 3.0f / 16.0f : 3.0f / (8.0f * n);

          Vector3D sum(0.0, 0.0, 0.0);
          HalfedgeIter h = v->halfedge();
          HalfedgeIter h_orig = h;

          do {
              VertexIter neighbor = h->twin()->vertex();
              sum += neighbor->position;
              h = h->twin()->next();
          } while (h != h_orig);

          v->newPosition = (1.0f - n * u) * v->position + (u * sum);
      }

      // 2. Compute the updated vertex positions associated with edges, and store it in Edge::newPosition.
      for (EdgeIter e = mesh.edgesBegin(); e != mesh.edgesEnd(); ++e) {

          HalfedgeIter h = e->halfedge();
          HalfedgeIter hTwin = h->twin();

          VertexIter v0 = h->vertex();
          VertexIter v2 = hTwin->vertex();

          VertexIter v1 = h->next()->next()->vertex();
          VertexIter v3 = hTwin->next()->next()->vertex();

          e->newPosition = (3.0 / 8.0) * (v0->position + v2->position) +
                           (1.0 / 8.0) * (v1->position + v3->position);

          e->isNew = false;
      }

    // 3. Split every edge in the mesh, in any order. For future reference, we're also going to store some
    // information about which subdivide edges come from splitting an edge in the original mesh, and which edges
    // are new, by setting the flat Edge::isNew. Note that in this loop, we only want to iterate over edges of
    // the original mesh---otherwise, we'll end up splitting edges that we just split (and the loop will never end!)
      std::vector<EdgeIter> edgesToSplit;
      for (EdgeIter e = mesh.edgesBegin(); e != mesh.edgesEnd(); ++e) {
          edgesToSplit.push_back(e);
      }

      for (EdgeIter& e : edgesToSplit) {
          VertexIter v = mesh.splitEdge(e);
          v->isNew = true;
          v->newPosition = e->newPosition;

          HalfedgeIter h = v->halfedge();
          do {
              EdgeIter edge = h->edge();
              edge->isNew = !edge->isNew;
              h = h->twin()->next();
          } while (h != v->halfedge());


          h->edge()->isNew = false;
          h->twin()->next()->twin()->next()->edge()->isNew = false;
          h->twin()->next()->edge()->isNew = true;
          h->next()->next()->edge()->isNew = true;
      }

    // 4. Flip any new edge that connects an old and new vertex.
      for (EdgeIter e = mesh.edgesBegin(); e != mesh.edgesEnd(); ++e) {
          if (e->isNew) {
              HalfedgeIter h = e->halfedge();
              VertexIter v1 = h->vertex();
              VertexIter v2 = h->twin()->vertex();

              if ((v1->isNew && !v2->isNew) || (!v1->isNew && v2->isNew)) {
                  mesh.flipEdge(e);
              }
          }
      }

      // 5. Copy the new vertex positions into final Vertex::position.
      for (VertexIter v = mesh.verticesBegin(); v != mesh.verticesEnd(); v++) {
          v->position = v->newPosition;
      }
  }
}
