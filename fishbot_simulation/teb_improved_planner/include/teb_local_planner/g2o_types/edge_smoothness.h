/*********************************************************************
 *
 * Software License Agreement (BSD License)
 *
 *  Copyright (c) 2016,
 *  TU Dortmund - Institute of Control Theory and Systems Engineering.
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions
 *  are met:
 *
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above
 *     copyright notice, this list of conditions and the following
 *     disclaimer in the documentation and/or other materials provided
 *     with the distribution.
 *   * Neither the name of the institute nor the names of its
 *     contributors may be used to endorse or promote products derived
 *     from this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 *  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 *  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 *  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 *  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 *  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 *  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 * 
 * Notes:
 * The following class is derived from a class defined by the
 * g2o-framework. g2o is licensed under the terms of the BSD License.
 * Refer to the base class source for detailed licensing information.
 *
 * Author: Christoph Rösmann
 *********************************************************************/

#ifndef EDGE_SMOOTHNESS_H
#define EDGE_SMOOTHNESS_H

#include <teb_local_planner/g2o_types/vertex_pose.h>
#include <teb_local_planner/g2o_types/vertex_timediff.h>
#include <teb_local_planner/g2o_types/base_teb_edges.h>
#include <teb_local_planner/teb_config.h>

namespace teb_improved_planner
{

/**
 * @class EdgeSmoothness
 * @brief Edge defining the cost function for trajectory smoothness (penalizing rapid changes in orientation).
 * 
 * The edge depends on three vertices \f$ \mathbf{s}_{i-1}, \mathbf{s}_i, \mathbf{s}_{i+1} \f$ and minimizes: \n
 * \f$ \min (\Delta\theta_{i+1} - \Delta\theta_i)^2 \cdot weight \f$. \n
 * where \f$ \Delta\theta_i = \theta_i - \theta_{i-1} \f$ and \f$ \Delta\theta_{i+1} = \theta_{i+1} - \theta_i \f$. \n
 * \e weight can be set using setInformation(). \n
 * The dimension of the error / cost vector is 1.
 * @see TebOptimalPlanner::AddEdgesSmoothness
 * @remarks Do not forget to call setTebConfig()
 */    
class EdgeSmoothness : public BaseTebMultiEdge<1, double>
{
public:

  /**
   * @brief Construct edge.
   */    
  EdgeSmoothness()
  {
    this->resize(3);
  }
    
  /**
   * @brief Actual cost function
   */   
  void computeError()
  {
    const VertexPose* pose_im1 = static_cast<const VertexPose*>(_vertices[0]);
    const VertexPose* pose_i   = static_cast<const VertexPose*>(_vertices[1]);
    const VertexPose* pose_ip1 = static_cast<const VertexPose*>(_vertices[2]);

    // Calculate first differences: Δθ_i = θ_i - θ_{i-1}
    double delta_theta_i = g2o::normalize_theta(pose_i->theta() - pose_im1->theta());
    // Calculate first differences: Δθ_{i+1} = θ_{i+1} - θ_i
    double delta_theta_ip1 = g2o::normalize_theta(pose_ip1->theta() - pose_i->theta());
    
    // Second difference: Δθ_{i+1} - Δθ_i
    _error[0] = delta_theta_ip1 - delta_theta_i;
    
    if (!std::isfinite(_error[0])) {
      _error[0] = 0.0;  // fallback: skip penalty on NaN/inf
    }
  }
  
  /**
   * @brief Helper function to set the vertices of this edge
   */
  void setVertices(VertexPose* pose_im1, VertexPose* pose_i, VertexPose* pose_ip1)
  {
    _vertices[0] = pose_im1;
    _vertices[1] = pose_i;
    _vertices[2] = pose_ip1;
  }
      
public: 
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

} // end namespace

#endif /* EDGE_SMOOTHNESS_H */