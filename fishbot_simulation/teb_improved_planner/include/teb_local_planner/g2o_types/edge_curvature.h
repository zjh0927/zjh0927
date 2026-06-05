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

#ifndef EDGE_CURVATURE_H
#define EDGE_CURVATURE_H

#include <teb_local_planner/g2o_types/vertex_pose.h>
#include <teb_local_planner/g2o_types/vertex_timediff.h>
#include <teb_local_planner/g2o_types/base_teb_edges.h>
#include <teb_local_planner/teb_config.h>
#include <cmath>

namespace teb_improved_planner
{

/**
 * @class EdgeCurvature
 * @brief Edge defining the cost function for curvature constraint (penalizing excessive curvature).
 * 
 * The edge depends on three vertices \f$ \mathbf{s}_{i-1}, \mathbf{s}_i, \mathbf{s}_{i+1} \f$ and minimizes: \n
 * \f$ \min \max(0, |\kappa_i| - \kappa_{max})^2 \cdot weight \f$. \n
 * where \f$ \kappa_i \f$ is the curvature at point i, computed from three consecutive points, and 
 * \f$ \kappa_{max} = 1 / R_{min} \f$ is the maximum allowed curvature (inverse of minimum turning radius). \n
 * The dimension of the error / cost vector is 1.
 * @see TebOptimalPlanner::AddEdgesCurvature
 * @remarks Do not forget to call setTebConfig()
 */    
class EdgeCurvature : public BaseTebMultiEdge<1, double>
{
public:

  /**
   * @brief Construct edge.
   */    
  EdgeCurvature()
  {
    this->resize(3);
  }
    
  /**
   * @brief Set the TEB configuration parameters
   * @param cfg TebConfig class
   */
  void setTebConfig(const TebConfig& cfg)
  {
    cfg_ = &cfg;
    max_curvature_ = 1.0 / cfg_->robot.min_turning_radius; // κ_max = 1 / R_min
  }
    
  /**
   * @brief Actual cost function
   */   
  void computeError()
  {
    const VertexPose* pose_im1 = static_cast<const VertexPose*>(_vertices[0]);
    const VertexPose* pose_i   = static_cast<const VertexPose*>(_vertices[1]);
    const VertexPose* pose_ip1 = static_cast<const VertexPose*>(_vertices[2]);

    // 获取三个连续点的坐标
    double x_im1 = pose_im1->x();
    double y_im1 = pose_im1->y();
    double x_i   = pose_i->x();
    double y_i   = pose_i->y();
    double x_ip1 = pose_ip1->x();
    double y_ip1 = pose_ip1->y();
    
    // 计算向量 v1 = P_i - P_{i-1} 和 v2 = P_{i+1} - P_i
    double dx1 = x_i - x_im1;
    double dy1 = y_i - y_im1;
    double dx2 = x_ip1 - x_i;
    double dy2 = y_ip1 - y_i;
    
    // 计算向量长度
    double norm1 = sqrt(dx1*dx1 + dy1*dy1);
    double norm2 = sqrt(dx2*dx2 + dy2*dy2);
    
    // 如果任一向量长度为零，曲率为0（不会产生惩罚）
    if (norm1 < 1e-6 || norm2 < 1e-6)
    {
      _error[0] = 0.0;
      return;
    }
    
    // 计算单位向量
    double ux1 = dx1 / norm1;
    double uy1 = dy1 / norm1;
    double ux2 = dx2 / norm2;
    double uy2 = dy2 / norm2;
    
    // 计算转角 Δθ（通过叉积计算sin(Δθ)，通过点积计算cos(Δθ)）
    double sin_dtheta = ux1 * uy2 - uy1 * ux2;  // cross product (z-component)
    double cos_dtheta = ux1 * ux2 + uy1 * uy2;  // dot product
    
    // 计算转角 Δθ（带符号，表示转弯方向）
    double dtheta = atan2(sin_dtheta, cos_dtheta);
    
    // 计算弦长（连续三点形成的三角形的底边）
    double dx_chord = x_ip1 - x_im1;
    double dy_chord = y_ip1 - y_im1;
    double chord_length = sqrt(dx_chord*dx_chord + dy_chord*dy_chord);
    
    // 计算曲率 κ = Δθ / 弦长
    double curvature = 0.0;
    if (chord_length > 1e-6)
    {
      curvature = fabs(dtheta) / chord_length;
    }
    
    // 计算误差：如果曲率超过最大允许值，则产生惩罚
    // 使用软约束：max(0, |κ| - κ_max)
    if (curvature > max_curvature_)
    {
      _error[0] = curvature - max_curvature_;
    }
    else
    {
      _error[0] = 0.0;
    }

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
      
private:
  double max_curvature_;  // 最大允许曲率 κ_max = 1 / R_min
      
public: 
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

} // end namespace

#endif /* EDGE_CURVATURE_H */