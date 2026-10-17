#pragma once
#include <vector>
#include <memory>
#include <chrono>
#include <cmath>
#include <algorithm>
#include <iostream>

#if __has_include(<Eigen/Sparse>)
#include <Eigen/Sparse>
#include <Eigen/Dense>
#elif __has_include(<eigen3/Eigen/Sparse>)
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/Dense>
#else
#error "Eigen3 library not found! Please ensure Eigen headers are in the include path."
#endif

#include "Particle.h"
#include "TrussElement.h"

namespace PD {

struct SolverStats {
    int iterations = 0;
    double residualForce = 0.0;
    double maxDisplacement = 0.0;
    double elapsedSimTime = 0.0;
    bool converged = false;
};

//! Projective Dynamics (PD) 空间杆系结构隐式交替优化求解器
//! 理论参考: Bouaziz et al., "Projective Dynamics: Fusing Constraint Projections for Fast Simulation", ACM TOG (SIGGRAPH 2014)
class PDSolver {
public:
    std::vector<Particle> particles;
    std::vector<std::unique_ptr<Element>> elements;

    // 求解控制参数
    int maxIterations = 500;            //!< 最大交替优化迭代次数
    double tolDisplacement = 1e-5;      //!< 节点位移步进收敛容差 (m, 0.01 mm)
    double chebyshevRho = 0.96;         //!< 切比雪夫/动量加速谱半径
    bool enableChebyshev = true;        //!< 是否开启动量外推加速
    int printInterval = 50;             //!< 日志输出间隔

    PDSolver() = default;

    int addParticle(const Vector3& pos0) {
        int id = static_cast<int>(particles.size());
        particles.emplace_back(id, pos0);
        return id;
    }

    int addTrussElement(int n1, int n2, double E, double A, double N0 = 0.0, bool isCable = false) {
        int id = static_cast<int>(elements.size());
        elements.push_back(std::make_unique<TrussElement>(id, n1, n2, E, A, N0, isCable));
        return id;
    }

    void setParticleSpring(int pId, const Vector3& kTrans, const Vector3& kRot = {0.0, 0.0, 0.0}, bool compressionOnlyZ = false) {
        if (pId >= 0 && pId < static_cast<int>(particles.size())) {
            particles[pId].setSpring(kTrans, kRot, compressionOnlyZ);
        }
    }

    void initialize() {
        for (auto& elem : elements) {
            elem->initialize(particles);
        }
        for (auto& p : particles) {
            p.pos = p.pos0;
            p.posPrev = p.pos0;
            p.vel.setZero();
            p.acc.setZero();
        }
    }

    //! 核心交替优化求解循环
    SolverStats solve() {
        initialize();

        SolverStats stats;
        const int nNodes = static_cast<int>(particles.size());
        const int nElems = static_cast<int>(elements.size());
        if (nNodes == 0 || nElems == 0) {
            stats.converged = true;
            return stats;
        }

        auto tStart = std::chrono::high_resolution_clock::now();

        // 1. 预先计算各单元的等效刚度权重 w_i = E*A / L0
        std::vector<double> elemWeights(nElems, 0.0);
        std::vector<double> elemLength0(nElems, 0.0);
        std::vector<bool>   elemIsCable(nElems, false);
        std::vector<std::pair<int, int>> elemNodes(nElems);

        for (int i = 0; i < nElems; ++i) {
            auto* truss = dynamic_cast<TrussElement*>(elements[i].get());
            if (truss) {
                elemWeights[i] = (truss->E * truss->A) / truss->L0;
                elemLength0[i] = truss->L0;
                elemIsCable[i] = truss->isCable;
                elemNodes[i]   = {truss->node1, truss->node2};
            }
        }

        // 2. 组装全局常数稀疏矩阵 K_x, K_y, K_z
        // 采用 penalty 施加刚性约束，保证拓扑一致与严格对称正定 (SPD)
        const double kRigidFix = 1e12; // 刚性固定惩罚刚度 (N/m)

        std::vector<Eigen::Triplet<double>> tripletsX, tripletsY, tripletsZ;
        tripletsX.reserve(nElems * 4 + nNodes);
        tripletsY.reserve(nElems * 4 + nNodes);
        tripletsZ.reserve(nElems * 4 + nNodes);

        // 图拉普拉斯部分 (对于 X, Y, Z 完全相同)
        for (int i = 0; i < nElems; ++i) {
            double w = elemWeights[i];
            int n1 = elemNodes[i].first;
            int n2 = elemNodes[i].second;

            tripletsX.emplace_back(n1, n1, w);
            tripletsX.emplace_back(n2, n2, w);
            tripletsX.emplace_back(n1, n2, -w);
            tripletsX.emplace_back(n2, n1, -w);

            tripletsY.emplace_back(n1, n1, w);
            tripletsY.emplace_back(n2, n2, w);
            tripletsY.emplace_back(n1, n2, -w);
            tripletsY.emplace_back(n2, n1, -w);

            tripletsZ.emplace_back(n1, n1, w);
            tripletsZ.emplace_back(n2, n2, w);
            tripletsZ.emplace_back(n1, n2, -w);
            tripletsZ.emplace_back(n2, n1, -w);
        }

        // 支座与外载常数向量
        Eigen::VectorXd b0_x = Eigen::VectorXd::Zero(nNodes);
        Eigen::VectorXd b0_y = Eigen::VectorXd::Zero(nNodes);
        Eigen::VectorXd b0_z = Eigen::VectorXd::Zero(nNodes);

        for (int j = 0; j < nNodes; ++j) {
            const auto& p = particles[j];
            b0_x[j] += p.fExt.x;
            b0_y[j] += p.fExt.y;
            b0_z[j] += p.fExt.z;

            // X 边界
            if (p.fixX) {
                tripletsX.emplace_back(j, j, kRigidFix);
                b0_x[j] += kRigidFix * p.pos0.x;
            } else if (p.springK.x > 0.0) {
                tripletsX.emplace_back(j, j, p.springK.x);
                b0_x[j] += p.springK.x * p.pos0.x;
            }

            // Y 边界
            if (p.fixY) {
                tripletsY.emplace_back(j, j, kRigidFix);
                b0_y[j] += kRigidFix * p.pos0.y;
            } else if (p.springK.y > 0.0) {
                tripletsY.emplace_back(j, j, p.springK.y);
                b0_y[j] += p.springK.y * p.pos0.y;
            }

            // Z 边界
            if (p.fixZ) {
                tripletsZ.emplace_back(j, j, kRigidFix);
                b0_z[j] += kRigidFix * p.pos0.z;
            } else if (p.springK.z > 0.0) {
                tripletsZ.emplace_back(j, j, p.springK.z);
                b0_z[j] += p.springK.z * p.pos0.z;
            }
        }

        Eigen::SparseMatrix<double> Kx(nNodes, nNodes);
        Eigen::SparseMatrix<double> Ky(nNodes, nNodes);
        Eigen::SparseMatrix<double> Kz(nNodes, nNodes);

        Kx.setFromTriplets(tripletsX.begin(), tripletsX.end());
        Ky.setFromTriplets(tripletsY.begin(), tripletsY.end());
        Kz.setFromTriplets(tripletsZ.begin(), tripletsZ.end());

        // 3. 预分解系统矩阵 (稀疏 Cholesky 分解，整个求解过程仅需一次！)
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solverX;
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solverY;
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solverZ;

        solverX.compute(Kx);
        solverY.compute(Ky);
        solverZ.compute(Kz);

        if (solverX.info() != Eigen::Success || solverY.info() != Eigen::Success || solverZ.info() != Eigen::Success) {
            std::cerr << "[PDSolver] Error: Cholesky decomposition failed!\n";
            stats.converged = false;
            return stats;
        }

        // 4. 坐标向量初始化与历史迭代记录
        Eigen::VectorXd qX(nNodes), qY(nNodes), qZ(nNodes);
        for (int j = 0; j < nNodes; ++j) {
            qX[j] = particles[j].pos0.x;
            qY[j] = particles[j].pos0.y;
            qZ[j] = particles[j].pos0.z;
        }

        Eigen::VectorXd qX_prev = qX, qY_prev = qY, qZ_prev = qZ;
        std::vector<Vector3> projP(nElems);

        // 5. 交替优化迭代循环
        int iter = 0;
        for (; iter < maxIterations; ++iter) {
            // ---- 动量外推加速 (Polyak Momentum) ----
            Eigen::VectorXd yX = qX, yY = qY, yZ = qZ;
            if (enableChebyshev && iter > 0) {
                double beta = chebyshevRho; // e.g. 0.90 ~ 0.96
                yX = qX + beta * (qX - qX_prev);
                yY = qY + beta * (qY - qY_prev);
                yZ = qZ + beta * (qZ - qZ_prev);
            }

            // ---- 局部步 (Local Projection Step) ----
            for (int i = 0; i < nElems; ++i) {
                int n1 = elemNodes[i].first;
                int n2 = elemNodes[i].second;

                double dx = yX[n2] - yX[n1];
                double dy = yY[n2] - yY[n1];
                double dz = yZ[n2] - yZ[n1];
                double len = std::sqrt(dx * dx + dy * dy + dz * dz);
                double L0 = elemLength0[i];

                if (elemIsCable[i] && len < L0) {
                    // 索构件受压松弛: 许可流形内部，投影等于当前向量 (无约束拉力)
                    projP[i] = Vector3(dx, dy, dz);
                } else {
                    if (len < 1e-12) {
                        projP[i] = Vector3(L0, 0.0, 0.0);
                    } else {
                        double scale = L0 / len;
                        projP[i] = Vector3(dx * scale, dy * scale, dz * scale);
                    }
                }
            }

            // ---- 全局步 (Global Solve Step) ----
            Eigen::VectorXd bx = b0_x;
            Eigen::VectorXd by = b0_y;
            Eigen::VectorXd bz = b0_z;

            for (int i = 0; i < nElems; ++i) {
                double w = elemWeights[i];
                int n1 = elemNodes[i].first;
                int n2 = elemNodes[i].second;
                const auto& p = projP[i];

                bx[n1] -= w * p.x;
                bx[n2] += w * p.x;

                by[n1] -= w * p.y;
                by[n2] += w * p.y;

                bz[n1] -= w * p.z;
                bz[n2] += w * p.z;
            }

            // 极速回代求解 (仅利用预分解的因子做前向/后向代入)
            Eigen::VectorXd qX_new = solverX.solve(bx);
            Eigen::VectorXd qY_new = solverY.solve(by);
            Eigen::VectorXd qZ_new = solverZ.solve(bz);

            // ---- 收敛判别 ----
            double maxDelta = 0.0;
            for (int j = 0; j < nNodes; ++j) {
                double diffX = qX_new[j] - qX[j];
                double diffY = qY_new[j] - qY[j];
                double diffZ = qZ_new[j] - qZ[j];
                double d = std::sqrt(diffX * diffX + diffY * diffY + diffZ * diffZ);
                if (d > maxDelta) maxDelta = d;
            }

            qX_prev = qX;
            qY_prev = qY;
            qZ_prev = qZ;

            qX = qX_new;
            qY = qY_new;
            qZ = qZ_new;

            if (maxDelta < tolDisplacement && iter >= 5) {
                stats.converged = true;
                break;
            }
        }

        auto tEnd = std::chrono::high_resolution_clock::now();
        stats.elapsedSimTime = std::chrono::duration<double>(tEnd - tStart).count();
        stats.iterations = iter + 1;

        // 6. 将最终坐标回填到质点，并计算最大位移
        double maxDisp = 0.0;
        for (int j = 0; j < nNodes; ++j) {
            particles[j].pos = Vector3(qX[j], qY[j], qZ[j]);
            double d = particles[j].displacement().norm();
            if (d > maxDisp) maxDisp = d;
        }
        stats.maxDisplacement = maxDisp;

        // 7. 计算构件内力与不平衡残差
        for (auto& p : particles) p.resetInternalForces();
        for (auto& elem : elements) elem->computeInternalForces(particles);

        double maxResidual = 0.0;
        for (int j = 0; j < nNodes; ++j) {
            if (particles[j].fixX && particles[j].fixY && particles[j].fixZ) continue;
            Vector3 res = particles[j].fExt - particles[j].fInt;
            if (!particles[j].fixX) res.x -= particles[j].springK.x * (particles[j].pos.x - particles[j].pos0.x);
            if (!particles[j].fixY) res.y -= particles[j].springK.y * (particles[j].pos.y - particles[j].pos0.y);
            if (!particles[j].fixZ) res.z -= particles[j].springK.z * (particles[j].pos.z - particles[j].pos0.z);
            double rNorm = res.norm();
            if (rNorm > maxResidual) maxResidual = rNorm;
        }
        stats.residualForce = maxResidual;

        return stats;
    }
};

} // namespace PD

namespace LiangYiTest {
    using PDSolver = PD::PDSolver;
    using SolverStats = PD::SolverStats;
}
