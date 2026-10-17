#include <print>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <chrono>
#include <algorithm>

#include "include/Vector3.h"
#include "include/Matrix3x3.h"
#include "include/Particle.h"
#include "include/Element.h"
#include "include/TrussElement.h"
#include "include/PDSolver.h"
#include "include/MstModel.h"
#include "include/MstParser.h"

// 结构支座反力结果项
struct ReactionResult {
    int nodeId = 0;
    std::string supportType;
    double x = 0.0, y = 0.0, z = 0.0;
    double rx = 0.0, ry = 0.0, rz = 0.0; // kN
    double norm = 0.0;                   // kN
};

// 工况描述
struct LoadCaseOption {
    std::string name;
    int category = 0; // 0=Dead, 1=Live, 2=Wind
    int caseId = 0;
    size_t count = 0;
    double sumFx = 0.0;
    double sumFy = 0.0;
    double sumFz = 0.0;
};

int main(int argc, char* argv[]) {
    std::println("=================================================================================");
    std::println("         Projective Dynamics (PD) 空间杆系结构求解器开源算例示例程序             ");
    std::println("  [理论论文]: Bouaziz et al., Projective Dynamics, ACM TOG (SIGGRAPH 2014)       ");
    std::println("  [重要声明]: 本程序仅供学术研究、教学与力学交流使用，请勿用于真实工程设计！     ");
    std::println("=================================================================================\n");

    std::string modelFile = "test.log";
    if (argc >= 2) {
        modelFile = argv[1];
    }

    std::println(">>> 正在加载 MST 空间结构模型文件: {} ...", modelFile);
    MST::Model model;
    if (!MST::Parser::parseFile(modelFile, model)) {
        std::println(stderr, "[-] 错误: 无法解析模型文件: {}", modelFile);
        return 1;
    }

    std::println("[+] 模型加载成功!");
    std::println("  • 节点总数: {} 个 (上弦: {} 个, 下弦: {} 个)",
                 model.nodes.size(),
                 std::count_if(model.nodes.begin(), model.nodes.end(), [](const auto& n) { return n.layer == 1; }),
                 std::count_if(model.nodes.begin(), model.nodes.end(), [](const auto& n) { return n.layer == 2; }));
    std::println("  • 杆件总数: {} 根", model.elements.size());
    std::println("  • 截面类型: {} 种", model.sectionProperties.size());
    std::println("  • 边界支座: {} 个\n",
                 !model.supportInfos.empty() ? model.supportInfos.size() : model.springSupports.size());

    // 1. 扫描可用荷载工况
    std::vector<LoadCaseOption> options;

    if (!model.deadLoads.empty()) {
        LoadCaseOption opt;
        opt.name = "恒载 (Dead Load)";
        opt.category = 0;
        opt.caseId = 1;
        opt.count = model.deadLoads.size();
        for (const auto& dl : model.deadLoads) {
            opt.sumFx += dl.fx; opt.sumFy += dl.fy; opt.sumFz += dl.fz;
        }
        options.push_back(opt);
    }

    std::map<int, std::vector<MST::CaseNodeLoad>> liveGroups;
    for (const auto& ll : model.liveLoads) liveGroups[ll.caseId].push_back(ll);
    for (const auto& [cid, llist] : liveGroups) {
        LoadCaseOption opt;
        opt.name = "活载工况 " + std::to_string(cid) + " (Live Load " + std::to_string(cid) + ")";
        opt.category = 1;
        opt.caseId = cid;
        opt.count = llist.size();
        for (const auto& ll : llist) {
            opt.sumFx += ll.fx; opt.sumFy += ll.fy; opt.sumFz += ll.fz;
        }
        options.push_back(opt);
    }

    std::map<int, std::vector<MST::CaseNodeLoad>> windGroups;
    for (const auto& wl : model.windLoads) windGroups[wl.caseId].push_back(wl);
    for (const auto& [cid, wlist] : windGroups) {
        LoadCaseOption opt;
        opt.name = "风载工况 " + std::to_string(cid) + " (Wind Load " + std::to_string(cid) + ")";
        opt.category = 2;
        opt.caseId = cid;
        opt.count = wlist.size();
        for (const auto& wl : wlist) {
            opt.sumFx += wl.fx; opt.sumFy += wl.fy; opt.sumFz += wl.fz;
        }
        options.push_back(opt);
    }

    if (options.empty()) {
        std::println(stderr, "[-] 错误: 模型中未发现任何有效节点荷载!");
        return 1;
    }

    std::println("---------------------------------------------------------------------------------");
    std::println("【模型荷载工况列表】");
    std::println("---------------------------------------------------------------------------------");
    std::println(" 编号 | {:<32} | 节点数 | 总荷载 ΣFx (kN) | 总荷载 ΣFy (kN) | 总荷载 ΣFz (kN)", "工况名称");
    std::println("------+----------------------------------+--------+-----------------+-----------------+-----------------");
    for (size_t i = 0; i < options.size(); ++i) {
        const auto& opt = options[i];
        std::println("  [{:>2}] | {:<32} | {:>6} | {:>15.2f} | {:>15.2f} | {:>15.2f}",
                     i + 1, opt.name, opt.count, opt.sumFx, opt.sumFy, opt.sumFz);
    }
    std::println("---------------------------------------------------------------------------------");

    int chosenIdx = 0;
    if (argc >= 3) {
        chosenIdx = std::stoi(argv[2]) - 1;
        if (chosenIdx < 0 || chosenIdx >= (int)options.size()) chosenIdx = 0;
    } else {
        // 默认推荐风荷载工况（体现空间拉索/受拉抗拔与三维反力特征），如无风载则选恒载
        for (size_t i = 0; i < options.size(); ++i) {
            if (options[i].category == 2) {
                chosenIdx = static_cast<int>(i);
                break;
            }
        }
        std::println("\n[提示] 默认选择工况 [{}] ({}) 进行求解...", chosenIdx + 1, options[chosenIdx].name);
        std::println("       如需切换工况，可在命令行运行: ./PDSolverDemo <model.log> <工况编号(1~{})>", options.size());
    }

    const auto& curOpt = options[chosenIdx];
    std::println("\n>>> 已选定工况: {} (作用节点数: {}, 合力: [{:.2f}, {:.2f}, {:.2f}] kN)",
                 curOpt.name, curOpt.count, curOpt.sumFx, curOpt.sumFy, curOpt.sumFz);

    // 2. 构建 Projective Dynamics 求解模型
    PD::PDSolver solver;
    solver.maxIterations = 500;
    solver.tolDisplacement = 1e-5; // 0.01 mm 位移收敛容差
    solver.enableChebyshev = true; // 开启动量外推加速
    solver.chebyshevRho = 0.94;    // 切比雪夫谱半径

    std::map<int, int> nodeToPIdx;
    for (size_t i = 0; i < model.nodes.size(); ++i) {
        const auto& n = model.nodes[i];
        int pId = solver.addParticle(PD::Vector3(n.x, n.y, n.z));
        nodeToPIdx[n.id] = pId;
    }

    for (const auto& elem : model.elements) {
        int p1 = nodeToPIdx[elem.n1];
        int p2 = nodeToPIdx[elem.n2];
        const auto* ai = model.getAreaInfo(elem.secId);
        double A = (ai ? ai->A : 621.25) * 1e-6; // mm^2 -> m^2
        double E = (ai ? ai->E : 206000.0) * 1e6; // MPa -> Pa (N/m^2)
        solver.addTrussElement(p1, p2, E, A);
    }

    // 3. 施加边界支座约束与外荷载
    if (!model.supportInfos.empty()) {
        for (const auto& si : model.supportInfos) {
            auto it = nodeToPIdx.find(si.nodeId);
            if (it == nodeToPIdx.end()) continue;
            int pId = it->second;
            auto& p = solver.particles[pId];

            const auto* sp = model.getSpringSupport(si.nodeId);
            double kxVal = (sp ? sp->kx : 3.0) * 1e6; // 3.0 kN/mm = 3.0e6 N/m
            double kyVal = (sp ? sp->ky : 3.0) * 1e6;

            if (si.dx == 1) { p.fixX = true; kxVal = 0.0; }
            else if (si.dx == 2) { p.fixX = false; }
            else { p.fixX = false; kxVal = 0.0; }

            if (si.dy == 1) { p.fixY = true; kyVal = 0.0; }
            else if (si.dy == 2) { p.fixY = false; }
            else { p.fixY = false; kyVal = 0.0; }

            if (si.dz == 1 || si.dz == 2) { p.fixZ = true; }
            else { p.fixZ = false; }

            p.setSpring(PD::Vector3(kxVal, kyVal, 0.0));
        }
    } else {
        for (const auto& sp : model.springSupports) {
            auto it = nodeToPIdx.find(sp.nodeId);
            if (it == nodeToPIdx.end()) continue;
            int pId = it->second;
            solver.particles[pId].fixZ = true;
            solver.particles[pId].setSpring(PD::Vector3(sp.kx * 1e6, sp.ky * 1e6, 0.0));
        }
    }

    // 施加外荷载 (kN -> N)
    if (curOpt.category == 0) {
        for (const auto& dl : model.deadLoads) {
            int pId = nodeToPIdx[dl.nodeId];
            solver.particles[pId].fExt = PD::Vector3(dl.fx * 1e3, dl.fy * 1e3, dl.fz * 1e3);
        }
    } else if (curOpt.category == 1) {
        for (const auto& ll : model.liveLoads) {
            if (ll.caseId == curOpt.caseId) {
                int pId = nodeToPIdx[ll.nodeId];
                solver.particles[pId].fExt = PD::Vector3(ll.fx * 1e3, ll.fy * 1e3, ll.fz * 1e3);
            }
        }
    } else if (curOpt.category == 2) {
        for (const auto& wl : model.windLoads) {
            if (wl.caseId == curOpt.caseId) {
                int pId = nodeToPIdx[wl.nodeId];
                solver.particles[pId].fExt = PD::Vector3(wl.fx * 1e3, wl.fy * 1e3, wl.fz * 1e3);
            }
        }
    }

    // 4. 调用 Projective Dynamics 求解器
    std::println(">>> 正在启动 Projective Dynamics 交替优化求解...");
    auto stats = solver.solve();

    std::println("\n=================================================================================");
    std::println("【PD 求解结果统计】");
    std::println("=================================================================================");
    std::println("  • 收敛状态: {}", stats.converged ? "已收敛 (Converged)" : "达到最大迭代步 (Max Steps Reached)");
    std::println("  • 迭代步数: {} 步", stats.iterations);
    std::println("  • 求解耗时: {:.2f} 毫秒 ({:.4f} 秒)", stats.elapsedSimTime * 1000.0, stats.elapsedSimTime);
    std::println("  • 最大位移: {:.3f} mm", stats.maxDisplacement * 1000.0);
    std::println("  • 最大内力残差: {:.4e} N", stats.residualForce);

    // 5. 提取构件轴力
    double maxTension = 0.0, maxCompression = 0.0;
    int maxTensionElemId = 0, maxCompElemId = 0;
    for (size_t i = 0; i < solver.elements.size(); ++i) {
        auto* truss = dynamic_cast<PD::TrussElement*>(solver.elements[i].get());
        if (!truss) continue;
        double f_kN = truss->axialForce / 1e3;
        if (f_kN > maxTension) {
            maxTension = f_kN;
            maxTensionElemId = model.elements[i].id;
        }
        if (-f_kN > maxCompression) {
            maxCompression = -f_kN;
            maxCompElemId = model.elements[i].id;
        }
    }
    std::println("  • 杆件最大拉力: {:>10.2f} kN (杆件 #{})", maxTension, maxTensionElemId);
    std::println("  • 杆件最大压力: {:>10.2f} kN (杆件 #{})", -maxCompression, maxCompElemId);

    // 6. 提取支座反力
    std::vector<ReactionResult> reactions;
    double sumRx = 0.0, sumRy = 0.0, sumRz = 0.0;

    if (!model.supportInfos.empty()) {
        for (const auto& si : model.supportInfos) {
            auto it = nodeToPIdx.find(si.nodeId);
            if (it == nodeToPIdx.end()) continue;
            int pId = it->second;
            const auto& p = solver.particles[pId];

            ReactionResult r;
            r.nodeId = si.nodeId;
            r.x = p.pos0.x; r.y = p.pos0.y; r.z = p.pos0.z;
            r.rx = (si.dx != 0) ? (p.fInt.x - p.fExt.x) / 1e3 : 0.0;
            r.ry = (si.dy != 0) ? (p.fInt.y - p.fExt.y) / 1e3 : 0.0;
            r.rz = (si.dz != 0) ? (p.fInt.z - p.fExt.z) / 1e3 : 0.0;
            r.norm = std::sqrt(r.rx * r.rx + r.ry * r.ry + r.rz * r.rz);

            if (si.dx == 1 && si.dy == 1) r.supportType = "固定铰支座";
            else if (si.dx == 0 && si.dy == 1) r.supportType = "单向滑移(X)";
            else if (si.dx == 1 && si.dy == 0) r.supportType = "单向滑移(Y)";
            else if (si.dx == 0 && si.dy == 0) r.supportType = "双向滑移(XY)";
            else r.supportType = "弹性支座";

            reactions.push_back(r);
            sumRx += r.rx; sumRy += r.ry; sumRz += r.rz;
        }
    } else {
        for (const auto& sp : model.springSupports) {
            auto it = nodeToPIdx.find(sp.nodeId);
            if (it == nodeToPIdx.end()) continue;
            int pId = it->second;
            const auto& p = solver.particles[pId];

            ReactionResult r;
            r.nodeId = sp.nodeId;
            r.supportType = "弹性支座";
            r.x = p.pos0.x; r.y = p.pos0.y; r.z = p.pos0.z;
            r.rx = (p.fInt.x - p.fExt.x) / 1e3;
            r.ry = (p.fInt.y - p.fExt.y) / 1e3;
            r.rz = (p.fInt.z - p.fExt.z) / 1e3;
            r.norm = std::sqrt(r.rx * r.rx + r.ry * r.ry + r.rz * r.rz);

            reactions.push_back(r);
            sumRx += r.rx; sumRy += r.ry; sumRz += r.rz;
        }
    }

    std::println("\n---------------------------------------------------------------------------------");
    std::println("【边界支座反力明细表】");
    std::println("---------------------------------------------------------------------------------");
    std::println(" 序号 | 节点号 | 支座类型     | 位置 (X, Y, Z) m     | 反力 Rx(kN) | 反力 Ry(kN) | 反力 Rz(kN) | 空间合力(kN)");
    std::println("------+--------+--------------+----------------------+-------------+-------------+-------------+-------------");
    for (size_t i = 0; i < reactions.size(); ++i) {
        const auto& r = reactions[i];
        std::println("  [{:>2}] | {:>6} | {:<12} | ({:>4.1f}, {:>4.1f}, {:>4.1f}) | {:>11.2f} | {:>11.2f} | {:>11.2f} | {:>12.2f}",
                     i + 1, r.nodeId, r.supportType, r.x, r.y, r.z, r.rx, r.ry, r.rz, r.norm);
    }
    std::println("------+--------+--------------+----------------------+-------------+-------------+-------------+-------------");
    std::println(" 合计 | 支座总反力向量 ΣR                    | {:>11.2f} | {:>11.2f} | {:>11.2f} |",
                 sumRx, sumRy, sumRz);
    std::println("---------------------------------------------------------------------------------");

    // 7. 全局静力平衡自校验
    double errX = std::abs(curOpt.sumFx + sumRx);
    double errY = std::abs(curOpt.sumFy + sumRy);
    double errZ = std::abs(curOpt.sumFz + sumRz);
    double totalApplied = std::sqrt(curOpt.sumFx * curOpt.sumFx + curOpt.sumFy * curOpt.sumFy + curOpt.sumFz * curOpt.sumFz);
    double totalResidual = std::sqrt(errX * errX + errY * errY + errZ * errZ);
    double relErrorPct = (totalApplied > 1e-3) ? (totalResidual / totalApplied * 100.0) : 0.0;

    std::println("\n=================================================================================");
    std::println("【结构全局静态力学自平衡校验】");
    std::println("=================================================================================");
    std::println("  • X 向平衡: ΣFx ({:+.2f}) + ΣRx ({:+.2f}) = {:+.4f} kN (残差: {:.4f} kN)", curOpt.sumFx, sumRx, curOpt.sumFx + sumRx, errX);
    std::println("  • Y 向平衡: ΣFy ({:+.2f}) + ΣRy ({:+.2f}) = {:+.4f} kN (残差: {:.4f} kN)", curOpt.sumFy, sumRy, curOpt.sumFy + sumRy, errY);
    std::println("  • Z 向平衡: ΣFz ({:+.2f}) + ΣRz ({:+.2f}) = {:+.4f} kN (残差: {:.4f} kN)", curOpt.sumFz, sumRz, curOpt.sumFz + sumRz, errZ);
    std::println("  • 空间矢量总合力残差: {:.4f} kN (相对误差: {:.4f} %)", totalResidual, relErrorPct);
    if (relErrorPct < 0.1) {
        std::println("  • 自平衡检验结论: [优秀] 全局平衡高精度闭合！");
    } else if (relErrorPct < 1.0) {
        std::println("  • 自平衡检验结论: [良好] 满足工程精度要求。");
    } else {
        std::println("  • 自平衡检验结论: [需关注] 残差偏大，建议增加最大迭代步数。");
    }
    std::println("=================================================================================\n");

    return 0;
}
